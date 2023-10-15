#!/bin/bash
#
# Analyze piRNA ping-pong at selected regions
#

refdir="$(pwd)/data/mm10" # Reference genome/annotation directory; full path
sampledir=$(pwd) # Directory with preprocessed and mapped samples

threads=12
rnd=$RANDOM

mkdir -p pingpong/run
cd pingpong/

################################################################################
samples=(
    "mmu.RIPSeq.MiwiIP.MiwiHet.P24.2"
    "mmu.RIPSeq.MiwiIP.MiwiRK.P24.2"
    "mmu.RIPSeq.MiwiIP.MiwiHet.P24.3"
    "mmu.RIPSeq.MiwiIP.MiwiRK.P24.3"
    "mmu.RIPSeq.MiliIP.MiwiHet.P24.1"
    "mmu.RIPSeq.MiliIP.MiwiRK.P24.1"
    "mmu.RIPSeq.MiliIP.MiwiHet.P24.2"
    "mmu.RIPSeq.MiliIP.MiwiRK.P24.2"
)

# BAM to BED
for name in "${samples[@]}"; do
    sdir=$name
    mkdir -p $sdir

    echo -n "
    samtools view -@ 1 -b -F 4 -F 256 -F 1024 -F 2048 $sampledir/$name/align/reads.1.Aligned.final.bam \
        | bedtools bamtobed -i stdin \
        | awk 'BEGIN{FS=\"\t\";OFS=\"\t\"}{print \$1,\$2,\$3,\$4,\$5,\$6,\$3-\$2}' \
        | gzip -c > $sdir/reads.1.genome.bed.gz"
done > run/cmds.$rnd.txt
cat run/cmds.$rnd.txt | rush '{}' -j $threads --verbose && rm run/cmds.$rnd.txt 

# split +/- mapped; get only starts
for name in "${samples[@]}"; do
    sdir=$name
    mkdir -p $sdir

    echo -n "
    zcat $sdir/reads.1.genome.bed.gz | grep -P \"\t\+\t\" | awk 'BEGIN{FS=\"\t\";OFS=\"\t\"}{print \$1,\$2,\$2+1,\$4,\$5,\$6,\$7}' \
        | gzip -c > $sdir/reads.1.genome-plus.start.bed.gz"
    echo -n "
    zcat $sdir/reads.1.genome.bed.gz | grep -P \"\t-\t\"| awk 'BEGIN{FS=\"\t\";OFS=\"\t\"}{print \$1,\$3-1,\$3,\$4,\$5,\$6,\$7}' \
        | gzip -c > $sdir/reads.1.genome-minus.start.bed.gz"
done > run/cmds.$rnd.txt
cat run/cmds.$rnd.txt | rush '{}' -j $threads --verbose && rm run/cmds.$rnd.txt 

# bedtools window (30 nt)
# should be reciprocal since we have -sw and the pairs are "unique" = there is only one match of + and - reads in the ping-pong pair
wind=30
for name in "${samples[@]}"; do
    sdir=$name
    mkdir -p $sdir

    echo -n "
    bedtools window -l 0 -r $wind \
        -Sm -sw \
        -a $sdir/reads.1.genome-plus.start.bed.gz -b $sdir/reads.1.genome-minus.start.bed.gz \
        | gzip -c > $sdir/reads.1.genome-plusminus.start.window$wind.bed.gz"
done > run/cmds.$rnd.txt
cat run/cmds.$rnd.txt | rush '{}' -j $threads --verbose && rm run/cmds.$rnd.txt 

# overlap with piRNA clusters
clusters=$refdir/pirna-clusters.bed

for name in "${samples[@]}"; do
    sdir=$name
    mkdir $sdir

    echo -n "
    bedtools intersect \
        -wa \
        -f 1.0 \
        -u \
        -a $sdir/reads.1.genome-plusminus.start.window$wind.bed.gz -b $clusters \
        | grep -P -v \"\\t0$\" \
        | gzip -c \
        > $sdir/reads.1.genome-plusminus.start.window$wind.pirna_clusters.bed.gz" # Unstranded because the annotation is unstranded
done > run/cmds.$rnd.txt
cat run/cmds.$rnd.txt | rush '{}' -j $threads --verbose && rm run/cmds.$rnd.txt 

# remove repeats from overlaps
repeats==$refdir/rmsk.bed

for name in "${samples[@]}"; do
    sdir=$name
    mkdir $sdir

    echo -n "
    bedtools subtract \
        -A \
        -a $sdir/reads.1.genome-plusminus.start.window$wind.pirna_clusters.bed.gz -b $repeats \
        | gzip -c \
        > $sdir/reads.1.genome-plusminus.start.window$wind.pirna_clusters.x_rpmk.bed.gz"
done > run/cmds.$rnd.txt
cat run/cmds.$rnd.txt | rush '{}' -j $threads --verbose && rm run/cmds.$rnd.txt 

# add sample name
header="chr\tstart\tend\tname\tscore\tstrand\tlength\tchr_ovl\tstart_ovl\tend_ovl\tname_ovl\tscore_ovl\tstrand_ovl\tlength_ovl"

for name in "${samples[@]}"; do
    sdir=$name
    mkdir $sdir

    for bed in $(find $sdir -regextype posix-extended -regex '.*(pirna_clusters).*.bed.gz'); do
        echo -e $header > ${bed%.gz}.tmp
        zcat $bed >> ${bed%.gz}.tmp
        cat ${bed%.gz}.tmp \
            | table-paste-col \
                --table - \
                --col-name Sample --col-val $name \
                --at-end \
            | gzip -c \
            > $bed && rm ${bed%.gz}.tmp &
    done
    wait
done

# make plots
for name in "${samples[@]}"; do
    sdir=$name
    mkdir $sdir

    # Only uniquely mapped reads
    zcat $sdir/reads.1.genome-plusminus.start.window30.pirna_clusters.x_rpmk.bed.gz \
        | src/R/pingpong.R --input stdin --output "$sdir/reads.1.genome-plusminus.start.window30.pirna_clusters.x_rpmk.pp.pdf" &
    wait

    # # Including multimapped
    # zcat $sdir/reads.1.genome-plusminus.start.window30.pirna_clusters.x_rpmk.bed.gz \
    #     | src/R/pingpong.R --input stdin --multi --output "$sdir/reads.1.genome-plusminus.start.window30.pirna_clusters.x_rpmk.includingMultimapped.pp.pdf" &
    # wait
done
