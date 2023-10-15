#!/bin/bash
#
# piRNA length distribution
#

refdir="$(pwd)/data/mm10" # Reference genome/annotation directory; full path
sampledir=$(pwd) # Directory with preprocessed and mapped samples

threads=12
rnd=$RANDOM

mkdir -p pirna-length/run
cd pirna-length/

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

# get piRNA mappings
for name in "${samples[@]}"; do
    sdir=$name
    mkdir -p $sdir

    # All alignments
    echo -n "
    samtools view -@ 1 -b -F 4 -F 256 -F 1024 -F 2048 $sampledir/$name/align/reads.1.Aligned.final.bam \
        | bedtools bamtobed -i stdin \
        | sed '1i chr\tstart\tend\tname\tscore\tstrand' | gzip -c > $sdir/reads.1.genome.bed.gz"
done > run/cmds.$rnd.txt
cat run/cmds.$rnd.txt | rush '{}' -j $threads --verbose && rm run/cmds.$rnd.txt

# get  overlaps
clusters="$refdir/pirna-clusters.bed"
repeats="$refdir/rmsk.bed"

for name in "${samples[@]}"; do
    sdir=$name
    mkdir $sdir

    # clusters
    echo -n "
    bedtools intersect \
        -wao \
        -f 1.0 \
        -a $sdir/reads.1.genome.bed.gz -b $clusters \
        | grep -P -v \"\\t0$\" \
        | gzip -c \
        > $sdir/reads.1.genome.pirna_clusters.bed.gz" # Unstranded because the annotation is unstranded
done > run/cmds.$rnd.txt
cat run/cmds.$rnd.txt | rush '{}' -j $threads --verbose && rm run/cmds.$rnd.txt

# remove repeats overlaps
for name in "${samples[@]}"; do
    sdir=$name
    mkdir $sdir

    echo -n "
    bedtools subtract \
        -A \
        -a $sdir/reads.1.genome.pirna_clusters.bed.gz -b $repeats \
        | gzip -c \
        > $sdir/reads.1.genome.pirna_clusters.x_rpmk.bed.gz"
done > run/cmds.$rnd.txt
cat run/cmds.$rnd.txt | rush '{}' -j $threads --verbose && rm run/cmds.$rnd.txt

# add sample name
header="chr\tstart\tend\tname\tscore\tstrand\tchr_ovl\tstart_ovl\tend_ovl\tname_ovl\tscore_ovl\tstrand_ovl\toverlap_length"

for name in "${samples[@]}"; do
    sdir=$name
    mkdir $sdir

    for bed in $sdir/reads.1.genome.*.*.bed.gz; do
        echo -e $header > ${bed%.gz}.tmp
        zcat $bed >> ${bed%.gz}.tmp
        cat ${bed%.gz}.tmp \
            | table-paste-col \
                --table - \
                --col-name Sample --col-val $name \
            | gzip -c \
            > $bed && rm ${bed%.gz}.tmp &
    done
    wait
done

# plot lengths - all mappings
table-cat \
    <(zcat mmu.RIPSeq.MiwiIP.MiwiHet.P24.3/reads.1.genome.bed.gz | table-paste-col \
                --table - \
                --col-name Sample --col-val "mmu.RIPSeq.MiwiIP.MiwiHet.P24.2") \
    <(zcat mmu.RIPSeq.MiwiIP.MiwiRK.P24.3/reads.1.genome.bed.gz | table-paste-col \
                --table - \
                --col-name Sample --col-val "mmu.RIPSeq.MiwiIP.MiwiRK.P24.2") \
    <(zcat mmu.RIPSeq.MiliIP.MiwiHet.P24.1/reads.1.genome.bed.gz | table-paste-col \
                --table - \
                --col-name Sample --col-val "mmu.RIPSeq.MiliIP.MiwiHet.P24.1") \
    <(zcat mmu.RIPSeq.MiliIP.MiwiRK.P24.1/reads.1.genome.bed.gz | table-paste-col \
                --table - \
                --col-name Sample --col-val "mmu.RIPSeq.MiliIP.MiwiRK.P24.1") \
| src/R/length-distro.R --input stdin --ofile length-distro.pdf --from 24 --to 32 --by 1 &

# plot lengths - piRNA cluster mappings
table-cat \
    <(zcat mmu.RIPSeq.MiwiIP.MiwiHet.P24.3/reads.1.genome.pirna_clusters.x_rpmk.bed.gz) \
    <(zcat mmu.RIPSeq.MiwiIP.MiwiRK.P24.3/reads.1.genome.pirna_clusters.x_rpmk.bed.gz) \
    <(zcat mmu.RIPSeq.MiliIP.MiwiHet.P24.1/reads.1.genome.pirna_clusters.x_rpmk.bed.gz) \
    <(zcat mmu.RIPSeq.MiliIP.MiwiRK.P24.1/reads.1.genome.pirna_clusters.x_rpmk.bed.gz) \
| src/R/length-distro.R --input stdin --ofile length-distro-pirna_clusters.x_rpmk.pdf --from 24 --to 32 --by 1 --unstranded &
wait