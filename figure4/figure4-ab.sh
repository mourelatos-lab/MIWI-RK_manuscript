#!/bin/bash
#
# Calculate theoretical piRNA targeting
#
# IMPORTANT: piRNA targeting MUST be run after RNA expression. It uses RNA expression to normalize degradome "occupancy" 
#

refdir="$(pwd)/data/mm10" # Reference genome/annotation directory; full path
sampledir=$(pwd) # Directory with preprocessed and mapped samples

threads=12
rnd=$RANDOM

mkdir -p pirna-target/run
cd pirna-target/

################################################################################
samples=(
    "mmu.PARESeq.polya.MiwiHet.P24.1"
    "mmu.PARESeq.polya.MiwiHet.P24.2"
    "mmu.PARESeq.polya.MiwiRK.P24.1"
    "mmu.PARESeq.polya.MiwiRK.P24.2"
)

# extract degradome 5' positions
# Important: Assumes the data are PE and sense (forward) strand specific
for name in "${samples[@]}"; do
    sdir=$name
    mkdir -p $sdir

    echo -n "
    samtools view -@ 1 -b -F 4 -F 256 -F 2048 -f 1 -f 64 -F 16 -f 32 $sampledir/$name/align/reads.1.Aligned.final.bam \
        | bedtools bamtobed -i stdin \
        | sed '1i chr\tstart\tend\tname\tscore\tstrand' \
        | gzip -c > $sdir/reads.1.genome.bed.gz && \
    samtools view -@ 1 -b -F 4 -F 256 -F 2048 -f 1 -f 64 -f 16 -F 32 $sampledir/$name/align/reads.1.Aligned.final.bam \
        | bedtools bamtobed -i stdin \
        | gzip -c >> $sdir/reads.1.genome.bed.gz"
done > run/cmds.$rnd.txt
cat run/cmds.$rnd.txt | rush '{}' -j $threads --verbose && rm run/cmds.$rnd.txt

# extend 5' position for GTBuster 
extend_degra=50

for name in "${samples[@]}"; do
    sdir=$name
    mkdir -p $sdir

    echo -n "
    zcat $sdir/reads.1.genome.bed.gz \
        | awk -v ext=\"$extend_degra\" 'BEGIN{FS=\"\t\"; OFS=\"\t\"}{if(\$6==\"+\") print \$1,\$2-ext,\$2+1+ext,\$4,\$5,\$6; else if(\$6==\"-\") print \$1,\$3-1-ext,\$3+ext,\$4,\$5,\$6}' \
        | awk 'BEGIN{FS=\"\t\"; OFS=\"\t\"}{if (\$2<0) {\$2=0} print \$0}' \
        | gzip -c > $sdir/reads.1.5p-ext${extend_degra}.bed.gz"
done > run/cmds.$rnd.txt
cat run/cmds.$rnd.txt | rush '{}' -j $threads --verbose && rm run/cmds.$rnd.txt

# extract unique mapping sequences, summarize positions (collapse/count) into score column, sort
for name in "${samples[@]}"; do
    sdir=$name
    mkdir -p $sdir

    echo -n "
    zcat $sdir/reads.1.5p-ext${extend_degra}.bed.gz \
        | ./src/R/count-bed.R | sort -T $sdir --parallel=2 -k 1,1 -k2,2n \
        | gzip -c > $sdir/reads.1.5p-ext${extend_degra}.collapsed.bed.gz"
done > run/cmds.$rnd.txt
cat run/cmds.$rnd.txt | rush '{}' -j $threads --verbose && rm run/cmds

# get FASTA tab
for name in "${samples[@]}"; do
    sdir=$name
    mkdir -p $sdir

    echo -n "
    bedtools getfasta -s -nameOnly -tab -fi $refdir/genome.fa -bed <(zcat $sdir/reads.1.5p-ext${extend_degra}.collapsed.bed.gz) \
        | sed 's/(+)//g; s/(-)//g' \
        | python3 src/Python/dedup-tab.py "stdin" "_" \
        | sed '1i name\tseq_genome' \
        | gzip -c > $sdir/reads.1.5p-ext${extend_degra}.collapsed.txt.gz"
done > run/cmds.$rnd.txt
cat run/cmds.$rnd.txt | rush '{}' -j $threads --verbose && rm run/cmds.$rnd.txt

# find common sequences between replicates
pairs=(
    "mmu.PARESeq.polya.MiwiHet.P24.1"   "mmu.PARESeq.polya.MiwiHet.P24.2"
    "mmu.PARESeq.polya.MiwiRK.P24.1"    "mmu.PARESeq.polya.MiwiRK.P24.2"
)
tLen=${#pairs[@]}

for (( i=0; i<${tLen}; i++ )); do
    n1=${pairs[$i]}
    ((i++))
    n2=${pairs[$i]}

    name=${n1%.*}
    sdir=$name
    mkdir -p $sdir

    col2=`zcat $n1/reads.1.5p-ext${extend_degra}.collapsed.txt.gz | head -1 | cut -f2` #
    echo -e "${n1}|${n2}\t${col2}" | gzip -c > $sdir/reads.1.5p-ext${extend_degra}.collapsed.txt.gz # Copy header

    zcat $n1/reads.1.5p-ext${extend_degra}.collapsed.txt.gz $n2/reads.1.5p-ext${extend_degra}.collapsed.txt.gz \
        | grep -v -P "name\tseq_genome" \
        | python3 src/Python/dedup-tab.py "stdin" "|" \
        | grep "|" \
        | gzip -c >> $sdir/reads.1.5p-ext${extend_degra}.collapsed.txt.gz &
done
wait

samples=(
    "mmu.PARESeq.polya.MiwiHet.P24"
    "mmu.PARESeq.polya.MiwiRK.P24"
)

# convert TAB to FASTA
for name in "${samples[@]}"; do
    sdir=$name
    mkdir -p $sdir/ref

    # IMPORTANT - We remove sites with less than 1 RPM (grep -v ",0.")
    zcat $sdir/reads.1.5p-ext${extend_degra}.collapsed.txt.gz | tail -n+2 \
        | grep -v ",0." \
        | sed 's/^/>/g' \
        | tr '\t' '\n' | gzip -c > $sdir/ref/reads.1.5p-ext${extend_degra}.collapsed.fa.gz &
done
wait

# chunk FASTA
lines=8000 # 4000 fasta sequences (header+sequence)

for name in "${samples[@]}"; do
    echo "Working on $name"
    sdir=$name
    mkdir -p $sdir/ref

    split -l $lines --additional-suffix=".fa" <(zcat $sdir/ref/reads.1.5p-ext${extend_degra}.collapsed.fa.gz) $name/ref/p
done

## get piRNA distro
trim_len=25 # How much to trim piRNA sequences for counting/merging/...

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

# calculate piRNA coverage
for name in "${samples[@]}"; do
    sdir=$name
    mkdir -p $sdir

    echo -e "chr\tstart\tend\tname\tscore\tstrand" | gzip -c > $sdir/reads.1.genome.bed.gz

    echo -n "
    samtools view -@ 1 -b -F 4 -F 256 -F 1024 -F 2048 $SAMPLE_DIR/$name/align/reads.1.Aligned.final.bam \
        | bedtools bamtobed -i stdin \
        | gzip -c >> $sdir/reads.1.genome.bed.gz" # sed '1i chr\tstart\tend\tname\tscore\tstrand' |
done > run/cmds.$rnd.txt
cat run/cmds.$rnd.txt | rush '{}' -j $threads --verbose && rm run/cmds.$rnd.txt

# extract piRNA sequences from BAM WITHOUT softclipping
# IMPORTANT: REMOVE softclipping and then get the sequence - BAM has the whole sequence inclucing the softclipped part
# IMPORTANT: Sequence of reads mapping to the - strand has to be reverse-complemented - BAM has matching sequence as they look on + strand
# IMPORTANT: Assuming the sequencing is sense strand/forward strand specific
for name in "${samples[@]}"; do
    sdir=$name
    mkdir -p $sdir

    echo -e "name\tseq_reads" | gzip -c > $sdir/reads.1.genome.sequence-reads.txt.gz

    echo -n "
    java -jar remove-softlip.jar --samoutputformat BAM \
        $SAMPLE_DIR/$name/align/reads.1.Aligned.final.bam > $sdir/tmp.sf.bam; \
    samtools view -@ 1 -F 4 -F 256 -F 1024 -F 2048 -F 16 $sdir/tmp.sf.bam \
        | cut -f1,10  \
        | gzip -c >> $sdir/reads.1.genome.sequence-reads.txt.gz; \
    samtools view -@ 1 -F 4 -F 256 -F 1024 -F 2048 -f 16 $sdir/tmp.sf.bam \
        | cut -f1,10  \
        | tr '\t' '\n' | while read L; do  echo -n \$L; echo -ne \"\t\"; read L; echo \"\$L\" | rev | tr \"ATGCNatcgn\" \"TACGNtacgn\"; done \
        | gzip -c >> $sdir/reads.1.genome.sequence-reads.txt.gz && rm $sdir/tmp.sf.bam"
done > run/cmds.$rnd.txt
cat run/cmds.$rnd.txt | rush '{}' -j 4 --verbose && rm run/cmds.$rnd.txt

# get sequences from the reference genome
for name in "${samples[@]}"; do
    sdir=$name
    mkdir -p $sdir

    echo -n "
    bedtools getfasta -s -name -tab -fi $refdir/genome.fa -bed <(zcat $sdir/reads.1.genome.bed.gz | tail -n+2) \
    | sed 's/(+)//g; s/(-)//g' \
    | sed '1i name\tseq_genome' | gzip -c > $sdir/reads.1.genome.sequence-genome.txt.gz"
done > run/cmds.$rnd.txt
cat run/cmds.$rnd.txt | rush '{}' -j $threads --verbose && rm run/cmds.$rnd.txt

echo "> ADD SEQUENCES TO POSITIONS - READS W/O SOFTCLIP <"

for name in "${samples[@]}"; do
    sdir=$name
    mkdir -p $sdir

    echo -n "
    table-join --left \
        --table1 <(zcat $sdir/reads.1.genome.bed.gz) \
        --key1 'name' \
        --table2 <(zcat $sdir/reads.1.genome.sequence-reads.txt.gz) \
        | gzip -c > $sdir/reads.1.genome.bed.tmp.gz && mv $sdir/reads.1.genome.bed.tmp.gz $sdir/reads.1.genome.bed.gz"
done > run/cmds.$rnd.txt
cat run/cmds.$rnd.txt | rush '{}' -j 2 --verbose && rm run/cmds.$rnd.txt

# add sequences to positions
for name in "${samples[@]}"; do
    sdir=$name
    mkdir -p $sdir

    echo -n "
    table-join --left \
        --table1 <(zcat $sdir/reads.1.genome.bed.gz) \
        --key1 'name' \
        --table2 <(zcat $sdir/reads.1.genome.sequence-genome.txt.gz) \
        | gzip -c > $sdir/reads.1.genome.bed.tmp.gz && mv $sdir/reads.1.genome.bed.tmp.gz $sdir/reads.1.genome.bed.gz"
done > run/cmds.$rnd.txt
cat run/cmds.$rnd.txt | rush '{}' -j 2 --verbose && rm run/cmds.$rnd.txt

# remove tRNA and rRNA
cat $refdir/rmsk_categ.tab \
    | src/misc/fix-ucsc-rmsk-categ-bed.pl \
    > $refdir/rmsk_categ.bed
grep -P "\ttRNA\t" $refdir/rmsk_categ.bed > $refdir/tRNA.bed
cat $refdir/tRNA.bed $refdir/rRNA.bed > $refdir/rRNA_tRNA.bed

for name in "${samples[@]}"; do
    sdir=$name
    mkdir -p $sdir

    echo -e "chr\tstart\tend\tname\tscore\tstrand\tseq_reads\tseq_genome" > $sdir/reads.1.genome.clean.bed
    echo -n "
    zcat $sdir/reads.1.genome.bed.gz | tail -n+2 \
        | awk 'BEGIN{FS=\"\t\";OFS=\"\t\"}{print \$2,\$3,\$4,\$1,\$5,\$6,\$7,\$8}' \
        | bedtools subtract -s -a stdin -b $refdir/rRNA_tRNA.bed \
        >> $sdir/reads.1.genome.clean.bed && gzip $sdir/reads.1.genome.clean.bed"
done  > run/cmds.$rnd.txt
cat run/cmds.$rnd.txt | rush '{}' -j $threads --verbose && rm run/cmds.$rnd.txt

# summarize positions
for name in "${samples[@]}"; do
    sdir=$name
    mkdir -p $sdir

    # Counts by position
    echo -n "
    zcat $sdir/reads.1.genome.clean.bed.gz \
        | ./src/R/pirna-distro-counts.R --ifile stdin --unique --ofile $sdir/reads.1.genome"

    # Counts by sequence
    echo -n "
    zcat $sdir/reads.1.genome.clean.bed.gz \
        | ./src/R/pirna-distro-counts.R --ifile stdin --unique --type seq --genome --trimseq 25 --ofile $sdir/reads.1.genome.seqs-genome"
done > run/cmds.$rnd.txt
cat run/cmds.$rnd.txt | rush '{}' -j 6 --verbose -c -C run/cmds.$rnd.txt.done && rm run/cmds.$rnd.txt run/cmds.$rnd.txt.done

# add sample names
for name in "${samples[@]}"; do
    sdir=$name
    mkdir -p $sdir

    cat $sdir/reads.1.genome.counts-starts.uniq.bed | \
        paste-table-col --ifile - --col-name sample --col-val $name \
        > $sdir/reads.1.genome.counts-starts.uniq.bed.tmp && mv $sdir/reads.1.genome.counts-starts.uniq.bed.tmp $sdir/reads.1.genome.counts-starts.uniq.bed &

    cat $sdir/reads.1.genome.seqs-genome.counts.trim25.uniq.tsv | \
        paste-table-col --ifile - --col-name sample --col-val $name \
        > $sdir/reads.1.genome.seqs-genome.counts.trim25.uniq.tsv.tmp && mv $sdir/reads.1.genome.seqs-genome.counts.trim25.uniq.tsv.tmp $sdir/reads.1.genome.seqs-genome.counts.trim25.uniq.tsv &
done
wait

# filter results
for name in "${samples[@]}"; do
    echo "Working on $name"
    sdir=$name
    mkdir -p $sdir

    ./src/R/pirna-distro-filter.R \
        --ifile1 $sdir/reads.1.genome.seqs-genome.counts.trim25.uniq.tsv \
        --type "seq" \
        --filtcommon FALSE \
        --odir $sdir \
        &> $sdir/reads.1.genome.seq-genome.pirna-filt.txt &

    mkdir -p $sdir
    ./src/R/pirna-distro-filter.R \
        --ifile1 $sdir/reads.1.genome.counts-starts.uniq.bed \
        --type "pos" \
        --filtcommon FALSE \
        --odir $sdir \
        &> $sdir/reads.1.genome.pos.pirna-filt.txt &
done
wait

# simple sample comparison
# IMPORTANT: Biological replicates should be compared together
# Biological replicates BUT we don't require piRNA to be present in both - we need all reads for RPM normalization later on! Don't discard ANY piRNAs
pairs=(
    "mmu.RIPSeq.MiwiIP.MiwiHet.P24.3"    "mmu.RIPSeq.MiwiIP.MiwiHet.P24.2"
    "mmu.RIPSeq.MiwiIP.MiwiRK.P24.3"    "mmu.RIPSeq.MiwiIP.MiwiRK.P24.2"
    "mmu.RIPSeq.MiliIP.MiwiRK.P24.1"    "mmu.RIPSeq.MiliIP.MiwiRK.P24.2"
    "mmu.RIPSeq.MiliIP.MiwiHet.P24.1"   "mmu.RIPSeq.MiliIP.MiwiHet.P24.2"
)

tLen=${#pairs[@]}

for (( i=0; i<${tLen}; i++ )); do
    n1=${pairs[$i]}
    ((i++))
    n2=${pairs[$i]}

    # Dual mode
    echo "Working on $n1 vs $n2"

    sdir1=${n1}
    sdir2=${n2}

    sdir=compare/sequence-genome/uniq
    mkdir -p $sdir
    ./src/R/pirna-distro-filter.R \
        --ifile1 $sdir1/reads.1.genome.seqs-genome.counts.trim25.uniq.tsv \
        --ifile2 $sdir2/reads.1.genome.seqs-genome.counts.trim25.uniq.tsv \
        --type "seq" \
        --filtcommon FALSE \
        --odir $sdir \
        &> $sdir/${n1}-${n2}.seq-genome.pirna-compare.txt &

    sdir=compare/starts/uniq
    mkdir -p $sdir
    ./src/R/pirna-distro-filter.R \
        --ifile1 $sdir1/reads.1.genome.counts-starts.uniq.bed \
        --ifile2 $sdir2/reads.1.genome.counts-starts.uniq.bed \
        --type "pos" \
        --filtcommon FALSE \
        --odir $sdir \
        &> $sdir/${n1}-${n2}.pos.pirna-compare.txt &
done
wait

# comparisons between conditions - not requiring common sequences between the samples (="replicates")
pairs=(
    "mmu.RIPSeq.MiwiIP.MiwiHet.P24.2"    "mmu.RIPSeq.MiwiIP.MiwiRK.P24.2"
    "mmu.RIPSeq.MiwiIP.MiwiHet.P24.3"    "mmu.RIPSeq.MiwiIP.MiwiRK.P24.3"
    "mmu.RIPSeq.MiliIP.MiwiHet.P24.1"    "mmu.RIPSeq.MiliIP.MiwiRK.P24.1"
    "mmu.RIPSeq.MiliIP.MiwiHet.P24.2"    "mmu.RIPSeq.MiliIP.MiwiRK.P24.2"
)

tLen=${#pairs[@]}

for (( i=0; i<${tLen}; i++ )); do
    n1=${pairs[$i]}
    ((i++))
    n2=${pairs[$i]}

    # Dual mode
    echo "Working on $n1 vs $n2"

    sdir1=${n1}
    sdir2=${n2}

    sdir=compare/sequence-genome/uniq/between-conds
    mkdir -p $sdir
    ./src/R/pirna-distro-filter.R \
        --ifile1 $sdir1/reads.1.genome.seqs-genome.counts.trim25.uniq.tsv \
        --ifile2 $sdir2/reads.1.genome.seqs-genome.counts.trim25.uniq.tsv \
        --type "seq" \
        --filtcommon FALSE \
        --odir $sdir &> $sdir/${n1}-${n2}.seq-genome.pirna-compare.txt &

    sdir=compare/starts/uniq/between-conds
    mkdir -p $sdir
    ./src/R/pirna-distro-filter.R \
        --ifile1 $sdir1/reads.1.genome.counts-starts.uniq.bed \
        --ifile2 $sdir2/reads.1.genome.counts-starts.uniq.bed \
        --type "pos" \
        --filtcommon FALSE \
        --odir $sdir &> $sdir/${n1}-${n2}.pos.pirna-compare.txt &
done
wait

# join tables - starts
samples=(
    "mmu.RIPSeq.MiwiIP.MiwiHet.P24.2"
    "mmu.RIPSeq.MiwiIP.MiwiHet.P24.3"
    "mmu.RIPSeq.MiwiIP.MiwiRK.P24.2"
    "mmu.RIPSeq.MiwiIP.MiwiRK.P24.3"
    "mmu.RIPSeq.MiliIP.MiwiHet.P24.1"
    "mmu.RIPSeq.MiliIP.MiwiHet.P24.2"
    "mmu.RIPSeq.MiliIP.MiwiRK.P24.1"
    "mmu.RIPSeq.MiliIP.MiwiRK.P24.2"
)

# extract positions - starts
sdir=compare/starts/uniq

for i in */pirna-counts.filt.pos.tsv; do
    tail -n+2 $i | cut -f 1,3-6 > $i.tmp
done

echo -e "name\tchr\tstart\tend\tstrand" > $sdir/merged.pirna-counts.filt.coords.txt
cat */pirna-counts.filt.pos.tsv.tmp | sort -T . --parallel=$threads | uniq \
    >> $sdir/merged.pirna-counts.filt.coords.txt && rm */pirna-counts.filt.pos.tsv.tmp

# merge counts - starts
for name in "${samples[@]}"; do
    echo "Working on $name"
    sdir=$name
    mkdir -p $sdir

    cat $sdir/pirna-counts.filt.pos.tsv | cut -f1-2 > $sdir/pirna-counts.filt.pos.tsv.tmp
done

sdir=compare/starts/uniq

cond="Miwi"
table-join \
    --full \
    --table1 mmu.RIPSeq.MiwiIP.MiwiHet.P24.1/pirna-counts.filt.pos.tsv.tmp  \
    --key1 'name' \
    --table2 mmu.RIPSeq.MiwiIP.MiwiHet.P24.2/pirna-counts.filt.pos.tsv.tmp \
    | table-join \
        --full \
        --table1 - \
        --key1 'name' \
        --table2 mmu.RIPSeq.MiwiIP.MiwiHet.P24.3/pirna-counts.filt.pos.tsv.tmp \
    | table-join \
        --full \
        --table1 - \
        --key1 'name' \
        --table2 mmu.RIPSeq.MiwiIP.MiwiRK.P24.1/pirna-counts.filt.pos.tsv.tmp \
    | table-join \
        --full \
        --table1 - \
        --key1 'name' \
        --table2 mmu.RIPSeq.MiwiIP.MiwiRK.P24.2/pirna-counts.filt.pos.tsv.tmp \
    | table-join \
        --full \
        --table1 - \
        --key1 'name' \
        --table2 mmu.RIPSeq.MiwiIP.MiwiRK.P24.3/pirna-counts.filt.pos.tsv.tmp \
    > $sdir/merged-${cond}.pirna-counts.filt.pos.tsv && rm mmu.RIPSeq.MiwiIP.Miwi*.P24.[1,2,3]/pirna-counts.filt.pos.tsv.tmp

table-join --left \
    --table1 $sdir/merged-${cond}.pirna-counts.filt.pos.tsv \
    --key1 'name' \
    --table2 $sdir/merged.pirna-counts.filt.coords.txt \
    > $sdir/merged-${cond}.pirna-counts.filt.pos.tsv.tmp && mv $sdir/merged-${cond}.pirna-counts.filt.pos.tsv.tmp $sdir/merged-${cond}.pirna-counts.filt.pos.tsv

cond="Mili"
table-join \
    --full \
    --table1 mmu.RIPSeq.MiliIP.MiwiHet.P24.1/pirna-counts.filt.pos.tsv.tmp  \
    --key1 'name' \
    --table2 mmu.RIPSeq.MiliIP.MiwiHet.P24.2/pirna-counts.filt.pos.tsv.tmp \
    | table-join \
        --full \
        --table1 - \
        --key1 'name' \
        --table2 mmu.RIPSeq.MiliIP.MiwiRK.P24.1/pirna-counts.filt.pos.tsv.tmp \
    | table-join \
        --full \
        --table1 - \
        --key1 'name' \
        --table2 mmu.RIPSeq.MiliIP.MiwiRK.P24.2/pirna-counts.filt.pos.tsv.tmp \
    > $sdir/merged-${cond}.pirna-counts.filt.pos.tsv && rm mmu.RIPSeq.MiliIP.Miwi*.P24.[1,2]/pirna-counts.filt.pos.tsv.tmp

table-join --left \
    --table1 $sdir/merged-${cond}.pirna-counts.filt.pos.tsv \
    --key1 'name' \
    --table2 $sdir/merged.pirna-counts.filt.coords.txt \
    > $sdir/merged-${cond}.pirna-counts.filt.pos.tsv.tmp && mv $sdir/merged-${cond}.pirna-counts.filt.pos.tsv.tmp $sdir/merged-${cond}.pirna-counts.filt.pos.tsv

# get differential expression - starts
compar="P24"
mkdir -p ${compar}/de

cond="Miwi"
# By feature
./src/R/pirna-distro-de.R \
    --ifile "compare/starts/uniq/merged-${cond}.pirna-counts.filt.pos.tsv" \
    --design "data/design-pirna.txt" \
    --odir "${compar}/de/starts/uniq/${cond}" \
    --samples "mmu.RIPSeq.MiwiIP.MiwiHet.P24.2,mmu.RIPSeq.MiwiIP.MiwiHet.P24.3,mmu.RIPSeq.MiwiIP.MiwiRK.P24.2,mmu.RIPSeq.MiwiIP.MiwiRK.P24.3" \
    --comparcond1 "Het" \
    --comparcond2 "RK" \
    --refcond "Het" \
    --pval 0.001 \
    --fc 2 \
    --byfeature &
# By coordinate
./src/R/pirna-distro-de.R \
    --ifile "compare/starts/uniq/merged-${cond}.pirna-counts.filt.pos.tsv" \
    --design "data/design-pirna.txt" \
    --odir "${compar}/de/starts/uniq/${cond}" \
    --samples "mmu.RIPSeq.MiwiIP.MiwiHet.P24.2,mmu.RIPSeq.MiwiIP.MiwiHet.P24.3,mmu.RIPSeq.MiwiIP.MiwiRK.P24.2,mmu.RIPSeq.MiwiIP.MiwiRK.P24.3" \
    --comparcond1 "Het" \
    --comparcond2 "RK" \
    --refcond "Het" \
    --pval 0.001 \
    --fc 2 &

cond="Mili"
# By feature
./src/R/pirna-distro-de.R \
    --ifile "compare/starts/uniq/merged-${cond}.pirna-counts.filt.pos.tsv" \
    --design "data/design-pirna.txt" \
    --odir "${compar}/de/starts/uniq/${cond}" \
    --samples "mmu.RIPSeq.MiliIP.MiwiHet.P24.1,mmu.RIPSeq.MiliIP.MiwiHet.P24.2,mmu.RIPSeq.MiliIP.MiwiRK.P24.1,mmu.RIPSeq.MiliIP.MiwiRK.P24.2" \
    --comparcond1 "Het" \
    --comparcond2 "RK" \
    --refcond "Het" \
    --pval 0.001 \
    --fc 2 \
    --byfeature &
# By coordinate
./src/R/pirna-distro-de.R \
    --ifile "compare/starts/uniq/merged-${cond}.pirna-counts.filt.pos.tsv" \
    --design "data/design-pirna.txt" \
    --odir "${compar}/de/starts/uniq/${cond}" \
    --samples "mmu.RIPSeq.MiliIP.MiwiHet.P24.1,mmu.RIPSeq.MiliIP.MiwiHet.P24.2,mmu.RIPSeq.MiliIP.MiwiRK.P24.1,mmu.RIPSeq.MiliIP.MiwiRK.P24.2" \
    --comparcond1 "Het" \
    --comparcond2 "RK" \
    --refcond "Het" \
    --pval 0.001 \
    --fc 2
wait

# join tables - sequences
for i in */pirna-counts.filt.seq.tsv; do
    cat $i | cut -f 1,3 > $i.tmp
done

sdir=compare/sequence-genome/uniq
mkdir -p $sdir

echo -e "name" > $sdir/merged.pirna-counts.filt.coords.txt
for i in */pirna-counts.filt.seq.tsv.tmp; do
    table-join --full \
        --table1 $sdir/merged.pirna-counts.filt.coords.txt \
        --key1 'name' \
        --table2 $i > $sdir/merged.pirna-counts.filt.coords.txt.tmp && mv $sdir/merged.pirna-counts.filt.coords.txt.tmp $sdir/merged.pirna-counts.filt.coords.txt
done

# merge counts - sequences
for name in "${samples[@]}"; do
    echo "Working on $name"
    sdir=$name
    mkdir -p $sdir

    cat $sdir/pirna-counts.filt.seq.tsv | cut -f1-2 > $sdir/pirna-counts.filt.seq.tsv.tmp
done

sdir=compare/sequence-genome/uniq

cond="Miwi"
table-join --full \
    --table1 mmu.RIPSeq.MiwiIP.MiwiHet.P24.1/pirna-counts.filt.seq.tsv.tmp \
    --key1 'name' \
    --table2 mmu.RIPSeq.MiwiIP.MiwiHet.P24.2/pirna-counts.filt.seq.tsv.tmp \
    | table-join --full \
        --table1 - \
        --key1 'name' \
        --table2 mmu.RIPSeq.MiwiIP.MiwiHet.P24.3/pirna-counts.filt.seq.tsv.tmp \
    | table-join --full \
            --table1 - \
            --key1 'name' \
            --table2 mmu.RIPSeq.MiwiIP.MiwiRK.P24.1/pirna-counts.filt.seq.tsv.tmp \
    | table-join --full \
            --table1 - \
            --key1 'name' \
            --table2 mmu.RIPSeq.MiwiIP.MiwiRK.P24.2/pirna-counts.filt.seq.tsv.tmp \
    | table-join --full \
            --table1 - \
            --key1 'name' \
            --table2 mmu.RIPSeq.MiwiIP.MiwiRK.P24.3/pirna-counts.filt.seq.tsv.tmp \
    > $sdir/merged-${cond}.pirna-counts.filt.seq.tsv && rm mmu.RIPSeq.MiwiIP.Miwi*.P24.[1,2,3]/pirna-counts.filt.seq.tsv.tmp

table-join --left \
    --table1 $sdir/merged-${cond}.pirna-counts.filt.seq.tsv \
    --key1 'name' \
    --table2 $sdir/merged.pirna-counts.filt.coords.txt \
    > $sdir/merged-${cond}.pirna-counts.filt.seq.tsv.tmp && mv $sdir/merged-${cond}.pirna-counts.filt.seq.tsv.tmp $sdir/merged-${cond}.pirna-counts.filt.seq.tsv

cond="Mili"
table-join --full \
    --table1 mmu.RIPSeq.MiliIP.MiwiHet.P24.1/pirna-counts.filt.seq.tsv.tmp \
    --key1 'name' \
    --table2 mmu.RIPSeq.MiliIP.MiwiHet.P24.2/pirna-counts.filt.seq.tsv.tmp \
    | table-join --full \
            --table1 - \
            --key1 'name' \
            --table2 mmu.RIPSeq.MiliIP.MiwiRK.P24.1/pirna-counts.filt.seq.tsv.tmp \
    | table-join --full \
            --table1 - \
            --key1 'name' \
            --table2 mmu.RIPSeq.MiliIP.MiwiRK.P24.2/pirna-counts.filt.seq.tsv.tmp \
    > $sdir/merged-${cond}.pirna-counts.filt.seq.tsv && rm mmu.RIPSeq.MiliIP.Miwi*.P24.[1,2]/pirna-counts.filt.seq.tsv.tmp

table-join --left \
    --table1 $sdir/merged-${cond}.pirna-counts.filt.seq.tsv \
    --key1 'name' \
    --table2 $sdir/merged.pirna-counts.filt.coords.txt \
    > $sdir/merged-${cond}.pirna-counts.filt.seq.tsv.tmp && mv $sdir/merged-${cond}.pirna-counts.filt.seq.tsv.tmp $sdir/merged-${cond}.pirna-counts.filt.seq.tsv

# get piRNA counts
name=`cat compare/sequence-genome/uniq/mmu.RIPSeq.MiwiIP.MiwiHet.P24.3-mmu.RIPSeq.MiwiIP.MiwiHet.P24.2.pirna-counts.filt.seq.tsv | head -1 | cut -f2`
sdir=${name}
mkdir $sdir
cat compare/sequence-genome/uniq/mmu.RIPSeq.MiwiIP.MiwiHet.P24.3-mmu.RIPSeq.MiwiIP.MiwiHet.P24.2.pirna-counts.filt.seq.tsv | tail -n+2 | cut -f1,2 > $sdir/pirna.txt

name=`cat compare/sequence-genome/uniq/mmu.RIPSeq.MiwiIP.MiwiHet.P24.3-mmu.RIPSeq.MiwiIP.MiwiHet.P24.2.pirna-counts.filt.seq.tsv | head -1 | cut -f3`
sdir=${name}
mkdir $sdir
cat compare/sequence-genome/uniq/mmu.RIPSeq.MiwiIP.MiwiHet.P24.3-mmu.RIPSeq.MiwiIP.MiwiHet.P24.2.pirna-counts.filt.seq.tsv | tail -n+2 | cut -f1,3 > $sdir/pirna.txt

name=`cat compare/sequence-genome/uniq/mmu.RIPSeq.MiwiIP.MiwiRK.P24.3-mmu.RIPSeq.MiwiIP.MiwiRK.P24.2.pirna-counts.filt.seq.tsv | head -1 | cut -f2`
sdir=${name}
mkdir $sdir
cat compare/sequence-genome/uniq/mmu.RIPSeq.MiwiIP.MiwiRK.P24.3-mmu.RIPSeq.MiwiIP.MiwiRK.P24.2.pirna-counts.filt.seq.tsv | tail -n+2 | cut -f1,2 > $sdir/pirna.txt

name=`cat compare/sequence-genome/uniq/mmu.RIPSeq.MiwiIP.MiwiRK.P24.3-mmu.RIPSeq.MiwiIP.MiwiRK.P24.2.pirna-counts.filt.seq.tsv | head -1 | cut -f3`
sdir=${name}
mkdir $sdir
cat compare/sequence-genome/uniq/mmu.RIPSeq.MiwiIP.MiwiRK.P24.3-mmu.RIPSeq.MiwiIP.MiwiRK.P24.2.pirna-counts.filt.seq.tsv | tail -n+2 | cut -f1,3 > $sdir/pirna.txt

name=`cat compare/sequence-genome/uniq/mmu.RIPSeq.MiliIP.MiwiHet.P24.1-mmu.RIPSeq.MiliIP.MiwiHet.P24.2.pirna-counts.filt.seq.tsv | head -1 | cut -f2`
sdir=${name}
mkdir $sdir
cat compare/sequence-genome/uniq/mmu.RIPSeq.MiliIP.MiwiHet.P24.1-mmu.RIPSeq.MiliIP.MiwiHet.P24.2.pirna-counts.filt.seq.tsv | tail -n+2 | cut -f1,2 > $sdir/pirna.txt

name=`cat compare/sequence-genome/uniq/mmu.RIPSeq.MiliIP.MiwiHet.P24.1-mmu.RIPSeq.MiliIP.MiwiHet.P24.2.pirna-counts.filt.seq.tsv | head -1 | cut -f3`
sdir=${name}
mkdir $sdir
cat compare/sequence-genome/uniq/mmu.RIPSeq.MiliIP.MiwiHet.P24.1-mmu.RIPSeq.MiliIP.MiwiHet.P24.2.pirna-counts.filt.seq.tsv | tail -n+2 | cut -f1,3 > $sdir/pirna.txt

name=`cat compare/sequence-genome/uniq/mmu.RIPSeq.MiliIP.MiwiRK.P24.1-mmu.RIPSeq.MiliIP.MiwiRK.P24.2.pirna-counts.filt.seq.tsv | head -1 | cut -f2`
sdir=${name}
mkdir $sdir
cat compare/sequence-genome/uniq/mmu.RIPSeq.MiliIP.MiwiRK.P24.1-mmu.RIPSeq.MiliIP.MiwiRK.P24.2.pirna-counts.filt.seq.tsv | tail -n+2 | cut -f1,2 > $sdir/pirna.txt

name=`cat compare/sequence-genome/uniq/mmu.RIPSeq.MiliIP.MiwiRK.P24.1-mmu.RIPSeq.MiliIP.MiwiRK.P24.2.pirna-counts.filt.seq.tsv | head -1 | cut -f3`
sdir=${name}
mkdir $sdir
cat compare/sequence-genome/uniq/mmu.RIPSeq.MiliIP.MiwiRK.P24.1-mmu.RIPSeq.MiliIP.MiwiRK.P24.2.pirna-counts.filt.seq.tsv | tail -n+2 | cut -f1,3 > $sdir/pirna.txt

# get piRNA RPM
samples=(
    "mmu.RIPSeq.MiwiIP.MiwiHet.P24.2"
    "mmu.RIPSeq.MiwiIP.MiwiHet.P24.3"
    "mmu.RIPSeq.MiwiIP.MiwiRK.P24.2"
    "mmu.RIPSeq.MiwiIP.MiwiRK.P24.3"
    "mmu.RIPSeq.MiliIP.MiwiHet.P24.1"
    "mmu.RIPSeq.MiliIP.MiwiHet.P24.2"
    "mmu.RIPSeq.MiliIP.MiwiRK.P24.1"
    "mmu.RIPSeq.MiliIP.MiwiRK.P24.2"
)

for name in "${samples[@]}"; do
    echo "Working on $name"
    sdir=$name

    DIV=`awk -F'\t' '{sum+=$2;} END{print sum;}' $sdir/pirna.txt` # get lib. size
    awk -v divisor=${DIV} 'BEGIN {FS="\t";OFS="\t"};{$3=$2/(divisor/1000000); print $0;}' $sdir/pirna.txt > tmp; mv tmp $sdir/pirna.txt # get RPMs
    cat $sdir/pirna.txt | grep -v -P "\t0" > $sdir/pirna.1rpm.txt # simply removing starting with 0.xxx (<1 rpm)
done

# get only piRNAs in both replicates
# Note: Assuming we have 2 replicates!
cond="Miwi"
cat mmu.RIPSeq.MiwiIP.MiwiHet.P24.2/pirna.1rpm.txt mmu.RIPSeq.MiwiIP.MiwiHet.P24.3/pirna.1rpm.txt \
    | cut -f1 | sort -T . --parallel $threads | uniq -c | tr -s ' ' | grep "^ 2 " | cut -d ' ' -f3 \
    > comparison/mmu.RIPSeq.${cond}IP.MiwiHet.P24.pirna.1rpm.txt
cat mmu.RIPSeq.MiwiIP.MiwiRK.P24.2/pirna.1rpm.txt mmu.RIPSeq.MiwiIP.MiwiRK.P24.3/pirna.1rpm.txt \
    | cut -f1 | sort -T . --parallel $threads | uniq -c | tr -s ' ' | grep "^ 2 " | cut -d ' ' -f3 \
    > comparison/mmu.RIPSeq.${cond}IP.MiwiRK.P24.pirna.1rpm.txt
cat comparison/mmu.RIPSeq.${cond}IP.MiwiHet.P24.pirna.1rpm.txt comparison/mmu.RIPSeq.${cond}IP.MiwiRK.P24.pirna.1rpm.txt \
    | sort -T . --parallel=$threads | uniq > comparison/mmu.RIPSeq.${cond}IP.MiwiHet_MiwiRK.P24.pirna.1rpm.txt

cond="Mili"
cat mmu.RIPSeq.MiliIP.MiwiHet.P24.1/pirna.1rpm.txt mmu.RIPSeq.MiliIP.MiwiHet.P24.2/pirna.1rpm.txt \
    | cut -f1 | sort -T . --parallel $threads | uniq -c | tr -s ' ' | grep "^ 2 " | cut -d ' ' -f3 \
    > comparison/mmu.RIPSeq.${cond}IP.MiwiHet.P24.pirna.1rpm.txt
cat mmu.RIPSeq.MiliIP.MiwiRK.P24.1/pirna.1rpm.txt mmu.RIPSeq.MiliIP.MiwiRK.P24.2/pirna.1rpm.txt \
    | cut -f1 | sort -T . --parallel $threads | uniq -c | tr -s ' ' | grep "^ 2 " | cut -d ' ' -f3 \
    > comparison/mmu.RIPSeq.${cond}IP.MiwiRK.P24.pirna.1rpm.txt
cat comparison/mmu.RIPSeq.${cond}IP.MiwiHet.P24.pirna.1rpm.txt comparison/mmu.RIPSeq.${cond}IP.MiwiRK.P24.pirna.1rpm.txt \
    | sort -T . --parallel=$threads | uniq > comparison/mmu.RIPSeq.${cond}IP.MiwiHet_MiwiRK.P24.pirna.1rpm.txt

# run GTBuster
#   Zamore piRNA targeting from Wu PH, Fu Y, Cecchini K, Özata DM, Arif A, Yu T, Colpan C, Gainetdinov I, Weng Z, Zamore PD. 
#   The evolutionarily conserved piRNA-producing locus pi6 is required for male mouse fertility. 
#   Nature genetics. 2020 Jul;52(7):728-39.; https://github.com/weng-lab/GTBuster

pairs=(
    "mmu.RIPSeq.MiwiIP.MiwiHet.P24.2"   "mmu.PARESeq.polya.MiwiHet.P24"
    "mmu.RIPSeq.MiwiIP.MiwiHet.P24.3"   "mmu.PARESeq.polya.MiwiHet.P24"
    "mmu.RIPSeq.MiwiIP.MiwiRK.P24.2"    "mmu.PARESeq.polya.MiwiRK.P24"
    "mmu.RIPSeq.MiwiIP.MiwiRK.P24.3"    "mmu.PARESeq.polya.MiwiRK.P24"
    "mmu.RIPSeq.MiliIP.MiwiHet.P24.1"   "mmu.PARESeq.polya.MiwiHet.P24"
    "mmu.RIPSeq.MiliIP.MiwiHet.P24.2"   "mmu.PARESeq.polya.MiwiHet.P24"
    "mmu.RIPSeq.MiliIP.MiwiRK.P24.1"    "mmu.PARESeq.polya.MiwiRK.P24"
    "mmu.RIPSeq.MiliIP.MiwiRK.P24.2"    "mmu.PARESeq.polya.MiwiRK.P24"
)
tLen=${#pairs[@]}

for (( i=0; i<${tLen}; i++ )); do
    n1=${pairs[$i]}
    ((i++))
    n2=${pairs[$i]}

    name=`echo "$n1" | sed 's/mmu.RIPSeq.//'`
    name2=`echo "$n2" | sed 's/mmu.PARESeq.polya.//'`
    echo "$name and $name2"
    sdir=$name
    mkdir -p $sdir

    cut -f1,3 $n1/pirna.1rpm.txt > $n1/pirna.1rpm.txt.tmp  # get only piRNA\tRPM for GTBuster

    bash src/GTBuster/pirna_target_finder.wrapper.v3.sh $n1/pirna.1rpm.txt.tmp $sdir $n2/ref/p $threads \
        && rm $n1/pirna.1rpm.txt.tmp

    cut -f1-14 $sdir/gt2 > $sdir/gt2.tab
done

# filter common positions for both replicates
pairs=(
    "MiwiIP.MiwiHet.P24.2"  "MiwiIP.MiwiHet.P24.3"
    "MiwiIP.MiwiRK.P24.2"   "MiwiIP.MiwiRK.P24.3"
    "MiliIP.MiwiHet.P24.1"  "MiliIP.MiwiHet.P24.2"
    "MiliIP.MiwiRK.P24.1"   "MiliIP.MiwiRK.P24.2"
)
tLen=${#pairs[@]}

for (( i=0; i<${tLen}; i++ )); do
    n1=${pairs[$i]}
    ((i++))
    n2=${pairs[$i]}

    name=`echo "${n1%.*}"`
    echo $name
    sdir=comparison/$name
    mkdir -p $sdir

    ./src/R/filter-gtbuster.R $n1/gt2.tab $n2/gt2.tab $sdir/gt2.filt.tab \
        &> $sdir/gt2.filt-stats.txt &
done
wait

# annotate results with genes
samples=(
    "MiwiIP.MiwiHet.P24"
    "MiwiIP.MiwiRK.P24"
    "MiliIP.MiwiHet.P24"
    "MiliIP.MiwiRK.P24"
)

for name in "${samples[@]}"; do
    echo "Working on $name"
    sdir=comparison/$name
    mkdir -p $sdir

    # Keep only annotated
    ./src/R/annotate-gtbuster.R --ifile $sdir/gt2.filt.tab --annot $refdir/ensembl_genes.gtf \
        --ofile $sdir/gt2.filt.annot.tab --region "three_prime_utr,five_prime_utr,CDS" \
        &> $sdir/gt2.filt.annot.stats.txt &
done
wait

# prepare piRNA differential expression GTBuster input
conds=(
    "Miwi"
    "Mili"
)

sdir=comparison

for cond in "${conds[@]}"; do
  echo $cond

    cut -f8,10 $sdir/${cond}IP.MiwiHet.P24/gt2.filt.annot.tab $sdir/${cond}IP.MiwiRK.P24/gt2.filt.annot.tab \
        | sort -T $sdir | uniq > $sdir/RKvsHet/${cond}/pirna-gene.gtbuster.txt
done

# link piRNA sequence counts
for cond in "${conds[@]}"; do
  echo $cond

    ln -sf compare/sequence-genome/uniq/merged-${cond}.pirna-counts.filt.seq.tsv \
        $sdir/RKvsHet/${cond}/merged.pirna-counts.filt.seq.tsv
done

# get piRNA differential expression by seed
# IMPORTANT: Some part of the following script are hardcoded! Please check before you use it

cond="Miwi"
./src/R/pirna-seed-de.R \
    --ifile "$sdir/RKvsHet/${cond}/merged.pirna-counts.filt.seq.tsv" \
    --pirna_list "$sdir/RKvsHet/${cond}/pirna-gene.gtbuster.txt" \
    --design "data/design-pirna.txt" \
    --odir "comparison/RKvsHet/${cond}/de/sequence/uniq" \
    --samples "mmu.RIPSeq.MiwiIP.MiwiHet.P24.2,mmu.RIPSeq.MiwiIP.MiwiHet.P24.3,mmu.RIPSeq.MiwiIP.MiwiRK.P24.2,mmu.RIPSeq.MiwiIP.MiwiRK.P24.3" \
    --comparcond1 "Het" \
    --comparcond2 "RK" \
    --refcond "Het" \
    --pval 0.001 \
    --fc 2 \
    --type "seq" \
    && pigz -p 4 -f comparison/RKvsHet/${cond}/de/sequence/uniq/de.all-bySeq.inclCoords.tsv

cond="Mili"
./src/R/pirna-seed-de.R \
    --ifile "$sdir/RKvsHet/${cond}/merged.pirna-counts.filt.seq.tsv" \
    --pirna_list "$sdir/RKvsHet/${cond}/pirna-gene.gtbuster.txt" \
    --design "data/design-pirna.txt" \
    --odir "comparison/RKvsHet/${cond}/de/sequence/uniq" \
    --samples "mmu.RIPSeq.MiliIP.MiwiHet.P24.1,mmu.RIPSeq.MiliIP.MiwiHet.P24.2,mmu.RIPSeq.MiliIP.MiwiRK.P24.1,mmu.RIPSeq.MiliIP.MiwiRK.P24.2" \
    --comparcond1 "Het" \
    --comparcond2 "RK" \
    --refcond "Het" \
    --pval 0.001 \
    --fc 2 \
    --type "seq" \
    && pigz -p 4 -f comparison/RKvsHet/${cond}/de/sequence/uniq/de.all-bySeq.inclCoords.tsv
wait

conda activate pirna

# make table for degradome cuts
samples=(
    "mmu.PARESeq.polya.MiwiHet.P24.1"
    "mmu.PARESeq.polya.MiwiHet.P24.2"
    "mmu.PARESeq.polya.MiwiRK.P24.1"
    "mmu.PARESeq.polya.MiwiRK.P24.2"
)

# extract degradome table with names
for name in "${samples[@]}"; do
    sdir=$name

    echo -e "name\t${name}_rawCounts_degra\t${name}_rpm_degra" > $sdir/reads.1.5p-ext${extend_degra}.collapsed.tab
    zcat $sdir/reads.1.5p-ext50.collapsed.bed.gz \
        | awk 'BEGIN{FS="\t";OFS="\t"}{print $1","$2","$3","$6,$5,$7}' >> $sdir/reads.1.5p-ext${extend_degra}.collapsed.tab &
done
wait

table-join --full \
    --table1 mmu.PARESeq.polya.MiwiHet.P24.1/reads.1.5p-ext${extend_degra}.collapsed.tab \
    --key1 "name" \
    --table2 mmu.PARESeq.polya.MiwiHet.P24.2/reads.1.5p-ext${extend_degra}.collapsed.tab \
    | table-join --full \
    --table1 - \
    --key1 "name" \
    --table2 mmu.PARESeq.polya.MiwiRK.P24.1/reads.1.5p-ext${extend_degra}.collapsed.tab \
    | table-join --full \
    --table1 - \
    --key1 "name" \
    --table2 mmu.PARESeq.polya.MiwiRK.P24.2/reads.1.5p-ext${extend_degra}.collapsed.tab \
    | sed "s/\bNA\b/0/g" \
    | pigz -c -p $threads > comparison/reads.1.5p-ext${extend_degra}.collapsed.all.tab.gz

## run degradome occupancy - adjusted to RNA-Seq (similar to Ribo-Seq occupancy)
samples=(
    "mmu.RNASeq.total.MiwiHet.P24.1"
    "mmu.RNASeq.total.MiwiHet.P24.2"
    "mmu.RNASeq.total.MiwiHet.P24.3"
    "mmu.RNASeq.total.MiwiRK.P24.1"
    "mmu.RNASeq.total.MiwiRK.P24.2"
    "mmu.RNASeq.total.MiwiRK.P24.3"
    "mmu.PARESeq.polya.MiwiHet.P24.1"
    "mmu.PARESeq.polya.MiwiHet.P24.2"
    "mmu.PARESeq.polya.MiwiRK.P24.1"
    "mmu.PARESeq.polya.MiwiRK.P24.2"
)

mkdir degradome

# subset RNA-Seq only for matching PARE-Seq references
name=`echo ${samples[@]} | tr ' ' '\n' | grep PARESeq | head -1`
tail -n+2 $sampledir/$name/counts/quant.sf | cut -f1 > list.tmp

samples=(
    "mmu.RNASeq.total.MiwiHet.P24.1"
    "mmu.RNASeq.total.MiwiHet.P24.2"
    "mmu.RNASeq.total.MiwiHet.P24.3"
    "mmu.RNASeq.total.MiwiRK.P24.1"
    "mmu.RNASeq.total.MiwiRK.P24.2"
    "mmu.RNASeq.total.MiwiRK.P24.3"
)

for name in ${samples[@]}; do
    sdir=$name/counts
    mkdir -p $sdir

    subset-table -i $sampledir/$name/counts/quant.sf -c Name -l list.tmp -o $sdir/quant.sf &
done
wait

name=`echo ${samples[@]} | tr ' ' '\n' | grep RNASeq | head -1`
tail -n+2 $name/counts/quant.sf | cut -f1 > list.tmp

samples=(
    "mmu.PARESeq.polya.MiwiHet.P24.1"
    "mmu.PARESeq.polya.MiwiHet.P24.2"
    "mmu.PARESeq.polya.MiwiRK.P24.1"
    "mmu.PARESeq.polya.MiwiRK.P24.2"
)

for name in ${samples[@]}; do
  sdir=$name/counts
  mkdir -p $sdir

  subset-table -i $sampledir/$name/counts/quant.sf -c Name -l list.tmp -o $sdir/quant.sf &
done
wait

rm list.tmp

#  calculate differential degradome "occupancy"
samples=(
    "mmu.RNASeq.total.MiwiHet.P24.1"
    "mmu.RNASeq.total.MiwiHet.P24.2"
    "mmu.RNASeq.total.MiwiHet.P24.3"
    "mmu.RNASeq.total.MiwiRK.P24.1"
    "mmu.RNASeq.total.MiwiRK.P24.2"
    "mmu.RNASeq.total.MiwiRK.P24.3"
    "mmu.PARESeq.polya.MiwiHet.P24.1"
    "mmu.PARESeq.polya.MiwiHet.P24.2"
    "mmu.PARESeq.polya.MiwiRK.P24.1"
    "mmu.PARESeq.polya.MiwiRK.P24.2"
)

sdir=degradome
mkdir -p $sdir
genes_to_viz=""

./src/R/de-occupancy.R \
    --samples ${samples[@]} \
    --idir . \
    --gtf $refdir/Mus_musculus.GRCm38.99.gtf \
    --odir $sdir \
    --design data/design-degra.txt \
    --pval 0.1 \
    --fc 2 \
    --reftype RNA \
    --compcond1 Het \
    --compcond2 RK \
    --refcond Het \
    --replicates TRUE \
    --counts "salmon" \
    --genes_to_viz $genes_to_viz \
    --rnaseq "../expression/rna/de/RKvsHet/de-table.xls" \
    &> $sdir/de.log
wait

# make piRNA targeted subsets
samples=(
    "mmu.PARESeq.polya.MiwiHet.P24.1"
    "mmu.PARESeq.polya.MiwiHet.P24.2"
    "mmu.PARESeq.polya.MiwiRK.P24.1"
    "mmu.PARESeq.polya.MiwiRK.P24.2"
)

conds=(
    "Miwi"
    "Mili"
)

for cond in "${conds[@]}"; do
    echo $cond

    mkdir comparison/${cond}
    mkdir -p comparison/${cond}IP.MiwiRK.P24

    cut -f1 ../pirna-targeting/results/comparison/${cond}IP.MiwiRK.P24/gt2.filt.tab \
        | sed -e 's/|/\n/g' -e 's/_/\n/g' | sed 's/,[^,]*$//' \
        | sort -T comparison/${cond}IP.MiwiRK.P24 --parallel=6 | uniq \
        > comparison/${cond}IP.MiwiRK.P24/gt2.filt.names.txt &
    mkdir -p comparison/${cond}IP.MiwiHet.P24
    cut -f1 ../pirna-targeting/results/comparison/${cond}IP.MiwiHet.P24/gt2.filt.tab \
        | sed -e 's/|/\n/g' -e 's/_/\n/g'| sed 's/,[^,]*$//' \
        | sort -T comparison/${cond}IP.MiwiHet.P24 --parallel=6 | uniq \
        > comparison/${cond}IP.MiwiHet.P24/gt2.filt.names.txt &
    wait
    cat comparison/${cond}IP.MiwiRK.P24/gt2.filt.names.txt comparison/${cond}IP.MiwiHet.P24/gt2.filt.names.txt \
        | sort -T comparison --parallel=6 | uniq \
        > comparison/${cond}/gt2.filt.names.txt &
done
wait

# make degraome count tables
for name in ${samples[@]}; do
    sdir=$name
    mkdir -p $sdir

    echo -e "name\t${name}" > $sdir/degra-counts.tab
    zcat ../pirna-targeting/results/$name/reads.1.5p-ext50.collapsed.bed.gz | cut -f 4 | cut -d ',' -f1-4 > $sdir/tmp
    paste $sdir/tmp <(zcat ../pirna-targeting/results/$name/reads.1.5p-ext50.collapsed.bed.gz | cut -f 5) >> $sdir/degra-counts.tab && rm $sdir/tmp

    for cond in "${conds[@]}"; do
        subset-table -i $sdir/degra-counts.tab -l comparison/${cond}/gt2.filt.names.txt -c name > $sdir/degra-counts.filt.${cond}.tab &
    done
done
wait

for cond in "${conds[@]}"; do
    echo $cond

    table-join --full \
        --table1 mmu.PARESeq.polya.MiwiHet.P24.1/degra-counts.filt.${cond}.tab \
        --key1 'name' \
        --table2 mmu.PARESeq.polya.MiwiHet.P24.2/degra-counts.filt.${cond}.tab \
    | table-join --full \
        --table1 - \
        --key1 'name' \
        --table2 mmu.PARESeq.polya.MiwiRK.P24.1/degra-counts.filt.${cond}.tab \
        | table-join --full \
            --table1 - \
            --key1 'name' \
            --table2 mmu.PARESeq.polya.MiwiRK.P24.2/degra-counts.filt.${cond}.tab \
    > comparison/degra-counts.filt.${cond}.tab &
done
wait

# merge piRNA, RNA-Seq, and degradome results
# Note: Degradome cuts for all the samples; use ONLY those containing piRNA targeted cuts, not overall degradome cuts (they'll have mostly Xrn1 and other stuff)
#           with updates from targeting (degradome/comparison/degra-counts.filt.tab file)
for cond in "${conds[@]}"; do
    echo $cond

    ./src/R/merge-annot-gtbuster-de.R \
        --gene_de "../expression/rna/de/RKvsHet/de-table.xls" \
        --gtbuster_het "comparison/${cond}IP.MiwiHet.P24/gt2.filt.annot.tab" \
        --gtbuster_rk "comparison/${cond}IP.MiwiRK.P24/gt2.filt.annot.tab" \
        --pirna_de "comparison/RKvsHet/${cond}/de/sequence/uniq/de.all-bySeq.tsv" \
        --adjpval_pirna 0.25 \
        --degra_cuts "degradome/comparison/degra-counts.filt.${cond}.tab" \
        --ofile "comparison/RKvsHet/${cond}/targeting/pirna_adjpval_025/pirna_target.tsv" &
done

# predict piRNA binding on whole transcripts WITHOUT degradome information - theoretical cuts
# IMPORTANT: This will produce huge files! Up to 1TB in total
mkdir common

samples=(
    "MiwiIP.MiwiHet.P24"
    "MiwiIP.MiwiRK.P24"
    "MiliIP.MiwiHet.P24"
    "MiliIP.MiwiRK.P24"
)

for name in "${samples[@]}"; do
    echo "Working on $name"
    sdir=comparison/$name
    mkdir -p $sdir

    ./src/R/annotate-gtbuster.R --ifile $sdir/gt2.filt.tab --annot $refdir/ensembl_genes.gtf \
        --ofile $sdir/gt2.filt.annot-transcript.tab --region "three_prime_utr,five_prime_utr,CDS" \
        --feature "transcript" &> $sdir/gt2.filt.annot-transcript.stats.txt &
done
wait

for cond in "${conds[@]}"; do
    echo $cond

    cut -f10 comparison/${cond}IP.Miwi*.P24/gt2.filt.annot-transcript.tab | sort -T common --parallel=$threads \
        | uniq | sed '/transcript_id/d'> common/targeted_transcripts.${cond}.txt &
done
wait

# get fasta
gffread -w $refdir/transcripts-total.fa -g $refdir/genome.fa $refdir/transcripts-total.gtf

for cond in "${conds[@]}"; do
    echo $cond

    seqtk subseq $refdir/transcripts-total.fa common/targeted_transcripts.${cond}.txt \
        | gzip -c > common/targeted_transcripts.${cond}.fa.gz &
done
wait

# merge piRNA sequences
for cond in "${conds[@]}"; do
    echo $cond

    cut -f1 mmu.RIPSeq.${cond}IP.Miwi*.P24.*/pirna.1rpm.txt | sort -T . --parallel=$threads | uniq \
        | awk 'BEGIN {OFS="\t"}{print $1,"999"}' > common/pirna.1rpm.${cond}.txt &
done
wait

# chunk FASTA
lines=50 # 4000 fasta sequences (header+sequence)

for cond in "${conds[@]}"; do
    echo $cond

    sdir=common/whole_transcript_target_pred-targeted/${cond}
    mkdir -p $sdir/ref

    split -l $lines --additional-suffix=".fa" <(zcat common/targeted_transcripts.${cond}.fa.gz) $sdir/ref/p
done

# run GTBuster
for cond in "${conds[@]}"; do
    echo $cond

    sdir=common/whole_transcript_target_pred-targeted/${cond}

    export TMPDIR=$sdir # export temporary temporary Python dir; you can also try TEMP and TMP if TMPDIR doesn't work
    bash src/GTBuster/pirna_target_finder.wrapper.v3-wholetranscript.sh common/pirna.1rpm.${cond}.txt $sdir $sdir/ref/p $threads

    # extract only first 14 columns
    cut -f1-14 $sdir/gt2 | gzip -c > $sdir/gt2.tab.gz && rm $sdir/gt2

    # prefilter to min number of full matches (same as in filter-gtbuster.R) because the file is way to big
    zcat $sdir/gt2.tab.gz | awk 'BEGIN{FS="\t"; OFS="\t"}{if($11 >= 14) print $0}' \
        | gzip -c > $sdir/gt2.tab.tmp.gz
    ./src/R/filter-gtbuster-nopos.R $sdir/gt2.tab.tmp.gz $sdir/gt2.filt.tab &> $sdir/gt2.filt-stats.txt && rm $sdir/gt2.tab.tmp.gz
done
wait

# get only expressed transcripts - transcripts with >= 100 reads in all samples together
samples=(
    "mmu.RNASeq.total.MiwiHet.P24.1"
    "mmu.RNASeq.total.MiwiHet.P24.2"
    "mmu.RNASeq.total.MiwiHet.P24.3"
    "mmu.RNASeq.total.MiwiRK.P24.1"
    "mmu.RNASeq.total.MiwiRK.P24.2"
    "mmu.RNASeq.total.MiwiRK.P24.3"
)

mkdir common

table-join --full \
    --table1 <(cut -f1,5 $sampledir/mmu.RNASeq.total.MiwiHet.P24.1/counts/quant.sf) \
    --key1 "Name" \
    --value1-as 'MiwiHet.P24.1' \
    --table2 <(cut -f1,5 $sampledir/mmu.RNASeq.total.MiwiHet.P24.2/counts/quant.sf) \
    --value2-as 'MiwiHet.P24.2' \
    | table-join --full \
        --table1 - \
        --key1 "Name" \
        --table2 <(cut -f1,5 $sampledir/mmu.RNASeq.total.MiwiHet.P24.3/counts/quant.sf) \
        --value2-as 'MiwiHet.P24.3' \
    | table-join --full \
        --table1 - \
        --key1 "Name" \
        --table2 <(cut -f1,5 $sampledir/mmu.RNASeq.total.MiwiRK.P24.1/counts/quant.sf) \
        --value2-as 'MiwiRK.P24.1' \
    | table-join --full \
        --table1 - \
        --key1 "Name" \
        --table2 <(cut -f1,5 $sampledir/mmu.RNASeq.total.MiwiRK.P24.2/counts/quant.sf) \
        --value2-as 'MiwiRK.P24.2' \
    | table-join --full \
        --table1 - \
        --key1 "Name" \
        --table2 <(cut -f1,5 $sampledir/mmu.RNASeq.total.MiwiRK.P24.3/counts/quant.sf) \
        --value2-as 'MiwiRK.P24.3' \
    > common/quant-all.sf

# rowsums with transcript in first column and filtering only >=100 reads in all the samples
tail -n+2 common/quant-all.sf \
    | awk 'BEGIN{OFS="\t"}{sum=0; for (i=2; i<=NF; i++) { sum+= $i } if (sum >= 100) print $1,sum}' \
    | cut -f1 \
    > common/expressed_transcripts.txt

# get fasta
seqtk subseq $refdir/transcripts-total.fa common/expressed_transcripts.txt | gzip -c > common/expressed_transcripts.fa.gz

# chunk FASTA
lines=1000

sdir=common/whole_transcript_target_pred-expressed
mkdir -p $sdir/ref

split -l $lines --additional-suffix=".fa" <(zcat common/expressed_transcripts.fa.gz) $sdir/ref/p

# run GTBuster
refdir_tmp=$sdir # temporary save reference fasta directory

for cond in "${conds[@]}"; do
    echo $cond

    sdir=common/whole_transcript_target_pred-expressed/${cond}
    mkdir -p $sdir
    export TMPDIR=$sdir # export temporary temporary Python dir; you can also try TEMP and TMP if TMPDIR doesn't work

    bash src/GTBuster/pirna_target_finder.wrapper.v3-wholetranscript.sh common/pirna.1rpm.${cond}.txt $sdir $refdir_tmp/ref/p $threads

    cut -f1-14 $sdir/gt2 | pigz -f -p $threads -c > $sdir/gt2.tab.gz && rm $sdir/gt2

    # prefilter to min number of full matches (same as in filter-gtbuster.R) because the file is way to big
    unpigz -p $threads -c $sdir/gt2.tab.gz | awk 'BEGIN{FS="\t"; OFS="\t"}{if($11 >= 14) print $0}' \
        | gzip -c > $sdir/gt2.tab.tmp.gz
    ./src/R/filter-gtbuster-nopos.R $sdir/gt2.tab.tmp.gz $sdir/gt2.filt.tab \
        &> $sdir/gt2.filt-stats.txt && rm $sdir/gt2.tab.tmp.gz # we pre-filter so this should not do anything

    echo -e "V1\tV2\tV3\tV4\tV5\tV6\tV7\tV8\tV9\tV10\tV11\tV12\tV13" > $sdir/gt2.filt.tmp.tab # add dummy header
    cat $sdir/gt2.filt.tab >> $sdir/gt2.filt.tmp.tab
    cat $sdir/gt2.filt.tmp.tab \
        | paste-table-col --ifile - --col-name "sample" --col-val "transcript" \
        > $sdir/gt2.filt.samp.tab && rm $sdir/gt2.filt.tmp.tab

    tail -n+2 $sdir/gt2.filt.samp.tab | cut -f2 | sort -T $sdir --parallel=$threads | uniq \
        > $sdir/targeted_transcripts.txt
done

# convert and filter cuts from piRNA+degradome prediction
for cond in "${conds[@]}"; do
    echo $cond

    cat comparison/${cond}IP.MiwiHet.P24/gt2.filt.annot-transcript.tab \
        | ./src/R/genomic-to-transcriptomic-ensembl.R --ifile stdin --gtf $refdir/Mus_musculus.GRCm38.99.gtf \
        | paste-table-col --ifile - --col-name "sample" --col-val "MiwiHet"  \
        > comparison/${cond}IP.MiwiHet.P24/gt2.filt.annot-transcript.samp.bed &

    cat comparison/${cond}IP.MiwiRK.P24/gt2.filt.annot-transcript.tab \
        | ./src/R/genomic-to-transcriptomic-ensembl.R --ifile stdin --gtf $refdir/Mus_musculus.GRCm38.99.gtf \
        | paste-table-col --ifile - --col-name "sample" --col-val "MiwiRK" \
        > comparison/${cond}IP.MiwiRK.P24/gt2.filt.annot-transcript.samp.bed &
    wait

    mkdir -p comparison/${cond}
    table-cat \
        comparison/${cond}IP.MiwiHet.P24/gt2.filt.annot-transcript.samp.bed comparison/${cond}IP.MiwiRK.P24/gt2.filt.annot-transcript.samp.bed \
        > comparison/${cond}/gt2.filt.annot-transcript.samp.all.bed &
done
wait

# split gene annotation into genic elements
src/R/gff-to-genic-elements-bed.R \
	--input $refdir/ensembl_genes.gtf \
	> $refdir/genic_elements.bed

cat \
    $refdir/genic_elements.cds.bed \
    $refdir/genic_elements.utr5.bed \
    $refdir/genic_elements.utr3.bed \
    > common/genic_elements.bed

cat $refdir/ensembl_genes.gtf \
    | awk 'BEGIN{FS="\t"}{split($9, a, "; "); if($3~"transcript") print a[1]"\t"a[3]"\t"$1":"$4"-"$5"\t"a[5]"\t"$7"\t"a[7]"\t"a[10]}' \
    | sed 's/gene_id "//' | sed 's/transcript_id "//' | sed 's/gene_name "//' | sed 's/gene_biotype "//' | \
    sed 's/transcript_biotype "//' | sed 's/"//g' | sed '1igene_id\ttranscript_id\tloci\tgene_name\tstrand\tgene_biotype\ttranscript_biotype' \
    > $refdir/ensembl_transcripts-gene.table.txt

for cond in "${conds[@]}"; do
    echo $cond

    ./src/R/summarize-cuts-transcripts-gtbuster.R \
        --predict common/whole_transcript_target_pred-expressed/${cond}/gt2.filt.samp.tab \
        --real comparison/${cond}/gt2.filt.annot-transcript.samp.all.bed \
        --pirna comparison/mmu.RIPSeq.${cond}IP.MiwiHet_MiwiRK.P24.pirna.1rpm.txt \
        --bed common/genic_elements.bed \
        --genetotrans $refdir/ensembl_transcripts-gene.table.txt \
        --expr common/expressed_transcripts.txt \
        --target_theor common/whole_transcript_target_pred-expressed/${cond}/gt2.filt.tab \
        --trans_target comparison/RKvsHet/${cond}/targeting/pirna_adjpval_025/pirna_target.tsv \
        --trans_de "../expression/results/rna/de/RKvsHet/de-table.xls" \
        --de_target comparison/RKvsHet/${cond}/targeting/pirna_adjpval_025/pirna_target-selected.tsv \
        --pval_mrna 0.1 \
        --pval_pirna 0.25 \
        --odir comparison/${cond} &
done
wait