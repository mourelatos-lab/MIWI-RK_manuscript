#!/bin/bash
#
# Genomic piRNA distribution
#

refdir="$(pwd)/data/mm10" # Reference genome/annotation directory; full path
sampledir=$(pwd) # Directory with preprocessed and mapped samples

threads=12
rnd=$RANDOM

mkdir -p pirna-distro/run
cd pirna-distro/

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

# populate piPipes references
# Run piPipes installer and install mm10 with following parameters
#piPipes-c93bde3/piPipes.sh  install -g mm10
#export rRNA_MM=2
#export hairpin_MM=1
#export genome_MM=2
#export transposon_MM=2
#export siRNA_bot=20
#export siRNA_top=25
#export piRNA_bot=24
#export piRNA_top=32

# run piPipes
assembly="mm10"

for name in ${samples[@]}; do
    sdir=$name/piPipes
    mkdir -p $sdir

    name1=`echo $name | sed 's/mmu.RIPSeq.//' | sed 's/\./-/g'`

    cp $sampledir/$name/fastq/reads.1.adtrim.fastq.gz $sdir/${name1}.fastq.gz

    piPipes small \
        -i $sdir/${name1}.fastq.gz \
        -g $assembly \
        -c $threads \
        -N allXmiRNA \
        -F $refdir/rRNA.fa \
        -o $sdir \
        1> $sdir/piPipes.stdout \
        2> $sdir/piPipes.stderr && rm $sdir/${name1}.fastq.gz
done