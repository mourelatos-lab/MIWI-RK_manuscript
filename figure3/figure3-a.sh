#!/bin/bash
#
# Calculate differential expression
#

refdir="$(pwd)/data/mm10" # Reference genome/annotation directory; full path
sampledir=$(pwd) # Directory with preprocessed and mapped samples

threads=12
rnd=$RANDOM

mkdir -p expression/run
cd expression

################################################################################

samples=(
    "mmu.RNASeq.total.MiwiHet.P24.1"
    "mmu.RNASeq.total.MiwiHet.P24.2"
    "mmu.RNASeq.total.MiwiHet.P24.3"
    "mmu.RNASeq.total.MiwiRK.P24.1"
    "mmu.RNASeq.total.MiwiRK.P24.2"
    "mmu.RNASeq.total.MiwiRK.P24.3"
)

# calculate RNA differential expression
sdir=rna
mkdir -p $sdir

./src/R/de.R \
    --samples "${samples[@]}" \
    --idir $sampledir \
    --gtf $refdir/Mus_musculus.GRCm38.99.gtf \
    --odir $sdir \
    --pval 0.05 \
    --fc 2 \
    --compcond1 Het \
    --compcond2 RK \
    --refcond Het \
    --design "data/design-rna.txt" \
    --replicates TRUE \
    &> $sdir/de.log