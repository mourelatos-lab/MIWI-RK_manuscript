#!/bin/bash
#
# Calculate differential ribosome occupancy
#
# IMPORTANT: Ribosome occupancy MUST be run after RNA expression. It uses RNA expression to normalize ribosome occupancy.
#

refdir="$(pwd)/data/mm10" # Reference genome/annotation directory; full path
sampledir=$(pwd) # Directory with preprocessed and mapped samples

threads=12
rnd=$RANDOM

mkdir -p ribosome/run
cd ribosome/

################################################################################
samples=(
    "mmu.RNASeq.total.MiwiHet.P24.1"
    "mmu.RNASeq.total.MiwiHet.P24.2"
    "mmu.RNASeq.total.MiwiHet.P24.3"
    "mmu.RNASeq.total.MiwiRK.P24.1"
    "mmu.RNASeq.total.MiwiRK.P24.2"
    "mmu.RNASeq.total.MiwiRK.P24.3"
    "mmu.RIBOSeq.total.MiwiHet.P24.1"
    "mmu.RIBOSeq.total.MiwiHet.P24.2"
    "mmu.RIBOSeq.total.MiwiHet.P24.3"
    "mmu.RIBOSeq.total.MiwiRK.P24.1"
    "mmu.RIBOSeq.total.MiwiRK.P24.2"
    "mmu.RIBOSeq.total.MiwiRK.P24.3"
)

sdir=ribo
mkdir -p $sdir

## Run my ribo-occupancy and rna-seq DE with edgeR
./src/R/de-occupancy.R \
    --samples ${samples[@]} \
    --idir $sampledir \
    --gtf $refdir/Mus_musculus.GRCm38.99.gtf \
    --odir $$sdir \
    --design $(pwd)/data/design-ribo.txt \
    --pval 0.05 \
    --fc 2 \
    --reftype RNA \
    --compcond1 Het \
    --compcond2 RK \
    --refcond Het \
    --replicates TRUE \
    --rnaseq "../expression/rna/de/RKvsHet/de-table.xls" \
    &> $sdir/de.log