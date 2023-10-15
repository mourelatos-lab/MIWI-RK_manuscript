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

# set temp dirs for Python (one of them should work)
mkdir $(pwd)/tmp
export TMPDIR=$(pwd)/tmp
export TMP=$(pwd)/tmp
export TEMP=$(pwd)/tmp

# variables for SQuIRE

read_len=47
strand="forward" # [forward|reverse|none]
threads=12
assembly="mm10"

# SQuIRE references
squire Fetch -b $assembly -o squire -f -c -r -g -x -p $threads
squire Clean -r squire/${assembly}_rmsk.txt -b $assembly -o squire/clean

# Set strandedness: '0' if unstranded eg Standard Illumina, 1 if first- strand eg Illumina Truseq, dUTP, NSR, NNSR, 2 if second-strand, eg Ligation, Standard SOLiD (optional,default=0)
if [ "$strand" == "forward" ]; then
    echo "Using $strand strand specificity."
    strand=2 # doesn't make sense but that's how squire wants it
elif [ "$strand" == "reverse" ]; then
    echo "Using $strand strand specificity."
    strand=1 # doesn't make sense but that's how squire wants it
else
    echo "Using $strand strand specificity."
    strand=0
fi

# SQuIRE map
for name in ${samples[@]}; do
    sdir=$name/squire
    mkdir -p $sdir

    squire Map \
      -1 $sampledir/$name/fastq/reads.1.adtrim.fastq.gz \
      -o $sdir -f squire \
      -r $read_len -n "reads.1.adtrim" -b $assembly --gtf squire/${assembly}_refGene.gtf -p $threads \
      &> $sdir/squire-map.log
    rm -r $sdir/reads.1.adtrim_STARgenome $sdir/reads.1.adtrim_STARpass1
done

# SQuIRE count
for name in ${samples[@]}; do
    sdir=$name/squire
    mkdir -p $sdir/tmp

    squire Count -m $sdir -c squire/clean -o $sdir -t $sdir/tmp -f squire \
        -r $read_len -n "reads.1.adtrim" -b $assembly -p $threads -s $strand -e 2 &> $sdir/squire-count.log
    rm -r $sdir/tmp
done

# SQuIRE differential expression
sdir=rna/squire_repeats/de
mkdir -p $sdir/subfamily

for name in ${samples[@]}; do
    cat ${name}/squire/reads.1.adtrim_subFcounts.txt | sed "s/reads\.1\.adtrim/${name}/g" > $sdir/${name}_subFcounts.txt &
    cat ${name}/squire/reads.1.adtrim_TEcounts.txt | sed "s/reads\.1\.adtrim/${name}/g" > $sdir/${name}_TEcounts.txt &
    cat ${name}/squire/reads.1.adtrim_refGenecounts.txt | sed "s/reads\.1\.adtrim/${name}/g" > $sdir/${name}_refGenecounts.txt &
    wait
done

group1="*MiwiRK.*"
group2="*MiwiHet.*"
group1_name="RK"
group2_name="Het"

squire Call -1 $group1 -2 $group2 -A $group1_name -B $group2_name -i $sdir -p $(echo $threads / 2 | bc) -N "squire_repeats" \
    -f "pdf" -o $sdir & # Call DE based on loci level
squire Call -1 $group1 -2 $group2 -A $group1_name -B $group2_name -i $sdir -p $(echo $threads / 2 | bc) -N "squire_repeats" \
    -f "pdf" -o $sdir/subfamily --subfamily & # Call DE based on subfamily counts
wait
done
