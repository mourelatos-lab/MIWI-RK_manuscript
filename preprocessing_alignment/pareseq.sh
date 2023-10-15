#!/bin/bash

sdir="PARESeq_sample" # Main sample directory
refdir="data/mm10" # Reference genome/annotation directory
threads=12

mkdir $sdir/logfiles
mkdir $sdir/align
mkdir $sdir/counts

### Adapter trimming
#adapters.fa
#>adapter/1
#AGATCGGAAGAGCACACGTC
#>adapter/2
#AGATCGGAAGAGCGTCGTGT
minlen=25
trimmomatic PE -threads $threads -phred33 \
    $sdir/fastq/reads.1.fastq.gz $sdir/fastq/reads.2.fastq.gz \
    $sdir/fastq/reads.1.adtrim.tmp.fastq.gz $sdir/fastq/reads.1.adtrim.unpaired.fastq.gz \
    $sdir/fastq/reads.2.adtrim.fastq.gz $sdir/fastq/reads.2.adtrim.unpaired.fastq.gz \
    $sdir/adapters.fa:2:30:5:3:true TRAILING:3 SLIDINGWINDOW:4:5 MINLEN:$minlen \
    &> $sdir/logfiles/trimmomatic.log

### Trim leading three nt from R1 (library prep. artifact) 
seqtk trimfq -b 3 $sdir/fastq/reads.1.adtrim.tmp.fastq.gz \
    | gzip -c > $sdir/fastq/reads.1.adtrim.fastq.gz

### Reference index - STAR
mkdir $refdir/STAR-2.7.2b-annot

STAR \
    --runMode genomeGenerate \
    --runThreadN $threads \
    --genomeDir $refdir/STAR-2.7.2b-annot/ \
    --genomeFastaFiles $refdir/genome.fa \
	--sjdbGTFfile $refdir/ensembl_genes.gtf \
	--sjdbOverhang 100 \
    &> $refdir/STAR-2.7.2b-annot/star-index.log

### Alignment
STAR --runThreadN $threads --genomeDir $refdir/STAR-2.7.2b-annot/ \
    --readFilesIn $sdir/fastq/reads.1.adtrim.fastq.gz $sdir/fastq/reads.2.adtrim.fastq.gz --readFilesCommand zcat \
    --limitBAMsortRAM 40000000000 \
    --outFileNamePrefix $sdir/align/reads.1. \
    --outMultimapperOrder Random --outSAMtype BAM SortedByCoordinate --outSAMattributes All --outSAMunmapped Within \
    --outSAMattrRGline ID:$(basename $sdir) PL:Illumina PU:$(basename $sdir) SM:$(basename $sdir) \
    --outFilterType BySJout --outFilterMultimapNmax 20 --outFilterMultimapScoreRange 1 --outFilterScoreMinOverLread 0.66 \
    --outFilterMatchNmin 15 --outFilterMatchNminOverLread 0.66 --outFilterMismatchNmax 999 --outFilterMismatchNoverLmax 0.05 \
    --outFilterMismatchNoverReadLmax 1 --outFilterIntronMotifs RemoveNoncanonicalUnannotated --seedSearchStartLmax 25 \
    --alignIntronMin 20 --alignIntronMax 1000000 --alignMatesGapMax 1000000 --alignSJoverhangMin 5 --alignSJDBoverhangMin 3 \
    --alignEndsType Extend5pOfRead1 --sjdbGTFfile $refdir/ensembl_genes.gtf --sjdbOverhang 100 --quantMode GeneCounts TranscriptomeSAM \
    --twopassMode Basic

### Sort bam
samtools sort -@ $threads -T $sdir/align/deleteme \
    $sdir/align/reads.1.Aligned.toTranscriptome.out.bam \
    > $sdir/align/reads.1.Aligned.toTranscriptome.sortedByCoord.out.bam
samtools index -@ $threads $sdir/align/reads.1.Aligned.toTranscriptome.sortedByCoord.out.bam

### Remove rRNA - genome
bedtools_strand="-s" # -s - same strand (sense); -S opposite strand (antisense)
bedtools intersect $bedtools_strand -v \
    -abam $sdir/align/reads.1.Aligned.sortedByCoord.out.bam \
    -b $refdir/rRNA.bed \
    > $sdir/align/reads.1.Aligned.sortedByCoord.clean.out.bam 
samtools index $sdir/align/reads.1.Aligned.sortedByCoord.clean.out.bam

### Remove rRNA - transcriptome
samtools view -@ $threads -bh -L $refdir/rRNA.trans.bed \
    -U $sdir/align/reads.1.Aligned.toTranscriptome.sortedByCoord.clean.out.bam \
    $sdir/align/reads.1.Aligned.toTranscriptome.sortedByCoord.out.bam
samtools index $sdir/align/reads.1.Aligned.toTranscriptome.sortedByCoord.clean.out.bam

### Mark duplicates - genome
samtools sort -n $sdir/align/reads.1.Aligned.sortedByCoord.clean.out.bam | samtools view -h - \
    | samblaster --ignoreUnmated --addMateTags \
    2>$sdir/logfiles/samblaster-genome.log \
    | samtools view -bh - | samtools sort - \
    > $sdir/align/reads.1.Aligned.sortedByCoord.final.bam
samtools index $sdir/align/reads.1.Aligned.sortedByCoord.final.bam

### Mark duplicates - transcriptome
samtools sort -n $sdir/align/reads.1.Aligned.toTranscriptome.sortedByCoord.clean.out.bam | samtools view -h - \
    | samblaster --ignoreUnmated --addMateTags \
    2>$sdir/logfiles/samblaster-transcriptome.log \
    | samtools view -bh - | samtools sort - \
    > $sdir/align/reads.1.Aligned.toTranscriptome.sortedByCoord.final.bam
samtools index $sdir/align/reads.1.Aligned.toTranscriptome.sortedByCoord.final.bam

### Extract poly(A) genes
cat $refdir/ensembl_genes.gtf \
    | grep -v "tag \"mRNA_start_NF\"" \
    | grep -v "tag \"mRNA_end_NF\"" \
    | grep -v "tag \"cds_start_NF\"" \
    | grep -v "tag \"cds_end_NF\"" \
    | grep -E 'transcript_biotype \"protein_coding\"|transcript_biotype \"lncRNA\"|transcript_biotype \"retained_intron\"|transcript_biotype \"processed_transcript\"|transcript_biotype \"non_stop_decay\"' \
    > $refdir/transcripts-polya.gtf

gffread -w $refdir/transcripts-polya.fa -g $refdir/genome.fa $refdir/transcripts-polya.gtf

### Reference index - salmon
library="polya"
salmon_ref="transcripts_index"

mkdir -p $refdir/salmon/$library

generateDecoyTranscriptome.sh \
    -j $threads \
    -b bedtools -m mashmap \
    -a $refdir/transcripts-polya.gtf \
    -g $refdir/genome.fa \
    -t $refdir/transcripts-polya.fa \
    -o $refdir/salmon/$library/$salmon_ref

salmon index \
    -p $threads \
    -t $refdir/salmon/$library/$salmon_ref/gentrome.fa \
    -i $refdir/salmon/$library/$salmon_ref \
    --decoys $refdir/salmon/$library/$salmon_ref/decoys.txt \
    -k 31

### Gene expression estimates
libtype="ISF"
salmon quant \
    -p $threads \
    -i $refdir/salmon/$library/$salmon_ref \
    --libType $libtype \
    -1 $sdir/fastq/reads.1.adtrim.fastq.gz -2 $sdir/fastq/reads.1.adtrim.fastq.gz \
    --validateMappings \
    -o $sdir/counts \
    1> $sdir/logfiles/salmon.log 2>&1