#!/bin/bash

sdir="RNASeq_sample" # Main sample directory
refdir="data/mm10" # Reference genome/annotation directory
threads=12

mkdir $sdir/logfiles
mkdir $sdir/align
mkdir $sdir/counts

### Trim three non-random nt from R1 5' end (library prep. artifact) 
minlen=15
cutadapt \
        -g ^GGG \
        -g ^GGC \
        -g ^GCC \
        -g ^CCC \
        -g ^CCG \
        -g ^CGG \
        -g ^CGC \
        -g ^GCG \
        --times 1 \
        --overlap 3 \
        --minimum-length $minlen \
        --error-rate 0.0 \
        --discard-untrimmed \
        -o $sdir/fastq/reads.1.trinucltrim.fastq.gz \
        $sdir/fastq/reads.1.fastq.gz \
        &> $sdir/logfiles/cutadapt.trinucltrim.log

### Trim poly(A) and adapter
# Note: Adapter in this lib. prep. is ligated to the poly(A) -> this step also removes adapter
adapter="AAAAAAAAAA"
cutadapt \
    -a $adapter \
    --overlap 5 \
    --minimum-length $minlen \
    --error-rate 0.1 \
    --discard-untrimmed \
    -o $sdir/fastq/reads.1.adtrim.fastq.gz \
    $sdir/fastq/reads.1.trinucltrim.fastq.gz \
    &> $sdir/logfiles/cutadapt.adtrim.log

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
    --readFilesIn $sdir/fastq/reads.1.adtrim.fastq.gz --readFilesCommand zcat \
    --limitBAMsortRAM 40000000000 \
    --outFileNamePrefix $sdir/align/reads.1. \
    --outMultimapperOrder Random --outSAMtype BAM SortedByCoordinate --outSAMattributes All --outSAMunmapped Within \
    --outSAMattrRGline ID:$(basename $sdir) PL:Illumina PU:$(basename $sdir) SM:$(basename $sdir) --outFilterType BySJout \
    --outFilterMultimapNmax 20 --outFilterMultimapScoreRange 1 --outFilterScoreMinOverLread 0.66 --outFilterMatchNmin 15 \
    --outFilterMatchNminOverLread 0.66 --outFilterMismatchNmax 999 --outFilterMismatchNoverLmax 0.05 --outFilterMismatchNoverReadLmax 1 \
    --outFilterIntronMotifs RemoveNoncanonicalUnannotated --alignIntronMin 20 --alignIntronMax 1000000 --alignMatesGapMax 1000000 \
    --alignSJoverhangMin 5 --alignSJDBoverhangMin 3 --alignEndsType Local --sjdbGTFfile $refdir/ensembl_genes.gtf --sjdbOverhang 100 \
    --quantMode GeneCounts TranscriptomeSAM --twopassMode None

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

### Extract total genes
cat $refdir/ensembl_genes.gtf \
    | grep -E 'transcript_biotype \"IG_C_gene\"|transcript_biotype \"IG_D_gene\"|transcript_biotype \"IG_J_gene\"|transcript_biotype \"IG_V_gene\"|transcript_biotype \"lncRNA\"|transcript_biotype \"miRNA\"|transcript_biotype \"misc_RNA\"|transcript_biotype \"nonsense_mediated_decay\"|transcript_biotype \"non_stop_decay\"|transcript_biotype \"processed_transcript\"|transcript_biotype \"protein_coding\"|transcript_biotype \"retained_intron\"|transcript_biotype \"ribozyme\"|transcript_biotype \"scaRNA\"|transcript_biotype \"scRNA\"|transcript_biotype \"snoRNA\"|transcript_biotype \"snRNA\"|transcript_biotype \"sRNA\"|transcript_biotype \"TEC\"|transcript_biotype \"TR_C_gene\"|transcript_biotype \"TR_D_gene\"|transcript_biotype \"TR_J_gene\"|transcript_biotype \"TR_V_gene\"|transcript_biotype \"IG_LV_gene\"' \
    > $refdir/transcripts-total.gtf

gffread -w $refdir/transcripts-total.fa -g $refdir/genome.fa $refdir/transcripts-total.gtf

### Reference index - salmon
library="total"
salmon_ref="transcripts_index-k15"

mkdir -p $refdir/salmon/$library

generateDecoyTranscriptome.sh \
    -j $threads \
    -b bedtools -m mashmap \
    -a $refdir/transcripts-total.gtf \
    -g $refdir/genome.fa \
    -t $refdir/transcripts-total.fa \
    -o $refdir/salmon/$library/$salmon_ref

salmon index \
    -p $threads \
    -t $refdir/salmon/$library/$salmon_ref/gentrome.fa \
    -i $refdir/salmon/$library/$salmon_ref \
    --decoys $refdir/salmon/$library/$salmon_ref/decoys.txt \
    -k 15

### Gene expression estimates
libtype="SF"
salmon quant \
    -p $threads \
    -i $refdir/salmon/$library/$salmon_ref \
    --libType $libtype \
    -r $sdir/fastq/reads.1.adtrim.fastq.gz \
    --validateMappings \
    -o $sdir/counts \
    1> $sdir/logfiles/salmon.log 2>&1