#!/bin/bash

sdir="RIPSeq_sample" # Main sample directory
refdir="data/mm10" # Reference genome/annotation directory
threads=12

mkdir $sdir/logfiles
mkdir $sdir/align
mkdir $sdir/counts

### Subsampling
# Note: Subsampling was applied only for samples with >85M raw reads
subsample=85000000
seqtk sample \
    -s100 $sdir/fastq/reads.1.fastq.gz $subsample \
    | gzip -c > $sdir/fastq/reads.1.subsample.fastq.gz

### Adapter trimming
adapter="CTGTAGGCACCATCAATAGA"
minlen=15
cutadapt -a $adapter \
    --overlap 10 --minimum-length $minlen --error-rate 0.10 --times 2 --discard-untrimmed \
    -o $sdir/fastq/reads.1.ad3trim.fastq.gz $sdir/fastq/reads.1.subsample.fastq.gz \
    &> $sdir/logfiles/cutadapt.adtrim.log

### Read collapsing (including UMI)
fqtrim -l 15 \
    -C $sdir/fastq/reads.1.3adtrim.fastq.gz \
    2> $sdir/logfiles/collapsing.log \
    | gzip -c > $sdir/fastq/reads.1.3adtrim.collapsed.fastq.gz

### Extract and remove UMI
umi_form="NNNNNNNN"
umi_tools extract \
    --stdin=$sdir/fastq/reads.1.3adtrim.collapsed.fastq.gz \
    --stdout=$sdir/fastq/reads.1.3adtrim.umi.fastq.gz \
    --3prime \
    --bc-pattern=$umi_form \
    --log=$sdir/logfiles/umiextract.log

### Reference index - STAR
mkdir $refdir/STAR-2.7.2b-noannot

STAR \
    --runMode genomeGenerate \
    --runThreadN $threads \
    --genomeDir STAR-2.7.2b-noannot/ \
    --genomeFastaFiles $refdir/genome.fa \
    &> $refdir/STAR-2.7.2b-noannot/star-index.log

### Alignment
STAR --runThreadN $threads --genomeDir $refdir/STAR-2.7.2b-noannot/ \
    --readFilesIn $sdir/fastq/reads.1.adtrim.fastq.gz --readFilesCommand zcat \
    --limitBAMsortRAM 40000000000 \
    --outFileNamePrefix $sdir/align/reads.1. \
    --outMultimapperOrder Random --outSAMtype BAM SortedByCoordinate --outSAMattributes All --outSAMunmapped Within \
    --outSAMattrRGline ID:$(basename $sdir) PL:Illumina PU:$(basename $sdir) SM:$(basename $sdir) \
    --outFilterType BySJout --outFilterMultimapNmax 1000000 --outFilterMultimapScoreRange 1 --outFilterScoreMinOverLread 0.66 \
    --outFilterMatchNmin 15 --outFilterMatchNminOverLread 0.66 --outFilterMismatchNmax 2 --outFilterMismatchNoverLmax 1 \
    --outFilterMismatchNoverReadLmax 1 --outFilterIntronMotifs RemoveNoncanonicalUnannotated --seedSearchStartLmax 10 \
    --alignIntronMin 20 --alignIntronMax 1 --alignSJoverhangMin 5 --alignSJDBoverhangMin 999 --alignEndsType EndToEnd \
    --sjdbGTFfile $refdir/ensembl_genes.gtf --sjdbOverhang 100 --quantMode GeneCounts TranscriptomeSAM --twopassMode None

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
umi_mm=0
umi_tools dedup \
    -I $sdir/align/reads.1.Aligned.sortedByCoord.clean.out.bam \
    -S $sdir/align/reads.1.Aligned.final.bam \
    -L $sdir/logfiles/umidedup.genome.log \
    --output-stats=$sdir/logfiles/umidedup.genome --temp-dir=$sdir \
    --buffer-whole-contig \
    --extract-umi-method=read_id --umi-separator="_" --method=directional --edit-distance-threshold=$umi_mm \
    --spliced-is-unique --multimapping-detection-method=NH
samtools index $sdir/align/reads.1.Aligned.toTranscriptome.final.bam

### Mark duplicates - transcriptome
umi_tools dedup \
    -I $sdir/align/reads.1.Aligned.toTranscriptome.sortedByCoord.clean.out.bam \
    -S $sdir/align/reads.1.Aligned.toTranscriptome.final.bam \
    -L $sdir/logfiles/umidedup.trans.log \
    --output-stats=$sdir/logfiles/umidedup.trans --temp-dir=$sdir \
    --buffer-whole-contig \
    --extract-umi-method=read_id --umi-separator="_" --method=directional --edit-distance-threshold=$umi_mm \
    --multimapping-detection-method=NH
samtools index $sdir/align/reads.1.Aligned.final.bam

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