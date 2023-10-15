# Miwi-RK_manuscript
Code for main findings in **N-terminal arginine methylation enables diversification of PIWI protein function** publication by Nicholas Vrettos`*`, Jan Oppelt`*`, Ansgar Zoch`*`, Paraskevi Sgourdou, Haruka Yoshida, Brian
Song, Ryan Fink, Dónal O’Carroll `‡`, and Zissimos Mourelatos `‡`. 

`*` These authors contributed equally to this work: NV, JO, AZ. `‡` Correspondence authors.

Publication details will be updated once the publication has been published.


## Datasets
All datasets have been deposited in the Sequence Read Archive under the BioProject accession [PRJNA977257](https://www.ncbi.nlm.nih.gov/bioproject/?term=PRJNA977257).


## References
- Mouse genome reference: [mm10; Ensembl release 99](ftp://ftp.ensembl.org/pub/release-99/fasta/mus_musculus/dna/Mus_musculus.GRCm38.dna.primary_assembly.fa.gz)
- Mouse gene annotation: [mm10; Ensembl release 99](ftp://ftp.ensembl.org/pub/release-99/gtf/mus_musculus/Mus_musculus.GRCm38.99.gtf.gz)
- rRNA annotation: [SILVA; accessed 2020-02-21](arb-silva.de); [mm10; Ensembl release 99](ftp://ftp.ensembl.org/pub/release-99/gtf/mus_musculus/Mus_musculus.GRCm38.99.gtf.gz)
  - rRNA is removed from the base gene annotation for `salmon` count estimates
- Repeat annotation: [mm10; UCSC](http://hgdownload.cse.ucsc.edu/goldenPath/mm10/database/rmsk.txt.gz)

## Code
The following section summarizes main parts of the code needed to reproduce publication findings. 


### Tools and software
Majority of tools were installed using `conda 23.3.1`.

Used software (alphabetically):
- `bedtools 2.29.0`
- `cutadapt 2.5`
- `FastQC 0.11.9`
- `fqtrim 0.9.7`
- `gffread 0.12.1`
- `piPipes c93bde3`
- `R 3.5.1`
- `salmon 1.2.1`
- `samblaster 0.1.26`
- `samtools 1.10`
- `seqtk 1.3`
- `SQuIRE 0.9.9.92`
- `STAR-2.7.2b`
- `trimmomatic 0.39`
- `umi_tools 1.0.1`

### Read preprocessing and alignment
- [PARE-Seq](preprocessing_alignment/pareseq.sh)
- [Ribo-Seq](preprocessing_alignment/riboseq.sh)
- [RIP-Seq](preprocessing_alignment/ripseq.sh)
- [RNA-Seq](preprocessing_alignment/rnaseq.sh)

Assuming the input reads are saved in `sampleName/fastq/reads.1.fastq.gz` and `sampleName/fastq/reads.2.fastq.gz` (for paired-end); references are saved in `data/assemblyName`.

Sample directory structure:
- `sampleName/fastq`: FASTQ files
- `sampleName/alignment`: BAM files
- `sampleName/counts`: Counts
- `sampleName/logfiles`: Logfiles

Data directory structure:
- `data/assemblyName/genome.fa`: reference genome sequence
- `data/assemblyName/ensembl_genes.gtf`: Ensembl gene annotation
- `data/assemblyName/rRNA.bed`: rRNA genomic locations

#### Preprocessing
All sRNA (RIP-Seq) datasets were subsampled to 85,000,000 raw reads to have approximately the same depth as the sample with the lowest number of reads.

RIP-Seq datasets include UMI. 

#### Alignment



### Figure 2 MIWI-NTRs sustain pachytene piRNA amplification of MIWI bound piRNAs and transposon control
E, F. 5'-5' distance (ping-pong) analyses and Z-scores of MILI-bound (e) and MIWI-bound (F) piRNAs mapping to pachytene 
clusters in indicated genotypes. G. Differential expression analysis of transposons between Miwi+/RK and MiwiRK/RK calculated 
from RNA-Seq libraries of P24 testes; top, MA plot (average expression relative to fold-change); red, adjusted p-value < 0.05 
and log2 fold-change >= 1; blue, adjusted p-value < 0.05 and log2 fold-change <=-1; grey, adjusted p-value >= 0.05 and/or 
log2 foldchange > -1 < 1.
#### E, F. 5'-5' distance (ping-pong)
#### G. Differential transposons expression

### Figure 3 MIWIRK impacts the transcriptome but not the translatome
Differential expression analysis calculated from RNA-Seq (A) and differential ribosome occupancy calculated from Ribo-seq 
(B) between Miwi+/RK and MiwiRK/RK P24 testes, visualized by Volcano (left) and MA (right) plots; red, adjusted p-value 
< 0.05 and log2 fold-change >= 1; blue, adjusted p-value < 0.05 and log2 fold-change <=-1; grey, adjusted p-value >= 0.05 
and/or log2 foldchange > -1 < 1.
#### A. Differential RNA expression
#### B. Differential ribosome occupancy

### Figure 4 MIWI-NTRs sustain piRNAs that cleave and destabilize select mRNAs essential for spermiogenesis
Metatranscript distribution of predicted, piRNA-mediated cleavages on mRNA targets, without (“theoretical”, A), or with 
degradome-Seq support (B).
#### A. Metatranscript distribution of predicted cleavages on mRNA targets
#### B. Metatranscript distribution of piRNA-mediated cleavages on mRNA targets

### Supplementary Figure 2 Characteristics of MILI- and MIWI- bound piRNAs in Miwi+/RK and MiwiRK/RK
A. Genomic distribution of MILI and MIWI piRNAs from Miwi+/RK and MiwiRK/RK P24 testes. B. Base composition of piRNAs mapping 
to pachytene clusters. C. Length distribution off all piRNAs (left) or piRNA derived exclusively from pachytene clusters 
(right). Lengths between 24 to 32 nucleotides are highlighted with dotted lines.
#### A. Genomic distribution of MILI and MIWI piRNAs
#### B. Base composition of piRNAs
#### C. piRNA length distribution
