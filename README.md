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
- piRNA clusters: [piRNAclusterDB; accessed 2020-07-30](https://www.smallrnagroup.uni-mainz.de/piRNAclusterDB/)


## Code
The following section summarizes main parts of the code needed to reproduce publication findings. 

### Tools and software
Majority of tools were installed using `conda 23.3.1`.

Some analyses used different software versions compared to the *main* environment. This is noted in the brackes.

Used software (alphabetically):
- `bedtools 2.29.0`
- `cutadapt 2.5`
- `FastQC 0.11.9`
- `fqtrim 0.9.7`
- `gffread 0.12.1`
- `jvarkit ebcbaba`
- `mashmap 2.0`
- `piPipes c93bde3`
- `R 3.5.1`
- `R 3.6.3` (differential expression/occupancy)
- `R 4.0` (gff-to-genic-elements-bed.R)
- `rush 0.4.2`
- `salmon 1.2.1`
- `SalmonTools 23eac84`
- `samblaster 0.1.26`
- `samtools 1.10`
- `seqtk 1.3`
- `SQuIRE 0.9.9.92`
- `STAR-2.7.2b`
- `trimmomatic 0.39`
- `umi_tools 1.0.1`

Used R packages (alphabetically):
- `biomaRt 2.42.0` (differential expression/occupancy)
- `dplyr 0.8.5`
- `dplyr 1.0.4` (differential expression/occupancy)
- `edgeR 3.28.0` (differential expression/occupancy)
- `GenomicFeatures 1.40.0` (gff-to-genic-elements-bed.R)
- `ggplot2 3.3.0`
- `ggplot2 3.3.3` (differential expression/occupancy)
- `gplots 3.1.1` (differential expression/occupancy)
- `grid 2.3` (differential expression/occupancy)
- `optparse 1.6.2`
- `optparse 1.6.6` (gff-to-genic-elements-bed.R)
- `RColorBrewer 1.1_2` (differential expression/occupancy)
- `reshape2 1.4.4` (differential expression/occupancy)
- `rio 0.5.16`
- `rtracklayer 1.46.0` (differential expression/occupancy)
- `tibble 3.0.6` (differential expression/occupancy)
- `tximport 1.14.0` (differential expression/occupancy)

Some general custom scripts are provided in [src](src).

### Read preprocessing and mapping
Assuming the input reads are saved in `sampleName/fastq/reads.1.fastq.gz` and `sampleName/fastq/reads.2.fastq.gz` (for paired-end); references are saved in `data/assemblyName`.

Sample directory structure:
- `sampleName/fastq`: FASTQ files
- `sampleName/alignment`: BAM files
- `sampleName/counts`: Counts
- `sampleName/logfiles`: Logfiles

Data directory structure:
- `data/assemblyName/genome.fa`: reference genome sequence
- `data/assemblyName/Mus_musculus.GRCm38.99.gtf`: Ensembl gene annotation; Ensembl release version necessary for differential expression calculation
- `data/assemblyName/ensembl_genes.gtf`: Ensembl gene annotation; Can be a link to full-name Ensembl annotation
- `data/assemblyName/rRNA.bed`: rRNA genomic locations
- `data/assemblyName/pirna-clusters.bed`: piRNA clusters
- `data/assemblyName/rmsk.bed`: RepeatMasker repeats
- `data/assemblyName/rmsk_categ.tab`: RepeatMasker repeats with categories -> Download UCSC RMSK table with columns swScore, genoName, genoStart, genoEnd, strand, repClass

We use `mm10` as `assemblyName`.

Code for preprocessing and mapping: 
- [PARE-Seq](preprocessing_mapping/pareseq.sh)
- [Ribo-Seq](preprocessing_mapping/riboseq.sh)
- [RIP-Seq](preprocessing_mapping/ripseq.sh)
- [RNA-Seq](preprocessing_mapping/rnaseq.sh)

#### Preprocessing notes
piRNA (**RIP-Seq**) datasets were subsampled to 85,000,000 raw reads to have approximately the same depth as the sample with the lowest number of reads.

**RIP-Seq** datasets include UMI. **Ribo-Seq** and **RNA-Seq** datasets include random 3 nt at the R1 5' end and synthetically added poly(A) tail.

### Analysis

#### Figure 2 MIWI-NTRs sustain pachytene piRNA amplification of MIWI bound piRNAs and transposon control
E, F. 5'-5' distance (ping-pong) analyses and Z-scores of MILI-bound (e) and MIWI-bound (F) piRNAs mapping to pachytene 
clusters in indicated genotypes. G. Differential expression analysis of transposons between Miwi+/RK and MiwiRK/RK calculated 
from RNA-Seq libraries of P24 testes; top, MA plot (average expression relative to fold-change); red, adjusted p-value < 0.05 
and log2 fold-change >= 1; blue, adjusted p-value < 0.05 and log2 fold-change <=-1; grey, adjusted p-value >= 0.05 and/or 
log2 foldchange > -1 < 1.
##### E, F. 5'-5' distance (ping-pong)
- [Figure 2 E,F](figure2/figure2-ef.sh)
##### G. Differential transposons expression
- [Figure 2 G](figure2/figure2-g.sh)

#### Figure 3 MIWIRK impacts the transcriptome but not the translatome
Differential expression analysis calculated from RNA-Seq (A) and differential ribosome occupancy calculated from Ribo-seq 
(B) between Miwi+/RK and MiwiRK/RK P24 testes, visualized by Volcano (left) and MA (right) plots; red, adjusted p-value 
< 0.05 and log2 fold-change >= 1; blue, adjusted p-value < 0.05 and log2 fold-change <=-1; grey, adjusted p-value >= 0.05 
and/or log2 foldchange > -1 < 1.
##### A. Differential RNA expression
- [Figure 3 A](figure3/figure3-a.sh)
##### B. Differential ribosome occupancy
- [Figure 3 B](figure3/figure3-b.sh)

#### Figure 4 MIWI-NTRs sustain piRNAs that cleave and destabilize select mRNAs essential for spermiogenesis
Metatranscript distribution of predicted, piRNA-mediated cleavages on mRNA targets, without ("theoretical", A), or with 
degradome-Seq support (B).
##### A. Metatranscript distribution of predicted cleavages on mRNA targets
- [Figure 4 A](figure4/figure4-ab.sh)
##### B. Metatranscript distribution of piRNA-mediated cleavages on mRNA targets
- [Figure 4 B](figure4/figure4-ab.sh)

#### Supplementary Figure 2 Characteristics of MILI- and MIWI- bound piRNAs in Miwi+/RK and MiwiRK/RK
A. Genomic distribution of MILI and MIWI piRNAs from Miwi+/RK and MiwiRK/RK P24 testes. B. Base composition of piRNAs mapping 
to pachytene clusters. C. Length distribution off all piRNAs (left) or piRNA derived exclusively from pachytene clusters 
(right). Lengths between 24 to 32 nucleotides are highlighted with dotted lines.
##### A. Genomic distribution of MILI and MIWI piRNAs
- [Supplementary Figure 2](suppfigure2/suppfigure2-a.sh)
##### B. Base composition of piRNAs
- [Supplementary Figure 2](suppfigure2/suppfigure2-b.sh)
##### C. piRNA length distribution
- [Supplementary Figure 2](suppfigure2/suppfigure2-c.sh)