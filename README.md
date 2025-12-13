# Differential Expression Profiling of Circular RNAs in OSCC Using  Bioinformatics Pipeline

## Project Overview

This repository corresponds to my MSc Bioinformatics thesis titled "Differential Expression Profiling of Circular RNAs in OSCC Using  Bioinformatics Pipeline". 
The study aims to identify and quantify circular RNAs (circRNAs) from RNA-Seq data derived from Oral Squamous Cell Carcinoma (OSCC) samples and matched controls. 
CircRNA detection was performed using *CIRI2*, followed by quantification and differential expression analysis using *CIRIquant*. 
The pipeline emphasizes reproducibility and is suitable for disease-focused circRNA transcriptomic studies.

## Background

Oral squamous cell carcinoma (OSCC) is one of the most common and aggressive
malignancies of the oral cavity, often diagnosed at advanced stages due to
the lack of reliable biomarkers for early detection and prognosis.
Recent studies have highlighted the regulatory roles of circular RNAs
(circRNAs) in cancer biology, including tumor progression, metastasis, and
therapy resistance. However, the functional relevance and expression
patterns of circRNAs in OSCC remain poorly understood, necessitating
systematic transcriptomic investigations.


## Tools Used
*BWA MEM-Read alignment 
*CIRI2-circRNA detection from RNA-Seq alignments
*CIRIquant-circRNA quantification and differential expression
*R progrraming - for volcano plots 
*Galaxy for -For functional enrichment

## Project Workflow

## Workflow Overview

1. Quality-controlled RNA-Seq reads were aligned to the reference genome
   using **BWA-MEM**, which is recommended for circRNA detection with CIRI2.
2. The resulting aligned SAM files were used as input for **CIRI2** to
   identify circRNA candidates.
3. **CIRIquant** was applied to quantify circRNA expression levels across
   OSCC and matched normal samples.
4. Differential expression analysis was performed between experimental
   conditions to identify significantly dysregulated circRNAs.
5. Significant circRNAs were selected for downstream functional analysis.
6.** Volcano plots** were generated for visualization of differentially
   expressed circRNAs.
7. **Gene Ontology (GO) and KEGG pathway** analyses were performed for functional
   enrichment analysis.


## Commands and arguments
----Align the reads to reference to generate SAM file using BWA-MEM tool 
bwa index -a bwtsw hg19.fa
bwa mem –T 19 hg19.fa read1.fq read2.fq > SRR_name.sam (paired-end reads)
bwa mem –T 19 hg19.fa read1.fq read2.fq 1> SRR_name.sam 2> SRR_name.log (to include for log files)

**run CIRI2**
perl CIRI2.pl -I SRR_name.sam -O outfile -F hg19.fa -A hg19.gtf
**Preparation for CIRIquant **
install anaconda 
# Download packed package
wget https://github.com/bioinfo-biols/CIRIquant/releases/download/v1.1.3/CIRIquant_v1.1.3.tar.gz
mkdir -p CIRIquant_env
tar zxvf CIRIquant_v1.1.3.tar.gz -C CIRIquant_env

# Configuration environments (required)
conda activate ./CIRIquant_env
cd CIRIquant_env
make
conda activate /path/to/CIRIquant_env
which CIRIquant
CIRIquant -t 4 \
          -1 ./test_1.fq.gz \
          -2 ./test_2.fq.gz \
          --config ./chr1.yml \
          -o ./test \
          -p test
prepare YAML-formated config file
reference:
  fasta: /home/zhangjy/Data/database/hg19.fa
  gtf: /home/zhangjy/Data/database/gencode.v19.annotation.gtf
  bwa_index: /home/zhangjy/Data/database/hg19/_BWAtmp/hg19
  hisat_index: /home/zhangjy/Data/database/hg19/_HISATtmp/hg19
// Example of config file
name: hg19
tools:
  bwa: /home/zhangjy/bin/bwa
  hisat2: /home/zhangjy/bin/hisat2
  stringtie: /home/zhangjy/bin/stringtie
  samtools: /home/zhangjy/bin/samtools

reference:
  fasta: /home/zhangjy/Data/database/hg19.fa
  gtf: /home/zhangjy/Data/database/gencode.v19.annotation.gtf
  bwa_index: /home/zhangjy/Data/database/hg19/_BWAtmp/hg19
  hisat_index: /home/zhangjy/Data/database/hg19/_HISATtmp/hg19
**Study with biological replicates**
Step1: Prepare CIRIquant output files
CONTROL1 ./c1/c1.gtf C 1
CONTROL2 ./c2/c2.gtf C 2
CONTROL3 ./c3/c3.gtf C 3
CASE1 ./t1/t1.gtf T 1
CASE2 ./t2/t2.gtf T 2
CASE3 ./t3/t3.gtf T 3
Then, run prep_CIRIquant to summarize the circRNA expression profile in all samples
 prep_CIRIquant -i sample.lst \
                 --lib library_info.csv \
                 --circ circRNA_info.csv \
                 --bsj circRNA_bsj.csv \
                 --ratio circRNA_ratio.csv
These count matrices (CSV files) can then be imported into R for use by DESeq2 and edgeR (using the DESeqDataSetFromMatrix and DGEList functions, respectively).
Step2: Prepare StringTie output
need to use prepDE.py from stringTie to generate the gene count matrix for normalization.
provide a text file sample_gene.lst containing sample IDs and path to StringTie outputs:
CONTROL1 ./c1/gene/c1_out.gtf
CONTROL2 ./c2/gene/c2_out.gtf
CONTROL3 ./c3/gene/c3_out.gtf
CASE1 ./t1/gene/t1_out.gtf
CASE2 ./t2/gene/t2_out.gtf
CASE3 ./t3/gene/t3_out.gtf
Then, run prepDE.py -i sample_gene.lst and use gene_count_matrix.csv generated under current working directory for further analysis.
Step3: Differential expression analysis
For differential analysis using CIRI_DE_replicate,needed to install a R environment and edgeR package from Bioconductor.
 CIRI_DE_replicate \
          --lib  library_info.csv \
          --bsj  circRNA_bsj.csv \
          --gene gene_count_matrix.csv \
          --out  circRNA_de.tsv \
          --out2 gene_de.tsv
##  Downstream Analysis and Visualization

After circRNA quantification using **CIRIquant**, downstream analyses were performed to interpret the biological significance of differentially expressed circRNAs.

###  Differential Expression Visualization

* Differentially expressed circRNAs were visualized using a **volcano plot** generated in **R**.
* The volcano plot represents:

  * **log2 fold change** on the x-axis
  * **−log10 adjusted p-value** on the y-axis
* Significantly upregulated and downregulated circRNAs were highlighted based on statistical cutoffs.

### Functional Enrichment Analysis

* Host genes of significantly differentially expressed circRNAs were subjected to:

  * **Gene Ontology (GO) enrichment analysis**
  * **KEGG pathway enrichment analysis**
* Enrichment analyses were performed using the **Galaxy web platform**, enabling reproducible and user-friendly functional annotation.
* Enriched biological processes and pathways relevant to OSCC were identified and interpreted.

## Applications

* circRNA expression profiling
* Differential circRNA analysis
* Candidate selection for ceRNA network construction
* Disease-specific circRNA biomarker discovery

## Notes

* Raw sequencing data are not shared due to data privacy constraints.
* This repository contains scripts and processed results for reproducibility.

## Author

**Jyoti Tiwari**
MSc Bioinformatics
Interests: circRNA analysis, RNA-Seq, cancer genomics, multi-omics integration

## Citation

If you use this workflow, please cite:

* Gao et al., *Bioinformatics* (CIRI2)
* Zhang et al., *Genome Biology* (CIRIquant)
