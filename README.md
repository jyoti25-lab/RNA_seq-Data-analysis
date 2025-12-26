Dessertation project
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
1. BWA MEM-Read alignment 
2. CIRI2-circRNA detection from RNA-Seq alignments
3. CIRIquant-circRNA quantification and differential expression
4. R progrraming - for volcano plots
5. Galaxy for -For functional enrichment

## Project Workflow
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
present in code.sh
------------------------------------
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

## Author
**Jyoti Tiwari**
MSc Bioinformatics
Interests: circRNA analysis, RNA-Seq, cancer genomics, multi-omics integration

## Citation
If you use this workflow, please cite:
* Gao et al., *Bioinformatics* (CIRI2)
* Zhang et al., *Genome Biology* (CIRIquant)
