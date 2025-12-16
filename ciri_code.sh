#RNA-seq analysis using CIRI2 and CIRIquant 
#index Reference genome
bwa index -a bwtsw hg19.fa
# Align reads using BWA-MEM
# BWA-MEM is recommended by CIRI2 for accurate back-splice junction detection
bwa mem -T 19 hg19.fa sample_R1.fastq sample_R2.fastq 1> sample.sam 2> sample.log
# install CIRI2 and run this script for circular rna detection 
perl CIRI2.pl -I sample.sam -O outfile -F hg19.fa -A hg19.gtf
#install ciriquant for quantification and differential expression 
#Install CIRIquant from source code
# create and activate virtual env
pip install virtualenv
virtualenv -p /path/to/your/python2/executable venv
source ./venv/bin/activate

# Install CIRIquant and its requirement automatically
tar zxvf CIRIquant.tar.gz
cd CIRIquant
python setup.py install

# Manual installation of required pacakges is also supported
pip install -r requirements.txt
#A YAML-formated config file is needed for CIRIquant to find software and reference needed
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
 #Predict circRNAs using CIRI2
  CIRIquant -t 4 \
          -1 ./test_1.fq.gz \
          -2 ./test_2.fq.gz \
          --config ./chr1.yml \
          -o ./test \
          -p test
#Study with biological replicates
#Requires a customed analysis pipeline of edgeR and prep_CIRIquant is used to generate matrix of circRNA expression level / junction ratio and CIRI_DE_replicate for DE analysis
Step1: Prepare CIRIquant output files
CONTROL1 ./c1/c1.gtf C 1
CONTROL2 ./c2/c2.gtf C 2
CONTROL3 ./c3/c3.gtf C 3
CASE1 ./t1/t1.gtf T 1
CASE2 ./t2/t2.gtf T 2
CASE3 ./t3/t3.gtf T 3
Then, run prep_CIRIquant to summarize the circRNA expression profile in all samples
#These count matrices (CSV files) can then be imported into R for use by DESeq2 and edgeR (using the DESeqDataSetFromMatrix and DGEList functions, respectively).
Step2: Prepare StringTie output
#The output of StringTie should locate under output_dir/gene/prefix_out.gtf. You need to use prepDE.py from stringTie to generate the gene count matrix for normalization.
CONTROL1 ./c1/gene/c1_out.gtf
CONTROL2 ./c2/gene/c2_out.gtf
CONTROL3 ./c3/gene/c3_out.gtf
CASE1 ./t1/gene/t1_out.gtf
CASE2 ./t2/gene/t2_out.gtf
CASE3 ./t3/gene/t3_out.gtf
#Then, run prepDE.py -i sample_gene.lst and use gene_count_matrix.csv generated under current working directory for further analysis.
Step3: Differential expression analysis
For differential analysis using CIRI_DE_replicate, you need to install a R environment and edgeR package from Bioconductor.
 CIRI_DE_replicate \
          --lib  library_info.csv \
          --bsj  circRNA_bsj.csv \
          --gene gene_count_matrix.csv \
          --out  circRNA_de.tsv \
          --out2 gene_de.tsv
  
