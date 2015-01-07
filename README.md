# hiCLIP reveals the atlas of mRNA secondary structures recognized by Staufen 1 #


**This package performs analysis of high-throughput sequencing data produced by hiCLIP, mRNA-Seq and ribosome profiling.** 



### 1. Data analysis summary

a. **High-throughput sequence data** The fastq files are available from [ArrayExpress](http://www.ebi.ac.uk/arrayexpress/) with the following access codes: iCLIP and hiCLIP: E-MTAB-2937, mRNA-Seq: E-MTAB-2940, ribosome profiling: E-MTAB-2941. 

b. **Demultiplexing libraries** Fastq files were demultiplexed. 



c. **Data pre-processing** Data such as gtf and fasta files were downloaded and parsed for the analysis. Part c and e are implemented as `Sweave` scripts.

 

d. **Mapping of high-throughput sequencing reads** hiCLIP, ribosome profiling and mRNA-Seq data were mapped.



e. **Analyzing high-throughput sequencing data.** The analysis produces plots for the major figures used in the manuscript and performs the associated statistical analyses.



f. **Data visulaization of hybrid reads with IGV** The hybrid reads were visualized with [IGV-2.1.20](http://www.broadinstitute.org/software/igv/download_1.5). The following files were visualized with IGV:



```

    Genome: data/processed/fasta/longest_mRNA_and_ncRNA.fa

    Transcript coordinate: data/processed/gtf/gr_tc_longest_gene_IGV.gtf

    IGV compatible sam files for hiCLIP data: results/manuscript/IGV/

```
g. **Data visualization of hybrid reads with Circos plot**
The hybrid reads were visualized using Circos plot, and the scripts and the data for the visualization are located in the "R/visulalization_circos_plot folder"


### 2. ENVIRONMENT

Tested on Linux.





### 3. USAGE

`python run_analysis.py -d /path/to/STAU1_hiCLIP -p /path/to/Python-2.7.1 -r /path/to/R-2.15.1`


Below is the example code to run the analysis.



```

qsub -b y ~/Python-2.7.1/python /user/STAU1_hiCLIP/exec/run_analysis.py -d /user/STAU1_hiCLIP -p ~/Python-2.7.1 -r ~/R-2.15.1
```



### 4. PREREQUISITE

1. Python-2.7.1, R-2.15.1 and pdflatex have to be installed.

2. R packages described in the methods section should be installed in `R-2.15.1/library/`

3. All other binary packages described in the methods section have to be compiled and located under `bin/` directory.  





### 5. ACKNOWLEDGEMENT

The package was developed by Yoichiro Sugimoto with the help of Kathi Zarnak and Nejc Haberman by referring to Python Cookbook, R Graphics Cookbook and online discussions, and the scripts in R/visulalization_circos_plot folder were written by Elodie Darbo. If any acknowledgement is missing, please let us know.



