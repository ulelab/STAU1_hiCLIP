## 1. Description
In order to estimate the importance of Alu repeat elements, the enrichment of STAU1 non-hybrid reads (STAU1 target in the cytoplasm), total STAU1 iCLIP reads (STAU1 target in the total fraction of the cell), hnRNP C iCLIP (known to interact with Alu elements), compared to mRNA-Seq reads were examined using [Repeat Enrichment Estimator](http://compbio.med.harvard.edu/repeats/login.jsp).

## 2. The raw data used for the analysis
The fasta file prepared in the STAU1_hiCLIP pipeline was examined.
In detail, the following files were submitted to the Repeat Enrichment Estimator web server (for the analysis using hg18), and the output files () were examined.

The files submitted as "ChIP Files".

a. hiCLIP_nonHybrid_merged.26.fa.gz
b. maerged_STAU1_iCLIP.fa.gz
c. id1_hnRNPc_RUTT.fa.gz

The files submitted as "Control Files".

d. mRNASeq_merged.26.fa.gz

## 3. Output
The output files (standard repeat enrichment table) from the Repeat Enrichment Estimator web server were saved with the follwing name in this directory.

Comparison between a and d: STAU1_nonhybrid_vs_mRNASeq_repenr.std.csv
Comparison between b and d: STAU1_totaliCLIP_vs_mRNASeq_repenr.std.csv
Comparison between c and d: hnRNPc_iCLIP_vs_mRNASeq_repenr.std.csv


