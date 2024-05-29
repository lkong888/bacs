# Absolute quantitative and base-resolution sequencing reveals comprehensive landscape of pseudouridine across the human transcriptome

Authors: Haiqi Xu<sup>1,2,\*</sup>, Linzhen Kong<sup>1,2,\*</sup>, Jingfei Cheng<sup>1,2</sup>, Khatoun Al Moussawi<sup>1</sup>, Xiufei Chen<sup>1,2</sup>, Aleema Iqbal<sup>1,2</sup>, Peter A. C. Wing<sup>3</sup>, James M. Harris<sup>4</sup>, Senko Tsukuda<sup>4</sup>, Azman Embarc-Buh<sup>5</sup>, Guifeng Wei<sup>6</sup>, Alfredo Castello<sup>5</sup>, Skirmantas Kriaucionis<sup>1</sup>, Jane A. McKeating<sup>4</sup>, Xin Lu<sup>1</sup>, and Chun-Xiao Song<sup>1,2,†</sup>

# Introduction
These scripts are for BACS data analysis. The following steps are included:
1. Data preprocessing
2. Alignment and filtering
3. Mutation counts and site calling

# Spike in
NNUNN1(fully modified NNΨNN)
GCTTCAAGTTGANNTNNCATCGCAAGTGCA

NNUNN2(unmodified NNUNN)
ATGTCTCGACGTNNTNNGTTACAGTACCGT

# Reference
Hg38 for human samples; 
For RNA viruses, the following reference genomes were used: Severe acute respiratory syndrome coronavirus 2 isolate Wuhan-Hu-1 (NC_045512.2), Recombinant Hepatitis C virus J6(5’UTR-NS2)/JFH1 (JF343782.1), Zika virus isolate ZIKV/H.sapiens/Brazil/Natal/2015 (NC_035889.1), Hepatitis Delta Virus sequence from the pSVL(D3) plasmid99 (Addgene plasmid #29335) (https://www.addgene.org/29335/), and Sindbis virus (NC_001547.1). For EBV samples, reads were aligned to Epstein-Barr virus (EBV) genome, strain B95-8 (V01555.2)


# 1. Data preprocessing, alignments, filtering, and mutation counts

For Ribo-depletion libraries: bacs_smallrna.smk, bacs_smallrna_callsite.smk 

For polyA libraries: bacs_alignment1.smk bacs_mRNA_part2.smk and bacs_mRNA_site.smk

# 2. Raw signals correlation between CAPSplus, CAPS, TAB-seq and ACE-seq(Figure 4c, Figure S13)
caps_filtering.r

raw_signal.sh

figure_mESC.r

# 3. Raw signal correlation between mESCs replicates (Figure S12)
mESC_replicate.r

# 4. Statistical testing of high confidence 5hmC sites and genomic annotations(Figure S16)
 annotation.sh
 
 CAPS_annotation.r
 
 # 5. Subsampling and coverage analysis(Figure S15)
  SI12.sh
  
  figure_mESC.r
  
 # 6. Downstream analysis for normal brain and glioblastoma (Figure 4d,f, Figure S17)
 brain_work.sh
 
 figure_brain.r
 
 # 8. Per-base quality plot(Figure S14)
 phred_score.sh
 
 # 9. IGV visualisation(Figure 4e)
 bedtobw.sh
 
 
 ####Note: the data were co-analysed by Jingfei Cheng and Linzhen Kong
  
  
  
  
  
  
