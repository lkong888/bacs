##module load R/4.0.3-foss-2020b
library(tidyverse)
setwd("/gpfs3/well/ludwig/users/ebu571/27Oct2023")

####
######obtain high confident list
allsite1 <- read.table("site/treat1_input1_smallrna_allsites.txt", header=T, sep=" ", colClasses=c("ref"="character") ) 
allsite2 <- read.table("site/treat2_input2_smallrna_allsites.txt", header=T, sep=" ", colClasses=c("ref"="character") ) 

rrnamis1 <- read.table("align/H9_2023Oct26_S9.rrna.snoRNA.tRNA.other.mismatch.txt", header=F, sep="\t", colClasses=c("V3"="character") ) 
colnames(rrnamis1) <- c("chr", "end", "ref", "mismatch_count")
rrnamis2 <- read.table("align/H10_2023Oct26_S10.rrna.snoRNA.tRNA.other.mismatch.txt", header=F, sep="\t", colClasses=c("V3"="character") ) 
colnames(rrnamis2) <- c("chr", "end", "ref", "mismatch_count")

allsite1 <- merge(allsite1,rrnamis1, by=c("chr", "end", "ref"))
allsite2 <- merge(allsite2,rrnamis2, by=c("chr", "end", "ref"))

merge1 <- list(allsite1, allsite2)
merge1 <- merge1 %>% reduce(full_join, by=c("chr", "start", "end", "ref", "motif")) %>% na.omit()
#colnames(merge1) <- c("motif","chr", "start", "end", "ref","treat1_depth","treat1_T","treat1_C", "treat1_conversion","treat1_gap","treat1_gap_rate","ctrl1_depth","ctrl1_T","ctrl1_C","ctrl1_conversion","ctrl1_gap","ctrl1_gap_rate","conversion1", "treat1_fp", "treat1_level","treat1_p", "treat1_fdr","treat2_depth","treat2_T","treat2_C", "treat2_conversion","treat2_gap","treat2_gap_rate","ctrl2_depth","ctrl2_T","ctrl2_C","ctrl2_conversion","ctrl2_gap","ctrl2_gap_rate","conversion2", "treat2_fp", "treat2_level","treat2_p", "treat2_fdr")
colnames(merge1) <- c("chr", "end", "ref","motif","start","treat1_depth","treat1_T","treat1_C", "treat1_conversion","treat1_gap","treat1_gap_rate","ctrl1_depth","ctrl1_T","ctrl1_C","ctrl1_conversion","ctrl1_gap","ctrl1_gap_rate","conversion1", "treat1_fp", "treat1_level","treat1_p", "treat1_fdr","mis_count1","treat2_depth","treat2_T","treat2_C", "treat2_conversion","treat2_gap","treat2_gap_rate","ctrl2_depth","ctrl2_T","ctrl2_C","ctrl2_conversion","ctrl2_gap","ctrl2_gap_rate","conversion2", "treat2_fp", "treat2_level","treat2_p", "treat2_fdr","mis_count2")

nrow(merge1) #14967

high_confident_site <- merge1[(merge1$treat1_T+merge1$treat1_C)>=20 &(merge1$treat2_T+merge1$treat2_C)>=20 & merge1$treat1_level>=0.05 & merge1$treat2_level>=0.05 & merge1$treat1_fdr<0.001 & merge1$treat2_fdr<0.001, ]
high_confident_site <- high_confident_site[(high_confident_site$ctrl1_conversion<=0.01 |high_confident_site$ctrl1_C<=2) &(high_confident_site$ctrl2_conversion<=0.01 |high_confident_site$ctrl2_C<=2), ]
high_confident_site$mismatch_ratio1 <- high_confident_site$mis_count1/(high_confident_site$ctrl1_C+high_confident_site$ctrl1_T+high_confident_site$mis_count1)
high_confident_site$mismatch_ratio2 <- high_confident_site$mis_count2/(high_confident_site$ctrl2_C+high_confident_site$ctrl2_T+high_confident_site$mis_count2)
high_confident_site <- high_confident_site[high_confident_site$mismatch_ratio1<=0.10 & high_confident_site$mismatch_ratio2<=0.10,]

nrow(high_confident_site) ##1093

write.table(merge1, "site/snoRNA_tRNA_all_sites.txt", quote = FALSE, col.names = T, row.names = F)
write.table(high_confident_site, "site/snoRNA_tRNA_called_sites.txt", quote = FALSE, col.names = T, row.names = F)

######For tRNA

#################################get the high-expressed isodecoder list
######raw read counts > 1% of total reads in the same anticodon catogary
###cat featureCounts/H10_2023Oct26_S10.counts | grep -v ^# |cut -f2,7,10 | grep tRNA | grep -v mt >featureCounts/tRNA.counts.txt
counts <- read.table("featureCounts/tRNA.counts.txt",header=T, sep="\t")
colnames(counts) <- c("chr", "input2", "input1")
counts$tmp <- sub("tRNA_(\\w+_\\w+)_(\\w+_\\w+)", "\\1", counts$chr)
counts$anticodon <- sub("tRNA_(\\w+_\\w+)_\\w+", "\\1", counts$tmp)


counts_agg1 <- aggregate(input1~anticodon, counts, sum)
colnames(counts_agg1) <- c("anticodon", "input1_sum")
counts_agg2 <- aggregate(input2~anticodon, counts, sum)
colnames(counts_agg2) <- c("anticodon", "input2_sum")
counts_merge <- merge(counts, counts_agg1, by=c("anticodon"))
counts_merge <- merge(counts_merge, counts_agg2, by=c("anticodon"))

excludelist <- counts_merge[counts_merge$input1 <= (counts_merge$input1_sum*0.02) &counts_merge$input2 <= (counts_merge$input2_sum*0.02), ] ##67
filteredlist <- counts_merge[counts_merge$input1 > (counts_merge$input1_sum*0.02) | counts_merge$input2 > (counts_merge$input2_sum*0.02), ] ##193


write.table(excludelist,"featureCounts/tRNA_counts_lessthan1per.txt",sep='\t',quote=F,row.names=F)
write.table(filteredlist,"featureCounts/tRNA_counts_filtered.txt",sep='\t',quote=F,row.names=F)

#############################tRNA filtered####################################
###fdr filtered and counts filtered 
###10% level in at least two replicates
site <- read.table("site/snoRNA_tRNA_called_sites.txt",header=T, sep=" ",colClasses=c("ref"="character")) ###1093

tRNA_site <- merge(site, filteredlist[,c("chr", "input1","input2", "input1_sum", "input2_sum")], by=c("chr")) ##626
tRNA_site1 <- tRNA_site[tRNA_site$treat1_level >= 0.10 | tRNA_site$treat2_level >= 0.10 , ]
##tRNA_site1 <- tRNA_site[(tRNA_site$treat1_level >= 0.10 & tRNA_site$treat2_level >= 0.10 &tRNA_site$treat3_level < 0.10) | (tRNA_site$treat1_level >= 0.10 & tRNA_site$treat3_level >= 0.10 & tRNA_site$treat2_level < 0.10) | (tRNA_site$treat2_level >= 0.10 & tRNA_site$treat3_level >= 0.10 & tRNA_site$treat1_level < 0.10) |  (tRNA_site$treat2_level >= 0.10 & tRNA_site$treat3_level >= 0.10 & tRNA_site$treat1_level >= 0.10), ]
nrow(tRNA_site1) ##609
write.table(tRNA_site1,"site/snoRNA_tRNA_called_sites_fdrpassed_tRNA_filtered.txt",sep='\t',quote=F,row.names=F)


