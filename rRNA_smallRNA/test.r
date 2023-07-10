##module load R/4.0.3-foss-2020b
library(dplyr)
library(reshape2)
library(tidyr)
library(ggplot2)
library(ggpubr)
library(cowplot)
library(RColorBrewer)
library(tidyverse)
setwd("/users/ludwig/ebu571/ebu571/20May2023_final")

read_sta <- read.table("determined_ambiguous_reads.txt", header=F, sep="\t")
colnames(read_sta)<- c("chr", "determined_reads", "ambiguous_reads")
read_sta$known_ratio <- read_sta$determined_reads/(read_sta$determined_reads+read_sta$ambiguous_reads)
sum(read_sta$determined_reads+read_sta$ambiguous_reads) ##23666714
sum(read_sta$determined_reads) ##14990816
sum(read_sta$ambiguous_reads) ##8675898

group1 <- read_sta[as.numeric(read_sta$known_ratio)<0.5 & read_sta$determined_reads >20,] ###186

snoRNA <- group1[!grepl("tRNA", group1$chr),] ##33
snoRNA$smp <- gsub(".*_SNO", "SNO", snoRNA$chr)
tRNA <- read_sta[grepl("tRNA", read_sta$chr) & !grepl("mt_tRNA", read_sta$chr),] ##260
mtRNA <- read_sta[grep("mt_tRNA", read_sta$chr),] ##22

melt_sno <- melt(snoRNA[,c("smp", "determined_reads", "ambiguous_reads")], id.vars="smp")

pdf("processed_plot/input1_snoRNA_read.pdf", width=60, height=4)
ggplot(melt_sno, aes(x=smp, y=value, fill=variable))+
geom_bar(stat = "identity", width=0.55,position = position_dodge(), alpha = 0.75)+  ####identity equals to sum; needs to make sure is the mean value; otherwise, use geom_bar(stat = "summary", fun = "mean", width=0.55)
theme_light() + ylab("reads")+xlab("")+
 theme(axis.text = element_text(color = "black")) +theme(legend.position="none") 
dev.off()

tRNA$smp <- gsub("_..$", "", tRNA$chr)
tRNA$smp <- gsub("_.$", "", tRNA$smp)

known_agg <- aggregate(determined_reads ~ smp, tRNA, sum)
amb_agg <- aggregate(ambiguous_reads ~ smp, tRNA, sum)

tRNA_merge <- merge(known_agg, amb_agg, by=c("smp"))
tRNA_melt <- tRNA_merge %>% melt(id.vars="smp")
tRNA_melt1 <- tRNA_melt[grepl("Pro",tRNA_melt$smp) |grepl("Tyr",tRNA_melt$smp), ]

pdf("processed_plot/input1_tRNA_read.pdf", width=60, height=4)
ggplot(tRNA_melt, aes(x=smp, y=value, fill=variable))+
geom_bar(stat = "identity", width=0.55,position = position_dodge(), alpha = 0.75)+  ####identity equals to sum; needs to make sure is the mean value; otherwise, use geom_bar(stat = "summary", fun = "mean", width=0.55)
theme_light() + ylab("reads")+xlab("")+
 theme(axis.text = element_text(color = "black")) +theme(legend.position="none") 
dev.off()

pdf("processed_plot/input1_tRNA_read_Pro_Tyr.pdf", width=10, height=5)
ggplot(tRNA_melt1, aes(x=smp, y=value, fill=variable))+
geom_bar(stat = "identity", width=0.55,position = position_dodge(), alpha = 0.75)+  ####identity equals to sum; needs to make sure is the mean value; otherwise, use geom_bar(stat = "summary", fun = "mean", width=0.55)
theme_light() + ylab("reads")+xlab("")+
 theme(axis.text = element_text(color = "black"))
dev.off()


tRNA1 <- tRNA[tRNA$known_ratio<0.35,] ##116
tRNA1 %>%  group_by(smp) %>% summarize(n = n())

############plot pairwise distance
aligndat <- read.table("resource/tRNA_high_confident_pairwise.txt", stringsAsFactors = FALSE)
colnames(aligndat) <- c("chr", "seq1", "seq2_name", "seq2", "nedit")

aligndat$smp <- gsub("_..$", "", aligndat$chr)
aligndat$smp <- gsub("_.$", "", aligndat$smp)

pdf("processed_plot/tRNA_pairwisedistance.pdf", width=60, height=4)
ggplot(aligndat, aes(x=smp, y=nedit, fill=smp))+
geom_dotplot(aes(colour = smp),binaxis = "y", stackdir = "center", position = "dodge")+  ####identity equals to sum; needs to make sure is the mean value; otherwise, use geom_bar(stat = "summary", fun = "mean", width=0.55)
theme_light() + ylab("edit_distance")+xlab("")+
 theme(axis.text = element_text(color = "black")) +theme(legend.position="none") 
dev.off()

aligndat1 <- aligndat[grepl("Pro",aligndat$smp) |grepl("Tyr",aligndat$smp), ]
pdf("processed_plot/tRNA_pairwisedistance_Pro_Tyr.pdf", width=10, height=5)
ggplot(aligndat1, aes(x=smp, y=nedit, fill=smp))+
geom_dotplot(aes(colour = smp),binaxis = "y", stackdir = "center", position = "dodge")+  ####identity equals to sum; needs to make sure is the mean value; otherwise, use geom_bar(stat = "summary", fun = "mean", width=0.55)
theme_light() + ylab("edit_distance")+xlab("")+
 theme(axis.text = element_text(color = "black")) +theme(legend.position="none") 
dev.off()

##########################################################################Some test code for tRNA and snoRNA calling
##module load R/4.0.3-foss-2020b
library(dplyr)
library(reshape2)
library(tidyr)
library(ggplot2)
library(cowplot)
library(tidyverse)
library(parallel)
setwd("/users/ludwig/ebu571/ebu571/20May2023_final")

###binomial test function
binom_t <- function(depth, mod, p) { # false positive rate/probability
  p_value <- mcmapply(as.numeric(mod), as.numeric(depth), 
                      FUN = function(x, y) {
                        if(!is.na(x) & !is.na(y)) {
                          if(y > 0) {
                            return(binom.test(x = x, n = y, p = p, alternative = "greater")$p.value)
                          } else {
                            return(1)
                          }
                        } else {
                          return(1)
                        }
                        
                      }, 
                      mc.cores = ncores)
  return(p_value)
}
############################################################### mutation calling on rRNA #########################################################
##Criteria
#####
## TPM >1 in control
##1. coverage over 20 in both treated and control 
##2. conversion rate <1% in input(background list) and the number of T-C mutation <=2
##3. Actual conversion rate above 5% (actual conversion rate=(y-fp)/(conversion-fp))
##4. Overlap in all replicates

##Sample infor: H7 vs H9; H8vs H10; H3vs H16; H23 vs H26
##input colnumn names: chr, start, end, ref, depth, T, C, mod_level, gap, gap_ratio, motif, ctrl_depth, ctrl_T, ctrl_C, ctrl_ratio, ctrl_gap, ctrl_gap_ratio
##corrected curve: y=(conversion-fp)/(1-0)*x+fp  x=(y-fp)/(conversion-fp)
##corrected curve for each motif

###read false positives
call_sites <- function(input, output1, output2, output3){
dat <- read.table(input, header=F, sep="\t", colClasses=c("V4"="character"))
###read modified and unmotified NNUNN results from figure.R
nnunn1 <- read.table("align/modified_nnunn.context.table", header=T, sep=" ")
nnunn1$motif <- paste(nnunn1$first,nnunn1$last, sep="T")

nnunn2 <- read.table("align/unmodified_nnunn.context.table", header=T, sep=" ")
nnunn2$motif <- paste(nnunn2$first,nnunn2$last, sep="T")

if(grepl("align/H7_H9", input) | grepl("align/treat2_input2", input)){
    spikein <- merge(nnunn1[,c("motif", "h7_ratio")], nnunn2[,c("motif", "h7_ratio")], by=c("motif"))
    colnames(spikein) <- c("motif", "conversion", "fp")
    spikein <- rbind(c(".", mean(spikein$conversion),mean(spikein$fp)), spikein)
} else if(grepl("align/H8_H10", input) | grepl("align/treat3_input2", input)){
    spikein <- merge(nnunn1[,c("motif", "h8_ratio")], nnunn2[,c("motif", "h8_ratio")], by=c("motif"))
    colnames(spikein) <- c("motif", "conversion", "fp")
    spikein <- rbind(c(".", mean(spikein$conversion),mean(spikein$fp)), spikein)
} else if(grepl("align/H3_H16", input) | grepl("align/H23_H26", input) | grepl("align/treat1_input1", input)){
    spikein <- merge(nnunn1[,c("motif", "h3_ratio")], nnunn2[,c("motif", "h23_ratio")], by=c("motif"))
    colnames(spikein) <- c("motif", "conversion", "fp")
    spikein <- rbind(c(".", mean(spikein$conversion),mean(spikein$fp)), spikein)
} else {print ("smaple infor not found")}


####append actual modification level and fp in given motif
colnames(dat) <- c("chr", "start", "end", "ref", "treat_depth", "treat_T", "treat_C", "treat_conversion", "treat_gap", "treat_gap_rate", "motif", "ctrl_depth", "ctrl_T", "ctrl_C", "ctrl_conversion", "ctrl_gap", "ctrl_gap_rate")
dat <- merge(dat, spikein, by=c("motif"),all.x=TRUE)
dat$treat_level <- (dat$treat_conversion-as.numeric(dat$fp))/(as.numeric(dat$conversion)-as.numeric(dat$fp))

###append TPM
if(grepl("align/treat2_input2_motif.bed", input) | grepl("align/treat3_input2_motif.bed", input)){
 tpm <- read.table("featureCounts/input2-chrM.snorna.trna.tpm", header=T, sep="\t")
 colnames(tpm) <- c("Geneid", "chr", "Start", "End", "Strand", "Length", "TPM")
 dat <- merge(dat, tpm, by=c("chr"),all.x=TRUE)

} else if(grepl("align/treat1_input1_motif.bed", input)){
  tpm <- read.table("featureCounts/input1-chrM.snorna.trna.tpm", header=T, sep="\t")
  colnames(tpm) <- c("Geneid", "chr", "Start", "End", "Strand", "Length", "TPM")
 dat <- merge(dat, tpm, by=c("chr"),all.x=TRUE)
} else if(grepl("align/treat2_input2.uniq.best_motif.bed", input) | grepl("align/treat3_input2.uniq.best_motif.bed", input)){
  tpm <- read.table("featureCounts/input2-chrM.snorna.trna.uniq.best.tpm", header=T, sep="\t")
 colnames(tpm) <- c("Geneid", "chr", "Start", "End", "Strand", "Length", "TPM")
 dat <- merge(dat, tpm, by=c("chr"),all.x=TRUE)
} else if(grepl("align/treat2_input2.mapq10_motif.bed", input) | grepl("align/treat3_input2.mapq10_motif.bed", input)){
  tpm <- read.table("featureCounts/input2-chrM.snorna.trna.mapq10.tpm", header=T, sep="\t")
 colnames(tpm) <- c("Geneid", "chr", "Start", "End", "Strand", "Length", "TPM")
 dat <- merge(dat, tpm, by=c("chr"),all.x=TRUE)
} else if(grepl("align/treat2_input2.mapq1_motif.bed", input) | grepl("align/treat3_input2.mapq1_motif.bed", input)){
  tpm <- read.table("featureCounts/input2-chrM.snorna.trna.mapq1.tpm", header=T, sep="\t")
 colnames(tpm) <- c("Geneid", "chr", "Start", "End", "Strand", "Length", "TPM")
 dat <- merge(dat, tpm, by=c("chr"),all.x=TRUE)
} else{
  dat <- dat
}


###binomial test
for (i in 1:nrow(dat)){
   fp <- dat[i,] %>% dplyr::select(fp) %>% as.numeric()
   dat$p[i] <- binom_t(depth = (dat$treat_T[i]+dat$treat_C[i]), mod = dat$treat_C[i], p = fp)
}

dat$fdr <- p.adjust(dat$p, method = "BH")
##called sites
#site <- dat[(dat$treat_T+dat$treat_C)>=20 & (dat$ctrl_T+dat$ctrl_C)>=20 & dat$TPM >1 & as.numeric(dat$ctrl_conversion)<0.01 & as.numeric(dat$treat_level)>=0.05,]
site <- dat[(dat$treat_T+dat$treat_C)>=20 & (dat$ctrl_T+dat$ctrl_C)>=20 & as.numeric(dat$treat_level)>=0.05,] %>% na.omit()
site1 <- site[as.numeric(site$ctrl_conversion)<0.01 | as.numeric(site$ctrl_C)<2,] %>% na.omit()
background <- dat[as.numeric(dat$ctrl_conversion)>=0.01 & as.numeric(dat$ctrl_C)>=2,]

write.table(site1, output1, quote = FALSE, col.names = T, row.names = F)
write.table(dat[(dat$treat_T+dat$treat_C)>=20 & (dat$ctrl_T+dat$ctrl_C)>=20,], output2, quote = FALSE, col.names = T, row.names = F)
write.table(background, output3, quote = FALSE, col.names = T, row.names = F)
}

#######includes ambiguous alignment
call_sites("align/treat1_input1_motif.bed", "site/treat1_input1_called_sites.txt", "site/treat1_input1_all_sites.txt", "site/treat1_input1_background_sites.txt")
call_sites("align/treat2_input2_motif.bed", "site/treat2_input2_called_sites.txt", "site/treat2_input2_all_sites.txt", "site/treat2_input2_background_sites.txt")
call_sites("align/treat3_input2_motif.bed", "site/treat3_input2_called_sites.txt", "site/treat3_input2_all_sites.txt", "site/treat3_input2_background_sites.txt")

site1 <- read.table("site/treat1_input1_called_sites.txt", header=T, sep=" ", colClasses="character" ) ### 1251
site2 <- read.table("site/treat2_input2_called_sites.txt", header=T, sep=" ", colClasses="character" ) ###1279
site3 <- read.table("site/treat3_input2_called_sites.txt", header=T, sep=" ", colClasses="character" ) ###1274

merge1 <- list(site1[,c("chr", "start", "end", "ref", "motif", "fp", "treat_conversion", "treat_level", "ctrl_conversion","p", "fdr", "TPM")], site2[,c("chr", "start", "end", "ref", "motif", "fp", "treat_conversion", "treat_level","p", "fdr")], site3[,c("chr", "start", "end", "ref", "motif", "fp", "treat_conversion", "treat_level", "ctrl_conversion","p", "fdr", "TPM")])
merge1 <- merge1 %>% reduce(full_join, by=c("chr", "start", "end", "ref", "motif")) %>% na.omit() ###1169
colnames(merge1) <- c("chr", "start", "end", "ref", "motif", "treat1_fp", "treat1_conversion", "treat1_level", "ctrl1_conversion","treat1_p", "treat1_fdr", "input1_TPM", "treat2_fp", "treat2_conversion", "treat2_level","treat2_p", "treat2_fdr", "treat3_fp", "treat3_conversion", "treat3_level", "ctrl2_conversion","treat3_p", "treat3_fdr", "input2_TPM")
nrow(merge1) #1169
write.table(merge1, "site/snoRNA_tRNA_called_sites.txt", quote = FALSE, col.names = T, row.names = F)


merge2 <- merge1 %>% reduce(full_join, by=c("chr", "start", "end", "ref", "motif")) 
colnames(merge2) <- c("chr", "start", "end", "ref", "motif", "treat1_fp", "treat1_conversion", "treat1_level", "ctrl1_conversion","treat1_p", "treat1_fdr", "input1_TPM", "treat2_fp", "treat2_conversion", "treat2_level","treat2_p", "treat2_fdr", "treat3_fp", "treat3_conversion", "treat3_level", "ctrl2_conversion","treat3_p", "treat3_fdr", "input2_TPM")
nrow(merge2) ## 1386
write.table(merge2, "site/snoRNA_tRNA_called_sites_na.txt", quote = FALSE, col.names = T, row.names = F)


##############background#####################################
input1 <- read.table("site/treat1_input1_background_sites.txt", header=T, sep=" ",colClasses="character")
input2 <- read.table("site/treat2_input2_background_sites.txt", header=T, sep=" ",colClasses="character")
input3 <- read.table("site/treat3_input2_background_sites.txt", header=T, sep=" ",colClasses="character")

merge1 <- list(input1[,c("chr", "start", "end", "ref", "motif", "fp", "treat_conversion", "treat_level", "ctrl_conversion","p", "fdr", "TPM")], input2[,c("chr", "start", "end", "ref", "motif", "fp", "treat_conversion", "treat_level","p", "fdr")], input3[,c("chr", "start", "end", "ref", "motif", "fp", "treat_conversion", "treat_level", "ctrl_conversion","p", "fdr", "TPM")])
merge1 <- merge1 %>% reduce(full_join, by=c("chr", "start", "end", "ref", "motif"))
colnames(merge1) <- c("chr", "start", "end", "ref", "motif", "treat1_fp", "treat1_conversion", "treat1_level", "ctrl1_conversion","treat1_p", "treat1_fdr", "input1_TPM", "treat2_fp", "treat2_conversion", "treat2_level","treat2_p", "treat2_fdr", "treat3_fp", "treat3_conversion", "treat3_level", "ctrl2_conversion","treat3_p", "treat3_fdr", "input2_TPM")
nrow(merge1) #233
write.table(merge1, "site/snoRNA_tRNA_background_sites.txt", quote = FALSE, col.names = T, row.names = F)


######all sites
input1 <- read.table("site/treat1_input1_all_sites.txt", header=T, sep=" ",colClasses="character")
input2 <- read.table("site/treat2_input2_all_sites.txt", header=T, sep=" ",colClasses="character")
input3 <- read.table("site/treat3_input2_all_sites.txt", header=T, sep=" ",colClasses="character")

merge1 <- list(input1[,c("chr", "start", "end", "ref", "motif", "fp", "treat_conversion", "treat_level", "ctrl_conversion","p", "fdr", "TPM")], input2[,c("chr", "start", "end", "ref", "motif", "fp", "treat_conversion", "treat_level","p", "fdr")], input3[,c("chr", "start", "end", "ref", "motif", "fp", "treat_conversion", "treat_level", "ctrl_conversion","p", "fdr", "TPM")])
merge1 <- merge1 %>% reduce(full_join, by=c("chr", "start", "end", "ref", "motif")) %>% na.omit()
colnames(merge1) <- c("chr", "start", "end", "ref", "motif", "treat1_fp", "treat1_conversion", "treat1_level", "ctrl1_conversion","treat1_p", "treat1_fdr", "input1_TPM", "treat2_fp", "treat2_conversion", "treat2_level","treat2_p", "treat2_fdr", "treat3_fp", "treat3_conversion", "treat3_level", "ctrl2_conversion","treat3_p", "treat3_fdr", "input2_TPM")
nrow(merge1) #15413
write.table(merge1, "site/snoRNA_tRNA_all_sites.txt", quote = FALSE, col.names = T, row.names = F)



################################################### for uniq and best sitement ,only for treat2 and treat3############################################
call_sites("align/treat1_input1.uniq.best_motif.bed", "site/treat1_input1.uniq.best_called_sites.txt", "site/treat1_input1.uniq.best_all_sites.txt", "site/treat1_input1.uniq.best_background_sites.txt")
call_sites("align/treat2_input2.uniq.best_motif.bed", "site/treat2_input2.uniq.best_called_sites.txt", "site/treat2_input2.uniq.best_all_sites.txt", "site/treat2_input2.uniq.best_background_sites.txt")
call_sites("align/treat3_input2.uniq.best_motif.bed", "site/treat3_input2.uniq.best_called_sites.txt", "site/treat3_input2.uniq.best_all_sites.txt", "site/treat3_input2.uniq.best_background_sites.txt")

site1 <- read.table("site/treat1_input1.uniq.best_called_sites.txt", header=T, sep=" ", colClasses="character" ) ### 1251
site2 <- read.table("site/treat2_input2.uniq.best_called_sites.txt", header=T, sep=" ", colClasses="character" ) ###1279
site3 <- read.table("site/treat3_input2.uniq.best_called_sites.txt", header=T, sep=" ", colClasses="character" ) ###1274


merge1 <- list(site2[,c("chr", "start", "end", "ref", "motif", "fp", "treat_conversion", "treat_level","p", "fdr")], site3[,c("chr", "start", "end", "ref", "motif", "fp", "treat_conversion", "treat_level", "ctrl_conversion","p", "fdr", "TPM")])
merge1 <- merge1 %>% reduce(full_join, by=c("chr", "start", "end", "ref", "motif")) %>% na.omit() ###1169
colnames(merge1) <- c("chr", "start", "end", "ref", "motif", "treat2_fp", "treat2_conversion", "treat2_level","treat2_p", "treat2_fdr", "treat3_fp", "treat3_conversion", "treat3_level", "ctrl2_conversion","treat3_p", "treat3_fdr", "input2_TPM")
nrow(merge1) #1088
write.table(merge1, "site/snoRNA_tRNA.uniq.best_called_sites.txt", quote = FALSE, col.names = T, row.names = F)

merge2 <- merge1 %>% reduce(full_join, by=c("chr", "start", "end", "ref", "motif")) 
colnames(merge2) <- c("chr", "start", "end", "ref", "motif", "treat2_fp", "treat2_conversion", "treat2_level","treat2_p", "treat2_fdr", "treat3_fp", "treat3_conversion", "treat3_level", "ctrl2_conversion","treat3_p", "treat3_fdr", "input2_TPM")
nrow(merge2) ## 1227
write.table(merge2, "site/snoRNA_tRNA.uniq.best_called_sites_na.txt", quote = FALSE, col.names = T, row.names = F)


##############background#####################################
input2 <- read.table("site/treat2_input2.uniq.best_background_sites.txt", header=T, sep=" ",colClasses="character")
input3 <- read.table("site/treat3_input2.uniq.best_background_sites.txt", header=T, sep=" ",colClasses="character")

merge1 <- list(input2[,c("chr", "start", "end", "ref", "motif", "fp", "treat_conversion", "treat_level","p", "fdr")], input3[,c("chr", "start", "end", "ref", "motif", "fp", "treat_conversion", "treat_level", "ctrl_conversion","p", "fdr", "TPM")])
merge1 <- merge1 %>% reduce(full_join, by=c("chr", "start", "end", "ref", "motif"))
colnames(merge1) <- c("chr", "start", "end", "ref", "motif","treat2_fp", "treat2_conversion", "treat2_level","treat2_p", "treat2_fdr", "treat3_fp", "treat3_conversion", "treat3_level", "ctrl2_conversion","treat3_p", "treat3_fdr", "input2_TPM")
nrow(merge1) #156
write.table(merge1, "site/snoRNA_tRNA.uniq.best_background_sites.txt", quote = FALSE, col.names = T, row.names = F)


######all sites
input2 <- read.table("site/treat2_input2.uniq.best_all_sites.txt", header=T, sep=" ",colClasses="character")
input3 <- read.table("site/treat3_input2.uniq.best_all_sites.txt", header=T, sep=" ",colClasses="character")

merge1 <- list(input2[,c("chr", "start", "end", "ref", "motif", "fp", "treat_conversion", "treat_level","p", "fdr")], input3[,c("chr", "start", "end", "ref", "motif", "fp", "treat_conversion", "treat_level", "ctrl_conversion","p", "fdr", "TPM")])
merge1 <- merge1 %>% reduce(full_join, by=c("chr", "start", "end", "ref", "motif")) %>% na.omit()
colnames(merge1) <- c("chr", "start", "end", "ref", "motif", "treat2_fp", "treat2_conversion", "treat2_level","treat2_p", "treat2_fdr", "treat3_fp", "treat3_conversion", "treat3_level", "ctrl2_conversion","treat3_p", "treat3_fdr", "input2_TPM")
nrow(merge1) #14627
write.table(merge1, "site/snoRNA_tRNA.uniq.best_all_sites.txt", quote = FALSE, col.names = T, row.names = F)



######################## useless
exclude1 <- read.table("align/treat2_snoRNA_tRNA_exclude.txt", header=F, sep="\t", colClasses="character" )  ##37
colnames(exclude1) <- c("chr", "pos", "ref", "A", "C", "G", "A_ratio", "C_ratio", "G_ratio")
exclude2 <- read.table("align/treat3_snoRNA_tRNA_exclude.txt", header=F, sep="\t", colClasses="character" ) ###24
colnames(exclude2) <- c("chr", "pos", "ref", "A", "C", "G", "A_ratio", "C_ratio", "G_ratio")
exclude <- merge(exclude1, exclude2, by=c("chr", "pos"), all.x=TRUE, all.y=TRUE) ##52
exclude$start <- as.numeric(exclude$pos) -1
exclude$end <- exclude$po
merge(exclude, merge1, by=c("chr", "start", "end"))




#########################################
################################################### for MAPQ 10############################################
call_sites("align/treat1_input1.mapq10_motif.bed", "site/treat1_input1.mapq10_called_sites.txt", "site/treat1_input1.mapq10_all_sites.txt", "site/treat1_input1.mapq10_background_sites.txt")
call_sites("align/treat2_input2.mapq10_motif.bed", "site/treat2_input2.mapq10_called_sites.txt", "site/treat2_input2.mapq10_all_sites.txt", "site/treat2_input2.mapq10_background_sites.txt")
call_sites("align/treat3_input2.mapq10_motif.bed", "site/treat3_input2.mapq10_called_sites.txt", "site/treat3_input2.mapq10_all_sites.txt", "site/treat3_input2.mapq10_background_sites.txt")

site1 <- read.table("site/treat1_input1.mapq10_called_sites.txt", header=T, sep=" ", colClasses="character" ) ### 1035
site2 <- read.table("site/treat2_input2.mapq10_called_sites.txt", header=T, sep=" ", colClasses="character" ) ###1107
site3 <- read.table("site/treat3_input2.mapq10_called_sites.txt", header=T, sep=" ", colClasses="character" ) ###1104


merge1 <- list(site2[,c("chr", "start", "end", "ref", "motif", "fp", "treat_conversion", "treat_level","p", "fdr")], site3[,c("chr", "start", "end", "ref", "motif", "fp", "treat_conversion", "treat_level", "ctrl_conversion","p", "fdr", "TPM")])
merge1 <- merge1 %>% reduce(full_join, by=c("chr", "start", "end", "ref", "motif")) %>% na.omit() ###1169
colnames(merge1) <- c("chr", "start", "end", "ref", "motif", "treat2_fp", "treat2_conversion", "treat2_level","treat2_p", "treat2_fdr", "treat3_fp", "treat3_conversion", "treat3_level", "ctrl2_conversion","treat3_p", "treat3_fdr", "input2_TPM")
nrow(merge1) #1031
write.table(merge1, "site/snoRNA_tRNA.mapq10_called_sites.txt", quote = FALSE, col.names = T, row.names = F)

################################################### for MAPQ 1############################################
call_sites("align/treat1_input1.mapq1_motif.bed", "site/treat1_input1.mapq1_called_sites.txt", "site/treat1_input1.mapq1_all_sites.txt", "site/treat1_input1.mapq1_background_sites.txt")
call_sites("align/treat2_input2.mapq1_motif.bed", "site/treat2_input2.mapq1_called_sites.txt", "site/treat2_input2.mapq1_all_sites.txt", "site/treat2_input2.mapq1_background_sites.txt")
call_sites("align/treat3_input2.mapq1_motif.bed", "site/treat3_input2.mapq1_called_sites.txt", "site/treat3_input2.mapq1_all_sites.txt", "site/treat3_input2.mapq1_background_sites.txt")

site1 <- read.table("site/treat1_input1.mapq1_called_sites.txt", header=T, sep=" ", colClasses="character" ) ### 1206
site2 <- read.table("site/treat2_input2.mapq1_called_sites.txt", header=T, sep=" ", colClasses="character" ) ###1255
site3 <- read.table("site/treat3_input2.mapq1_called_sites.txt", header=T, sep=" ", colClasses="character" ) ###1244


merge1 <- list(site2[,c("chr", "start", "end", "ref", "motif", "fp", "treat_conversion", "treat_level","p", "fdr")], site3[,c("chr", "start", "end", "ref", "motif", "fp", "treat_conversion", "treat_level", "ctrl_conversion","p", "fdr", "TPM")])
merge1 <- merge1 %>% reduce(full_join, by=c("chr", "start", "end", "ref", "motif")) %>% na.omit() ###1169
colnames(merge1) <- c("chr", "start", "end", "ref", "motif", "treat2_fp", "treat2_conversion", "treat2_level","treat2_p", "treat2_fdr", "treat3_fp", "treat3_conversion", "treat3_level", "ctrl2_conversion","treat3_p", "treat3_fdr", "input2_TPM")
nrow(merge1) #1193
write.table(merge1, "site/snoRNA_tRNA.mapq1_called_sites.txt", quote = FALSE, col.names = T, row.names = F)




##############################plot TPM of two input
tpm1 <- read.table("featureCounts/input1.tpm", header=T, sep="\t")
colnames(tpm1) <- c("geneid", "chr", "start","end","strand","length","input1", "input2", "treat1", "treat2", "treat3")
####get tpm only for tRNA
tpm1 <- tpm1[grepl("tRNA", tpm1$chr) & !grepl("mt_tRNA", tpm1$chr),] %>% dplyr::select(chr,input1, input2, treat1, treat2, treat3)

meltdat <- melt(tpm1[,c("chr","input1", "input2")], id.vars="chr")

pdf("processed_plot/tRNA_TPM.pdf")
ggplot(meltdat, aes(x=variable, y=value,fill=variable))+
geom_boxplot()+  ####identity equals to sum; needs to make sure is the mean value; otherwise, use geom_bar(stat = "summary", fun = "mean", width=0.55)
theme_light() +ylim(0,7500)
 theme(axis.text = element_text(color = "black")) +theme(legend.position="none") 
dev.off()

###############plot coverage over anticodon
meltdat$smp <- sub("tRNA_(\\w+_\\w+)_\\w+", "\\1", meltdat$chr)
tpm1$smp <- sub("tRNA_(\\w+_\\w+)_\\w+", "\\1", tpm1$chr)
pdf("processed_plot/tRNA_TPM_anticodon.pdf", width=80, height=4)
ggplot(tpm1, aes(x=smp, y=log2(input1)))+
geom_point()+  ####identity equals to sum; needs to make sure is the mean value; otherwise, use geom_bar(stat = "summary", fun = "mean", width=0.55)
theme_light() +
 theme(axis.text = element_text(color = "black")) +theme(legend.position="none") 
dev.off()









################# output anticodons whose TPM is below the 1st quantile of all anticodons
summary(tpm1$input1)
 #   Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
 #   0.09   158.58   747.00  3127.18  2723.67 86790.31 

summary(tpm1$input2)
  #  Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
  #  0.29   172.30   830.42  3116.93  2871.31 61897.74 

nrow(tpm1[tpm1$input1<=158.58,]) ##65
nrow(tpm1[tpm1$input2<=172.30,]) ##65
write.table(tpm1[tpm1$input2<=172.30,], "site/input2_TPM_lower_quantile.txt", quote = FALSE, col.names = T, row.names = F)

merge(tpm1[tpm1$input1<=158.58,], tpm1[tpm1$input2<=172.30,], by=c("chr"),all.x=TRUE, all.y=TRUE) ###63


######fdr filtered

cat snoRNA_tRNA_called_sites.txt | grep -v ^chr | sed 's/ /\t/g' | grep tRNA | grep -v mt_tRNA | awk '$11<0.001 && $18<0.001 && $25<0.001' > snoRNA_tRNA_called_sites_fdr_filtered.txt ###756
cat snoRNA_tRNA_called_sites.txt | grep -v ^chr | sed 's/ /\t/g' | grep tRNA | grep -v mt_tRNA | awk '$11>=0.001 || $18>=0.001 || $25>=0.001' > snoRNA_tRNA_called_sites_fdr_unpassed.txt ###26

site1 <- read.table("site/snoRNA_tRNA_called_sites.txt", , header=T, sep=" ")
trna <- site1[grepl("tRNA", site1$chr) & !grepl("mt_tRNA", site1$chr),] ###782

site <- trna[trna$treat1_fdr<0.001 & trna$treat2_fdr<0.001 & trna$treat3_fdr<0.001,] ##768
write.table(site, "site/snoRNA_tRNA_called_sites_fdr_filtered.txt", quote = FALSE, col.names = F, row.names = F)

exclude <- trna[trna$treat1_fdr>=0.001 | trna$treat2_fdr>=0.001 | trna$treat3_fdr>=0.001,] ###14
write.table(exclude, "site/snoRNA_tRNA_called_sites_fdr_unpassed.txt", quote = FALSE, col.names = F, row.names = F)
