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

if(grepl("align/H7_H9_allT_motif.bed", input) | grepl("align/treat2_input2_motif.bed", input)){
    spikein <- merge(nnunn1[,c("motif", "h7_ratio")], nnunn2[,c("motif", "h7_ratio")], by=c("motif"))
    colnames(spikein) <- c("motif", "conversion", "fp")
    spikein <- rbind(c(".", mean(spikein$conversion),mean(spikein$fp)), spikein)
} else if(grepl("align/H8_H10_allT_motif.bed", input) | grepl("align/treat3_input2_motif.bed", input)){
    spikein <- merge(nnunn1[,c("motif", "h8_ratio")], nnunn2[,c("motif", "h8_ratio")], by=c("motif"))
    colnames(spikein) <- c("motif", "conversion", "fp")
    spikein <- rbind(c(".", mean(spikein$conversion),mean(spikein$fp)), spikein)
} else if(grepl("align/H3_H16_allT_motif.bed", input) | grepl("align/H23_H26_allT_motif.bed", input) | grepl("align/treat1_input1_motif.bed", input)){
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
 tpm1 <- read.table("featureCounts/input2-chrM.snorna.tpm", header=T, sep="\t")
 colnames(tpm1) <- c("Geneid", "chr", "Start", "End", "Strand", "Length", "TPM")
 tpm2 <- read.table("featureCounts/input2-chrM.tRNA.tpm", header=T, sep="\t")
 colnames(tpm2) <- c("Geneid", "chr", "Start", "End", "Strand", "Length", "TPM")
 tpm <- rbind(tpm1[,c("chr","TPM")], tpm2[,c("chr","TPM")])
 dat <- merge(dat, tpm, by=c("chr"),all.x=TRUE)

} else if(grepl("align/treat1_input1_motif.bed", input)){
  tpm1 <- read.table("featureCounts/input1-chrM.snorna.tpm", header=T, sep="\t")
  colnames(tpm1) <- c("Geneid", "chr", "Start", "End", "Strand", "Length", "TPM")
  tpm2 <- read.table("featureCounts/input1-chrM.tRNA.tpm", header=T, sep="\t")
  colnames(tpm2) <- c("Geneid", "chr", "Start", "End", "Strand", "Length", "TPM")
 tpm <- rbind(tpm1[,c("chr","TPM")], tpm2[,c("chr","TPM")])
 dat <- merge(dat, tpm, by=c("chr"),all.x=TRUE)

} else {
  dat <- dat
}


###binomial test
#for (i in 1:nrow(dat)){
 #  fp <- dat[i,] %>% dplyr::select(fp) %>% as.numeric()
 #  dat$p[i] <- binom_t(depth = (dat$treat_T[i]+dat$treat_C[i]), mod = dat$treat_C[i], p = fp)
#}
##& as.numeric(dat$ctrl_C)<=2
#dat$fdr <- p.adjust(dat$p, method = "BH")
###called sites
#site <- dat[(dat$treat_T+dat$treat_C)>=20 & (dat$ctrl_T+dat$ctrl_C)>=20 & dat$TPM >1 & as.numeric(dat$ctrl_conversion)<0.01 & as.numeric(dat$treat_level)>=0.05,]
site <- dat[(dat$treat_T+dat$treat_C)>=20 & (dat$ctrl_T+dat$ctrl_C)>=20 & as.numeric(dat$treat_level)>=0.05,] %>% na.omit()
site1 <- site[as.numeric(site$ctrl_conversion)<0.01 | as.numeric(site$ctrl_C)<2,] %>% na.omit()
background <- dat[as.numeric(dat$ctrl_conversion)>=0.01 & as.numeric(dat$ctrl_C)>=2,]

write.table(site1, output1, quote = FALSE, col.names = T, row.names = F)
write.table(site, output2, quote = FALSE, col.names = T, row.names = F)
write.table(background, output3, quote = FALSE, col.names = T, row.names = F)
}


call_sites("align/H7_H9_allT_motif.bed", "align/H7_H9_called_sites.txt","align/H7_H9_all_sites.txt", "align/H7_H9_background.txt")
call_sites("align/H8_H10_allT_motif.bed", "align/H8_H10_called_sites.txt","align/H8_H10_all_sites.txt", "align/H8_H10_background.txt")
call_sites("align/H3_H16_allT_motif.bed", "align/H3_H16_called_sites.txt","align/H3_H16_all_sites.txt", "align/H3_H16_background.txt")
call_sites("align/H23_H26_allT_motif.bed", "align/H23_H26_called_sites.txt","align/H23_H26_all_sites.txt","align/H23_H26_background.txt")


########merge sites on rRNA
site1 <- read.table("align/H7_H9_called_sites.txt", header=T, sep=" ", colClasses="character" ) ##162
site2 <- read.table("align/H8_H10_called_sites.txt", header=T, sep=" ", colClasses="character" ) ##161
site3 <- read.table("align/H3_H16_called_sites.txt", header=T, sep=" ", colClasses="character" ) ##158
site4 <- read.table("align/H23_H26_called_sites.txt", header=T, sep=" ",colClasses="character") ##156

tmp1 <- merge(site1[,c("chr", "start", "end", "ref", "motif", "fp", "treat_conversion", "treat_level", "ctrl_conversion")], site2[,c("chr", "start", "end", "ref", "motif", "fp", "treat_conversion", "treat_level", "ctrl_conversion")], by=c("chr", "start", "end", "ref", "motif")) %>% na.omit()
colnames(tmp1) <- c("chr", "start", "end", "ref", "motif", "h7_fp", "h7_conversion", "h7_level", "h9_conversion", "h8_fp","h8_conversion", "h8_level", "h10_conversion")

tmp2 <- merge(tmp1, site3[,c("chr", "start", "end", "ref", "motif", "fp", "treat_conversion", "treat_level", "ctrl_conversion")], by=c("chr", "start", "end", "ref", "motif")) %>% na.omit()
colnames(tmp2) <- c("chr", "start", "end", "ref", "motif", "h7_fp", "h7_conversion", "h7_level", "h9_conversion", "h8_fp", "h8_conversion", "h8_level", "h10_conversion", "h3_fp","h3_conversion", "h3_level", "h16_conversion")

merge_dat <- merge(tmp2, site4[,c("chr", "start", "end", "ref", "motif", "fp", "treat_conversion", "treat_level", "ctrl_conversion")], by=c("chr", "start", "end", "ref", "motif")) %>% na.omit()
colnames(merge_dat) <- c("chr", "start", "end", "ref", "motif", "h7_fp", "h7_conversion", "h7_level", "h9_conversion", "h8_fp", "h8_conversion", "h8_level", "h10_conversion", "h3_fp","h3_conversion", "h3_level", "h16_conversion", "h23_fp","h23_conversion", "h23_level", "h26_conversion")

nrow(merge_dat)###152
write.table(merge_dat, "align/rRNA_called_sites.txt", quote = FALSE, col.names = T, row.names = F)

#######Find background mutation
###coverage higher than 20
###mutation ratio higher than 1%
input1 <- read.table("align/H7_H9_background.txt", header=T, sep=" ",colClasses="character")
input2 <- read.table("align/H8_H10_background.txt", header=T, sep=" ",colClasses="character")
input3 <- read.table("align/H3_H16_background.txt", header=T, sep=" ",colClasses="character")
input4 <- read.table("align/H23_H26_background.txt", header=T, sep=" ",colClasses="character")

tmp1 <- merge(input1[,c("chr", "start", "end", "ref", "motif", "fp", "treat_conversion", "treat_level", "ctrl_conversion")], input2[,c("chr", "start", "end", "ref", "motif", "fp", "treat_conversion", "treat_level", "ctrl_conversion")], by=c("chr", "start", "end", "ref", "motif"), all.x=T, all.y=T) 
colnames(tmp1) <- c("chr", "start", "end", "ref", "motif", "h7_fp", "h7_conversion", "h7_level", "h9_conversion", "h8_fp","h8_conversion", "h8_level", "h10_conversion")

tmp2 <- merge(tmp1, input3[,c("chr", "start", "end", "ref", "motif", "fp", "treat_conversion", "treat_level", "ctrl_conversion")], by=c("chr", "start", "end", "ref", "motif"),all.x=T, all.y=T) 
colnames(tmp2) <- c("chr", "start", "end", "ref", "motif", "h7_fp", "h7_conversion", "h7_level", "h9_conversion", "h8_fp", "h8_conversion", "h8_level", "h10_conversion", "h3_fp","h3_conversion", "h3_level", "h16_conversion")

merge_dat <- merge(tmp2, input4[,c("chr", "start", "end", "ref", "motif", "fp", "treat_conversion", "treat_level", "ctrl_conversion")], by=c("chr", "start", "end", "ref", "motif"),all.x=T, all.y=T) 
colnames(merge_dat) <- c("chr", "start", "end", "ref", "motif", "h7_fp", "h7_conversion", "h7_level", "h9_conversion", "h8_fp", "h8_conversion", "h8_level", "h10_conversion", "h3_fp","h3_conversion", "h3_level", "h16_conversion", "h23_fp","h23_conversion", "h23_level", "h26_conversion")

nrow(merge_dat)###19
write.table(merge_dat, "align/rRNA_background.txt", quote = FALSE, col.names = T, row.names = F)






















h7 <- read.table("align/H7_H9_allT_motif.bed", header=F, sep="\t")
h8 <- read.table("align/H8_H10_allT_motif.bed", header=F, sep="\t")
h3 <- read.table("align/H3_H16_allT_motif.bed", header=F, sep="\t")
h23 <- read.table("align/H23_H26_allT_motif.bed", header=F, sep="\t")

####H7 vs H9
colnames(h7) <- c("chr", "start", "end", "ref", "h7_depth", "h7_T", "h7_C", "h7_conversion", "h7_gap", "h7_gap_rate", "motif", "h9_depth", "h9_T", "h9_C", "h9_conversion", "h9_gap", "h9_gap_rate")
h7$h7_level <- (h7$h7_conversion+0.0026)/0.8211
###Filtering: coverage over 20, fold change higher than 1.5, actual conversion rate > 5%
site1 <- h7[(h7$h7_T+h7$h7_C)>=20 & (h7$h9_T+h7$h9_C)>=20 & (h7$h7_conversion-h7$h9_conversion)/h7$h9_conversion>=1.5 & h7$h7_level>=0.05,]

####H8 vs H10
colnames(h8) <- c("chr", "start", "end", "ref", "h8_depth", "h8_T", "h8_C", "h8_conversion", "h8_gap", "h8_gap_rate", "motif", "h10_depth", "h10_T", "h10_C", "h10_conversion", "h10_gap", "h10_gap_rate")
h8$h8_level <- (h8$h8_conversion+0.0026)/0.8211
###Filtering: coverage over 20, fold change higher than 1.5, actual conversion rate > 5%
site2 <- h8[(h8$h8_T+h8$h8_C)>=20 & (h8$h10_T+h8$h10_C)>=20 & (h8$h8_conversion-h8$h10_conversion)/h8$h10_conversion>=1.5 & h8$h8_level>=0.05,]

###H3 vs H16
colnames(h3) <- c("chr", "start", "end", "ref", "h3_depth", "h3_T", "h3_C", "h3_conversion", "h3_gap", "h3_gap_rate", "motif", "h16_depth", "h16_T", "h16_C", "h16_conversion", "h16_gap", "h16_gap_rate")
h3$h3_level <- (h3$h3_conversion+0.0026)/0.8211
###Filtering: coverage over 20, fold change higher than 1.5, actual conversion rate > 5%
site3 <- h3[(h3$h3_T+h3$h3_C)>=20 & (h3$h16_T+h3$h16_C)>=20 & (h3$h3_conversion-h3$h16_conversion)/h3$h16_conversion>=1.5 & h3$h3_level>=0.05,]

###H23 vs H26
colnames(h23) <- c("chr", "start", "end", "ref", "h23_depth", "h23_T", "h23_C", "h23_conversion", "h23_gap", "h23_gap_rate", "motif", "h26_depth", "h26_T", "h26_C", "h26_conversion", "h26_gap", "h26_gap_rate")
h23$h23_level <- (h23$h23_conversion+0.0026)/0.8211
###Filtering: coverage over 20, fold change higher than 1.5, actual conversion rate > 5%
site4 <- h23[(h23$h23_T+h23$h23_C)>=20 & (h23$h26_T+h23$h26_C)>=20 & (h23$h23_conversion-h23$h26_conversion)/h23$h26_conversion>=1.5 & h23$h23_level>=0.05,]


tmp1 <- merge(site1, site2, by=c("chr", "start", "end", "ref", "motif"))
tmp2 <- merge(tmp1, site3, by=c("chr", "start", "end", "ref", "motif"))
merge_dat <- merge(tmp2, site4, by=c("chr", "start", "end", "ref", "motif"))
called_sites <- merge_dat
nrow(called_sites)
write.table(called_sites, "align/rRNA_called_sites.txt", quote = FALSE, col.names = T, row.names = F)
