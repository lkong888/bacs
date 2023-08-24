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
##1. coverage over 20 in both treated and control 
##2. conversion rate <1% in input(background list) or the number of T-C mutation <2
##3. Actual conversion rate above 5% (actual conversion rate=(y-fp)/(conversion-fp))
##4. Overlap in all replicates

##Sample infor: H7 vs H9; H8vs H10; H3vs H16; H23 vs H26
##input colnumn names: chr, start, end, ref, depth, T, C, mod_level, gap, gap_ratio, motif, ctrl_depth, ctrl_T, ctrl_C, ctrl_ratio, ctrl_gap, ctrl_gap_ratio
##corrected curve: y=(conversion-fp)/(1-0)*x+fp  x=(y-fp)/(conversion-fp)
##fp and conversion are motif-specific

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
if(nrow(dat[grepl("N",dat$motif),])>0){  ###mt-reference contains N
    for (i in 1:nrow(dat)){
      if(grepl("N",dat[i,]$motif)){
        dat[i,]$motif <- "."
      }
    }
}
dat <- merge(dat, spikein, by=c("motif"),all.x=TRUE)
dat$treat_level <- (dat$treat_conversion-as.numeric(dat$fp))/(as.numeric(dat$conversion)-as.numeric(dat$fp))

###append TPM
if(grepl("align/treat2_input2_motif.bed", input)){
 tpm1 <- read.table("featureCounts/input2-chrM.snorna.trna.tpm", header=T, sep="\t")
 colnames(tpm1) <- c("Geneid", "chr", "Start", "End", "Strand", "Length", "TPM_input2")
 dat <- merge(dat, tpm1, by=c("chr"),all.x=TRUE)
  tpm2 <- read.table("featureCounts/treat2-chrM.snorna.trna.tpm", header=T, sep="\t")
 colnames(tpm2) <- c("Geneid", "chr", "Start", "End", "Strand", "Length", "TPM_treat2")
 dat <- merge(dat, tpm2, by=c("chr"),all.x=TRUE)

} else if(grepl("align/treat1_input1_motif.bed", input)){
   tpm1 <- read.table("featureCounts/input1-chrM.snorna.trna.tpm", header=T, sep="\t")
 colnames(tpm1) <- c("Geneid", "chr", "Start", "End", "Strand", "Length", "TPM_input1")
 dat <- merge(dat, tpm1, by=c("chr"),all.x=TRUE)
  tpm2 <- read.table("featureCounts/treat1-chrM.snorna.trna.tpm", header=T, sep="\t")
 colnames(tpm2) <- c("Geneid", "chr", "Start", "End", "Strand", "Length", "TPM_treat1")
 dat <- merge(dat, tpm2, by=c("chr"),all.x=TRUE)
} else if(grepl("align/treat3_input2_motif.bed", input)){
 tpm1 <- read.table("featureCounts/input2-chrM.snorna.trna.tpm", header=T, sep="\t")
 colnames(tpm1) <- c("Geneid", "chr", "Start", "End", "Strand", "Length", "TPM_input2")
 dat <- merge(dat, tpm1, by=c("chr"),all.x=TRUE)
  tpm2 <- read.table("featureCounts/treat3-chrM.snorna.trna.tpm", header=T, sep="\t")
 colnames(tpm2) <- c("Geneid", "chr", "Start", "End", "Strand", "Length", "TPM_treat3")
 dat <- merge(dat, tpm2, by=c("chr"),all.x=TRUE)
} else {
  dat <- dat
}


###binomial test
for (i in 1:nrow(dat)){
   fp <- dat[i,] %>% dplyr::select(fp) %>% as.numeric()
   dat$p[i] <- binom_t(depth = (dat$treat_T[i]+dat$treat_C[i]), mod = dat$treat_C[i], p = fp)
}

dat$fdr <- p.adjust(dat$p, method = "BH")
##called sites
site <- dat[(dat$treat_T+dat$treat_C)>=20 & (dat$ctrl_T+dat$ctrl_C)>=20 & as.numeric(dat$treat_level)>=0.05,] %>% na.omit()
site1 <- site[as.numeric(site$ctrl_conversion)<0.01 | as.numeric(site$ctrl_C)<2,] %>% na.omit()
background <- dat[as.numeric(dat$ctrl_conversion)>=0.01 & as.numeric(dat$ctrl_C)>=2,]

write.table(site1, output1, quote = FALSE, col.names = T, row.names = F)
write.table(dat[(dat$treat_T+dat$treat_C)>=20 & (dat$ctrl_T+dat$ctrl_C)>=20,], output2, quote = FALSE, col.names = T, row.names = F)
write.table(background, output3, quote = FALSE, col.names = T, row.names = F)
}

call_sites("align/treat1_input1_motif.bed", "site/treat1_input1_called_sites.txt", "site/treat1_input1_all_sites.txt", "site/treat1_input1_background_sites.txt")
call_sites("align/treat2_input2_motif.bed", "site/treat2_input2_called_sites.txt", "site/treat2_input2_all_sites.txt", "site/treat2_input2_background_sites.txt")
call_sites("align/treat3_input2_motif.bed", "site/treat3_input2_called_sites.txt", "site/treat3_input2_all_sites.txt", "site/treat3_input2_background_sites.txt")

site1 <- read.table("site/treat1_input1_called_sites.txt", header=T, sep=" ", colClasses="character" ) 
site2 <- read.table("site/treat2_input2_called_sites.txt", header=T, sep=" ", colClasses="character" ) 
site3 <- read.table("site/treat3_input2_called_sites.txt", header=T, sep=" ", colClasses="character" ) 

merge1 <- list(site1[,c("chr", "start", "end", "ref", "motif", "fp", "treat_conversion", "treat_level", "ctrl_conversion","p", "fdr", "TPM_input1","TPM_treat1")], site2[,c("chr", "start", "end", "ref", "motif", "fp", "treat_conversion", "treat_level","p", "fdr","TPM_treat2")], site3[,c("chr", "start", "end", "ref", "motif", "fp", "treat_conversion", "treat_level", "ctrl_conversion","p", "fdr", "TPM_input2","TPM_treat3")])
merge1 <- merge1 %>% reduce(full_join, by=c("chr", "start", "end", "ref", "motif")) %>% na.omit() ###1169
colnames(merge1) <- c("chr", "start", "end", "ref", "motif", "treat1_fp", "treat1_conversion", "treat1_level", "ctrl1_conversion","treat1_p", "treat1_fdr", "TPM_input1","TPM_treat1", "treat2_fp", "treat2_conversion", "treat2_level","treat2_p", "treat2_fdr", "TPM_treat2", "treat3_fp", "treat3_conversion", "treat3_level", "ctrl2_conversion","treat3_p", "treat3_fdr", "TPM_input2","TPM_treat3")
nrow(merge1) #1131
write.table(merge1, "site/snoRNA_tRNA_called_sites.txt", quote = FALSE, col.names = T, row.names = F)
#####fdr filtered
site <- read.table("site/snoRNA_tRNA_called_sites.txt", header=T, sep=" ",colClasses=c("ref"="character"))
fdr_site <- site[site$treat1_fdr<=0.001 &site$treat2_fdr<=0.001 &site$treat3_fdr<=0.001, ] ###1100
fdr_unpassed <- site[site$treat1_fdr>0.001 | site$treat2_fdr>0.001 | site$treat3_fdr>0.001, ] ##31

write.table(fdr_site,"site/snoRNA_tRNA_called_sites_fdrpassed.txt",sep='\t',quote=F,row.names=F)
write.table(fdr_unpassed,"site/snoRNA_tRNA_called_sites_fdr_unpassed.txt",sep='\t',quote=F,row.names=F)



merge2 <- merge1 %>% reduce(full_join, by=c("chr", "start", "end", "ref", "motif")) 
colnames(merge2) <- c("chr", "start", "end", "ref", "motif", "treat1_fp", "treat1_conversion", "treat1_level", "ctrl1_conversion","treat1_p", "treat1_fdr", "TPM_input1","TPM_treat1", "treat2_fp", "treat2_conversion", "treat2_level","treat2_p", "treat2_fdr", "TPM_treat2", "treat3_fp", "treat3_conversion", "treat3_level", "ctrl2_conversion","treat3_p", "treat3_fdr", "TPM_input2","TPM_treat3")
nrow(merge2) ## 1353
write.table(merge2, "site/snoRNA_tRNA_called_sites_na.txt", quote = FALSE, col.names = T, row.names = F)


##############background#####################################
input1 <- read.table("site/treat1_input1_background_sites.txt", header=T, sep=" ",colClasses="character")
input2 <- read.table("site/treat2_input2_background_sites.txt", header=T, sep=" ",colClasses="character")
input3 <- read.table("site/treat3_input2_background_sites.txt", header=T, sep=" ",colClasses="character")

merge1 <- list(input1[,c("chr", "start", "end", "ref", "motif", "fp", "treat_conversion", "treat_level", "ctrl_conversion","p", "fdr", "TPM_input1","TPM_treat1")], input2[,c("chr", "start", "end", "ref", "motif", "fp", "treat_conversion", "treat_level","p", "fdr","TPM_treat2")], input3[,c("chr", "start", "end", "ref", "motif", "fp", "treat_conversion", "treat_level", "ctrl_conversion","p", "fdr", "TPM_input2","TPM_treat3")])
merge1 <- merge1 %>% reduce(full_join, by=c("chr", "start", "end", "ref", "motif"))
colnames(merge1) <- c("chr", "start", "end", "ref", "motif", "treat1_fp", "treat1_conversion", "treat1_level", "ctrl1_conversion","treat1_p", "treat1_fdr", "TPM_input1","TPM_treat1", "treat2_fp", "treat2_conversion", "treat2_level","treat2_p", "treat2_fdr", "TPM_treat2", "treat3_fp", "treat3_conversion", "treat3_level", "ctrl2_conversion","treat3_p", "treat3_fdr", "TPM_input2","TPM_treat3")
nrow(merge1) #219
write.table(merge1, "site/snoRNA_tRNA_background_sites.txt", quote = FALSE, col.names = T, row.names = F)


######all sites
input1 <- read.table("site/treat1_input1_all_sites.txt", header=T, sep=" ",colClasses="character")
input2 <- read.table("site/treat2_input2_all_sites.txt", header=T, sep=" ",colClasses="character")
input3 <- read.table("site/treat3_input2_all_sites.txt", header=T, sep=" ",colClasses="character")

merge1 <- list(input1[,c("chr", "start", "end", "ref", "motif", "fp", "treat_conversion", "treat_level", "ctrl_conversion","p", "fdr", "TPM_input1","TPM_treat1")], input2[,c("chr", "start", "end", "ref", "motif", "fp", "treat_conversion", "treat_level","p", "fdr","TPM_treat2")], input3[,c("chr", "start", "end", "ref", "motif", "fp", "treat_conversion", "treat_level", "ctrl_conversion","p", "fdr", "TPM_input2","TPM_treat3")])
merge1 <- merge1 %>% reduce(full_join, by=c("chr", "start", "end", "ref", "motif")) %>% na.omit()
colnames(merge1) <- c("chr", "start", "end", "ref", "motif", "treat1_fp", "treat1_conversion", "treat1_level", "ctrl1_conversion","treat1_p", "treat1_fdr", "TPM_input1","TPM_treat1", "treat2_fp", "treat2_conversion", "treat2_level","treat2_p", "treat2_fdr", "TPM_treat2", "treat3_fp", "treat3_conversion", "treat3_level", "ctrl2_conversion","treat3_p", "treat3_fdr", "TPM_input2","TPM_treat3")
nrow(merge1) #15233
write.table(merge1, "site/snoRNA_tRNA_all_sites.txt", quote = FALSE, col.names = T, row.names = F)


