##module load R/4.0.3-foss-2020b
library(dplyr)
library(reshape2)
library(tidyr)
library(ggplot2)
library(cowplot)
library(tidyverse)
library(parallel)
library(Biostrings)
setwd("/users/ludwig/ebu571/ebu571/23June2023")

options <- commandArgs(trailingOnly = TRUE)
input <- options[1]
output <- options[2]

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
############################################################### mutation calling on mRNA #########################################################
##Criteria
#####
## TPM >1 in control
##1. coverage over 20 in both treated and control 
##2. conversion rate <1% in input(background list) and the number of T-C mutation <=2
##3. Actual conversion rate above 5% (actual conversion rate=(y-fp)/(conversion-fp))
##4. Overlap in all replicates

##Sample infor: treat1 vs input1; treat2 vs input2
##input colnumn names: chr, start, end, ref, depth, T, C, mod_level, gap, gap_ratio, motif, ctrl_depth, ctrl_T, ctrl_C, ctrl_ratio, ctrl_gap, ctrl_gap_ratio
##corrected curve: y=(conversion-fp)/(1-0)*x+fp  x=(y-fp)/(conversion-fp)
##corrected curve for each motif

###read false positives
call_sites <- function(input, output){
dat <- read.table(input, header=F, sep="\t", colClasses=c("V5"="character"))
###read modified and unmotified NNUNN results from figure.R
nnunn1 <- read.table("mergebam/modified_nnunn1.context.table", header=T, sep=" ")
nnunn1$motif <- paste(nnunn1$first,nnunn1$last, sep="T")

nnunn2 <- read.table("mergebam/unmodified_nnunn2.context.table", header=T, sep=" ")
nnunn2$motif <- paste(nnunn2$first,nnunn2$last, sep="T")

if(grepl("treat1",input)){
    spikein <- merge(nnunn1[,c("motif", "treat1_ratio")], nnunn2[,c("motif", "treat1_ratio")], by=c("motif"))
    colnames(spikein) <- c("motif", "conversion", "fp")
    meanconversion <- mean(spikein$conversion)
    meanfp <- mean(spikein$fp)
    spikein <- rbind(c(".", meanconversion,meanfp), spikein)
    spikein <- rbind(c("NTTCA", meanconversion,meanfp), spikein)
} else if(grepl("treat2", input)){
    spikein <- merge(nnunn1[,c("motif", "treat2_ratio")], nnunn2[,c("motif", "treat2_ratio")], by=c("motif"))
    colnames(spikein) <- c("motif", "conversion", "fp")
    meanconversion <- mean(spikein$conversion)
    meanfp <- mean(spikein$fp)
    spikein <- rbind(c(".", meanconversion,meanfp), spikein)
    spikein <- rbind(c("NTTCA", meanconversion,meanfp), spikein)
} else {print ("smaple infor not found")}


####append actual modification level and fp in given motif
colnames(dat) <- c("chr", "start", "end", "strand", "ref", "treat_depth", "treat_T", "treat_C", "treat_conversion", "treat_gap", "treat_gap_rate", "motif", "control_ref","ctrl_depth", "ctrl_T", "ctrl_C", "ctrl_conversion", "ctrl_gap", "ctrl_gap_rate")

####transform motif for ref=A
reverseComplement <- function(sequence) {
  complement <- gsub("A", "t", sequence)
  complement <- gsub("T", "a", complement)
  complement <- gsub("G", "c", complement)
  complement <- gsub("C", "g", complement)
  complement <- paste(rev(strsplit(complement, "")[[1]]), collapse = "")
  return(toupper(complement))
}
datA <-dat[dat$ref=="A",]
for (i in 1:nrow(datA)){
  datA$motif[i] <- reverseComplement(datA$motif[i])
}
dat <- rbind(datA,dat[dat$ref=="T",])
#dat[dat$ref=="A",]$motif <- reverseComplement(dat[dat$ref=="A",]$motif)

dat <- merge(dat, spikein, by=c("motif"),all.x=TRUE)
dat$treat_level <- (dat$treat_conversion-as.numeric(dat$fp))/(as.numeric(dat$conversion)-as.numeric(dat$fp))
##############filter by depth 
dat <- dat[(dat$treat_T+dat$treat_C)>=20 & (dat$ctrl_T+dat$ctrl_C)>=20,]
print("start test")
###binomial test
for (i in 1:nrow(dat)){
   fp <- dat[i,] %>% dplyr::select(fp) %>% as.numeric()
   dat$p[i] <- binom_t(depth = (dat$treat_T[i]+dat$treat_C[i]), mod = dat$treat_C[i], p = fp)
}
print("finished test")
dat$fdr <- p.adjust(dat$p, method = "BH")
##called sites
#site <- dat[(dat$treat_T+dat$treat_C)>=20 & (dat$ctrl_T+dat$ctrl_C)>=20 & dat$TPM >1 & as.numeric(dat$ctrl_conversion)<0.01 & as.numeric(dat$treat_level)>=0.05,]
#site <- dat[(dat$treat_T+dat$treat_C)>=20 & (dat$ctrl_T+dat$ctrl_C)>=20 & as.numeric(dat$treat_level)>=0.05,] %>% na.omit()
#site1 <- site[as.numeric(site$ctrl_conversion)<0.01 | as.numeric(site$ctrl_C)<2,] %>% na.omit()
#background <- dat[as.numeric(dat$ctrl_conversion)>=0.01 & as.numeric(dat$ctrl_C)>=2,]

#write.table(site1, output1, quote = FALSE, col.names = T, row.names = F)
#write.table(dat[(dat$treat_T+dat$treat_C)>=20 & (dat$ctrl_T+dat$ctrl_C)>=20,], output2, quote = FALSE, col.names = T, row.names = F)
write.table(dat, output, quote = FALSE, col.names = T, row.names = F)
}

#call_sites("mergebam/treat1_input1_allT_motif.bed", "site/treat1_input1_sites.txt")
#call_sites("mergebam/treat2_input2_allT_motif.bed", "site/treat2_input2_sites.txt")
call_sites(input, output)