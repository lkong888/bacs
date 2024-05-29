##module load R/4.0.3-foss-2020b
library(dplyr)
library(reshape2)
library(tidyr)
library(ggplot2)
library(cowplot)
library(tidyverse)
library(parallel)
library(Biostrings)

options <- commandArgs(trailingOnly = TRUE)
input <- options[1]
output <- options[2]
ivt <- options[3]
spikein <- options[4]

spikein <- read.table(spikein, header=T, sep=" ")
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
colnames(dat) <- c("chr", "start", "end", "strand", "ref", "treat_depth", "treat_T", "treat_C", "treat_conversion", "treat_gap", "treat_gap_rate", "motif", "control_ref","ctrl_depth", "ctrl_T", "ctrl_C", "ctrl_conversion", "ctrl_gap", "ctrl_gap_rate")
if(nrow(dat[grepl("N",dat$motif),])>0){
    for (i in 1:nrow(dat)){
      if(grepl("N",dat[i,]$motif)){
        dat[i,]$motif <- "."
      }
    }
}
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
dat$treat_level1 <- (dat$treat_conversion-as.numeric(dat$fp))/(as.numeric(dat$conversion)-as.numeric(dat$fp))

###read ivt
ivt <- read.table(ivt, header=F, sep="\t")
colnames(ivt) <- c("chr", "start", "end", "strand", "ref", "ivt_depth", "ivt_T", "ivt_C", "gap")
ivt <- ivt[(ivt$ivt_T+ivt$ivt_C)>=20,]
ivt$ivt_fp <- ivt$ivt_C/(ivt$ivt_T+ivt$ivt_C)
dat <- merge(dat, ivt[,c("chr", "start", "end", "ivt_T","ivt_C","ivt_fp")])
dat$treat_level2 <- (dat$treat_conversion-as.numeric(dat$ivt_fp))/(as.numeric(dat$conversion)-as.numeric(dat$ivt_fp))

##############filter by depth 
dat <- dat[(dat$treat_T+dat$treat_C)>=20 & (dat$ctrl_T+dat$ctrl_C)>=20,]
print("start test")

for (i in 1:nrow(dat)){
   fp <- dat[i,] %>% dplyr::select(ivt_fp) %>% as.numeric()
   dat$p[i] <- binom_t(depth = (dat$treat_T[i]+dat$treat_C[i]), mod = dat$treat_C[i], p = fp)
   data <- matrix(c(dat$treat_T[i], dat$treat_C[i], dat$ivt_T[i], dat$ivt_C[i]), nrow = 2)
   test <- chisq.test(data)
   dat$p_con[i] <- test$p.value
}
print("finished test")








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
