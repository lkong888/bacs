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

#####################read conversion and false positives 
"""
h1 <- read.table("mergebam/H1.clean.nnunn.sort_NNUNN1.contex.table", header=T, sep=" ")
h2 <- read.table("mergebam/H2.clean.nnunn.sort_NNUNN1.contex.table", header=T, sep=" ")
h3 <- read.table("mergebam/H3.clean.nnunn.sort_NNUNN1.contex.table", header=T, sep=" ")
h4 <- read.table("mergebam/H4.clean.nnunn.sort_NNUNN1.contex.table", header=T, sep=" ")

mergemdf <- list(h1[,c("first", "last", "ratio")], h2[,c("first", "last", "ratio")], h3[,c("first", "last", "ratio")], h4[,c("first", "last", "ratio")])
merge1 <- mergemdf %>% reduce(full_join, by=c("first", "last"))
colnames(merge1) <- c("first", "last", "h1_ratio", "h2_ratio", "h3_ratio", "h4_ratio")
merge1$treat1_ratio <- (merge1$h1_ratio+merge1$h2_ratio)/2
merge1$treat2_ratio <- (merge1$h3_ratio+merge1$h4_ratio)/2
merge1$ratio <- rowSums(merge1[,3:6])/4
#write.table(merge1, "mergebam/modified_nnunn1.context.table", quote = FALSE, col.names = T, row.names = F)

####################
h1 <- read.table("mergebam/H1.clean.nnunn.sort_NNUNN2.contex.table", header=T, sep=" ")
h2 <- read.table("mergebam/H2.clean.nnunn.sort_NNUNN2.contex.table", header=T, sep=" ")
h3 <- read.table("mergebam/H3.clean.nnunn.sort_NNUNN2.contex.table", header=T, sep=" ")
h4 <- read.table("mergebam/H4.clean.nnunn.sort_NNUNN2.contex.table", header=T, sep=" ")

mergemdf2 <- list(h1[,c("first", "last", "ratio")], h2[,c("first", "last", "ratio")], h3[,c("first", "last", "ratio")], h4[,c("first", "last", "ratio")])
merge2 <- mergemdf2 %>% reduce(full_join, by=c("first", "last"))
colnames(merge2) <- c("first", "last", "h1_ratio", "h2_ratio", "h3_ratio", "h4_ratio")
merge2$treat1_ratio <- (merge2$h1_ratio+merge2$h2_ratio)/2
merge2$treat2_ratio <- (merge2$h3_ratio+merge2$h4_ratio)/2
merge2$ratio <- rowSums(merge2[,3:6])/4
#write.table(merge2, "mergebam/unmodified_nnunn2.context.table", quote = FALSE, col.names = T, row.names = F)
"""
###read false positives
call_sites <- function(input){
dat <- read.table(input, header=F, sep="\t", colClasses=c("V5"="character"))
###read modified and unmotified NNUNN results from figure.R
nnunn1 <- read.table("mergebam/modified_nnunn1.context.table", header=T, sep=" ")
nnunn1$motif <- paste(nnunn1$first,nnunn1$last, sep="T")

nnunn2 <- read.table("mergebam/unmodified_nnunn2.context.table", header=T, sep=" ")
nnunn2$motif <- paste(nnunn2$first,nnunn2$last, sep="T")

if(grepl("treat1",input)){
    spikein <- merge(nnunn1[,c("motif", "treat1_ratio")], nnunn2[,c("motif", "treat1_ratio")], by=c("motif"))
    colnames(spikein) <- c("motif", "conversion", "fp")
    spikein <- rbind(c(".", mean(spikein$conversion),mean(spikein$fp)), spikein)
} else if(grepl("treat2", input)){
    spikein <- merge(nnunn1[,c("motif", "treat2_ratio")], nnunn2[,c("motif", "treat2_ratio")], by=c("motif"))
    colnames(spikein) <- c("motif", "conversion", "fp")
    spikein <- rbind(c(".", mean(spikein$conversion),mean(spikein$fp)), spikein)
} else {print ("smaple infor not found")}


####append actual modification level and fp in given motif
colnames(dat) <- c("chr", "start", "end", "strand", "ref", "treat_depth", "treat_T", "treat_C","motif")
dat$treat_conversion <- dat$treat_C/(dat$treat_T+dat$treat_C)
dat <- merge(dat, spikein, by=c("motif"),all.x=TRUE)
dat$treat_level <- (dat$treat_conversion-as.numeric(dat$fp))/(as.numeric(dat$conversion)-as.numeric(dat$fp))
##############filter by depth 
dat <- dat[(dat$treat_T+dat$treat_C)>=20,]

####transform motif for ref=A
reverseComplement <- function(sequence) {
  complement <- gsub("A", "t", sequence)
  complement <- gsub("T", "a", complement)
  complement <- gsub("G", "c", complement)
  complement <- gsub("C", "g", complement)
  return(toupper(rev(complement)))
}
dat[dat$ref=="A",]$motif <- reverseComplement(dat[dat$ref=="A",]$motif)

###binomial test
for (i in 1:nrow(dat)){
   fp <- dat[i,] %>% dplyr::select(fp) %>% as.numeric()
   dat$p[i] <- binom_t(depth = (dat$treat_T[i]+dat$treat_C[i]), mod = dat$treat_C[i], p = fp)
}

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

call_sites("splitbam/treat1.sta.txt", "site/treat1_all_sites.txt")
call_sites("splitbam/treat1.sta.txt", "site/treat2_all_sites.txt")





#####################read conversion and false positives 

h1 <- read.table("align/H1_2023Jun19_S1.clean.nnunn.sort_NNUNN1.contex.table", header=T, sep=" ")
h2 <- read.table("align/H2_2023Jun19_S2.clean.nnunn.sort_NNUNN1.contex.table", header=T, sep=" ")
h3 <- read.table("align/H3_2023Jun19_S3.clean.nnunn.sort_NNUNN1.contex.table", header=T, sep=" ")
h4 <- read.table("align/H4_2023Jun19_S4.clean.nnunn.sort_NNUNN1.contex.table", header=T, sep=" ")

mergemdf <- list(h1[,c("first", "last", "ratio")], h2[,c("first", "last", "ratio")], h3[,c("first", "last", "ratio")], h4[,c("first", "last", "ratio")])
merge1 <- mergemdf %>% reduce(full_join, by=c("first", "last"))
colnames(merge1) <- c("first", "last", "h1_ratio", "h2_ratio", "h3_ratio", "h4_ratio")
merge1$treat1_ratio <- (merge1$h1_ratio+merge1$h2_ratio)/2
merge1$treat2_ratio <- (merge1$h3_ratio+merge1$h4_ratio)/2
merge1$ratio <- rowSums(merge1[,3:6])/4
write.table(merge1, "align/modified_nnunn1.context.table", quote = FALSE, col.names = T, row.names = F)

####################
h1 <- read.table("align/H1_2023Jun19_S1.clean.nnunn.sort_NNUNN2.contex.table", header=T, sep=" ")
h2 <- read.table("align/H2_2023Jun19_S2.clean.nnunn.sort_NNUNN2.contex.table", header=T, sep=" ")
h3 <- read.table("align/H3_2023Jun19_S3.clean.nnunn.sort_NNUNN2.contex.table", header=T, sep=" ")
h4 <- read.table("align/H4_2023Jun19_S4.clean.nnunn.sort_NNUNN2.contex.table", header=T, sep=" ")

mergemdf2 <- list(h1[,c("first", "last", "ratio")], h2[,c("first", "last", "ratio")], h3[,c("first", "last", "ratio")], h4[,c("first", "last", "ratio")])
merge2 <- mergemdf2 %>% reduce(full_join, by=c("first", "last"))
colnames(merge2) <- c("first", "last", "h1_ratio", "h2_ratio", "h3_ratio", "h4_ratio")
merge2$treat1_ratio <- (merge2$h1_ratio+merge2$h2_ratio)/2
merge2$treat2_ratio <- (merge2$h3_ratio+merge2$h4_ratio)/2
merge2$ratio <- rowSums(merge2[,3:6])/4
write.table(merge2, "align/unmodified_nnunn2.context.table", quote = FALSE, col.names = T, row.names = F)
