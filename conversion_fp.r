##module load R/4.0.3-foss-2020b
library(tidyverse)

options <- commandArgs(trailingOnly = TRUE)
ctrl_nnunn1 <- options[1]
treat_nnunn1 <- options[2]

ctrl_nnunn2 <- options[3]
treat_nnunn2 <- options[4]

out_file <- options[5]

ctrl1 <- read.table(ctrl_nnunn1, header=T, sep=" ")
treat1 <- read.table(treat_nnunn1, header=T, sep=" ")

mergemdf1 <- list(ctrl1[,c("first", "last", "ratio")], treat1[,c("first", "last", "ratio")])
merge1 <- mergemdf1 %>% reduce(full_join, by=c("first", "last"))
colnames(merge1) <- c("first", "last", "ctrl1_ratio", "treat1_ratio")
merge1$motif <- paste(merge1$first,merge1$last, sep="T")

####################
ctrl2 <- read.table(ctrl_nnunn2, header=T, sep=" ")
treat2 <- read.table(treat_nnunn2, header=T, sep=" ")

mergemdf2 <- list(ctrl2[,c("first", "last", "ratio")], treat2[,c("first", "last", "ratio")])
merge2 <- mergemdf2 %>% reduce(full_join, by=c("first", "last"))
colnames(merge2) <- c("first", "last","ctrl2_ratio", "treat2_ratio")
merge2$motif <- paste(merge2$first,merge2$last, sep="T")


spikein <- merge(merge1[,c("motif", "treat1_ratio")], merge2[,c("motif", "treat2_ratio")], by=c("motif"))
colnames(spikein) <- c("motif", "conversion", "fp")

#######read fp from nnunn1
nnunn1_fp <- read.table("fp.clean.nnunn.sort_NNUNN1.contex.table", header=T, sep=" ")
nnunn1_fp$motif <- paste(nnunn1_fp$first,nnunn1_fp$last, sep="T")

spikein <- merge(spikein, nnunn1_fp[,c("motif", "ratio")], by=c("motif"))
colnames(spikein) <- c("motif", "conversion","fp", "fp1")
spikein$fp <- pmax(spikein$fp, spikein$fp1)
spikein <- spikein %>% dplyr::select(motif,conversion,fp)
spikein <- rbind(c(".", mean(spikein$conversion),mean(spikein$fp)), spikein)

write.table(spikein, out_file, quote = FALSE, col.names = T, row.names = F)