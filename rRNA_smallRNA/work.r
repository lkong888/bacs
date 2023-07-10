##module load R/4.0.3-foss-2020b
library(dplyr)
library(reshape2)
library(tidyr)
library(ggplot2)
library(ggpubr)
setwd("/users/ludwig/ebu571/ebu571/20May2023_final")
############################70-mer spikein##################################
data_summary <- function(data, varname, groupnames){
  require(plyr)
  summary_func <- function(x, col){
    c(mean = mean(x[[col]], na.rm=TRUE),
      sd = sd(x[[col]], na.rm=TRUE))
  }
  data_sum<-ddply(data, groupnames, .fun=summary_func,
                  varname)
  data_sum <- rename(data_sum, c("mean" = varname))
 return(data_sum)
}


####compare modification before and after deduplication sample: H7 and H8
nondedup <- read.table("align/H7_2023May19_S7.rRNA.filter.sort_pseusite.mpile.txt", header=F, sep="\t")
colnames(nondedup) <- c("chr", "pos", "ref", "depth", "T", "C", "gap")
nondedup$ratio <- nondedup$C/(nondedup$C+nondedup$T)
dedup <- read.table("align/H7_2023May19_S7.rRNA.snoRNA.tRNA_merge.filter.sort_pseusite.mpile.txt", header=F, sep="\t")
colnames(dedup) <- c("chr", "pos", "ref", "depth", "T", "C", "gap")
dedup$ratio <- dedup$C/(dedup$C+dedup$T)

######check 18S
nondedup_18s <- nondedup[nondedup$chr=="NR_003286.4_RNA18SN5", ]
dedup_18s <- dedup[dedup$chr=="NR_003286.4_RNA18SN5", ]
merge_18s <- merge(nondedup_18s, dedup_18s, by=c("chr", "pos"))
seldat_18s <- merge_18s[,c("chr", "pos", "ref.x", "depth.x", "ratio.x", "depth.y", "ratio.y")]
colnames(seldat_18s) <- c("chr", "pos", "ref", "depth.nondedup", "ratio.nondedup", "depth.dedup", "ratio.dedup")
dat1 <- melt(seldat_18s[,c("pos", "ratio.nondedup", "ratio.dedup")], id.vars=c("pos"))
pdf("plot/H7_dedup_vs_nondedup_18Sratio.pdf", width=50,height = 10)
ggplot(dat1, aes(x=pos, y=round(value,4)*100, fill=variable))+ geom_bar(stat = "identity", position = position_dodge(), alpha = 0.75)
dev.off()

######check 28S
nondedup_28s <- nondedup[nondedup$chr=="NR_003287.4_RNA28SN5", ]
dedup_28s <- dedup[dedup$chr=="NR_003287.4_RNA28SN5", ]
merge_28s <- merge(nondedup_28s, dedup_28s, by=c("chr", "pos"))
seldat_28s <- merge_28s[,c("chr", "pos", "ref.x", "depth.x", "ratio.x", "depth.y", "ratio.y")]
colnames(seldat_28s) <- c("chr", "pos", "ref", "depth.nondedup", "ratio.nondedup", "depth.dedup", "ratio.dedup")
dat2 <- melt(seldat_28s[,c("pos", "ratio.nondedup", "ratio.dedup")], id.vars=c("pos"))
pdf("plot/H7_dedup_vs_nondedup_28Sratio.pdf", width=50,height = 10)
ggplot(dat2, aes(x=pos, y=round(value,4)*100, fill=variable))+ geom_bar(stat = "identity", position = position_dodge(), alpha = 0.75)
dev.off()

#####check 5.8S
nondedup_5.8s <- nondedup[nondedup$chr=="NR_003285.3_RNA5_8SN5", ]
dedup_5.8s <- dedup[dedup$chr=="NR_003285.3_RNA5_8SN5", ]
merge_5.8s <- merge(nondedup_5.8s, dedup_5.8s, by=c("chr", "pos"))
seldat_5.8s <- merge_5.8s[,c("chr", "pos", "ref.x", "depth.x", "ratio.x", "depth.y", "ratio.y")]
colnames(seldat_5.8s) <- c("chr", "pos", "ref", "depth.nondedup", "ratio.nondedup", "depth.dedup", "ratio.dedup")
dat3 <- melt(seldat_5.8s[,c("pos", "ratio.nondedup", "ratio.dedup")], id.vars=c("pos"))
pdf("plot/H7_dedup_vs_nondedup_5.8Sratio.pdf", width=50,height = 10)
ggplot(dat3, aes(x=pos, y=round(value,4)*100, fill=variable))+ geom_bar(stat = "identity", position = position_dodge(), alpha = 0.75)
dev.off()

########check depth for 18S 28S and 5.8S

nondedup_18s <- nondedup[nondedup$chr=="NR_003286.4_RNA18SN5", ]
dedup_18s <- dedup[dedup$chr=="NR_003286.4_RNA18SN5", ]
merge_18s <- merge(nondedup_18s, dedup_18s, by=c("chr", "pos"))
seldat_18s <- merge_18s[,c("chr", "pos", "ref.x", "depth.x", "ratio.x", "depth.y", "ratio.y")]
colnames(seldat_18s) <- c("chr", "pos", "ref", "depth.nondedup", "ratio.nondedup", "depth.dedup", "ratio.dedup")
dat1 <- melt(seldat_18s[,c("pos", "depth.nondedup", "depth.dedup")], id.vars=c("pos"))
pdf("plot/H7_dedup_vs_nondedup_18Sdepth.pdf", width=50,height = 10)
ggplot(dat1, aes(x=pos, y=round(value,4)*100, fill=variable))+ geom_bar(stat = "identity", position = position_dodge(), alpha = 0.75)
dev.off()

######check 28S
nondedup_28s <- nondedup[nondedup$chr=="NR_003287.4_RNA28SN5", ]
dedup_28s <- dedup[dedup$chr=="NR_003287.4_RNA28SN5", ]
merge_28s <- merge(nondedup_28s, dedup_28s, by=c("chr", "pos"))
seldat_28s <- merge_28s[,c("chr", "pos", "ref.x", "depth.x", "ratio.x", "depth.y", "ratio.y")]
colnames(seldat_28s) <- c("chr", "pos", "ref", "depth.nondedup", "ratio.nondedup", "depth.dedup", "ratio.dedup")
dat2 <- melt(seldat_28s[,c("pos", "depth.nondedup", "depth.dedup")], id.vars=c("pos"))
pdf("plot/H7_dedup_vs_nondedup_28Sdepth.pdf", width=50,height = 10)
ggplot(dat2, aes(x=pos, y=round(value,4)*100, fill=variable))+ geom_bar(stat = "identity", position = position_dodge(), alpha = 0.75)
dev.off()

#####check 5.8S
nondedup_5.8s <- nondedup[nondedup$chr=="NR_003285.3_RNA5_8SN5", ]
dedup_5.8s <- dedup[dedup$chr=="NR_003285.3_RNA5_8SN5", ]
merge_5.8s <- merge(nondedup_5.8s, dedup_5.8s, by=c("chr", "pos"))
seldat_5.8s <- merge_5.8s[,c("chr", "pos", "ref.x", "depth.x", "ratio.x", "depth.y", "ratio.y")]
colnames(seldat_5.8s) <- c("chr", "pos", "ref", "depth.nondedup", "ratio.nondedup", "depth.dedup", "ratio.dedup")
dat3 <- melt(seldat_5.8s[,c("pos", "depth.nondedup", "depth.dedup")], id.vars=c("pos"))
pdf("plot/H7_dedup_vs_nondedup_5.8Sdepth.pdf", width=50,height = 10)
ggplot(dat3, aes(x=pos, y=round(value,4)*100, fill=variable))+ geom_bar(stat = "identity", position = position_dodge(), alpha = 0.75)
dev.off()

##################################################################################################for H8
nondedup <- read.table("align/H8_2023May19_S8.rRNA.filter.sort_pseusite.mpile.txt", header=F, sep="\t")
colnames(nondedup) <- c("chr", "pos", "ref", "depth", "T", "C", "gap")
nondedup$ratio <- nondedup$C/(nondedup$C+nondedup$T)
dedup <- read.table("align/H8_2023May19_S8.rRNA.snoRNA.tRNA_merge.filter.sort_pseusite.mpile.txt", header=F, sep="\t")
colnames(dedup) <- c("chr", "pos", "ref", "depth", "T", "C", "gap")
dedup$ratio <- dedup$C/(dedup$C+dedup$T)

######check 18S
nondedup_18s <- nondedup[nondedup$chr=="NR_003286.4_RNA18SN5", ]
dedup_18s <- dedup[dedup$chr=="NR_003286.4_RNA18SN5", ]
merge_18s <- merge(nondedup_18s, dedup_18s, by=c("chr", "pos"))
seldat_18s <- merge_18s[,c("chr", "pos", "ref.x", "depth.x", "ratio.x", "depth.y", "ratio.y")]
colnames(seldat_18s) <- c("chr", "pos", "ref", "depth.nondedup", "ratio.nondedup", "depth.dedup", "ratio.dedup")
dat1 <- melt(seldat_18s[,c("pos", "ratio.nondedup", "ratio.dedup")], id.vars=c("pos"))
pdf("plot/H8_dedup_vs_nondedup_18Sratio.pdf", width=50,height = 10)
ggplot(dat1, aes(x=pos, y=round(value,4)*100, fill=variable))+ geom_bar(stat = "identity", position = position_dodge(), alpha = 0.75)
dev.off()

######check 28S
nondedup_28s <- nondedup[nondedup$chr=="NR_003287.4_RNA28SN5", ]
dedup_28s <- dedup[dedup$chr=="NR_003287.4_RNA28SN5", ]
merge_28s <- merge(nondedup_28s, dedup_28s, by=c("chr", "pos"))
seldat_28s <- merge_28s[,c("chr", "pos", "ref.x", "depth.x", "ratio.x", "depth.y", "ratio.y")]
colnames(seldat_28s) <- c("chr", "pos", "ref", "depth.nondedup", "ratio.nondedup", "depth.dedup", "ratio.dedup")
dat2 <- melt(seldat_28s[,c("pos", "ratio.nondedup", "ratio.dedup")], id.vars=c("pos"))
pdf("plot/H8_dedup_vs_nondedup_28Sratio.pdf", width=50,height = 10)
ggplot(dat2, aes(x=pos, y=round(value,4)*100, fill=variable))+ geom_bar(stat = "identity", position = position_dodge(), alpha = 0.75)
dev.off()

#####check 5.8S
nondedup_5.8s <- nondedup[nondedup$chr=="NR_003285.3_RNA5_8SN5", ]
dedup_5.8s <- dedup[dedup$chr=="NR_003285.3_RNA5_8SN5", ]
merge_5.8s <- merge(nondedup_5.8s, dedup_5.8s, by=c("chr", "pos"))
seldat_5.8s <- merge_5.8s[,c("chr", "pos", "ref.x", "depth.x", "ratio.x", "depth.y", "ratio.y")]
colnames(seldat_5.8s) <- c("chr", "pos", "ref", "depth.nondedup", "ratio.nondedup", "depth.dedup", "ratio.dedup")
dat3 <- melt(seldat_5.8s[,c("pos", "ratio.nondedup", "ratio.dedup")], id.vars=c("pos"))
pdf("plot/H8_dedup_vs_nondedup_5.8Sratio.pdf", width=50,height = 10)
ggplot(dat3, aes(x=pos, y=round(value,4)*100, fill=variable))+ geom_bar(stat = "identity", position = position_dodge(), alpha = 0.75)
dev.off()

########check depth for 18S 28S and 5.8S

nondedup_18s <- nondedup[nondedup$chr=="NR_003286.4_RNA18SN5", ]
dedup_18s <- dedup[dedup$chr=="NR_003286.4_RNA18SN5", ]
merge_18s <- merge(nondedup_18s, dedup_18s, by=c("chr", "pos"))
seldat_18s <- merge_18s[,c("chr", "pos", "ref.x", "depth.x", "ratio.x", "depth.y", "ratio.y")]
colnames(seldat_18s) <- c("chr", "pos", "ref", "depth.nondedup", "ratio.nondedup", "depth.dedup", "ratio.dedup")
dat1 <- melt(seldat_18s[,c("pos", "depth.nondedup", "depth.dedup")], id.vars=c("pos"))
pdf("plot/H8_dedup_vs_nondedup_18Sdepth.pdf", width=50,height = 10)
ggplot(dat1, aes(x=pos, y=round(value,4)*100, fill=variable))+ geom_bar(stat = "identity", position = position_dodge(), alpha = 0.75)
dev.off()

######check 28S
nondedup_28s <- nondedup[nondedup$chr=="NR_003287.4_RNA28SN5", ]
dedup_28s <- dedup[dedup$chr=="NR_003287.4_RNA28SN5", ]
merge_28s <- merge(nondedup_28s, dedup_28s, by=c("chr", "pos"))
seldat_28s <- merge_28s[,c("chr", "pos", "ref.x", "depth.x", "ratio.x", "depth.y", "ratio.y")]
colnames(seldat_28s) <- c("chr", "pos", "ref", "depth.nondedup", "ratio.nondedup", "depth.dedup", "ratio.dedup")
dat2 <- melt(seldat_28s[,c("pos", "depth.nondedup", "depth.dedup")], id.vars=c("pos"))
pdf("plot/H8_dedup_vs_nondedup_28Sdepth.pdf", width=50,height = 10)
ggplot(dat2, aes(x=pos, y=round(value,4)*100, fill=variable))+ geom_bar(stat = "identity", position = position_dodge(), alpha = 0.75)
dev.off()

#####check 5.8S
nondedup_5.8s <- nondedup[nondedup$chr=="NR_003285.3_RNA5_8SN5", ]
dedup_5.8s <- dedup[dedup$chr=="NR_003285.3_RNA5_8SN5", ]
merge_5.8s <- merge(nondedup_5.8s, dedup_5.8s, by=c("chr", "pos"))
seldat_5.8s <- merge_5.8s[,c("chr", "pos", "ref.x", "depth.x", "ratio.x", "depth.y", "ratio.y")]
colnames(seldat_5.8s) <- c("chr", "pos", "ref", "depth.nondedup", "ratio.nondedup", "depth.dedup", "ratio.dedup")
dat3 <- melt(seldat_5.8s[,c("pos", "depth.nondedup", "depth.dedup")], id.vars=c("pos"))
pdf("plot/H8_dedup_vs_nondedup_5.8Sdepth.pdf", width=50,height = 10)
ggplot(dat3, aes(x=pos, y=round(value,4)*100, fill=variable))+ geom_bar(stat = "identity", position = position_dodge(), alpha = 0.75)
dev.off()







##############################plot TPM of two input
tpm1 <- read.table("featureCounts/input1.tpm", header=F, sep="\t")
