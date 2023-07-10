##module load R/4.0.3-foss-2020b
####this script is for figures on spikein and internal control
library(dplyr)
library(reshape2)
library(tidyr)
library(ggplot2)
library(ggpubr)
library(cowplot)
library(RColorBrewer)
library(tidyverse)
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


dat <- read.table("align/spike_in.txt", header=T, sep="\t")
dat$sample <- gsub(".rRNA.filter.dedup.mpile.txt","",dat$smp)
dat$smp <- gsub("align/","",dat$sample)
dat <- dat[!grepl("H10_2023May19_S10", dat$smp) & !grepl("H9_2023May19_S9",dat$smp) & !grepl("H16_2023Apr16_S15",dat$smp) & !grepl("H26_2023Apr16_S22",dat$smp),]
meltdat <- melt(dat[,c("smp", "UtoC", "UtoR")], id.vars=c("smp"))
df2 <- data_summary(meltdat, varname="value", 
                    groupnames=c("variable"))
##pdf("plot/70mer_spikein.pdf", width=2,height = 3)
pdf("plot/70mer_spikein.pdf")
ggplot(df2, aes(x=variable, y=round(value*100,2), fill=variable))+
geom_bar(stat = "identity", width=0.55)+  ####identity equals to sum; needs to make sure is the mean value; otherwise, use geom_bar(stat = "summary", fun = "mean", width=0.55)
ylab("mutation ratio(%)")+xlab("")+
ylim(0,100)+
geom_errorbar(aes(ymin=(value-sd)*100, ymax=(value+sd)*100), width=.04,
                 position=position_dodge(.9)) +
geom_text(data = df2, aes(label=round(value*100, 2)), vjust=-1.5)+  theme_light() +
  theme(axis.text = element_text(color = "black")) +theme(legend.position="none")
dev.off()


#######plot conversion and false positives on NNUNN spikeins using H7 and H8
###CONVERSION
nnunn1 <- read.table("align/nnunn1.txt", header=T, sep="\t")
nnunn1$sample <- gsub(".clean.nnunn.sort_NNUNN1.contex.table","",nnunn1$smp)
nnunn1$smp <- gsub("align/","",nnunn1$sample)
nnunn1 <- nnunn1[grepl("H3_2023Apr16_S3", nnunn1$smp) | grepl("H7_2023May19_S7", nnunn1$smp) | grepl("H8_2023May19_S8", nnunn1$smp),] %>% dplyr::select(smp, NNUNN1_ratio)
colnames(nnunn1) <- c("smp", "ratio")

meltdf1 <- nnunn1 %>% melt(id.vars=c("smp"))
meandf1 <- data_summary(meltdf1, varname="value", 
                    groupnames=c("variable")) ###mean and sd

p1 <- ggplot(meandf1, aes(x=variable, y=round(value,2)*100, fill=variable))+
geom_bar(stat = "identity", width=0.55)+  ####identity equals to sum; needs to make sure is the mean value; otherwise, use geom_bar(stat = "summary", fun = "mean", width=0.55)
ylim(0,100)+
geom_errorbar(aes(ymin=(value-sd)*100, ymax=(value+sd)*100), width=.05,
                 position=position_dodge(.9)) +
geom_text(data = meandf1, aes(label=round(value, 2)*100), vjust=-1.5)+  theme_light() + ylab("conversion(%)")+xlab("modified NNUNN")+
  theme(axis.text = element_text(color = "black")) +theme(legend.position="none") 


###FALSE POSITIVES ON unmodified 2kb
dat <- read.table("align/spike_in.txt", header=T, sep="\t")
dat$sample <- gsub(".rRNA.filter.dedup.mpile.txt","",dat$smp)
dat$smp <- gsub("align/","",dat$sample)
dat <- dat[!grepl("H10_2023May19_S10", dat$smp) & !grepl("H9_2023May19_S9",dat$smp) & !grepl("H16_2023Apr16_S15",dat$smp) & !grepl("H26_2023Apr16_S22",dat$smp),]
meltdat <- melt(dat[,c("smp", "fp_ratio")], id.vars=c("smp"))
meandf2 <- data_summary(meltdat, varname="value", 
                    groupnames=c("variable"))

p2 <- ggplot(meandf2, aes(x=variable, y=round(value,4)*100, fill=variable))+
geom_bar(stat = "identity", width=0.55)+  ####identity equals to sum; needs to make sure is the mean value; otherwise, use geom_bar(stat = "summary", fun = "mean", width=0.55)
ylab("conversion(%)")+xlab("unmodified 2kb")+
ylim(0,5)+
geom_errorbar(aes(ymin=(value-sd)*100, ymax=(value+sd)*100), width=.05,
                 position=position_dodge(.9)) +
geom_text(data = meandf2, aes(label=round(value, 4)*100), vjust=-1.5)+  theme_light() +
  theme(axis.text = element_text(color = "black")) +theme(legend.position="none")

p <- plot_grid(p1, p2, nrow=1)
ggsave("plot/nnunn1_2kb.pdf",p,width=4,height = 3)


######conversion rate and false positives on different motis
###write intergrated output for nnunn1
h1 <- read.table("align/H1_2023May19_S1.clean.nnunn.sort_NNUNN1.contex.table", header=T, sep=" ")
h2 <- read.table("align/H2_2023May19_S2.clean.nnunn.sort_NNUNN1.contex.table", header=T, sep=" ")
h3 <- read.table("align/H3_2023May19_S3.clean.nnunn.sort_NNUNN1.contex.table", header=T, sep=" ")
h4 <- read.table("align/H4_2023May19_S4.clean.nnunn.sort_NNUNN1.contex.table", header=T, sep=" ")
h5 <- read.table("align/H5_2023May19_S5.clean.nnunn.sort_NNUNN1.contex.table", header=T, sep=" ")
h6 <- read.table("align/H6_2023May19_S6.clean.nnunn.sort_NNUNN1.contex.table", header=T, sep=" ")
h7 <- read.table("align/H7_2023May19_S7.clean.nnunn.sort_NNUNN1.contex.table", header=T, sep=" ")
h8 <- read.table("align/H8_2023May19_S8.clean.nnunn.sort_NNUNN1.contex.table", header=T, sep=" ")
merge1 <- list(h1[,c("first", "last", "ratio")], h2[,c("first", "last", "ratio")],h3[,c("first", "last", "ratio")],h4[,c("first", "last", "ratio")],h5[,c("first", "last", "ratio")],h6[,c("first", "last", "ratio")],h7[,c("first", "last", "ratio")],h8[,c("first", "last", "ratio")])
merge1 <- merge1 %>% reduce(full_join, by=c("first", "last"))
colnames(merge1) <- c("first", "last", "h1_ratio", "h2_ratio", "h3_ratio", "h4_ratio", "h5_ratio", "h6_ratio", "h7_ratio", "h8_ratio")
write.table(merge1, "align/modified_nnunn.context_all.table", quote = FALSE, col.names = T, row.names = F)


####plot conversion rate for NNUNN1
h7 <- read.table("align/H7_2023May19_S7.clean.nnunn.sort_NNUNN1.contex.table", header=T, sep=" ")
h8 <- read.table("align/H8_2023May19_S8.clean.nnunn.sort_NNUNN1.contex.table", header=T, sep=" ")
h3 <- read.table("align/H3_2023Apr16_S3.clean.nnunn.sort_NNUNN1.contex.table", header=T, sep=" ")

mergemdf <- merge(h7[,c("first", "last", "ratio")], h8[,c("first", "last", "ratio")], by=c("first", "last"))
colnames(mergemdf) <- c("first", "last", "h7_ratio", "h8_ratio")
mergemdf <- merge(mergemdf, h3[,c("first", "last", "ratio")], by=c("first", "last"))
colnames(mergemdf) <- c("first", "last", "h7_ratio", "h8_ratio", "h3_ratio")

mergemdf$ratio <- round((mergemdf$h7_ratio+mergemdf$h8_ratio+mergemdf$h3_ratio)/3,2)*100
mdf1 <- spread(mergemdf[,c("first", "last", "ratio")], key = last, value = ratio) %>%  as.data.frame()
final_matrix <- data.frame(mdf1[,-1], row.names=as.character(mdf1[,1])) %>% as.matrix(rownames = TRUE)

mergemdf$ratio <- (mergemdf$h7_ratio+mergemdf$h8_ratio+mergemdf$h3_ratio)/3
write.table(mergemdf, "align/modified_nnunn.context.table", quote = FALSE, col.names = T, row.names = F)

breaksList = seq(0, 100, by = 1)
library(pheatmap)
library(RColorBrewer)

pdf("plot/modified_nnunn.ratio.pdf")
pheatmap(final_matrix,cluster_cols = FALSE, cluster_rows = FALSE, display_numbers = TRUE, number_format="%.0i", fontsize_number = 12, color=colorRampPalette(c("white", "red"))(length(breaksList)), breaks = breaksList, legend_breaks = seq(0, 100, 50), legend_labels = c("0", "50", "100"),border_color = "grey70", number_color="grey0")
dev.off()

#########use ggplot heatmap########################################
pdf("plot/modified_nnunn.ratio.pdf")
ggplot(mergemdf[,c("first", "last", "ratio")], aes(x = first, y = last, fill = ratio)) + geom_tile(color = "gray", size=0.1)+ coord_fixed()+ ##to make Square tiles
  geom_text(aes(label = ratio), color = "black", size = 6) +
  scale_fill_gradient2(low = "white", high = "red", limits=c(0,100))
dev.off()




###for unmodified nnunn1 h1-h8(nnunn2)
h1 <- read.table("align/H1_2023May19_S1.clean.nnunn.sort_NNUNN2.contex.table", header=T, sep=" ")
h2 <- read.table("align/H2_2023May19_S2.clean.nnunn.sort_NNUNN2.contex.table", header=T, sep=" ")
h3 <- read.table("align/H3_2023May19_S3.clean.nnunn.sort_NNUNN2.contex.table", header=T, sep=" ")
h4 <- read.table("align/H4_2023May19_S4.clean.nnunn.sort_NNUNN2.contex.table", header=T, sep=" ")
h5 <- read.table("align/H5_2023May19_S5.clean.nnunn.sort_NNUNN2.contex.table", header=T, sep=" ")
h6 <- read.table("align/H6_2023May19_S6.clean.nnunn.sort_NNUNN2.contex.table", header=T, sep=" ")
h7 <- read.table("align/H7_2023May19_S7.clean.nnunn.sort_NNUNN2.contex.table", header=T, sep=" ")
h8 <- read.table("align/H8_2023May19_S8.clean.nnunn.sort_NNUNN2.contex.table", header=T, sep=" ")

merge1 <- list(h1[,c("first", "last", "ratio")], h2[,c("first", "last", "ratio")],h3[,c("first", "last", "ratio")],h4[,c("first", "last", "ratio")],h5[,c("first", "last", "ratio")],h6[,c("first", "last", "ratio")],h7[,c("first", "last", "ratio")],h8[,c("first", "last", "ratio")])
merge1 <- merge1 %>% reduce(full_join, by=c("first", "last"))
colnames(merge1) <- c("first", "last", "h1_ratio", "h2_ratio", "h3_ratio", "h4_ratio", "h5_ratio", "h6_ratio", "h7_ratio", "h8_ratio")

merge1$ratio <- round(rowSums(merge1[,3:10])/8, 2)*100
mdf2 <- spread(merge1[,c("first", "last", "ratio")], key = last, value = ratio) %>%  as.data.frame()
final_matrix1 <- data.frame(mdf2[,-1], row.names=as.character(mdf2[,1])) %>% as.matrix(rownames = TRUE)
merge1$ratio <- rowSums(merge1[,3:10])/8
write.table(merge1, "align/unmodified_nnunn2.context.table", quote = FALSE, col.names = T, row.names = F)


library(pheatmap)
pdf("plot/unmodified_nnunn2.ratio.pdf")
pheatmap(final_matrix1,cluster_cols = FALSE, cluster_rows = FALSE, display_numbers = TRUE, number_format="%.1f", fontsize_number = 12, color=colorRampPalette(c("white", "red"))(length(breaksList)), breaks = breaksList, legend_breaks = seq(0, 100, 50), legend_labels = c("0", "50", "100"), border_color = "grey", number_color="grey0")
dev.off()

pdf("plot/unmodified_nnunn2.ratio.pdf")
ggplot(merge1[,c("first", "last", "ratio")], aes(x = first, y = last, fill = ratio)) + geom_tile(color = "gray", size=0.1)+ coord_fixed()+ ##to make Square tiles
  geom_text(aes(label = ratio), color = "black", size = 6) +
  scale_fill_gradient2(low = "white", high = "red", limits=c(0,100))
dev.off()


##for unmodified nnunn1 h1 and h23 (nnunn1)
h1_1 <- read.table("align/H1_2023May19_S1.clean.nnunn.sort_NNUNN1.contex.table", header=T, sep=" ")
h23 <- read.table("align/H23_2023Apr16_S19.clean.nnunn.sort_NNUNN1.contex.table", header=T, sep=" ")
merge1 <- list(h1_1[,c("first", "last", "ratio")],h23[,c("first", "last", "ratio")])
merge1 <- merge1 %>% reduce(full_join, by=c("first", "last"))
colnames(merge1) <- c("first", "last", "h1_1_ratio", "h23_ratio")

merge1$ratio <- round(rowSums(merge1[,3:4])/2, 2)*100
mdf2 <- spread(merge1[,c("first", "last", "ratio")], key = last, value = ratio) %>%  as.data.frame()
final_matrix1 <- data.frame(mdf2[,-1], row.names=as.character(mdf2[,1])) %>% as.matrix(rownames = TRUE)
merge1$ratio <- rowSums(merge1[,3:4])/2
write.table(merge1, "align/unmodified_nnunn1.context.table", quote = FALSE, col.names = T, row.names = F)


library(pheatmap)
pdf("plot/unmodified_nnunn1.ratio.pdf")
pheatmap(final_matrix1,cluster_cols = FALSE, cluster_rows = FALSE, display_numbers = TRUE, number_format="%.1f", fontsize_number = 12, color=colorRampPalette(c("white", "red"))(length(breaksList)), breaks = breaksList, legend_breaks = seq(0, 100, 50), legend_labels = c("0", "50", "100"), border_color = "grey",number_color="grey0")
dev.off()

pdf("plot/unmodified_nnunn1.ratio.pdf")
ggplot(merge1[,c("first", "last", "ratio")], aes(x = first, y = last, fill = ratio)) + geom_tile(color = "gray", size=0.1)+ coord_fixed()+ ##to make Square tiles
  geom_text(aes(label = ratio), color = "black", size = 6) +
  scale_fill_gradient2(low = "white", high = "red", limits=c(0,100))
  dev.off()


##############output unmodified_nnunn.context.table for calling
merge1 <- list(h1[,c("first", "last", "ratio")], h2[,c("first", "last", "ratio")],h3[,c("first", "last", "ratio")],h4[,c("first", "last", "ratio")],h5[,c("first", "last", "ratio")],h6[,c("first", "last", "ratio")],h7[,c("first", "last", "ratio")],h8[,c("first", "last", "ratio")],h1_1[,c("first", "last", "ratio")],h23[,c("first", "last", "ratio")])
merge1 <- merge1 %>% reduce(full_join, by=c("first", "last"))
colnames(merge1) <- c("first", "last", "h1_ratio", "h2_ratio", "h3_ratio", "h4_ratio", "h5_ratio", "h6_ratio", "h7_ratio", "h8_ratio","h1_1_ratio", "h23_ratio")
merge1$ratio <- rowSums(merge1[,3:4])/2
write.table(merge1, "align/unmodified_nnunn.context.table", quote = FALSE, col.names = T, row.names = F)



########################################################################plot 18S 28S and 5.8S and U2##############################################
dat1 <- read.table("align/H7_H9_allT_motif.bed", header=F, sep="\t") 
dat2 <- read.table("align/H8_H10_allT_motif.bed", header=F, sep="\t") 
dat3 <- read.table("align/H3_H16_allT_motif.bed", header=F, sep="\t") 
dat4 <- read.table("align/H23_H26_allT_motif.bed", header=F, sep="\t") 
colnames(dat1) <- c("chr", "start", "end", "ref", "treat_depth", "treat_T", "treat_C", "treat_conversion", "treat_gap", "treat_gap_rate", "motif", "ctrl_depth", "ctrl_T", "ctrl_C", "ctrl_conversion", "ctrl_gap", "ctrl_gap_rate")
colnames(dat2) <- c("chr", "start", "end", "ref", "treat_depth", "treat_T", "treat_C", "treat_conversion", "treat_gap", "treat_gap_rate", "motif", "ctrl_depth", "ctrl_T", "ctrl_C", "ctrl_conversion", "ctrl_gap", "ctrl_gap_rate")
colnames(dat3) <- c("chr", "start", "end", "ref", "treat_depth", "treat_T", "treat_C", "treat_conversion", "treat_gap", "treat_gap_rate", "motif", "ctrl_depth", "ctrl_T", "ctrl_C", "ctrl_conversion", "ctrl_gap", "ctrl_gap_rate")
colnames(dat4) <- c("chr", "start", "end", "ref", "treat_depth", "treat_T", "treat_C", "treat_conversion", "treat_gap", "treat_gap_rate", "motif", "ctrl_depth", "ctrl_T", "ctrl_C", "ctrl_conversion", "ctrl_gap", "ctrl_gap_rate")
dat1 <- dat1[dat1$chr=="NR_003285.3_RNA5_8SN5" | dat1$chr=="NR_003286.4_RNA18SN5" | dat1$chr=="NR_003287.4_RNA28SN5",]
dat2 <- dat2[dat2$chr=="NR_003285.3_RNA5_8SN5" | dat2$chr=="NR_003286.4_RNA18SN5" | dat2$chr=="NR_003287.4_RNA28SN5",]
dat3 <- dat3[dat3$chr=="NR_003285.3_RNA5_8SN5" | dat3$chr=="NR_003286.4_RNA18SN5" | dat3$chr=="NR_003287.4_RNA28SN5",]
dat4 <- dat4[dat4$chr=="NR_003285.3_RNA5_8SN5" | dat4$chr=="NR_003286.4_RNA18SN5" | dat4$chr=="NR_003287.4_RNA28SN5",]

mergedat <- list(dat1[,c("chr", "start", "end", "ref","treat_conversion","ctrl_conversion")], dat2[,c("chr", "start", "end", "ref","treat_conversion","ctrl_conversion")], dat3[,c("chr", "start", "end", "ref","treat_conversion","ctrl_conversion")], dat4[,c("chr", "start", "end", "ref","treat_conversion","ctrl_conversion")])
mergedat <- mergedat %>% reduce(full_join, by=c("chr", "start", "end", "ref"))%>% na.omit()
colnames(mergedat) <- c("chr", "start", "end", "ref","h7_conversion","h9_conversion","h8_conversion","h10_conversion","h3_conversion","h16_conversion","h23_conversion","h26_conversion")

mergedat$treat_conversion<- (mergedat$h7_conversion+mergedat$h8_conversion+mergedat$h3_conversion+mergedat$h23_conversion)/4
mergedat$ctrl_conversion <- (mergedat$h9_conversion+mergedat$h10_conversion+mergedat$h16_conversion+mergedat$h26_conversion)/4

seldat <- mergedat[,c("chr", "start", "end", "ctrl_conversion", "treat_conversion")]

#####################generate pseudouridine and other modification list
##cat align/rRNA_called_sites.txt | sed 's/ /\t/g' | awk '$1~"18S" || $1 ~"28S" || $1~"5_8S"' | awk '{OFS="\t"}{print $1, $2, $3,"pseu"}' > align/pseudouridine_list.txt
####add type information
###known pseudoU
###U
###other modifications
###

type <- read.table("align/pseudouridine_list.txt", header=F, sep="\t") 
colnames(type) <- c("chr", "start", "end", "types")
type <- rbind(c("NR_003286.4_RNA18SN5", "35", "36", "pseu"), type)
###append other modification
type <- rbind(c("NR_003286.4_RNA18SN5", "1247", "1248", "other"), type)
type <- rbind(c("NR_003287.4_RNA28SN5", "4529", "4530", "other"), type)

final_dat <- merge(seldat, type, by=c("chr", "start", "end"), all.x=TRUE)
final_dat[is.na(final_dat)] <- c("U")

pdf("plot/rRNA_scatter_legend.pdf", width=5,height = 4)
ggplot(final_dat, aes(x=round(ctrl_conversion*100,2), y=round(treat_conversion*100,2), color=types))+
geom_point(size = 2)+
scale_color_manual(values=c('#ABABAB','#000000', '#FF0000'))+ #ABABAB #FF0000
ylab("conversion ratio(%)")+xlab("")+
theme_bw()+scale_x_continuous(limits = c(-1,100), expand = c(0, 0)) +
  scale_y_continuous(limits = c(-1,100), expand = c(0, 0))+theme(axis.text = element_text(color = "black")) +geom_abline(intercept = 0,slope = 1,linetype=2)
dev.off()

pdf("plot/rRNA_scatter_figure.pdf")
ggplot(final_dat, aes(x=round(ctrl_conversion*100,2), y=round(treat_conversion*100,2), color=types))+
geom_point(size = 3)+
scale_color_manual(values=c('#ABABAB','#000000', '#FF0000'))+ #ABABAB #FF0000
ylab("conversion ratio(%)")+xlab("")+
theme_bw()+scale_x_continuous(limits = c(-0.8,100), expand = c(0, 0)) +
  scale_y_continuous(limits = c(-0.8,100), expand = c(0, 0))+theme(axis.text = element_text(color = "black"))+theme(legend.position="none") +geom_abline(intercept = 0,slope = 1,linetype=2)
dev.off()


#######################################################box plot for 18S 28S and 5.8S
###5.8S
dat <- seldat[seldat$chr=="NR_003285.3_RNA5_8SN5",] %>% dplyr::select(end,ctrl_conversion, treat_conversion)

meldat <- melt(dat, id.var="end")
pdf("plot/5.8S_plot.pdf", width=6,height = 4)
ggplot(meldat, aes(group=variable, col=variable))+geom_segment(aes(x=end, y=round(value,4)*100,xend=end,yend=0), size=0.75)+geom_point(aes(x=end, y=round(value,4)*100),size=4) +ylab("conversion ratio(%)")+xlab("")+
theme_bw()+scale_x_continuous(limits = c(0,200), expand = c(0, 0)) + scale_y_continuous(limits = c(-2,102), expand = c(0, 0))+scale_color_manual(values=c('#ABABAB','#FE767C'))
dev.off()



###18S
dat <- seldat[seldat$chr=="NR_003286.4_RNA18SN5",] %>% dplyr::select(end,ctrl_conversion, treat_conversion)

meldat <- melt(dat, id.var="end")
pdf("plot/18S_plot.pdf", width=12,height = 4)
ggplot(meldat, aes(group=variable, col=variable))+geom_segment(aes(x=end, y=round(value,4)*100,xend=end,yend=0), size=0.75)+geom_point(aes(x=end, y=round(value,4)*100),size=4) +ylab("conversion ratio(%)")+xlab("")+
theme_bw()+scale_x_continuous(limits = c(0,2000), expand = c(0, 0)) + scale_y_continuous(limits = c(-2,102), expand = c(0, 0))+scale_color_manual(values=c('#ABABAB','#FE767C'))
dev.off()



###28S
dat <- seldat[seldat$chr=="NR_003287.4_RNA28SN5",] %>% dplyr::select(end,ctrl_conversion, treat_conversion)
meldat <- melt(dat, id.var="end")
pdf("plot/28S_plot.pdf",width=16,height = 4)
ggplot(meldat, aes(group=variable, col=variable))+geom_segment(aes(x=end, y=round(value,4)*100,xend=end,yend=0), size=0.75)+geom_point(aes(x=end, y=round(value,4)*100),size=4) +ylab("conversion ratio(%)")+xlab("")+
theme_bw()+scale_x_continuous(limits = c(0,5500), expand = c(0, 0)) + scale_y_continuous(limits = c(-2,102), expand = c(0, 0))+scale_color_manual(values=c('#ABABAB','#FE767C'))
dev.off()



dat <- seldat[seldat$chr=="NR_003287.4_RNA28SN5",] %>% dplyr::select(end,ctrl_conversion, treat_conversion)

meldat <- melt(dat, id.var="end")
pdf("plot/28S_3700_3800.pdf",width=16,height = 4)
ggplot(meldat, aes(x=end, y=round(value,4)*100, fill=variable, color=variable))+ geom_bar(stat = "identity", width=0.50,position = position_dodge(), alpha = 0.75)+geom_point(size=2)+
ylab("conversion ratio(%)")+xlab("")+
theme_bw()+scale_x_continuous(limits = c(3700,3800), expand = c(0, 0)) + 
  scale_y_continuous(limits = c(0,100), expand = c(0, 0))+scale_fill_manual(values=c('#ABABAB','#FE767C'))+scale_color_manual(values=c('#ABABAB','#FE767C'))
dev.off()


####################PIE chart
df <- data.frame(
  group = c("5.8S", "18S", "28S", "16S mt-rRNA", "12S mt-rRNA"),
  value = c(2, 40, 62, 1,6)
  )

pdf("plot/rRNA_pie.pdf")
ggplot(df, aes(x="", y=value, fill=group))+
geom_bar(width = 1, stat = "identity")+coord_polar("y") + scale_fill_brewer(palette="Set2") +  ###or RdYlBu
  geom_text(aes(label = value),position = position_stack(vjust = 0.5))+
  theme_void()
dev.off()

###"#B2182B" "#D6604D" "#F4A582" "#FDDBC7" "#D1E5F0" "#92C5DE" "#4393C3" "#2166AC"
###brewer.pal(n = 8, name = "RdBu")


########################################raw signal correlation scatter plot ##########################
###H7 vs H8 H3 vs H23
###H7 vs H3
###H7 vs H23

dat1 <- read.table("align/H7_H9_allT_motif.bed", header=F, sep="\t") 
dat2 <- read.table("align/H8_H10_allT_motif.bed", header=F, sep="\t") 
dat3 <- read.table("align/H3_H16_allT_motif.bed", header=F, sep="\t") 
dat4 <- read.table("align/H23_H26_allT_motif.bed", header=F, sep="\t") 
colnames(dat1) <- c("chr", "start", "end", "ref", "treat_depth", "treat_T", "treat_C", "treat_conversion", "treat_gap", "treat_gap_rate", "motif", "ctrl_depth", "ctrl_T", "ctrl_C", "ctrl_conversion", "ctrl_gap", "ctrl_gap_rate")
colnames(dat2) <- c("chr", "start", "end", "ref", "treat_depth", "treat_T", "treat_C", "treat_conversion", "treat_gap", "treat_gap_rate", "motif", "ctrl_depth", "ctrl_T", "ctrl_C", "ctrl_conversion", "ctrl_gap", "ctrl_gap_rate")
colnames(dat3) <- c("chr", "start", "end", "ref", "treat_depth", "treat_T", "treat_C", "treat_conversion", "treat_gap", "treat_gap_rate", "motif", "ctrl_depth", "ctrl_T", "ctrl_C", "ctrl_conversion", "ctrl_gap", "ctrl_gap_rate")
colnames(dat4) <- c("chr", "start", "end", "ref", "treat_depth", "treat_T", "treat_C", "treat_conversion", "treat_gap", "treat_gap_rate", "motif", "ctrl_depth", "ctrl_T", "ctrl_C", "ctrl_conversion", "ctrl_gap", "ctrl_gap_rate")
dat1 <- dat1[dat1$chr=="NR_003285.3_RNA5_8SN5" | dat1$chr=="NR_003286.4_RNA18SN5" | dat1$chr=="NR_003287.4_RNA28SN5" | dat1$chr=="NR_023363.1_RNA5S1",]
dat2 <- dat2[dat2$chr=="NR_003285.3_RNA5_8SN5" | dat2$chr=="NR_003286.4_RNA18SN5" | dat2$chr=="NR_003287.4_RNA28SN5" | dat2$chr=="NR_023363.1_RNA5S1",]
dat3 <- dat3[dat3$chr=="NR_003285.3_RNA5_8SN5" | dat3$chr=="NR_003286.4_RNA18SN5" | dat3$chr=="NR_003287.4_RNA28SN5"| dat3$chr=="NR_023363.1_RNA5S1",]
dat4 <- dat4[dat4$chr=="NR_003285.3_RNA5_8SN5" | dat4$chr=="NR_003286.4_RNA18SN5" | dat4$chr=="NR_003287.4_RNA28SN5"| dat4$chr=="NR_023363.1_RNA5S1",]

merge1 <- merge(dat1[,c("chr", "start", "end", "treat_conversion")], dat2[,c("chr", "start", "end", "treat_conversion")], by=c("chr", "start", "end"))
colnames(merge1) <- c("chr", "start", "end", "signal1", "signal2")
pdf("plot/rRNA_raw_signal_H7_vs_H8.pdf",width=4, height=4)
smoothScatter(merge1$signal1*100, merge1$signal2*100, xlab = "H7 raw signal for rRNA", ylab = "H8 raw signal for rRNA")
dev.off()
cor(merge1$signal1*100, merge1$signal2*100) ###0.999979


merge1 <- merge(dat3[,c("chr", "start", "end", "treat_conversion")], dat4[,c("chr", "start", "end", "treat_conversion")], by=c("chr", "start", "end"))
colnames(merge1) <- c("chr", "start", "end", "signal1", "signal2")
pdf("plot/rRNA_raw_signal_H3_vs_H23.pdf",width=4, height=4)
smoothScatter(merge1$signal1*100, merge1$signal2*100, xlab = "H3 raw signal for rRNA", ylab = "H23 raw signal for rRNA")
dev.off()
cor(merge1$signal1*100, merge1$signal2*100) ###0.9999201


merge1 <- merge(dat1[,c("chr", "start", "end", "treat_conversion")], dat3[,c("chr", "start", "end", "treat_conversion")], by=c("chr", "start", "end"))
colnames(merge1) <- c("chr", "start", "end", "signal1", "signal2")
pdf("plot/rRNA_raw_signal_H7_vs_H3.pdf",width=4, height=4)
smoothScatter(merge1$signal1*100, merge1$signal2*100, xlab = "H7 raw signal for rRNA", ylab = "H3 raw signal for rRNA")
dev.off()
cor(merge1$signal1*100, merge1$signal2*100) ###0.9998254


merge1 <- merge(dat1[,c("chr", "start", "end", "treat_conversion")], dat4[,c("chr", "start", "end", "treat_conversion")], by=c("chr", "start", "end"))
colnames(merge1) <- c("chr", "start", "end", "signal1", "signal2")
pdf("plot/rRNA_raw_signal_H7_vs_H23.pdf",width=4, height=4)
smoothScatter(merge1$signal1*100, merge1$signal2*100, xlab = "H7 raw signal for rRNA", ylab = "H23 raw signal for rRNA")
dev.off()
cor(merge1$signal1*100, merge1$signal2*100) ### 0.9998223

#################Venn gram for comparison with MS data
library(VennDiagram)##########to be revised, include pos 36 in 18S rRNA. the number in this plot and the pie chart needs to be revised
pdf("plot/bacs_ms_venn.pdf")
draw.pairwise.venn(area1=104 , area2=105,cross.area=103, 
                   category=c("bacs","SILNAS MS"),col=c("#FE767C","#666666"), cat.pos=c(-27, 27), cat.dist=c(0.03, 0.03), cat.col=c("#FE767C","#666666"), height = 300 , 
        width = 300)
dev.off()


########################
####a pie chart for major and minor splicesome RNAs
###a stacked plot for 18S, 28S and 5.8S rRNA

mergedat <- read.table("align/rRNA_called_sites.txt", header=T, sep=" ", colClasses="character") 
mergedat <- mergedat[mergedat$chr!="NR_003285.3_RNA5_8SN5" & mergedat$chr!="NR_003286.4_RNA18SN5" & mergedat$chr!="NR_003287.4_RNA28SN5" & mergedat$chr!="2kb_spikein"& mergedat$chr!="MT_RNR1_ENST00000389680.2_ENSE00001544499" & mergedat$chr!="MT_RNR2_ENST00000387347.2_ENSE00001544497" & mergedat$chr!="spike_in",]
site <- mergedat[mergedat$h7_level>=0.05 & mergedat$h8_level>=0.05 & mergedat$h3_level>=0.05 & mergedat$h23_level>=0.05,]
dat <- site %>% group_by(chr) %>% summarize(n = n()) 
dat1 <-dat[dat$chr!="NR_001566.1_TERC" & dat$chr!="NR_023317.1_RNU7_1" & dat$chr!="NR_001445.2_RN7SK",]
dat1$chr <- factor(dat1$chr, levels=c("NR_004430.2_RNU1_1","NR_002716.3_RNU2_1", "NR_003925.1_RNU4_1", "NR_002756.2_RNU5A_1", "NR_004394.1_RNU6_1", "NR_029422.2RNU12", "NR_023343.1_RNU4ATAC", "NR_023344.1_RNU6ATAC"))

  chr                      n
  <fct>                <int>
1 NR_002716.3_RNU2_1      14
2 NR_002756.2_RNU5A_1      4
3 NR_003925.1_RNU4_1       3
4 NR_004394.1_RNU6_1       4
5 NR_004430.2_RNU1_1       2
6 NR_023343.1_RNU4ATAC     2
7 NR_023344.1_RNU6ATAC     1
8 NR_029422.2RNU12         2

pdf("plot/spliceosomal_RNA_pie.pdf")
ggplot(dat1, aes(x="", y=n, fill=chr))+
geom_bar(width = 1, stat = "identity")+coord_polar("y") + scale_fill_brewer(palette="Set2") +  ###or RdYlBu
  geom_text(aes(label = n),position = position_stack(vjust = 0.5))+
  theme_void()
dev.off()

#################Venn gram for comparison with MS data on major spliceosomal RNA
library(VennDiagram)##########
pdf("plot/bacs_ms_venn_on_spliceosomal.pdf")
draw.pairwise.venn(area1=27 , area2=28,cross.area=27, 
                   category=c("bacs","SILNAS MS"),col=c("#FE767C","#666666"), cat.pos=c(-27, 27), cat.dist=c(0.03, 0.03), cat.col=c("#FE767C","#666666"), height = 300 , 
        width = 300)
dev.off()

################plot conversion/input on U2
dat1 <- read.table("align/H7_H9_allT_motif.bed", header=F, sep="\t") 
dat2 <- read.table("align/H8_H10_allT_motif.bed", header=F, sep="\t") 
dat3 <- read.table("align/H3_H16_allT_motif.bed", header=F, sep="\t") 
dat4 <- read.table("align/H23_H26_allT_motif.bed", header=F, sep="\t") 
colnames(dat1) <- c("chr", "start", "end", "ref", "treat_depth", "treat_T", "treat_C", "treat_conversion", "treat_gap", "treat_gap_rate", "motif", "ctrl_depth", "ctrl_T", "ctrl_C", "ctrl_conversion", "ctrl_gap", "ctrl_gap_rate")
colnames(dat2) <- c("chr", "start", "end", "ref", "treat_depth", "treat_T", "treat_C", "treat_conversion", "treat_gap", "treat_gap_rate", "motif", "ctrl_depth", "ctrl_T", "ctrl_C", "ctrl_conversion", "ctrl_gap", "ctrl_gap_rate")
colnames(dat3) <- c("chr", "start", "end", "ref", "treat_depth", "treat_T", "treat_C", "treat_conversion", "treat_gap", "treat_gap_rate", "motif", "ctrl_depth", "ctrl_T", "ctrl_C", "ctrl_conversion", "ctrl_gap", "ctrl_gap_rate")
colnames(dat4) <- c("chr", "start", "end", "ref", "treat_depth", "treat_T", "treat_C", "treat_conversion", "treat_gap", "treat_gap_rate", "motif", "ctrl_depth", "ctrl_T", "ctrl_C", "ctrl_conversion", "ctrl_gap", "ctrl_gap_rate")
dat1 <- dat1[dat1$chr=="NR_002716.3_RNU2_1",]
dat2 <- dat2[dat2$chr=="NR_002716.3_RNU2_1",]
dat3 <- dat3[dat3$chr=="NR_002716.3_RNU2_1",]
dat4 <- dat4[dat4$chr=="NR_002716.3_RNU2_1",]

mergedat <- list(dat1[,c("chr", "start", "end", "ref","treat_conversion","ctrl_conversion")], dat2[,c("chr", "start", "end", "ref","treat_conversion","ctrl_conversion")], dat3[,c("chr", "start", "end", "ref","treat_conversion","ctrl_conversion")], dat4[,c("chr", "start", "end", "ref","treat_conversion","ctrl_conversion")])
mergedat <- mergedat %>% reduce(full_join, by=c("chr", "start", "end", "ref"))%>% na.omit()
colnames(mergedat) <- c("chr", "start", "end", "ref","h7_conversion","h9_conversion","h8_conversion","h10_conversion","h3_conversion","h16_conversion","h23_conversion","h26_conversion")

mergedat$treat_conversion<- (mergedat$h7_conversion+mergedat$h8_conversion+mergedat$h3_conversion+mergedat$h23_conversion)/4
mergedat$ctrl_conversion <- (mergedat$h9_conversion+mergedat$h10_conversion+mergedat$h16_conversion+mergedat$h26_conversion)/4

seldat <- mergedat[,c("chr", "start", "end", "ctrl_conversion", "treat_conversion")]

dat <- seldat[seldat$chr=="NR_002716.3_RNU2_1",] %>% dplyr::select(end,ctrl_conversion, treat_conversion)

meldat <- melt(dat, id.var="end")
pdf("plot/U2_plot.pdf", width=12,height = 4)
ggplot(meldat, aes(group=variable, col=variable))+geom_segment(aes(x=end, y=round(value,4)*100,xend=end,yend=0), size=0.75)+geom_point(aes(x=end, y=round(value,4)*100),size=4) +ylab("conversion ratio(%)")+xlab("")+
theme_bw()+scale_x_continuous(limits = c(0,200), expand = c(0, 0)) + scale_y_continuous(limits = c(-2,102), expand = c(0, 0))+scale_color_manual(values=c('#ABABAB','#FE767C'))
dev.off()


################plot conversion/input on TERC and U4
dat1 <- read.table("align/H7_H9_allT_motif.bed", header=F, sep="\t") 
dat2 <- read.table("align/H8_H10_allT_motif.bed", header=F, sep="\t") 
dat3 <- read.table("align/H3_H16_allT_motif.bed", header=F, sep="\t") 
dat4 <- read.table("align/H23_H26_allT_motif.bed", header=F, sep="\t") 
colnames(dat1) <- c("chr", "start", "end", "ref", "treat_depth", "treat_T", "treat_C", "treat_conversion", "treat_gap", "treat_gap_rate", "motif", "ctrl_depth", "ctrl_T", "ctrl_C", "ctrl_conversion", "ctrl_gap", "ctrl_gap_rate")
colnames(dat2) <- c("chr", "start", "end", "ref", "treat_depth", "treat_T", "treat_C", "treat_conversion", "treat_gap", "treat_gap_rate", "motif", "ctrl_depth", "ctrl_T", "ctrl_C", "ctrl_conversion", "ctrl_gap", "ctrl_gap_rate")
colnames(dat3) <- c("chr", "start", "end", "ref", "treat_depth", "treat_T", "treat_C", "treat_conversion", "treat_gap", "treat_gap_rate", "motif", "ctrl_depth", "ctrl_T", "ctrl_C", "ctrl_conversion", "ctrl_gap", "ctrl_gap_rate")
colnames(dat4) <- c("chr", "start", "end", "ref", "treat_depth", "treat_T", "treat_C", "treat_conversion", "treat_gap", "treat_gap_rate", "motif", "ctrl_depth", "ctrl_T", "ctrl_C", "ctrl_conversion", "ctrl_gap", "ctrl_gap_rate")
dat1 <- dat1[dat1$chr=="NR_001566.1_TERC" | dat1$chr=="NR_023343.1_RNU4ATAC",]
dat2 <- dat2[dat2$chr=="NR_001566.1_TERC" | dat2$chr=="NR_023343.1_RNU4ATAC",]
dat3 <- dat3[dat3$chr=="NR_001566.1_TERC" | dat3$chr=="NR_023343.1_RNU4ATAC",]
dat4 <- dat4[dat4$chr=="NR_001566.1_TERC" | dat4$chr=="NR_023343.1_RNU4ATAC",]

mergedat <- list(dat1[,c("chr", "start", "end", "ref","treat_conversion","ctrl_conversion")], dat2[,c("chr", "start", "end", "ref","treat_conversion","ctrl_conversion")], dat3[,c("chr", "start", "end", "ref","treat_conversion","ctrl_conversion")], dat4[,c("chr", "start", "end", "ref","treat_conversion","ctrl_conversion")])
mergedat <- mergedat %>% reduce(full_join, by=c("chr", "start", "end", "ref"))%>% na.omit()
colnames(mergedat) <- c("chr", "start", "end", "ref","h7_conversion","h9_conversion","h8_conversion","h10_conversion","h3_conversion","h16_conversion","h23_conversion","h26_conversion")

mergedat$treat_conversion<- (mergedat$h7_conversion+mergedat$h8_conversion+mergedat$h3_conversion+mergedat$h23_conversion)/4
mergedat$ctrl_conversion <- (mergedat$h9_conversion+mergedat$h10_conversion+mergedat$h16_conversion+mergedat$h26_conversion)/4

seldat <- mergedat[,c("chr", "start", "end", "ctrl_conversion", "treat_conversion")]

terc <- seldat[seldat$chr=="NR_001566.1_TERC",] %>% dplyr::select(end,ctrl_conversion, treat_conversion)

melterc <- melt(terc, id.var="end")
pdf("plot/TERC_plot.pdf", width=18,height = 4)
ggplot(melterc, aes(group=variable, col=variable))+geom_segment(aes(x=end, y=round(value,4)*100,xend=end,yend=0), size=0.65)+geom_point(aes(x=end, y=round(value,4)*100),size=4) +ylab("conversion ratio(%)")+xlab("")+geom_hline(yintercept=5, linetype="dashed", color = "grey")+
theme_bw() +scale_x_continuous(limits = c(0,500), expand = c(0, 0))+ scale_y_continuous(limits = c(-2,102), expand = c(0, 0))+scale_color_manual(values=c('#ABABAB','#FE767C'))
dev.off()

u4atac <- seldat[seldat$chr=="NR_023343.1_RNU4ATAC",] %>% dplyr::select(end,ctrl_conversion, treat_conversion)

melu4atac <- melt(u4atac, id.var="end")
pdf("plot/U4ATAC_plot.pdf", width=12,height = 4)
ggplot(melu4atac, aes(group=variable, col=variable))+geom_segment(aes(x=end, y=round(value,4)*100,xend=end,yend=0), size=0.65)+geom_point(aes(x=end, y=round(value,4)*100),size=4) +ylab("conversion ratio(%)")+xlab("")+
theme_bw() +scale_x_continuous(limits = c(0,150), expand = c(0, 0))+ scale_y_continuous(limits = c(-2,102), expand = c(0, 0))+scale_color_manual(values=c('#ABABAB','#FE767C'))
dev.off()



######################################################plot for comparision of rRNA between cytosolic and rRNA
mergedat <- read.table("align/rRNA_called_sites.txt", header=T, sep=" ", colClasses="character") 
mergedat <- mergedat[mergedat$chr=="NR_003285.3_RNA5_8SN5" | mergedat$chr=="NR_003286.4_RNA18SN5" | mergedat$chr=="NR_003287.4_RNA28SN5" |  mergedat$chr=="MT_RNR1_ENST00000389680.2_ENSE00001544499" | mergedat$chr=="MT_RNR2_ENST00000387347.2_ENSE00001544497",]
site <- mergedat[mergedat$h7_level>=0.05 & mergedat$h8_level>=0.05 & mergedat$h3_level>=0.05 & mergedat$h23_level>=0.05,]
site$mean_level <- (as.numeric(site$h3_level)+as.numeric(site$h23_level)+as.numeric(site$h7_level)+as.numeric(site$h8_level))/4

site$name <- gsub("MT_RNR1_ENST00000389680.2_ENSE00001544499", "mt_rRNA",site$chr)
site$name <- gsub("MT_RNR2_ENST00000387347.2_ENSE00001544497", "mt_rRNA",site$name)
site$name <- gsub("NR_003285.3_RNA5_8SN5", "rRNA",site$name)
site$name <- gsub("NR_003286.4_RNA18SN5", "rRNA",site$name)
site$name <- gsub("NR_003287.4_RNA28SN5", "rRNA",site$name)

selsite <- site[,c("name", "mean_level")]
selsite$name <- factor(selsite$name, level=c("rRNA", "mt_rRNA"))
selsite[selsite$mean_level>=1,]$mean_level <- 1

p <- ggplot(selsite,aes(x=name,y=mean_level*100, fill=name))+ geom_boxplot(width=0.8, outlier.size=0.5)+ #+ geom_jitter(size=2,shape=16, colour = "grey",position=position_jitter(0.2))
  #stat_summary(fun = mean, geom = "point",color = "black") +
  theme_light() + xlab("") +  theme(axis.text = element_text(color = "black")) + ylim(0,100) +
  ylab("level") +
  scale_fill_manual(values=c("#006699","#F4A582"))

ggsave("plot/modification_level_cytorRNA_vsmtrRNA.pdf",p, width=3, height = 3)

mean(selsite[selsite$name=="rRNA",]$mean_level) ##0.8145089
mean(selsite[selsite$name=="mt_rRNA",]$mean_level) #0.2410905
