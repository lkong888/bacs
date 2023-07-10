##module load R/4.0.3-foss-2020b
####this script is for figures on snoRNA
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

####venn plot on site distribution of snoRNA
site <- read.table("site/snoRNA_tRNA_called_sites_fdrpassed.txt",header=T, sep="\t",colClasses=c("ref"="character")) 
site <- site[!grepl("tRNA", site$chr),] ###282

pattern <- "C/D_box|H/ACA_box|Cajal_body"
site$category <- regmatches(site$chr, regexpr(pattern, site$chr, ignore.case = TRUE))

#######read TERC site
terc <- read.table("align/rRNA_called_sites.txt", header=T, sep=" ", colClasses=c("ref"="character")) 
terc <- terc[terc$chr=="NR_001566.1_TERC",]
terc <- terc[terc$h7_level>=0.05 & terc$h8_level>=0.05 & terc$h3_level>=0.05 & terc$h23_level>=0.05,]
nrow(terc) ##7

count <- site %>% group_by(category) %>% summarize(n = n())
count <- rbind(count, c("TERC", "7"))

  category   n    
  <chr>      <chr>
1 C/D_box    192  
2 Cajal_body 28   
3 H/ACA_box  62   
4 TERC       7 

count$category <- factor(count$category, level=c("C/D_box", "H/ACA_box", "Cajal_body", "TERC"))
pdf("plot/snoRNA_pie.pdf")
ggplot(count, aes(x="", y=as.numeric(n), fill=category))+
geom_bar(width = 1, stat = "identity")+coord_polar("y") + scale_fill_brewer(palette="Set2") +  ###or RdYlBu
  geom_text(aes(label = as.numeric(n)),position = position_stack(vjust = 0.5))+
  theme_void()
dev.off()




#############################stack plot for snoRNA
snosite <- site[,c("chr", "start", "end", "motif", "treat1_level", "treat2_level", "treat3_level", "category")]
snosite$mean_level <- (snosite$treat1_level+snosite$treat2_level+snosite$treat3_level)/3
snosite$group <- cut(snosite$mean_level*100, breaks=c(5, 20,50,100))
snosite$group <-factor(snosite$group, level=c("(50,100]","(20,50]", "(5,20]"))
snosite$category <- factor(snosite$category, level=c("C/D_box", "H/ACA_box", "Cajal_body"))
write.table(snosite,"plot/snosite_motif_table.txt",sep='\t',quote=F,row.names=F)
###create stack plot
pdf("plot/snoRNA_stack_plot.pdf", width=7,height=6)
ggplot(snosite, aes(x = category, fill = group)) + 
  geom_bar(width=0.6) +scale_fill_manual(values=c('#f03e3e', '#fdb462', '#8bc53f')) +geom_text(stat = "count", aes(label = after_stat(count)),position = position_stack(vjust = 0.5))+theme_bw()+ylab("site counts")
dev.off()

snosite %>% group_by(category,group) %>% summarize(n = n())
  category   group        n
  <chr>      <fct>    <int>
1 C/D_box    (5,20]      83
2 C/D_box    (20,50]     42
3 C/D_box    (50,100]    67
4 Cajal_body (5,20]      12
5 Cajal_body (20,50]     11
6 Cajal_body (50,100]     5
7 H/ACA_box  (5,20]      30
8 H/ACA_box  (20,50]     12
9 H/ACA_box  (50,100]    20

########create dot plot for three different category ####to be revised
pdf("plot/snoRNA_point.pdf")
ggplot(snosite, aes(x = category, y=mean_level, fill = category)) + geom_jitter(aes(colour = category),width = 0.2)+scale_color_brewer(palette="Set2")+theme_bw()
dev.off()

pdf("plot/snoRNA_dot.pdf", width=7,height=6)
ggplot(snosite, aes(x = category, y=mean_level)) +geom_boxplot(fill="white",width=0.6)+geom_jitter(size=3,shape=16,aes(colour = category),width = 0.2)+scale_colour_brewer(palette="Set2")+theme_bw()
dev.off()


########plot motif frequency for snoRNA
level <- aggregate(mean_level ~ motif, snosite, mean)  ###nrow(level) 115
motif_freq <-snosite %>% group_by(motif) %>% summarize(n = n()) ###nrow 115
motif <- merge(level, motif_freq, by=c("motif")) ###115
motif$group <-cut(motif$n, breaks=c(0,5,30))

pdf("plot/snoRNA_motif_frequency.pdf",width=4, height=4)
ggplot(motif, aes(x = n, y=mean_level*100, fill = motif)) + geom_point(aes(color=group), size=3)+theme(legend.position = "none") +scale_color_manual(values=c("#808080","#e41a1c"))+
geom_text(aes(label=ifelse(n>5,as.character(motif),'')))+xlim(0,40)+ylim(0,100)+ylab("average modification level")+xlab("motif frequency")+theme_bw() +theme(legend.position = "none") 
dev.off()


################remove some isoforms
snosite1 <- read.table("plot/snosite_motif_table_modified.txt", header=T, sep="\t")
nrow(snosite1) ##196


level <- aggregate(mean_level ~ motif, snosite1, mean)  ###nrow(level) 115
motif_freq <-snosite1 %>% group_by(motif) %>% summarize(n = n()) ###nrow 115
motif <- merge(level, motif_freq, by=c("motif")) ###115
motif$group <-cut(motif$n, breaks=c(0,2,20))

pdf("plot/snoRNA_motif_frequency_removeisoform.pdf",width=4, height=4)
ggplot(motif, aes(x = n, y=mean_level*100, fill = motif)) + geom_point(aes(color=group), size=3)+theme(legend.position = "none") +scale_color_manual(values=c("#808080","#e41a1c"))+
geom_text(aes(label=ifelse(n>2,as.character(motif),'')))+xlim(0,10)+ylim(0,100)+ylab("average modification level")+xlab("motif frequency")+theme_bw() +theme(legend.position = "none") 
dev.off()




############################snoRNA annotation
# Extract the desired part before the second underscore for each value
snosite$refname <- apply(snosite, 1, function(x) {
  input <- x["chr"]
  result <- unlist(strsplit(input, "_"))[1:2]
  paste(result, collapse = "_")
})
write.table(snosite,"site/snosite_namechanged.txt",sep='\t',quote=F,row.names=F)
###read annotation file
snoanno <- read.table("resource/snoRNA_annptation.txt",header=F, sep=" ") 
colnames(snoanno) <- c("snoid", "refname", "start", "end", "structure", "seq")

snomerge <- merge(snosite, snoanno, by=c(refname))




#############################accumulation of pseudosites accross 5' to 3'
snobin <- read.table("site/snoRNA_density_40bins.txt", header=F, sep="\t")
colnames(snobin) <- c("chr", "start", "end", "binnum", "count")
aggsno <- aggregate(count ~ binnum, snobin, sum)
aggsno$density <- aggsno$count/sum(aggsno$count)

pdf("plot/snoRNA_density_40.pdf",width=4, height=4)
ggplot(aggsno, aes(x = binnum, y=density))+geom_line()
dev.off()