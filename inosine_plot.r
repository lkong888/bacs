##module load R/4.0.3-foss-2020b
library(dplyr)
library(reshape2)
library(tidyr)
library(ggplot2)
library(cowplot)
library(tidyverse)
library(parallel)
library(Biostrings)
library(RColorBrewer)
library(ggpubr)
library(stringi)
setwd("/users/ludwig/ebu571/ebu571/bacs")

#########################################################plot after annotation######################################################
annotated_site <- read.table("site/treat.unique.annotation.bed",  header=TRUE) ##1337
filtered_site <- annotated_site[annotated_site$strand_infor=="proper_strand"&annotated_site$transcript_type !="snRNA" &annotated_site$transcript_detail!="intronic", ]
nrow(filtered_site) ##1335
psi_gene_count <- filtered_site %>% group_by(gene_name) %>% summarize(n = n())
nrow(psi_gene_count)  ##1103

inosine <- read.table("inosine/inosine.unique.annotation.bed", header=TRUE)
filtered_inosine <- inosine[inosine$strand_infor=="proper_strand",] ##852

nb.cols <- 4
mycolors <- colorRampPalette(brewer.pal(4, "Set2"))(nb.cols)

feature_count <- inosine %>% group_by(transcript_detail) %>% summarize(n=n())
feature_count$transcript_detail <- factor(feature_count$transcript_detail, level=c("5UTR", "CDS", "3UTR", "exon"))

  transcript_detail     n
  <fct>             <int>
1 3UTR                711
2 5UTR                 12
3 CDS                  14
4 exon                115

pdf("plot/inosine.feature.distribution_pie.pdf")
ggplot(feature_count, aes(x="", y=n, fill=transcript_detail))+
geom_bar(width = 1, stat = "identity")+coord_polar("y") + scale_fill_manual(values = mycolors) + 
  geom_text(aes(label = n),position = position_stack(vjust = 0.5))+
  theme_void()
dev.off()

nrow(filtered_inosine) ###852
inosine_gene_count <- filtered_inosine %>% group_by(gene_name) %>% summarize(n = n())
nrow(inosine_gene_count) ##289

####overlapped gene
overlapgene <- merge(psi_gene_count,inosine_gene_count,by=c("gene_name")) #28

library(VennDiagram)##########
pdf("plot/bacs_psi_inosine_venn_on_gene.pdf")
draw.pairwise.venn(area1=1103 , area2=289,cross.area=28, 
                   category=c("psi","inosine"),col=c("#1b9e77","#404080"), fill=c("#1b9e77","#404080"),cat.pos=c(-27, 27), cat.dist=c(0.03, 0.03), cat.col=c("#1b9e77","#404080"), height = 300 , 
        width = 300)
dev.off()

psi_gene_count$type <- "psi"
inosine_gene_count$type <- "inosine"
gene_count <- rbind(psi_gene_count[,c("n", "type")],inosine_gene_count[,c("n", "type")])
pdf("plot/psi_inosine_gene_count_distribution_histogram.pdf",width=5,height=4)
ggplot(gene_count, aes(x=n, fill=type)) +  geom_histogram(alpha=0.6, position = "dodge",binwidth=1)+scale_fill_manual(values=c("#69b3a2", "#404080")) + theme_bw() +
    labs(fill="")
dev.off()


inosine_re <- read.table("inosine/inosine_site_repeat_annotated.bed", header=FALSE) ##1062

##rm duplicated annotations
inosine_re <- inosine_re[!duplicated(inosine_re[,c("V1","V2","V3")]),] ###1062
inosine_re_count <- inosine_re %>% group_by(V39) %>% summarize(n=n())##58
for (i in 1:nrow(inosine_re_count)){
  if(grepl("DNA",inosine_re_count[i,]$V39)){
    inosine_re_count$type[i] <- "DNA_transposons"
  }else if(grepl("LINE", inosine_re_count[i,]$V39)){
    inosine_re_count$type[i] <- "LINE"
  }else if(grepl("SINE", inosine_re_count[i,]$V39)){
    inosine_re_count$type[i] <- "SINE"
  }else if(grepl("LTR", inosine_re_count[i,]$V39)){
    inosine_re_count$type[i] <- "LTR"
  }else if(grepl("Satellite", inosine_re_count[i,]$V39)){
    inosine_re_count$type[i] <- "Satellite"
  }else if(grepl("Simple", inosine_re_count[i,]$V39)){
    inosine_re_count$type[i] <- "Simple_repeat"
  }else{
    inosine_re_count$type[i] <- "others"
  }
}
inosine_re_count_agg <- aggregate(n ~ type, inosine_re_count,sum)

df2 <- inosine_re_count_agg %>% 
  mutate(csum = rev(cumsum(rev(n))), 
         pos = n/2 + lead(csum, 1),
         pos = if_else(is.na(pos), n/2, pos))
library(ggrepel)
pdf("plot/inosine.repeat.distribution_pie.pdf")
ggplot(inosine_re_count_agg, aes(x="", y=n, fill=type))+
geom_bar(width = 1, stat = "identity")+coord_polar("y") +
geom_label_repel(data = df2,
                   aes(y = pos, label = paste0(n)),
                   size = 4.5, nudge_x = 1, show.legend = FALSE) +

  #geom_text(aes(label = n),position = position_stack(vjust = 0.5))+
  theme_void()
dev.off()


#################for
psi_re <- read.table("site/psi_repeat.bed", header=FALSE) ##46
psi_re_count <- psi_re %>% group_by(V8) %>% summarize(n=n())
for (i in 1:nrow(psi_re_count)){
  if(grepl("DNA",psi_re_count[i,]$V8)){
    psi_re_count$type[i] <- "DNA_transposons"
  }else if(grepl("LINE", psi_re_count[i,]$V8)){
    psi_re_count$type[i] <- "LINE"
  }else if(grepl("SINE", psi_re_count[i,]$V8)){
    psi_re_count$type[i] <- "SINE"
  }else if(grepl("LTR", psi_re_count[i,]$V8)){
    psi_re_count$type[i] <- "LTR"
  }else if(grepl("Satellite", psi_re_count[i,]$V8)){
    psi_re_count$type[i] <- "Satellite"
  }else if(grepl("Simple", psi_re_count[i,]$V8)){
    psi_re_count$type[i] <- "Simple_repeat"
  }else{
    psi_re_count$type[i] <- "others"
  }
}
psi_re_count_agg <- aggregate(n ~ type, psi_re_count,sum)
psi_re_count_agg <- psi_re_count_agg[psi_re_count_agg$type!="rRNA" & psi_re_count_agg$type!="snRNA" & psi_re_count_agg$type!="tRNA",]

pdf("plot/psi.repeat.distribution_pie.pdf")
ggplot(psi_re_count_agg, aes(x="", y=n, fill=type))+
geom_bar(width = 1, stat = "identity")+coord_polar("y") +
  geom_text(aes(label = n),position = position_stack(vjust = 0.5))+
  theme_void()
dev.off()


###organize inosine site
inosine <- read.table("inosine/input_treat.mRNA_inosine_site_deltaG.txt",  header=TRUE) ##1181
nrow(inosine[inosine$snp!="SNP",]) ##1121

inosine_annotation <- read.table("inosine/inosine.unique.annotation.bed", header=TRUE)
filtered_inosine <- inosine_annotation[inosine_annotation$strand_infor=="proper_strand",] ##852

inosine_re <- read.table("inosine/inosine_site_repeat_annotated.bed", header=FALSE) ##1062

final_inosine <- merge(inosine[inosine$snp!="SNP",],inosine_re[,c("V1", "V2", "V3", "V37")], by.x=c("chr", "start", "end"), by.y=c("V1", "V2", "V3"),all.x=TRUE)##1121
final_inosine <- merge(final_inosine, filtered_inosine, by=c("chr", "start", "end", "strand"), all.x=TRUE) ##1121
write.table(final_inosine,"site/hela_inosine_site.txt",sep='\t',quote=F,row.names=F)