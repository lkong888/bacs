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
library(stringr)
library(VennDiagram)
setwd("/users/ludwig/ebu571/ebu571/bacs")
##for tRNA and others: /users/ludwig/ebu571/ebu571/27Oct2023/site/snoRNA_tRNA_called_sites.txt and /users/ludwig/ebu571/ebu571/27Oct2023/site/bacs_rrna_all_sites.txt 

####venn plot on site distribution of snoRNA
site <- read.table("/users/ludwig/ebu571/ebu571/27Oct2023/site/snoRNA_tRNA_called_sites.txt", header=TRUE, sep=" ") 
site <- site[!grepl("tRNA", site$chr),] ###304

pattern <- "C/D_box|H/ACA_box|Cajal_body"
site$category <- regmatches(site$chr, regexpr(pattern, site$chr, ignore.case = TRUE))

#######read TERC site
terc <- read.table("/users/ludwig/ebu571/ebu571/27Oct2023/site/bacs_rrna_called_sites.txt", header=TRUE, sep=" ") 
terc <- terc[terc$chr=="NR_001566.1_TERC",] ###7

count <- site %>% group_by(category) %>% summarize(n = n())
count <- rbind(count, c("TERC", "7"))

  category   n    
  <chr>      <chr>
1 C/D_box    205  
2 Cajal_body 32   
3 H/ACA_box  67   
4 TERC       7  

count$category <- factor(count$category, level=c("C/D_box", "H/ACA_box", "Cajal_body", "TERC"))
pdf("plot/snoRNA_pie.pdf")
ggplot(count, aes(x="", y=as.numeric(n), fill=category))+
geom_bar(width = 1, stat = "identity")+coord_polar("y") + scale_fill_brewer(palette="Set2") +  ###or RdYlBu
  geom_text(aes(label = as.numeric(n)),position = position_stack(vjust = 0.5))+
  theme_void()
dev.off()


#############################stack plot for snoRNA
snosite <- site[,c("chr", "start", "end", "motif", "treat1_level", "treat2_level", "category")]
snosite$mean_level <- (snosite$treat1_level+snosite$treat2_level)/2
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
  <fct>      <fct>    <int>
1 C/D_box    (50,100]    76
2 C/D_box    (20,50]     35
3 C/D_box    (5,20]      94
4 H/ACA_box  (50,100]    22
5 H/ACA_box  (20,50]     12
6 H/ACA_box  (5,20]      33
7 Cajal_body (50,100]     9
8 Cajal_body (20,50]     10
9 Cajal_body (5,20]      13

########create dot plot for three different category ####to be revised

pdf("plot/snoRNA_dot.pdf", width=7,height=6)
ggplot(snosite, aes(x = category, y=mean_level)) +geom_boxplot(fill="white",width=0.6)+geom_jitter(size=3,shape=16,aes(colour = category),width = 0.2)+scale_colour_brewer(palette="Set2")+theme_bw()
dev.off()
   C/D_box  H/ACA_box Cajal_body 
       205         67         32 


############################snoRNA annotation
# Extract the desired part before the second underscore for each value
snosite$refname <- apply(snosite, 1, function(x) {
  input <- x["chr"]
  result <- unlist(strsplit(input, "_"))[1:2]
  paste(result, collapse = "_")
})
write.table(snosite,"site/snosite_namechanged.txt",sep='\t',quote=F,row.names=F)

######the fasta sequence from the database and our reference are not exactly the same. So this is to map the feature to our reference and obtain the corresbonding position
snofa <- read.table("/users/ludwig/ebu571/ebu571/20May2023_final/resource/snoRNA_reference_fullline.bed",header=F, sep="\t")
snofa$refname <- apply(snofa, 1, function(x) {
  input <- x["V1"]
  result <- unlist(strsplit(input, "_"))[1:2]
  paste(result, collapse = "_")
})
snoanno <- read.table("/users/ludwig/ebu571/ebu571/20May2023_final/resource/snoRNA_annptation.txt",header=F) 
colnames(snoanno) <- c("snoid", "refname", "start", "end", "structure", "seq")
snoanno_fa <- merge(snofa, snoanno, by=c("refname"))

for (i in 1:nrow(snoanno_fa)){
  string_infor <- str_locate_all(snoanno_fa[i,"V2"],snoanno_fa[i,"seq"])
  if(nrow(string_infor[[1]])==1){
    snoanno_fa[i,"start1"] <- string_infor[[1]][1,"start"]-1
    snoanno_fa[i,"end1"] <- string_infor[[1]][1,"end"]
  }else if(nrow(string_infor[[1]])>=1){
    oldend <- snoanno_fa[i,"end"]
    mindistance <- abs(string_infor[[1]][,"end"]-oldend)
    cloest_end_row <- which.min(mindistance)
    snoanno_fa[i,"end1"] <- string_infor[[1]][cloest_end_row,"end"]
    snoanno_fa[i,"start1"] <- string_infor[[1]][cloest_end_row,"start"]-1
  }else{
    print(i)
  }
}

#### """"disconcordant references in the following rows
[1] 413
[1] 709
[1] 1412
[1] 1425
[1] 1426
[1] 1430
[1] 1432

write.table(snoanno_fa,"resource/sno_annotation_fasta.txt",sep='\t',quote=F,row.names=F)

####################comparasion to other methods needs to be done
snobin <- read.table("/users/ludwig/ebu571/ebu571/20May2023_final/plot/sno_overlap.txt", header=T, sep="\t")

                   
pdf(file="plot/snosite_bacs_p_venn.pdf",width=4, height=4)
draw.pairwise.venn(area1 = 304,    # Draw pairwise venn diagram
                   area2 = 11,
                   cross.area = 9,category=c("bacs","p_seq"),col=c("#FE767C","#33D982"),fill=c("#FE767C","#33D982"), cat.cex=0,cat.pos=c(-27, 27), cat.dist=c(0.03, 0.03), cat.col=c("#FE767C","#666666"), height = 50 , 
        width = 50, lwd=0.5,ext.text=FALSE)
dev.off()

pdf(file="plot/snosite_bacs_bid_venn.pdf",width=4, height=4)
draw.pairwise.venn(area1 = 304,    # Draw pairwise venn diagram
                   area2 = 39,
                   cross.area = 35,category=c("bacs","bid_seq"),col=c("#FE767C","#53D2EE"),fill=c("#FE767C","#53D2EE"), cat.cex=0,cat.pos=c(-27, 27), cat.dist=c(0.03, 0.03), cat.col=c("#FE767C","#666666"), height = 50 , 
        width = 50, lwd=0.5, ext.text=FALSE)
dev.off()


####metagene plot for snoRNA
site <- read.table("site/snosite_namechanged.txt",header=T, sep="\t") ###304

#################after removing some isoforms
######For C/D box
snosite1 <- read.table("site/snosite_namechanged_modified.txt", header=T, sep="\t")
fai <- read.table("/users/ludwig/ebu571/ebu571/20May2023_final/resource/snoRNA_reference.fasta.fai", header=F, sep="\t")
nrow(snosite1) ##219
mergesite1 <- merge(snosite1[,c("chr", "start", "end","category")], fai[,c("V1","V2")], by.x=c("chr"), by.y=c("V1"), all.x=TRUE)
mergesite1$relative_pos <- mergesite1$end/mergesite1$V2
pdf("plot/snoRNA_metagene_plot_CD_box_rmisoform.pdf",width=5,height=2)
p1<- ggplot(mergesite1[mergesite1$category=="C/D_box",], aes(x=relative_pos))+##,y = ..scaled..
  geom_density(size=0.6)+theme_bw() +theme(axis.text.x=element_blank(),
      axis.ticks.x=element_blank(),
      axis.text.y=element_blank(),axis.title.x = element_blank(),axis.title.y = element_blank(),
      axis.ticks.y=element_blank(),panel.background = element_rect(colour = "black", size=0.5,fill=NA)) +scale_x_continuous(limits = c(0,1), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0,2.5),expand = c(0, 0))
dev.off()
###with legend
p1<- ggplot(mergesite1[mergesite1$category=="C/D_box",], aes(x=relative_pos))+##,y = ..scaled..
  geom_density(size=1)+theme_bw() +scale_x_continuous(limits = c(0,1), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0,2.5),expand = c(0, 0))

pdf("plot/snoRNA_metagene_plot_HA_box_rmisoform.pdf",width=5,height=2)
p2 <- ggplot(mergesite1[mergesite1$category=="H/ACA_box",], aes(x=relative_pos)) + ##,y = ..scaled..
  geom_density(size=0.6)+theme_bw() +theme(axis.text.x=element_blank(),
      axis.ticks.x=element_blank(),
      axis.text.y=element_blank(),axis.title.x = element_blank(),axis.title.y = element_blank(),
      axis.ticks.y=element_blank(),panel.background = element_rect(colour = "black", size=0.5,fill=NA))+scale_x_continuous(limits = c(0,1), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0,2.5),expand = c(0, 0))
dev.off()

##with legend
p2 <- ggplot(mergesite1[mergesite1$category=="H/ACA_box",], aes(x=relative_pos)) + ##,y = ..scaled..
  geom_density()+theme_bw()+scale_x_continuous(limits = c(0,1), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0,2.5),expand = c(0, 0))

pdf("plot/snoRNA_metagene_plot_full_rmisoform.pdf")
p3<- ggplot(mergesite1, aes(x=relative_pos,y = ..scaled..))+
  geom_density()+theme_bw() +theme(axis.text.x=element_blank(),
      axis.ticks.x=element_blank(),
      axis.text.y=element_blank(),
      axis.ticks.y=element_blank()) +scale_x_continuous(limits = c(0,1), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0,1.02), expand = c(0, 0))
dev.off()

pdf("plot/snoRNA_metagene_plot_Cajal_body.pdf")
p4 <- ggplot(mergesite1[mergesite1$category=="Cajal_body",], aes(x=relative_pos,y = ..scaled..)) + 
  geom_density()+theme_bw() +theme(axis.text.x=element_blank(),
      axis.ticks.x=element_blank(),
      axis.text.y=element_blank(),
      axis.ticks.y=element_blank()) +scale_x_continuous(limits = c(0,1), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0,1.02), expand = c(0, 0))
dev.off()

###############plot C/D box, HA box heatmap
#obtain refnames for sites after removing isoforms
rmso_site <- merge(snosite1[,c("chr","start","end","motif","category")], site[,c("chr","start","end","motif","category", "refname")], by=c("chr","start","end","motif","category")) ##197
CD_box_refname <- rmso_site[rmso_site$category=="C/D_box",] %>% dplyr::select("refname","chr") %>% unique() ##70
##read snobd annotation 
fai <- read.table("/users/ludwig/ebu571/ebu571/20May2023_final/resource/snoRNA_reference.fasta.fai", header=F, sep="\t")
snoanno <- read.table("/users/ludwig/ebu571/ebu571/20May2023_final/resource/snoRNA_annptation.txt",header=F) 
colnames(snoanno) <- c("snoid", "refname", "start", "end", "structure", "seq")
CD_box_list <- merge(CD_box_refname, snoanno, by.x=c("refname"), by.y=c("refname"))
CD_box_list <- merge(CD_box_list[,c("refname", "chr", "start", "end", "structure")], fai[,c("V1","V2")], by.x=c("chr"), by.y=c("V1"), all.x=TRUE)
nrow(CD_box_list) ##291
##remove reference that have much longer reference than ours
CD_box_list <- CD_box_list[CD_box_list$end<=CD_box_list$V2,] ##288
###transform region to position
CD_box_list_transformed <- CD_box_list %>%
  mutate(position = map2(start, end, seq)) %>%
  unnest(position) %>%
  select(-start, -end)

CD_box_list_transformed$relative_pos <- CD_box_list_transformed$position/CD_box_list_transformed$V2

seldat <- CD_box_list_transformed[CD_box_list_transformed$structure=="C_box" |CD_box_list_transformed$structure=="C'_box" |CD_box_list_transformed$structure=="D'_box"|CD_box_list_transformed$structure=="D_box", ] %>% dplyr::select(structure, relative_pos)

C_count <- as.data.frame(table(cut(seldat[seldat$structure=="C_box" | seldat$structure=="C'_box" ,]$relative_pos, breaks=seq(0,1, by=0.01)), dnn="Range"))
D_count <- as.data.frame(table(cut(seldat[seldat$structure=="D'_box" |seldat$structure== "D_box",]$relative_pos, breaks=seq(0,1, by=0.01)), dnn="Range"))

merge_CD <- merge(C_count,D_count, by=c("Range"))
colnames(merge_CD) <- c("range","Cbox", "Dbox")
merge_CD$Cbox <- -merge_CD$Cbox
melt_CD <- melt(merge_CD, id.vars=c("range"))
melt_CD$variable <- factor(melt_CD$variable, level=c("Dbox","Cbox"))
p5 <- ggplot(melt_CD, aes(x=range,y=variable,fill=value)) + geom_tile(size=0.01)+##to make Square tiles
  scale_fill_gradient2(low = "#00BA38", high = "#F8766D")+theme(legend.position = "none")+theme(axis.text.x=element_blank(),
      axis.ticks.x=element_blank(),
      axis.text.y=element_blank(),
      axis.ticks.y=element_blank(),axis.title.x = element_blank(),axis.title.y = element_blank(),
      panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
      panel.background = element_rect(colour = "black", size=0.5,fill=NA))


#with legend
p5 <- ggplot(melt_CD, aes(x=range,y=variable,fill=value)) + geom_tile(size=0.01)+##to make Square tiles
  scale_fill_gradient2(low = "#00BA38", high = "#F8766D")+theme_bw()

pdf("plot/snoRNA_CD_box_heatmap.pdf",width=5,height=1)
dev.off()
###############similarly, plot H box and ACA box
HA_box_refname <- rmso_site[rmso_site$category=="H/ACA_box",] %>% dplyr::select("refname","chr") %>% unique() ##45
HA_box_list <- merge(HA_box_refname, snoanno, by.x=c("refname"), by.y=c("refname"))
HA_box_list <- merge(HA_box_list[,c("refname", "chr", "start", "end", "structure")], fai[,c("V1","V2")], by.x=c("chr"), by.y=c("V1"), all.x=TRUE)
nrow(HA_box_list) ##162
##remove reference that have much longer reference than ours
HA_box_list <- HA_box_list[HA_box_list$end<=HA_box_list$V2,] ##162
###transform region to position
HA_box_list_transformed <- HA_box_list %>%
  mutate(position = map2(start, end, seq)) %>%
  unnest(position) %>%
  select(-start, -end)

HA_box_list_transformed$relative_pos <- HA_box_list_transformed$position/HA_box_list_transformed$V2

seldat_HA <- HA_box_list_transformed[HA_box_list_transformed$structure=="H_box" |HA_box_list_transformed$structure=="ACA_box", ] %>% dplyr::select(structure, relative_pos)

H_count <- as.data.frame(table(cut(seldat_HA[seldat_HA$structure=="H_box",]$relative_pos, breaks=seq(0,1, by=0.01)), dnn="Range"))
A_count <- as.data.frame(table(cut(seldat_HA[seldat_HA$structure=="ACA_box",]$relative_pos, breaks=seq(0,1, by=0.01)), dnn="Range"))
merge_HA <- merge(H_count, A_count, by=c("Range"))
colnames(merge_HA) <- c("range","Hbox", "Abox")
merge_HA$Hbox <- -merge_HA$Hbox
melt_HA <- melt(merge_HA, id.vars=c("range"))
melt_HA$variable <- factor(melt_HA$variable, level=c("Abox", "Hbox"))
pdf("plot/snoRNA_HA_box_heatmap.pdf",width=5,height=1)
p6 <- ggplot(melt_HA, aes(x=range,y=variable,fill=value)) + geom_tile(size=0.01)+##to make Square tiles
  scale_fill_gradient2(low = "#619CFF", high = "#ED8141")+theme(legend.position = "none")+theme(axis.text.x=element_blank(),
      axis.ticks.x=element_blank(),
      axis.text.y=element_blank(),
      axis.ticks.y=element_blank(),axis.title.x = element_blank(),axis.title.y = element_blank(),
      panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
      panel.background = element_rect(colour = "black", size=0.5,fill=NA))

dev.off()

###with legend
p6 <- ggplot(melt_HA, aes(x=range,y=variable,fill=value)) + geom_tile(size=0.01)+##to make Square tiles
  scale_fill_gradient2(low = "#619CFF", high = "#ED8141")+theme_bw()


pdf("plot/snosite_metagene_plot.pdf",width=10,height=3)
ggarrange(p1,p2,p5,p6,labels = c("A", "B", "C","D"),ncol = 2, nrow = 2,heights=c(2.5,0.5),widths=c(6.25,6.25))
dev.off()

###with legend 
pdf("plot/snosite_metagene_plot_legend.pdf",width=25,height=8)
ggarrange(p1,p2,p5,p6,labels = c("A", "B", "C","D"),ncol = 2, nrow = 2,heights=c(2.5,0.5),widths=c(6.25,6.25))
dev.off()

