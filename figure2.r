##module load R/4.0.3-foss-2020b
library(tidyverse)
library(dplyr)
library(reshape2)
library(tidyr)
library(ggplot2)
library(ggpubr)
library(cowplot)
library(RColorBrewer)
library(tidyverse)
setwd("/users/ludwig/ebu571/ebu571/bacs")

##data source: for rRNA sites: /users/ludwig/ebu571/ebu571/07Nov2023/site/bacs_rrna_called_sites.txt  and /users/ludwig/ebu571/ebu571/07Nov2023/site/bacs_rrna_all_sites check code: /users/ludwig/ebu571/ebu571/07Nov2023/code/process.r
##for snRNA: /users/ludwig/ebu571/ebu571/27Oct2023/site/bacs_rrna_called_sites.txt and /users/ludwig/ebu571/ebu571/27Oct2023/site/bacs_rrna_all_sites.txt 


####a pie chart for major and minor splicesome RNAs
###a stacked plot for 18S, 28S and 5.8S rRNA

mergedat <- read.table("/users/ludwig/ebu571/ebu571/27Oct2023/site/bacs_rrna_called_sites.txt", header=T, sep=" ", colClasses="character") 
mergedat <- mergedat[mergedat$chr!="NR_003285.3_RNA5_8SN5" & mergedat$chr!="NR_003286.4_RNA18SN5" & mergedat$chr!="NR_003287.4_RNA28SN5" & mergedat$chr!="2kb_spikein"& mergedat$chr!="MT_RNR1_ENST00000389680.2_ENSE00001544499" & mergedat$chr!="MT_RNR2_ENST00000387347.2_ENSE00001544497" & mergedat$chr!="spike_in",]

dat <- mergedat %>% group_by(chr) %>% summarize(n = n()) 
dat1 <-dat[dat$chr!="NR_001566.1_TERC" & dat$chr!="NR_023317.1_RNU7_1" & dat$chr!="NR_001445.2_RN7SK" & dat$chr!="NR_003051.3_RMRP",]
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

#####specific for Human snRNA, including C666.1, raji, elijah
all_sn <- read.table("site/bacs_snrna_allcelllines_merge_sites.txt", header=T, sep=" ", colClasses="character") 
all_sn <- all_sn[!grepl("MT",all_sn$chr) & !grepl("TERC",all_sn$chr),]
sn_chrcount <- all_sn %>% group_by(chr) %>% summarize(n=n())

   chr                      n
   <chr>                <int>
 1 NR_001445.2_RN7SK        2
 2 NR_002715.1_RN7SL1       1
 3 NR_002716.3_RNU2_1      14
 4 NR_002756.2_RNU5A_1      4
 5 NR_003925.1_RNU4_1       4
 6 NR_004394.1_RNU6_1       4
 7 NR_004430.2_RNU1_1       2
 8 NR_023343.1_RNU4ATAC     3
 9 NR_023344.1_RNU6ATAC     2
10 NR_029422.2RNU12         2

#################Venn gram for comparison with MS data on major spliceosomal RNA
library(VennDiagram)##########
pdf("plot/bacs_ms_venn_on_spliceosomal.pdf")
draw.pairwise.venn(area1=28 , area2=28,cross.area=28, 
                   category=c("bacs","SILNAS MS"),col=c("#FE767C","#666666"), cat.pos=c(-27, 27), cat.dist=c(0.03, 0.03), cat.col=c("#FE767C","#666666"), height = 300 , 
        width = 300)
dev.off()

################plot conversion/input on U2
allsite <- read.table("/users/ludwig/ebu571/ebu571/27Oct2023/site/bacs_rrna_all_sites.txt", header=T, sep=" ", colClasses="character") 
u2allsite <- allsite[allsite$chr=="NR_002716.3_RNU2_1",] %>% dplyr::select(chr,start,end,motif,treat1_conversion,ctrl1_conversion,treat2_conversion,ctrl2_conversion) %>% mutate_at(c('start','end','treat1_conversion', 'ctrl1_conversion', 'treat2_conversion','ctrl2_conversion'), as.numeric)

u2allsite$treat_conversion <- (u2allsite$treat1_conversion+u2allsite$treat2_conversion)/2
u2allsite$ctrl_conversion <- (u2allsite$ctrl1_conversion+u2allsite$ctrl2_conversion)/2

dat <- u2allsite %>% dplyr::select(end,ctrl_conversion, treat_conversion)

meldat <- melt(dat, id.var="end")
pdf("plot/U2_plot.pdf", width=12,height = 4)
ggplot(meldat, aes(group=variable, col=variable))+geom_segment(aes(x=end, y=round(value,4)*100,xend=end,yend=0), size=0.75)+geom_point(aes(x=end, y=round(value,4)*100),size=4) +ylab("conversion ratio(%)")+xlab("")+
theme_bw()+scale_x_continuous(limits = c(0,200), expand = c(0, 0)) + scale_y_continuous(limits = c(-2,102), expand = c(0, 0))+scale_color_manual(values=c('#ABABAB','#FE767C'))
dev.off()


################plot conversion/input on TERC and U4
selsite <- allsite[allsite$chr=="NR_001566.1_TERC" |allsite$chr=="NR_023343.1_RNU4ATAC" ,] %>% dplyr::select(chr,start,end,motif,treat1_conversion,ctrl1_conversion,treat2_conversion,ctrl2_conversion) %>% mutate_at(c('start','end','treat1_conversion', 'ctrl1_conversion', 'treat2_conversion','ctrl2_conversion'), as.numeric)

selsite$treat_conversion <- (selsite$treat1_conversion+selsite$treat2_conversion)/2
selsite$ctrl_conversion <- (selsite$ctrl1_conversion+selsite$ctrl2_conversion)/2

terc <- selsite[selsite$chr=="NR_001566.1_TERC",] %>% dplyr::select(end,ctrl_conversion, treat_conversion)

u4atac <- selsite[selsite$chr=="NR_023343.1_RNU4ATAC",] %>% dplyr::select(end,ctrl_conversion, treat_conversion)

melu4atac <- melt(u4atac, id.var="end")
pdf("plot/U4ATAC_plot.pdf", width=12,height = 4)
ggplot(melu4atac, aes(group=variable, col=variable))+geom_segment(aes(x=end, y=round(value,4)*100,xend=end,yend=0), size=0.65)+geom_point(aes(x=end, y=round(value,4)*100),size=4) +ylab("conversion ratio(%)")+xlab("")+
theme_bw() +scale_x_continuous(limits = c(0,150), expand = c(0, 0))+ scale_y_continuous(limits = c(-2,102), expand = c(0, 0))+scale_color_manual(values=c('#ABABAB','#FE767C'))
dev.off()


##############Plot treat_level for TERC
tercsite <- allsite[allsite$chr=="NR_001566.1_TERC",] %>% dplyr::select(chr,start,end,motif,treat1_level,treat2_level) %>% mutate_at(c('start','end','treat1_level', 'treat2_level'), as.numeric)
tercsite$level <- (tercsite$treat1_level+tercsite$treat2_level)/2
terc <- tercsite[tercsite$chr=="NR_001566.1_TERC",] %>% dplyr::select(end,level)


melterc <- melt(terc, id.var="end")
melterc[melterc$value<0,]$value <- 0
pdf("plot/TERC_plot.pdf", width=18,height = 4) 
ggplot(melterc, aes(group=variable, col=variable))+geom_segment(aes(x=end, y=round(value,4)*100,xend=end,yend=0), size=0.65)+geom_point(aes(x=end, y=round(value,4)*100),size=4) +ylab("ratio(%)")+xlab("")+geom_hline(yintercept=5, linetype="dashed", color = "black")+
theme_bw() +scale_x_continuous(limits = c(0,500), expand = c(0, 0))+ scale_y_continuous(limits = c(-3,103), expand = c(0, 0))+scale_color_manual(values=c('#77A9D7'))
dev.off()



###################plot for rRNA#################################
rrna_allsite <- read.table("/users/ludwig/ebu571/ebu571/07Nov2023/site/bacs_rrna_all_sites.txt", header=TRUE, sep=" ") 
rrna_allsite <- rrna_allsite[rrna_allsite$chr=="NR_003285.3_RNA5_8SN5" | rrna_allsite$chr=="NR_003286.4_RNA18SN5" | rrna_allsite$chr=="NR_003287.4_RNA28SN5",]
rrna_allsite$ctrl_conversion <- (rrna_allsite$ctrl1_conversion+rrna_allsite$ctrl2_conversion)/2
rrna_allsite$treat_conversion <- (rrna_allsite$treat1_conversion+rrna_allsite$treat2_conversion)/2
seldat <- rrna_allsite[,c("chr", "start", "end", "ctrl_conversion", "treat_conversion")]
#####################generate pseudouridine and other modification list from /users/ludwig/ebu571/ebu571/20May2023_final/code/figure1_2.r
##cat align/rRNA_called_sites.txt | sed 's/ /\t/g' | awk '$1~"18S" || $1 ~"28S" || $1~"5_8S"' | awk '{OFS="\t"}{print $1, $2, $3,"pseu"}' > align/pseudouridine_list.txt
####add type information
###known pseudoU
###U
###other modifications
###

type <- read.table("/users/ludwig/ebu571/ebu571/20May2023_final/align/pseudouridine_list.txt", header=F, sep="\t") 
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

###########segment plot
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

####################PIE chart
df <- data.frame(
  group = c("5.8S", "18S", "28S", "16S mt-rRNA", "12S mt-rRNA"),
  value = c(2, 40, 62, 1,8)
  )

pdf("plot/rRNA_pie.pdf")
ggplot(df, aes(x="", y=value, fill=group))+
geom_bar(width = 1, stat = "identity")+coord_polar("y") + scale_fill_brewer(palette="Set2") +  ###or RdYlBu
  geom_text(aes(label = value),position = position_stack(vjust = 0.5))+
  theme_void()
dev.off()


#####raw signal correlation
pdf("plot/rRNA_raw_signal_correlation.pdf")
smoothScatter(rrna_allsite$treat1_conversion*100, rrna_allsite$treat2_conversion*100, xlab = "rep1 raw signal for rRNA", ylab = "rep2 raw signal for rRNA")
dev.off()
cor(rrna_allsite$treat1_conversion*100, rrna_allsite$treat2_conversion*100) ###0.9999867

#################Venn gram for comparison with MS data
library(VennDiagram)##########to be revised, include pos 36 in 18S rRNA. the number in this plot and the pie chart needs to be revised
pdf("plot/bacs_rrna_ms_venn.pdf")
draw.pairwise.venn(area1=110 , area2=105,cross.area=105, 
                   category=c("bacs","SILNAS MS"),col=c("#FE767C","#666666"), cat.pos=c(-27, 27), cat.dist=c(0.03, 0.03), cat.col=c("#FE767C","#666666"), height = 300 , 
        width = 300)
dev.off()

######################################################plot for comparision of rRNA between cytosolic and rRNA
rrna_site <- read.table("/users/ludwig/ebu571/ebu571/07Nov2023/site/bacs_rrna_called_sites.txt", header=TRUE, sep=" ") 
site <- rrna_site[rrna_site$chr=="NR_003285.3_RNA5_8SN5" | rrna_site$chr=="NR_003286.4_RNA18SN5" | rrna_site$chr=="NR_003287.4_RNA28SN5" |  rrna_site$chr=="MT_RNR1_ENST00000389680.2_ENSE00001544499" | rrna_site$chr=="MT_RNR2_ENST00000387347.2_ENSE00001544497",]

site$mean_level <- (as.numeric(site$treat2_level)+as.numeric(site$treat1_level))/2

site$name <- gsub("MT_RNR1_ENST00000389680.2_ENSE00001544499", "mt_rRNA",site$chr)
site$name <- gsub("MT_RNR2_ENST00000387347.2_ENSE00001544497", "mt_rRNA",site$name)
site$name <- gsub("NR_003285.3_RNA5_8SN5", "rRNA",site$name)
site$name <- gsub("NR_003286.4_RNA18SN5", "rRNA",site$name)
site$name <- gsub("NR_003287.4_RNA28SN5", "rRNA",site$name)

selsite <- site[,c("name", "mean_level")]
selsite$name <- factor(selsite$name, level=c("rRNA", "mt_rRNA"))
selsite[selsite$mean_level>=1,]$mean_level <- 1

p <- ggplot(selsite,aes(x=name,y=mean_level*100, fill=name))+ geom_boxplot(width=0.6, outlier.size=0.5)+ #+ geom_jitter(size=2,shape=16, colour = "grey",position=position_jitter(0.2))
  #stat_summary(fun = mean, geom = "point",color = "black") +
  theme_light() + xlab("") +  theme(axis.text = element_text(color = "black")) + ylim(0,100) +
  ylab("level") +
  scale_fill_manual(values=c("#006699","#F4A582"))

ggsave("plot/modification_level_cytorRNA_vsmtrRNA.pdf",p, width=3, height = 3)

mean(selsite[selsite$name=="rRNA",]$mean_level) ##0.832353
mean(selsite[selsite$name=="mt_rRNA",]$mean_level) #0.2240638

mt_rRNA    rRNA 
      9     104 