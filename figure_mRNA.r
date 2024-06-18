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

##read ##annotated site 1464. 18 sites are on intergenic region
annotated_site <- read.table("site/treat.unique.annotation.bed",  header=TRUE) ##1337
site <- read.table("site/bacs_mRNA_site.txt",  header=TRUE) ###1376
write.table(merge(site[site$chr!="chrM",],annotated_site, by=c("chr", "start", "end", "strand")),"site/bacs.mRNA.annotated.txt",sep='\t',quote=F,row.names=F)
inputmis <- read.table("site/treat.unique.annotation.mismatch.bed", header=F, sep="\t") ##1337
test <- merge(site[site$chr!="chrM",],inputmis, by.x=c("chr", "start", "end", "strand"),by.y=c("V1", "V2", "V3", "V4"))
test$mis_ratio <- test$V40/test$ctrl_depth
max(test$mis_ratio) ###0.04081633
##check mismatch criteria

table(annotated_site$transcript_type, annotated_site$transcript_detail)
                                   3UTR 5UTR CDS exon
  lncRNA                              0    0   0   16
  nonsense_mediated_decay             7    1   1    0
  processed_pseudogene                0    0   0   16
  processed_transcript                0    0   0    2
  protein_coding                    409   68 808    0
  protein_coding_CDS_not_defined      0    0   0    1
  retained_intron                     0    0   0    4
  snRNA                               0    0   0    2
  transcribed_processed_pseudogene    0    0   0    2

filtered_site <- annotated_site[annotated_site$strand_infor=="proper_strand" &annotated_site$transcript_type !="snRNA" &annotated_site$transcript_detail!="intronic", ]
nrow(filtered_site) ##1335




table(filtered_site$transcript_detail,filtered_site$strand_infor)

       proper_strand
  3UTR           416
  5UTR            69
  CDS            809
  exon            41

###Finally, sites are divided into 4 catogaries. 5UTR,3UTR,CDS, ncRNA(exon)
feature_count <- data.frame(table(filtered_site$transcript_detail,filtered_site$strand_infor))

feature_count$Var1 <- factor(feature_count$Var1, level=c("5UTR", "CDS", "3UTR", "exon"))
feature_count$Freq <- as.numeric(feature_count$Freq)

###corrected stop codon information.
##there are two sites in stop codon(UAG). Regard it as CDS.
feature_count[feature_count$Var1=="3UTR",]$Freq <- feature_count[feature_count$Var1=="3UTR",]$Freq -2
feature_count[feature_count$Var1=="CDS",]$Freq <- feature_count[feature_count$Var1=="CDS",]$Freq +2

  Var1          Var2 Freq
1 3UTR proper_strand  414
2 5UTR proper_strand   69
3  CDS proper_strand  811
4 exon proper_strand   41

###Piechart
nb.cols <- 4
mycolors <- colorRampPalette(brewer.pal(4, "Set2"))(nb.cols)

pdf("plot/mRNA.feature.distribution_pie.pdf")
ggplot(feature_count, aes(x="", y=Freq, fill=Var1))+
geom_bar(width = 1, stat = "identity")+coord_polar("y") + scale_fill_manual(values = mycolors) +  ###or RdYlBu
  geom_text(aes(label = Freq),position = position_stack(vjust = 0.5))+
  theme_void()
dev.off()


############
####Read sites
site <- read.table("site/bacs_mRNA_site.txt",  header=TRUE) ###1378
chrsite <- merge(filtered_site[,c("chr","start","end","transcript_id","transcript_name","transcript_detail")],site[site$chr!="chrM",], by=c("chr", "start", "end")) ##1335
##remove improper strand; remove intergenic; remove intron
###Pie chart for different levels (except chrM)
chrsite$mean_level <- chrsite$treat_level2
chrsite[chrsite$mean_level>1,]$mean_level <- 1
chrsite$level_group <- cut(chrsite$mean_level,breaks=c(0,0.2,0.5,1))
level_count <- chrsite %>% group_by(level_group) %>% summarize(n = n())
level_count$level_group <- factor(level_count$level_group,level=c("(0.2,0.5]", "(0,0.2]","(0.5,1]"))
  level_group     n
  <fct>       <int>
1 (0,0.2]      1097
2 (0.2,0.5]     162
3 (0.5,1]        76

pdf("plot/mRNA_level_pie.pdf")
ggplot(level_count, aes(x="", y=as.numeric(n), fill=level_group))+
geom_bar(width = 1, stat = "identity")+coord_polar("y") + scale_fill_manual(values=c('#fdb462', '#8bc53f','#f03e3e')) +  ###or RdYlBu
  geom_text(aes(label = as.numeric(n)),position = position_stack(vjust = 0.5))+
  theme_void()
dev.off()

###read inosine
inosine <- read.table("inosine/input_treat.mRNA_inosine_site_deltaG.txt",  header=TRUE) ##1338
nrow(inosine[inosine$snp!="SNP",]) ##1121
inosine_site <- inosine[inosine$snp!="SNP",]
median(inosine_site$A2G_ratio_ctrl) ##0.20
mean(inosine_site$A2G_ratio_ctrl) ##0.2351807
inosine_site$level_group <- cut(inosine_site$A2G_ratio_ctrl,breaks=c(0,0.2,0.5,1))
level_count2 <- inosine_site %>% group_by(level_group) %>% summarize(n = n())
  level_group     n
  <fct>       <int>
1 (0,0.2]       587
2 (0.2,0.5]     451
3 (0.5,1]        83
inosine_site$category <- "inosine"
chrsite$category <- "psi"

level_dat <- rbind(chrsite[,c("level_group", "category")], inosine_site[,c("level_group", "category")])
level_dat$level_group <- factor(level_dat$level_group, level=c("(0.5,1]","(0.2,0.5]", "(0,0.2]"))
level_dat$category <- factor(level_dat$category, level=c("psi","inosine"))
pdf("plot/mRNA_level_stack_plot.pdf",width=2.5,height=4)
ggplot(level_dat, aes(x = category, fill = level_group))+geom_bar(width=0.6) +theme_bw() +scale_fill_manual(values=c('#f03e3e','#fdb462', '#8bc53f'))+theme(legend.position = "none")+theme(axis.title.x = element_blank(),axis.title.y = element_blank())+
geom_text(stat = "count", aes(label = after_stat(count)),position = position_stack(vjust = 0.5))+ylim(0,1500)
dev.off()

##
chrsite$category <- "Hela"
pdf("plot/mRNA_violin.pdf",width=2,height=4)
ggplot(chrsite, aes(x = category, y=mean_level*100),fill=category) + geom_violin(trim=TRUE, fill='#1B9E77')+stat_summary(fun.y=mean, geom="point", size=2, color="red")+theme_bw()+
  stat_summary(fun = mean, geom="text", aes(label = round(..y.., 3)), hjust = -0.2, vjust=-1.5)+scale_y_continuous(limits = c(0,100),expand = c(0, 1))
dev.off()

######################################### level distribution
inosine_site$mean_level <- inosine_site$A2G_ratio_ctrl
level_dat2 <- rbind(chrsite[,c("mean_level", "category")], inosine_site[,c("mean_level", "category")])
level_dat2$category <- factor(level_dat2$category, level=c("inosine","Hela"))
pdf("plot/psi_inosine_level_distribution_histogram.pdf",width=6,height=4)
ggplot(level_dat2, aes(x=mean_level*100, fill=category)) +  geom_histogram(alpha=0.6,binwidth=2,position = "identity")+scale_fill_manual(values=c("#404080","#1b9e77")) + theme_bw() +theme(axis.title.x = element_blank(),axis.title.y = element_blank())+scale_y_continuous(limits=c(0, 405)) + #scale_x_continuous(limits=c(0, 105),breaks=c(0,20,40,60,80,105))
labs(fill="")
dev.off()


############motif frequency 
##for this plot. use module load R/4.3.1-foss-2022b instead
library(ggplot2)
library(ggbreak) 
library(patchwork)

#########plot motif only on mRNA
mRNA_site <- chrsite[chrsite$transcript_detail!="exon",] ##1294
motifdat <- aggregate(mean_level ~ motif, mRNA_site, mean) ##222
motif_count <- mRNA_site %>% group_by(motif) %>% summarize(n=n())
motif <- merge(motifdat, motif_count, by=c("motif"))
motif <- motif[order(motif$n, decreasing=TRUE),]

motif$group <- stri_replace_all_regex(motif$motif,
                                  pattern=c('TGTAG', 'TCTAG', 'GTTCA',"GTTCC","GTTCG","GTTCT","CTTTG","ACTTT","CATTT","GTTTG","TTTTT","GTTAA"),
                                  replacement=c('USUAG', 'USUAG', 'GUUCN','GUUCN','GUUCN','GUUCN',"U_rich","U_rich","U_rich","U_rich","U_rich","GUUAA"),
                                  vectorize=FALSE)

motif[motif$group!="USUAG" & motif$group!="GUUCN" & motif$group !="U_rich" &motif$group!="GUUAA" &motif$group!="AATTG",]$group <- "others"
motif$group <- factor(motif$group, level=c("USUAG", "GUUCN", "GUUAA","U_rich", "AATTG","others"))
colnames(motif) <- c("motif", "mean_level", "n", "group")
write.table(motif,"plot/mRNA_motif_count.txt",sep='\t',quote=F,row.names=F)

p1 <- ggplot(motif, aes(x = n, y=mean_level*100, fill = group)) + geom_point(aes(color=group), size=2.5)+theme(legend.position = "none") +scale_color_manual(values=c("#CB0505","#1B9E77", "#984EA3","#E78AC3","#FC8D62","#808080"))+
geom_text(aes(label=ifelse(group!="others",as.character(motif),'')))+
scale_y_continuous(limits=c(0,60),expand = c(0, 0))+scale_x_continuous(limits = c(0, 40), breaks = seq(0,40,5),expand = c(0, 0))+ylab("average modification level")+xlab("motif frequency")+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
panel.background = element_blank(), axis.line = element_line(colour = "black")) +theme(legend.position = "none")

p2 <- ggplot(motif, aes(x = n, y=mean_level*100, fill = group)) + geom_point(aes(color=group), size=2.5)+theme(legend.position = "none") +scale_color_manual(values=c("#CB0505","#1B9E77", "#984EA3","#E78AC3","#FC8D62","#808080"))+
geom_text(aes(label=ifelse(group!="others",as.character(motif),'')))+
ylim(0,60)+scale_x_continuous(limits = c(40, 100), breaks = seq(40,100,30),expand = c(0, 0))+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
panel.background = element_blank(), axis.text.y=element_blank(),axis.ticks.y=element_blank(),axis.title.y = element_blank(),axis.line.x = element_line(colour = "black")) +theme(legend.position = "none")

pdf("plot/mRNA_motif_frequency_nolegend.pdf",width=5,height=3)
ggarrange(p1,p2,ncol = 2, nrow = 1,widths=c(3,1))
dev.off()


########violin plot for top 20 motifs
top <- merge(motif[motif$group!="others",],mRNA_site,by.x=c("motif"), by.y=c("motif")) ##401
top$motif <- factor(top$motif, level=c('TGTAG', 'TCTAG', 'GTTCA',"GTTCC","GTTCG","GTTCT","GTTAA","AATTG","ACTTT","CATTT","CTTTG","GTTTG","TTTTT"))
top[top$mean_level.y>1,]$mean_level.y <- 1
pdf("plot/mRNA_motif_boxplot.pdf",width=7, height=3)
ggplot(top, aes(x = motif, y=mean_level.y*100, fill = group)) + geom_boxplot(outlier.shape=NA,width=0.5, lwd=0.4) +scale_fill_manual(values=c("#CB0505","#1B9E77", "#984EA3","#FC8D62","#E78AC3"))+
stat_summary(fun = mean, geom="text", aes(label = round(..y.., 3)), hjust = -0.2, vjust=-1.5)+theme_bw()+theme(legend.position = "none")
dev.off()

pdf("plot/mRNA_motif_jitter.pdf",width=7, height=3)
ggplot(top, aes(x = motif, y=mean_level.y*100, color = group)) + geom_jitter(size=2,shape=16,alpha = 0.8,position=position_jitter(0.2))+scale_color_manual(values=c("#CB0505","#1B9E77", "#984EA3","#E78AC3","#FC8D62"))+
 stat_summary(fun = mean,geom = "crossbar", width = 0.3,colour = "gray40")+theme_bw()+theme(legend.position = "none")
dev.off()
#geom_boxplot(width=0.1,outlier.shape=NA,fill="white")

#######T, TT and TTT distribution
mRNA_site$numberU <- gsub("[A|C|G|T][A|C|G]T[A|C|G][A|C|G|T]","single",mRNA_site$motif)
mRNA_site$numberU <- gsub("[A|C|G]TT[A|C|G][A|C|G|T]","double",mRNA_site$numberU)
mRNA_site$numberU <- gsub("[A|C|G|T][A|C|G]TT[A|C|G]","double",mRNA_site$numberU)
mRNA_site$numberU <- gsub("[A|C|G|T][A|C|G|T]TTT","multiple",mRNA_site$numberU)
mRNA_site$numberU <- gsub("[A|C|G|T]TTT[A|C|G|T]","multiple",mRNA_site$numberU)
mRNA_site$numberU <- gsub("TTT[A|C|G|T][A|C|G|T]","multiple",mRNA_site$numberU)
countU <- mRNA_site %>% group_by(numberU) %>% summarize(n=n())
  numberU      n
  <chr>    <int>
1 double     452
2 multiple   264
3 single     578

nrow(mRNA_site[grepl("TT", mRNA_site$motif) & !grepl("TTT",mRNA_site$motif),]) # 452
nrow(mRNA_site[grepl("TTT",mRNA_site$motif),]) ## 264
nrow(mRNA_site[!grepl("TT", mRNA_site$motif) & !grepl("TTT",mRNA_site$motif),]) # 578
countU$numberU <- factor(countU$numberU,level=c("single", "double", "multiple"))
pdf("plot/mRNA_consecutiveU_pie.pdf",width=4,height=4)
ggplot(countU, aes(x="", y=as.numeric(n), fill=numberU))+
geom_bar(width = 1, stat = "identity")+coord_polar("y") + scale_fill_manual(values=c("#FEE8C8","#FDBB84","#EF6548")) +  ###or RdYlBu
  geom_text(aes(label = as.numeric(n)),position = position_stack(vjust = 0.5))+
  theme_void()
dev.off()

##violin plot for consecutive U
mRNA_site$numberU <- factor(mRNA_site$numberU, level=c("single", "double", "multiple"))
pdf("plot/mRNA_consecutiveU_violin.pdf",width=4, height=4)
ggplot(mRNA_site, aes(x = numberU, y=mean_level*100, fill = numberU)) + geom_violin(trim=TRUE)+theme(legend.position = "none") + geom_boxplot(width=0.1,outlier.shape=NA,fill="white")+#scale_fill_mannual(values=c(""))
stat_summary(fun = median, geom="text", aes(label = round(..y.., 3)), hjust = -0.2, vjust=-1.5)+scale_fill_manual(values=c("#FEE8C8","#FDBB84","#EF6548"))
dev.off()


#####plot for codon usage
site_codon <- merge(filtered_site,site[site$chr!="chrM",], by=c("chr", "start", "end")) ##1336

site_codon <- site_codon[site_codon$codon!="none",] ###811
count_codon <- site_codon %>% group_by(codon) %>% summarize(n=n())
count_codon <- count_codon[order(count_codon$n,decreasing=TRUE),]

site_codon$codon_pos <- as.numeric(site_codon$coding_pos)-as.numeric(site_codon$codon_start)+1

##translate codon into aa
genetic_code <- c(
    TTT = "Phe", TTC = "Phe", TTA = "Leu", TTG = "Leu",
    CTT = "Leu", CTC = "Leu", CTA = "Leu", CTG = "Leu",
    ATT = "Ile", ATC = "Ile", ATA = "Ile", ATG = "Met",
    GTT = "Val", GTC = "Val", GTA = "Val", GTG = "Val",
    TCT = "Ser", TCC = "Ser", TCA = "Ser", TCG = "Ser",
    CCT = "Pro", CCC = "Pro", CCA = "Pro", CCG = "Pro",
    ACT = "Thr", ACC = "Thr", ACA = "Thr", ACG = "Thr",
    GCT = "Ala", GCC = "Ala", GCA = "Ala", GCG = "Ala",
    TAT = "Tyr", TAC = "Tyr", TAA = "Stop", TAG = "Stop",
    CAT = "His", CAC = "His", CAA = "Gln", CAG = "Gln",
    AAT = "Asn", AAC = "Asn", AAA = "Lys", AAG = "Lys",
    GAT = "Asp", GAC = "Asp", GAA = "Glu", GAG = "Glu",
    TGT = "Cys", TGC = "Cys", TGA = "Stop", TGG = "Trp",
    CGT = "Arg", CGC = "Arg", CGA = "Arg", CGG = "Arg",
    AGT = "Ser", AGC = "Ser", AGA = "Arg", AGG = "Arg",
    GGT = "Gly", GGC = "Gly", GGA = "Gly", GGG = "Gly"
)
# Translate a single codon into an amino acid (three-letter abbreviation)
translate_codon <- function(codon) {
    return(genetic_code[[codon]])
}
site_codon$amino_acid <- sapply(site_codon$codon, translate_codon)

count_aa <- site_codon %>% group_by(amino_acid,codon_pos) %>% summarize(n=n())
count_aa <- count_aa[order(count_aa$n,decreasing=TRUE),]


site_codon$codon_pos <- as.factor(site_codon$codon_pos)
site_codon$amino_acid <- factor(site_codon$amino_acid,level=c("Phe", "Leu", "Val", "Ile", "Ser", "Tyr", "Asp", "Cys", "Thr","Met", "Ala","Trp","Gly", "Asn", "His", "Pro", "Arg", "Stop"))
pdf("plot/mRNA_codon_stack_plot.pdf",width=7,height=3)
ggplot(site_codon, aes(x = amino_acid, fill = codon_pos)) + 
  geom_bar(width=0.6) +theme_bw()+ylab("site counts") +scale_fill_manual(values=c("#FB6962", "#79DE79","#1ECBE1"))#scale_fill_brewer(palette = "Dark2")#scale_fill_manual(values=c("#E41A1C", "#377EB8","#4DAF4A"))
dev.off()
##+geom_text(stat = "count", aes(label = after_stat(count)),position = position_stack(vjust = 0.5))
count_aa <- site_codon %>% group_by(amino_acid) %>% summarize(n=n())
count_aa <- count_aa[order(count_aa$n,decreasing=TRUE),]

#######bar plot for codon_pos
count_pos <- site_codon %>% group_by(codon_pos) %>% summarize(n=n())

  codon_pos     n
  <fct>     <int>
1 1           226
2 2           359
3 3           226
pdf("plot/mRNA_codon_pos_bar_plot.pdf",width=3,height=3)
ggplot(count_pos, aes(x = codon_pos, y=n,fill = codon_pos)) + 
  geom_bar(width=0.6,stat='identity') +geom_text(aes(label=n),position = position_stack(vjust = 0.5)) +scale_fill_manual(values=c("#FB6962", "#79DE79","#1ECBE1")) +theme_bw()+ylab("site counts")+theme(legend.position = "none")
dev.off()

### stack plot for codon and pos

site_codon$codon_pos <- as.factor(site_codon$codon_pos)
site_codon$codon <- factor(site_codon$codon,level=c("TTT", "TTC", "TTG", "ATT", "GTT", "TAT", "ATC", "CTT", "GTG", "GAT", "GTA", "CTC","TCC", "TCT", "TGT", "CTG", "TAC", "ACT", "ATG", "GTC" ,"TCA" ,"GCT" ,"TTA" ,"TGG","GGT" ,"ATA" ,"AAT" ,"CAT" ,"CCT" ,"TCG", "AGT", "TGC", "CGT" ,"CTA" ,"TAG"))
pdf("plot/mRNA_codon_pos_stack_plot.pdf",width=12,height=3)
ggplot(site_codon, aes(x = codon, fill = codon_pos)) + 
  geom_bar(width=0.6) +scale_fill_manual(values=c("#FB6962", "#79DE79","#1ECBE1")) +theme_bw()+ylab("site counts")##+geom_text(stat = "count", aes(label = after_stat(count)),position = position_stack(vjust = 0.5)) 
dev.off()

##################TPM  correlation between input and treat sample 
expression <- read.table("/users/ludwig/ebu571/ebu571/27Oct2023/rpkm/H1_2023Oct26_S1.tpm",  header=T, sep="\t")
colnames(expression) <- c("geneid","chr", "start", "end", "strand","length", "input1", "input2","treat1", "treat2", "ivt1", "ivt2")
cor(expression$input1,expression$input2) ##0.9946491
cor(expression$input1,expression$treat1) ##0.9675942
cor(expression$input1,expression$treat2) ##0.9717791
expression$loginput1 <- log2(expression$input1+1)
expression$loginput2 <- log2(expression$input2+1)
expression$logtreat1 <- log2(expression$treat1+1)
expression$logtreat2 <- log2(expression$treat2+1)
expression$logivt1 <- log2(expression$ivt1+1)
expression$logivt2 <- log2(expression$ivt2+1)
cor(expression$loginput1,expression$loginput2) ##0.993334
cor(expression$loginput1,expression$logtreat1) ##0.9484423
cor(expression$loginput1,expression$logtreat2) ##0.9396075

pdf("plot/tpm_treat1_input1_correlation.pdf")
smoothScatter(expression$loginput1,expression$logtreat1)
dev.off()

pdf("plot/tpm_input1_input2_correlation.pdf")
smoothScatter(expression$loginput1,expression$loginput2)
dev.off()

exon_annotation <- read.table("/users/ludwig/ebu571/ebu571/resource/gencode.v43.annotation.bed", header=F) ##3424907
colnames(exon_annotation) <- c("chr", "start","end","strand","source","feature", "gene_id","transcript_id","gene_type","gene_name","transcript_type","transcript_name","exon_number")
gene_annotation <- exon_annotation[exon_annotation$feature=="gene",] ##62703
merge_exp <- merge(expression,unique(gene_annotation[,c("gene_name","gene_type")]),by.x=c("geneid"),by.y=c("gene_name"))
nrow(merge_exp) ##61432

##nrow(expression) 61180
merge_exp[duplicated(merge_exp$geneid),] %>% nrow() ##252
select_gene <- merge_exp[merge_exp$gene_type=="protein_coding",] ##20000
cor(select_gene$loginput1,select_gene$loginput2) ##0.9986752
cor(select_gene$loginput1,select_gene$logtreat1) ##0.9946678
cor(select_gene$loginput2,select_gene$logtreat2) ##0.99561

pdf("plot/tpm_pc_treat1_input1_correlation.pdf")
smoothScatter(select_gene$loginput1,select_gene$logtreat1,xlim=c(0,15),ylim=c(0,15))
dev.off()

pdf("plot/tpm_pc_input1_input2_correlation.pdf")
smoothScatter(select_gene$loginput1,select_gene$loginput2,xlim=c(0,15),ylim=c(0,15))
dev.off()

pdf("plot/tpm_pc_treat2_input2_correlation.pdf")
smoothScatter(select_gene$loginput2,select_gene$logtreat2,xlim=c(0,15),ylim=c(0,15))
dev.off()

select_dat <- select_gene[,c("geneid","loginput1", "loginput2", "logtreat1", "logtreat2", "logivt1", "logivt2")]
cor_table <- data.frame(matrix(nrow = 0, ncol = 3)) 
for (i in 2:7){
  for (j in 2:7){
    corvalue <- cor(select_dat[,i],select_dat[,j])
    test <-data.frame(i,j, corvalue)
    cor_table<- rbind(cor_table,test)
  }
  
}
colnames(cor_table) <- c("i","j", "cor")
write.table(cor_table,"plot/bacs_tpm_correlation.txt",sep='\t',quote=F,col.names = TRUE, row.names = FALSE)

cor_table$i <- factor(as.character(cor_table$i), level=c("2","3","4","5","6","7"))
cor_table$j <- factor(as.character(cor_table$j), level=c("2","3","4","5","6","7"))
cor_table <- read.table("plot/bacs_tpm_correlation.txt", header=T)

library(corrplot)
cot_value <- spread(cor_table,j,cor) %>% as.data.frame()
cot_value <- data.frame(cot_value[,-1], row.names=as.character(cot_value[,1])) %>% as.matrix(rownames = TRUE,colnames=TRUE)

pdf("plot/tpm_correlation_point.pdf")
corrplot(cot_value,type = 'lower',addCoef.col = 'white')
dev.off()
######################plot the RNA fractions
mapping_sta <- read.table("/gpfs3/well/ludwig/users/ebu571/27Oct2023/align/mrna_all.mapping.txt", header=TRUE)
mapping_sta$total_mapped_reads <- (mapping_sta$rrna_mapped+mapping_sta$snorna_map+mapping_sta$trna_map+mapping_sta$mRNA_map)

dat <- data.frame(matrix(NA,    # Create empty data frame
                          nrow = 0,
                          ncol = 5))
dat[1,] <- data.frame("input1", sum(mapping_sta[c(1),]$rrna_mapped),sum(mapping_sta[c(1),]$snorna_map),sum(mapping_sta[c(1),]$trna_map), sum(mapping_sta[c(1),]$mRNA_map))
dat[2,] <- data.frame("input2", sum(mapping_sta[c(2),]$rrna_mapped),sum(mapping_sta[c(2),]$snorna_map),sum(mapping_sta[c(2),]$trna_map), sum(mapping_sta[c(2),]$mRNA_map))
dat[3,] <- data.frame("treat1", sum(mapping_sta[c(3),]$rrna_mapped),sum(mapping_sta[c(3),]$snorna_map),sum(mapping_sta[c(3),]$trna_map), sum(mapping_sta[c(3),]$mRNA_map))
dat[4,] <- data.frame("treat2", sum(mapping_sta[c(4),]$rrna_mapped),sum(mapping_sta[c(4),]$snorna_map),sum(mapping_sta[c(4),]$trna_map), sum(mapping_sta[c(4),]$mRNA_map))
dat[5,] <- data.frame("ivt1", sum(mapping_sta[c(5),]$rrna_mapped),sum(mapping_sta[c(3),]$snorna_map),sum(mapping_sta[c(3),]$trna_map), sum(mapping_sta[c(3),]$mRNA_map))
dat[6,] <- data.frame("ivt2", sum(mapping_sta[c(6),]$rrna_mapped),sum(mapping_sta[c(4),]$snorna_map),sum(mapping_sta[c(4),]$trna_map), sum(mapping_sta[c(4),]$mRNA_map))

colnames(dat) <- c("smp","rrna_mapped", "snorna_map","trna_map","mrna_map")
dat$rrna_ratio <- dat$rrna_mapped/(dat$rrna_mapped+dat$snorna_map+dat$trna_map+dat$mrna_map)
dat$snorna_ratio <- dat$snorna_map/(dat$rrna_mapped+dat$snorna_map+dat$trna_map+dat$mrna_map)
dat$trna_ratio <- dat$trna_map/(dat$rrna_mapped+dat$snorna_map+dat$trna_map+dat$mrna_map)
dat$mrna_ratio <- dat$mrna_map/(dat$rrna_mapped+dat$snorna_map+dat$trna_map+dat$mrna_map)
meltdat <- melt(dat[,1:5],id.vars=c("smp"))
meltdat$variable <- factor(meltdat$variable, level=c("rrna_mapped","snorna_map", "trna_map","mrna_map"))

p1 <- ggplot(meltdat[meltdat$smp=="input1",], aes(x="", y=value, fill=variable))+
geom_bar(width = 1, stat = "identity")+coord_polar("y") + scale_fill_manual(values = c("#BEAED4","#FFFF99" , "#FDC086","#7FC97F")) +  ###or RdYlBu
  geom_text(aes(label = value),position = position_stack(vjust = 0.5))+
  theme_void()

p2 <- ggplot(meltdat[meltdat$smp=="input2",], aes(x="", y=value, fill=variable))+
geom_bar(width = 1, stat = "identity")+coord_polar("y") + scale_fill_manual(values = c("#BEAED4","#FFFF99" , "#FDC086","#7FC97F")) +  ###or RdYlBu
  geom_text(aes(label = value),position = position_stack(vjust = 0.5))+
  theme_void()

p3 <- ggplot(meltdat[meltdat$smp=="treat1",], aes(x="", y=value, fill=variable))+
geom_bar(width = 1, stat = "identity")+coord_polar("y") + scale_fill_manual(values = c("#BEAED4","#FFFF99" , "#FDC086","#7FC97F")) +  ###or RdYlBu
  geom_text(aes(label = value),position = position_stack(vjust = 0.5))+
  theme_void()


p4 <- ggplot(meltdat[meltdat$smp=="treat2",], aes(x="", y=value, fill=variable))+
geom_bar(width = 1, stat = "identity")+coord_polar("y") + scale_fill_manual(values = c("#BEAED4","#FFFF99" , "#FDC086","#7FC97F")) +  ###or RdYlBu
  geom_text(aes(label = value),position = position_stack(vjust = 0.5))+
  theme_void()

p5 <- ggplot(meltdat[meltdat$smp=="ivt1",], aes(x="", y=value, fill=variable))+
geom_bar(width = 1, stat = "identity")+coord_polar("y") + scale_fill_manual(values = c("#BEAED4","#FFFF99" , "#FDC086","#7FC97F")) +  ###or RdYlBu
  geom_text(aes(label = value),position = position_stack(vjust = 0.5))+
  theme_void()


p6 <- ggplot(meltdat[meltdat$smp=="ivt2",], aes(x="", y=value, fill=variable))+
geom_bar(width = 1, stat = "identity")+coord_polar("y") + scale_fill_manual(values = c("#BEAED4","#FFFF99" , "#FDC086","#7FC97F")) +  ###or RdYlBu
  geom_text(aes(label = value),position = position_stack(vjust = 0.5))+
  theme_void()

pdf("plot/mRNA_reads_pie_plot.pdf", width=6,height=6)
ggarrange(p1,p2,p3,p4,p5,p6,ncol = 2, nrow = 2,heights=c(2,2))
dev.off()

###rRNA site
mergedat <- read.table("/users/ludwig/ebu571/ebu571/07Nov2023/site/bacs_rrna_called_sites.txt", header=T, sep=" ", colClasses="character") 
mergedat <- mergedat[mergedat$chr=="NR_003285.3_RNA5_8SN5" | mergedat$chr=="NR_003286.4_RNA18SN5" | mergedat$chr=="NR_003287.4_RNA28SN5",] ##104
rrna_site <- mergedat
rrna_site$mean_level <- (as.numeric(rrna_site$treat1_level)+as.numeric(rrna_site$treat2_level))/2
rrna_site$type <-c("rrna")
###tRNA and snoRNA site
trna <- read.table("/users/ludwig/ebu571/ebu571/27Oct2023/site/snoRNA_tRNA_called_sites_fdrpassed_tRNA_filtered.txt",header=T, sep="\t",colClasses=c("ref"="character")) ###609
trna$mean_level <- (trna$treat1_level+trna$treat2_level)/2
trna$type <- c("trna")
snosite <-read.table("/users/ludwig/ebu571/ebu571/27Oct2023/site/snoRNA_tRNA_called_sites.txt",header=T, sep=" ",colClasses=c("ref"="character")) 
snosite <- snosite[!grepl("tRNA", snosite$chr),] ###304
snosite$mean_level <- (snosite$treat1_level+snosite$treat2_level)/2
snosite$type <- c("snorna")
##mRNA and ncRNA sites
annotated_site <- read.table("site/bacs.mRNA.annotated.txt",  header=TRUE) ##1337
filtered_site <- annotated_site[annotated_site$strand_infor=="proper_strand" &annotated_site$transcript_type !="snRNA" &annotated_site$transcript_detail!="intronic", ]
nrow(filtered_site) ##1335
chrsite <- filtered_site
chrsite$mean_level <- chrsite$treat_level2
chrsite$type <- c("transcriptome")

> summary(chrsite$mean_level)
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
0.05000 0.06608 0.09090 0.15328 0.15464 1.03986 

combinedat <- rbind(rrna_site[,c("type", "mean_level")],trna[,c("type", "mean_level")], snosite[,c("type","mean_level")], chrsite[,c("type", "mean_level")])
combinedat[combinedat$mean_level>1,]$mean_level <- 1
combinedat$type <- factor(combinedat$type, level=c("rrna","snorna","trna", "transcriptome"))

p <- ggplot(combinedat,aes(x=type,y=mean_level*100, fill=type))+ geom_boxplot(width=0.6, outlier.shape = NA)+ #+ geom_jitter(size=2,shape=16, colour = "grey",position=position_jitter(0.2))
  #stat_summary(fun = mean, geom = "point",color = "black") +
  theme_light() + xlab("") +  theme(axis.text = element_text(color = "black")) + ylim(0,100) +
  ylab("level") +scale_fill_manual(values=c("#BEAED4","#FFFF99" , "#FDC086","#7FC97F"))
ggsave("plot/modification_level_different_rna_type.pdf",p, width=4, height = 3)



##############################################compare pseudou sites on tRNA and rRNA
tRNA_index <- read.table("resource/tRNA_position_index.txt", header=T, sep="\t",colClasses="character")
colnames(tRNA_index) <- c("chr", "trna_pos", "ref", "pos")
tRNA_site <- read.table("/users/ludwig/ebu571/ebu571/27Oct2023/site/snoRNA_tRNA_called_sites_fdrpassed_tRNA_filtered.txt",header=T, sep="\t",colClasses=c("ref"="character")) 

trna <- merge(tRNA_index, tRNA_site, by.x=c("chr", "pos"), by.y=c("chr", "end"))
trna$mean_level <- (trna$treat1_level+trna$treat2_level)/2
trna <- trna[!grepl("mt",trna$chr),] ##609

###if mean level is over 100%, then corrected the level as 100%
trna[trna$mean_level>=1,]$mean_level <- 1

trna_pos55_guucn <- trna[trna$trna_pos=="55" &(trna$motif=="GTTCA" |trna$motif=="GTTCT"|trna$motif=="GTTCC"|trna$motif=="GTTCG"),] ##188
trna_pos55_guucn$type <- c("tRNA")
trna_pos55_guucn$motif <- c("GUUCN")
trna_pos13_usuag <- trna[trna$trna_pos=="13" &(trna$motif=="TCTAG" |trna$motif=="TGTAG"),] ##28
trna_pos13_usuag$type <- c("tRNA")
trna_pos13_usuag$motif <- c("USUAG")
####for motif dat in mRNA only
nrow(mRNA_site) ##1294
mRNA_site$group <- stri_replace_all_regex(mRNA_site$motif,
                                  pattern=c('TGTAG', 'TCTAG', 'GTTCA',"GTTCC","GTTCG","GTTCT","CTTTG","ACTTT","CATTT","GTTTC","TTTTT","GTTAA","GTTGA"),
                                  replacement=c('USUAG', 'USUAG', 'GUUCN','GUUCN','GUUCN','GUUCN',"U_rich","U_rich","U_rich","U_rich","U_rich","GUURA","GUURA"),
                                  vectorize=FALSE)
mrna_guucn <- mRNA_site[mRNA_site$motif=="GTTCA" | mRNA_site$motif=="GTTCC" |mRNA_site$motif=="GTTCG"|mRNA_site$motif=="GTTCT",] ##125
mrna_usuag <- mRNA_site[mRNA_site$motif=="TCTAG" | mRNA_site$motif=="TGTAG",] ##95
mrna_guucn$type <- c("mRNA")
mrna_guucn$motif <- c("GUUCN")
mrna_usuag$type <- c("mRNA")
mrna_usuag$motif <- c("USUAG")
comdat <- rbind(trna_pos55_guucn[,c("type","motif", "mean_level")],trna_pos13_usuag[,c("type","motif", "mean_level")],mrna_guucn[,c("type","motif", "mean_level")],mrna_usuag[,c("type","motif", "mean_level")])
nrow(comdat) ###429
comdat$type <- factor(comdat$type,level=c("tRNA","mRNA"))
p <- ggplot(comdat,aes(x=motif,y=mean_level*100, fill=type)) +geom_boxplot(width=0.6, outlier.shape = NA)+# geom_point(aes(fill = variable), size = 1, shape = 21, position = position_jitterdodge()) #geom_jitter(color="black", size=0.4) +##outlier.shape = NA, width=0.8 alpha=0.9
  theme_light() + xlab("") +  theme(axis.text = element_text(color = "black")) + ylim(0,100) +
  ylab("level") +
  scale_fill_manual(values=c("#FDC086","#7FC97F"))
#theme(legend.position = "none") 
ggsave("plot/mRNA_tRNA_compare_boxplot.pdf",p, width=4, height = 4) ##geom_boxplot(outlier.shape = NA, width=0.8)



###############gene list for GO annalysis
annotated_site <- read.table("site/bacs.mRNA.annotated.txt",  header=TRUE) 
filtered_site <- annotated_site[annotated_site$strand_infor=="proper_strand" &annotated_site$transcript_type !="snRNA" &annotated_site$transcript_detail!="intronic", ]
nrow(filtered_site) ##1335
mRNA_site <- filtered_site[filtered_site$transcript_detail!="exon",] ##1294
gene_name <- unique(mRNA_site[,c("gene_name")]) ##1070
write.table(data.frame(gene_name),"plot/mRNA_genename.txt",sep='\t',quote=F,row.names=F)



##library(enrichR)
dbs_go <- c("GO_Biological_Process_2023")
dbs_pw <- c("KEGG_2021_Human")
Enriched_go <- enrichr(genes = gene_name, databases = dbs_go)
Enriched_pw <- enrichr(genes = gene_name, databases = dbs_pw)
pdf("plot/mRNA_enrich_go_bp.pdf")
plotEnrich(Enriched_go[[1]], showTerms = 15, y = "Count", orderBy = "P.value")
dev.off()


 df <- transform(Enriched_go[[1]], count = colsplit(Overlap, "/", names = c('1', '2'))) %>% head(n=10)
 df <- df[order(df$count.1,decreasing=TRUE),]
 df$Term <- factor(df$Term,level=rev(df$Term))
 pdf("plot/mRNA_enrich_go_bp1.pdf",width=8,height=4)
ggplot(df,aes(x=count.1/1070,y=Term,size=count.1))+geom_point(aes(color = Adjusted.P.value), alpha = 1.0)+theme_bw()+scale_color_gradient(low="red", high="blue")
dev.off()

df1 <- transform(Enriched_pw[[1]], count = colsplit(Overlap, "/", names = c('1', '2'))) %>% head(n=15)
 df1$Term <- factor(df1$Term)
 df1$Term <- factor(df1$Term,level=rev(df1$Term))
 pdf("plot/mRNA_enrich_keggpw.pdf",width=10,height=4)
ggplot(df1,aes(x=count.1/1070,y=Term,size=count.1))+geom_point(aes(color = Adjusted.P.value), alpha = 1.0)+theme_bw()
dev.off()



#################metagene plot
meta_dat <- read.table("/users/ludwig/ebu571/ebu571/bacs/site/treat.metaplot.txt",  header=FALSE) ##1294
pdf("plot/psi_metagene_plot.pdf",width=5,height=4)
## A smooth density plot
plot(density(meta_dat$V1,bw = 0.08),xaxt='n',xlab='',ylab='',main='',tck=-0.02,ylim=c(0,0.8))
abline (v = 1, lty  = 1, col = "red")
abline (v = 2, lty  = 1, col = "red")
dev.off()

####metaplot for inosine
meta_dat_inosine <- read.table("/users/ludwig/ebu571/ebu571/bacs/inosine/inosine.metaplot.txt",  header=FALSE) ##737

pdf("plot/inosine_metagene_plot.pdf",width=5,height=4)
## A smooth density plot
plot(density(meta_dat_inosine$V1,bw = 0.2),xaxt='n',xlab='',ylab='',main='',tck=-0.02,ylim=c(0,0.4))
abline (v = 1, lty  = 1, col = "red")
abline (v = 2, lty  = 1, col = "red")
dev.off()