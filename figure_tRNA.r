##module load R/4.0.3-foss-2020b
####this script is for figures on tRNA and snoRNA
library(dplyr)
library(reshape2)
library(tidyr)
library(ggplot2)
library(ggpubr)
library(cowplot)
library(RColorBrewer)
library(tidyverse)
library(gridExtra)
setwd("/users/ludwig/ebu571/ebu571/bacs")
##for tRNA and others: /users/ludwig/ebu571/ebu571/27Oct2023/site/snoRNA_tRNA_called_sites_new.txt and /users/ludwig/ebu571/ebu571/27Oct2023/site/snoRNA_tRNA_called_sites_fdrpassed_tRNA_filtered.txt
##for mt-site: /users/ludwig/ebu571/ebu571/27Oct2023/site/snoRNA_tRNA_called_sites.txt
####venn plot on site distribution of cytosolic tRNA
###read called tables
tRNA_site <- read.table("/users/ludwig/ebu571/ebu571/27Oct2023/site/snoRNA_tRNA_called_sites_fdrpassed_tRNA_filtered.txt",header=T, sep="\t",colClasses=c("ref"="character")) ###609

tRNA_site$tmp <- gsub("^tRNA_", "", tRNA_site$chr)
tRNA_site$aa <- gsub("_.*", "", tRNA_site$tmp)

filteredlist <- read.table("/users/ludwig/ebu571/ebu571/27Oct2023/featureCounts/tRNA_counts_filtered.txt", header=T, sep="\t")

all_tRNA <- merge(filteredlist[,c("chr", "anticodon")], tRNA_site[,c("chr", "motif")],by=c("chr"), all.x=TRUE)
all_tRNA[is.na(all_tRNA$motif), ]
all_tRNA$tmp <- gsub("^tRNA_", "", all_tRNA$chr)
all_tRNA$aa <- gsub("_.*", "", all_tRNA$tmp)

count <- tRNA_site %>% group_by(aa) %>% summarize(n = n()) 
nb.cols <- 21
mycolors <- colorRampPalette(brewer.pal(8, "RdYlBu"))(nb.cols)

pdf("plot/tRNA_aa_pie.pdf")
ggplot(count, aes(x="", y=n, fill=aa))+
geom_bar(width = 1, stat = "identity")+coord_polar("y") + scale_fill_manual(values = mycolors) +  ###or RdYlBu
  geom_text(aes(label = n),position = position_stack(vjust = 0.5))+
  theme_void()
dev.off()

#######obtain isoform information
count1 <- tRNA_site %>% group_by(chr,aa) %>% summarize(n = n()) 
count2 <- merge(filteredlist[,c("chr", "anticodon")],count1,by=c("chr"), all.x=TRUE, all.y=TRUE)
count2[is.na(count2$n),4] <- 0
count2$tmp <- gsub("^tRNA_", "", count2$chr)
count2$aa <- gsub("_.*", "", count2$tmp)
count1 <- count2
isoform <- data.frame(matrix(nrow=0,ncol=2))
for (i in unique(count1$aa)){
  aainfor <- count1[count1$aa==i,]
  isoform_num <- median(aainfor$n) %>% as.numeric()
  isoform <- rbind(isoform,c(i,isoform_num))
}
colnames(isoform) <- c("aa", "median_count_per_aa")
#isoform[nrow(isoform)+1,]<- c("iMet", "0")

########mt_tRNA_sites 
site <- read.table("/users/ludwig/ebu571/ebu571/27Oct2023/site/snoRNA_tRNA_called_sites.txt",header=TRUE, sep=" ",colClasses=c("ref"="character")) 
mtsite <- site[grepl("mt", site$chr),] ##54
mtsite$aa <- gsub("mt_tRNA_", "", mtsite$chr)
#mtsite$aa <- gsub("_.*", "", mtsite$tmp)

mtcount <- mtsite %>% group_by(aa) %>% summarize(n = n()) 
nb.cols <- 19
mycolors <- colorRampPalette(brewer.pal(8, "RdYlBu"))(nb.cols)

pdf("plot/mtRNA_aa_pie.pdf")
ggplot(mtcount, aes(x="", y=n, fill=aa))+
geom_bar(width = 1, stat = "identity")+coord_polar("y") + scale_fill_manual(values = mycolors) +  ###or RdYlBu
  geom_text(aes(label = n),position = position_stack(vjust = 0.5))+
  theme_void()
dev.off()



mtcount <- rbind(mtcount,c("Leu", as.numeric(median(mtcount[mtcount$aa=="Leu_CUN" | mtcount$aa=="Leu_UUR",]$n))))
mtcount[grep("Ser", mtcount$aa),1]<- c("Ser")
mtcount[mtcount$aa=="Ser",2]<- as.character(median(c(0,2)))        ###two isoform
mtcount <- mtcount[mtcount$aa!="Leu_CUN" & mtcount$aa!="Leu_UUR", ]

aa <- c("Ala","Arg","Asn","Asp","Cys","Gln","Glu","Gly","His","Ile","Leu","Lys","Met","Phe","Pro","Ser","Tyr","Val","Thr","Trp", "SeC","iMet")
trna_list <- merge(data.frame(aa), isoform, by=c("aa"), all.x=TRUE,all.y=TRUE)
trna_mttrna <- merge(trna_list, mtcount, by=c("aa"), all.x=TRUE,all.y=TRUE) 
trna_mttrna[is.na(trna_mttrna$n),3] <- 0
colnames(trna_mttrna) <- c("aa","tRNA_median","mt_median")
trna_mttrna$tRNA_median <- as.numeric(trna_mttrna$tRNA_median)
trna_mttrna$mt_median <- as.numeric(trna_mttrna$mt_median)
trna_mttrna$aa <- factor(trna_mttrna$aa, level=rev(trna_mttrna$aa))

pdf("plot/mtRNA_aa_bar.pdf",width=5, height=10)
mtcount %>% 
ggplot(aes(x=reorder(aa,n), y=n, fill=aa))+
geom_bar(stat="identity", fill="#f68060", alpha=.6, width=.4)+
  geom_text(aes(label = n),position = position_stack(vjust = 0.5))+ coord_flip()+
  theme_bw()+ylab("count")+xlab("aa")
dev.off()

#######intergrative bar graph
p1 <- trna_mttrna %>% 
ggplot(aes(x=aa, y = -tRNA_median, fill=aa))+
geom_bar(stat="identity", fill="#f68060", alpha=.6, width=.4)+
  geom_text(aes(label = round(tRNA_median,1)),position = position_stack(vjust = 0))+ coord_flip()+
  theme_bw()+
  ylim(-6,0)+
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.y = element_blank(),
        plot.margin = margin(5.5, 0, 5.5, 5.5))
#  scale_x_discrete(limits = rev(levels(trna_mttrna$aa)))+
p2 <- trna_mttrna %>% 
ggplot(aes(x=aa, y=mt_median, fill=aa))+
geom_bar(stat="identity", fill="#f68060", alpha=.6, width=.4)+
  geom_text(aes(label = mt_median),position = position_stack(vjust = 1))+ coord_flip()+ylim(0,6)+
  theme_bw()+  theme(axis.title.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.ticks.length = unit(0, "pt"),
        plot.margin = margin(5.5, 5.5, 5.5, 0),
        axis.text.y.left = element_text(margin = margin(0, 5.5, 0, 5.5)))
#library(gridExtra)
pdf("plot/tRNA_mttRNA_intergrative_barchart_legend.pdf", width=8, height=13)
grid.arrange(p1,
    p2,
    widths=c(0.6,0.6),
    ncol=2
)
dev.off()



p1 <- trna_mttrna %>% 
ggplot(aes(x=aa, y = -tRNA_median, fill=aa))+
geom_bar(stat="identity", fill="#006699", alpha=.6, width=.5)+
  geom_text(aes(label = round(tRNA_median,1)),position = position_stack(vjust = 0))+ coord_flip()+
  theme_bw()+
  ylim(-6,0)+
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.y = element_blank(),
        plot.margin = margin(5.5, 0, 5.5, 5.5))

p2 <- trna_mttrna %>% 
ggplot(aes(x=aa, y=mt_median, fill=aa))+
geom_bar(stat="identity", fill="#EE7621", alpha=.6, width=.5)+ 
  geom_text(aes(label = mt_median),position = position_stack(vjust = 1))+ coord_flip()+ylim(0,6)+
  theme_bw()+  theme(axis.title.y = element_blank(),
        axis.ticks.y = element_blank(),
        plot.margin = margin(5.5, 5.5, 5.5, 0),
        axis.text.y=element_blank())
#library(gridExtra)
pdf("plot/tRNA_mttRNA_intergrative_barchart.pdf", width=4, height=6.5)
grid.arrange(p1,
    p2,
    widths=c(0.6,0.6),
    ncol=2
)
dev.off()


#######################################heatmap for tRNA#############################################
##################get tRNA table 
tRNA_index <- read.table("resource/tRNA_position_index.txt", header=T, sep="\t",colClasses="character")
colnames(tRNA_index) <- c("chr", "trna_pos", "ref", "pos")

tRNA_site <- tRNA_site[,c("chr", "end","treat1_level", "treat2_level")]


trna <- merge(tRNA_index, tRNA_site, by.x=c("chr", "pos"), by.y=c("chr", "end"), all.x=TRUE)
trna[is.na(trna)] <- 0
trna$mean <- (trna$treat1_level+trna$treat2_level)/2
trna <- trna[,c("chr", "trna_pos", "mean")]
###if mean level is over 100%, then corrected the level as 100%
trna[trna$mean>=1,]$mean <- 1
trna$trna_pos <- factor(trna$trna_pos, level=c("0","1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","17A","18","19","20","20A","20B","21","22","23","24","25","26","27","28","29","30","31","32","33","34","35","36","37","38","39","40","41","42","43","44","45","e11","e12","e13","e14","e15","e16","e17","e1","e2","e3","e4","e5","e27","e26","e25","e24","e23","e22","e21","46","47","48","49","50","51","52","53","54","55","56","57","58","59","60","61","62","63","64","65","66","67","68","69","70","71","72","73","74","75","76"))
tRNA <- trna[!grepl("mt",trna$chr),]

final <- spread(tRNA, key = trna_pos, value = mean) %>%  as.data.frame()
write.table(final,"plot/tRNA_heatmap_table.txt",sep='\t',quote=F,row.names=F)

pdf("plot/tRNA_heatmap.pdf", width=50, height=25)
ggplot(tRNA, aes(x = trna_pos, y = chr, fill = mean*100)) + geom_tile()+ coord_fixed()+ ##to make Square tiles color = "gray", size=0.1
scale_fill_gradient2(low = "white", high = "red", limits=c(0,100), guide = guide_colorbar(frame.colour = "black",ticks = TRUE,))+ theme(panel.border=element_rect(fill = NA, colour='black',size=1))
dev.off()

#####part of the anticodon
seltrna <- read.table("resource/representative_tRNA.txt", header=F, sep="\t",colClasses="character")
colnames(seltrna) <- c("chr")
selsite <- merge(trna,seltrna, by=c("chr"))  ##unique(selsite$chr) %>% length() 47
pdf("plot/tRNA_heatmap_representatives.pdf", width=20, height=12.5)
ggplot(selsite, aes(x = trna_pos, y = chr, fill = mean*100)) + geom_tile(color = "white", size=0.01)+coord_flip()+ #coord_fixed()+ ##to make Square tiles color = "gray", size=0.1
scale_fill_gradient2(low = "white", high = "red", limits=c(0,100), guide = guide_colorbar(frame.colour = "black",ticks = TRUE,))+ theme(panel.border=element_rect(fill = NA, colour='black',size=1))+theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()

#################heatmap for mtRNA
mtsite <- read.table("/users/ludwig/ebu571/ebu571/27Oct2023/site/snoRNA_tRNA_called_sites.txt",header=TRUE, sep=" ",colClasses=c("ref"="character")) 
mtsite <- mtsite[,c("chr", "end","treat1_level", "treat2_level")]

#mtindex <- tRNA_index[grepl("mt", tRNA_index$chr),]
mttrna <- merge(tRNA_index[grepl("mt", tRNA_index$chr) & !grepl("e",tRNA_index$trna_pos),], mtsite, by.x=c("chr", "pos"), by.y=c("chr", "end"), all.x=TRUE)
mttrna[is.na(mttrna)] <- 0
mttrna$mean <- (mttrna$treat1_level+mttrna$treat2_level)/2
mttrna <- mttrna[,c("chr", "trna_pos", "mean")]
mttrna$trna_pos <- factor(mttrna$trna_pos, level=c("0","1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","17A","18","19","20","20A","20B","21","22","23","24","25","26","27","28","29","30","31","32","33","34","35","36","37","38","39","40","41","42","43","44","45","46","47","48","49","50","51","52","53","54","55","56","57","58","59","60","61","62","63","64","65","66","67","68","69","70","71","72","73","74","75","76"))

final <- spread(mttrna, key = trna_pos, value = mean) %>%  as.data.frame()
write.table(final,"plot/mtRNA_heatmap_table.txt",sep='\t',quote=F,row.names=F)

pdf("plot/mtRNA_heatmap.pdf", width=16, height=8)
ggplot(mttrna, aes(x = trna_pos, y = chr, fill = mean*100)) + geom_tile()+ coord_fixed()+ ##to make Square tiles color = "gray", size=0.1
  #geom_text(aes(label = mean), color = "black", size = 6) +
  scale_fill_gradient2(low = "white", high = "red", limits=c(0,100))
dev.off()

###plot raw conversion for all U sites in mt-tRNA
allsite <- read.table("/users/ludwig/ebu571/ebu571/27Oct2023/site/snoRNA_tRNA_all_sites.txt",header=T, sep=" ",colClasses=c("ref"="character")) ##colClasses=c("ref"="character")
mtall <- allsite[grepl("mt", allsite$chr),] %>% dplyr::select(chr, end, treat1_conversion, treat2_conversion, ctrl1_conversion, ctrl2_conversion)
mtall$delta_conversion <- ((mtall$treat1_conversion-mtall$ctrl1_conversion) + (mtall$treat2_conversion - mtall$ctrl2_conversion)) /2

mttrna <- merge(tRNA_index[grepl("mt", tRNA_index$chr) & !grepl("e",tRNA_index$trna_pos),], mtall, by.x=c("chr", "pos"), by.y=c("chr", "end"), all.x=TRUE)
mttrna[is.na(mttrna)] <- 0

mttrna <- mttrna[,c("chr", "trna_pos", "delta_conversion")]
mttrna$trna_pos <- factor(mttrna$trna_pos, level=c("0","1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","17A","18","19","20","20A","20B","21","22","23","24","25","26","27","28","29","30","31","32","33","34","35","36","37","38","39","40","41","42","43","44","45","46","47","48","49","50","51","52","53","54","55","56","57","58","59","60","61","62","63","64","65","66","67","68","69","70","71","72","73","74","75","76"))

final <- spread(mttrna, key = trna_pos, value = delta_conversion) %>%  as.data.frame()
write.table(final,"plot/mtRNA_raw_conversion_heatmap_table.txt",sep='\t',quote=F,row.names=F)

pdf("plot/mtRNA_raw_conversion_heatmap.pdf", width=80, height=10)
ggplot(mttrna, aes(x = trna_pos, y = chr, fill = delta_conversion*100)) + geom_tile()+ coord_fixed()+ ##to make Square tiles 
  #geom_text(aes(label = delta_conversion), color = "black", size = 6) +
  scale_fill_gradient2(low = "white", high = "red", limits=c(-5,105))#+, limits=c(0,100)
dev.off()


###plot level for all U sites in mt-tRNA
#allsite <- read.table("site/snoRNA_tRNA_all_sites.txt",header=T, sep=" ",colClasses=c("ref"="character")) ##colClasses=c("ref"="character")
mtall <- allsite[grepl("mt", allsite$chr),] %>% dplyr::select(chr, end, treat1_level, treat2_level)
mtall$mean_level <- (mtall$treat1_level+mtall$treat2_level)/2

mttrna <- merge(tRNA_index[grepl("mt", tRNA_index$chr) & !grepl("e",tRNA_index$trna_pos),], mtall, by.x=c("chr", "pos"), by.y=c("chr", "end"), all.x=TRUE)
mttrna[is.na(mttrna)] <- 0

mttrna <- mttrna[,c("chr", "trna_pos", "mean_level")]
mttrna[mttrna$mean<0,]$mean_level <- 0
mttrna[mttrna$mean>1,]$mean_level <- 1
mttrna$trna_pos <- factor(mttrna$trna_pos, level=c("0","1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","17A","18","19","20","20A","20B","21","22","23","24","25","26","27","28","29","30","31","32","33","34","35","36","37","38","39","40","41","42","43","44","45","46","47","48","49","50","51","52","53","54","55","56","57","58","59","60","61","62","63","64","65","66","67","68","69","70","71","72","73","74","75","76"))

final <- spread(mttrna, key = trna_pos, value = mean_level) %>%  as.data.frame()
write.table(final,"plot/mtRNA_mean_level_all_site_heatmap_table.txt",sep='\t',quote=F,row.names=F)

pdf("plot/mtRNA_mean_level_all_site_heatmap.pdf", width=20, height=10)
ggplot(mttrna, aes(x = trna_pos, y = chr, fill = mean_level*100)) + geom_tile(color = "white", size=0.01)+ coord_flip()+#coord_fixed()+ ##to make Square tiles 
  #geom_text(aes(label = delta_conversion), color = "black", size = 6) +
  scale_fill_gradient2(low = "white", high = "red", limits=c(0,100), guide = guide_colorbar(frame.colour = "black",ticks = TRUE,))+ theme(panel.border=element_rect(fill = NA, colour='black',size=1))+theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()

#################violin plot at each position  27-28, 38,39,40,54 and 55
tRNA_site <- read.table("/users/ludwig/ebu571/ebu571/27Oct2023/site/snoRNA_tRNA_called_sites_fdrpassed_tRNA_filtered.txt",header=T, sep="\t",colClasses=c("ref"="character")) ##609

tRNA_index <- read.table("resource/tRNA_position_index.txt", header=T, sep="\t",colClasses="character")
colnames(tRNA_index) <- c("chr", "trna_pos", "ref", "pos")

trna <- merge(tRNA_index, tRNA_site, by.x=c("chr", "pos"), by.y=c("chr", "end"), all.x=TRUE)
trna[is.na(trna)] <- 0
trna$mean <- (trna$treat1_level+trna$treat2_level)/2
trna <- trna[,c("chr", "trna_pos", "mean")]
###if mean level is over 100%, then corrected the level as 100%
trna[trna$mean>=1,]$mean <- 1
trna$trna_pos <- factor(trna$trna_pos, level=c("0","1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","17A","18","19","20","20A","20B","21","22","23","24","25","26","27","28","29","30","31","32","33","34","35","36","37","38","39","40","41","42","43","44","45","e11","e12","e13","e14","e15","e16","e17","e1","e2","e3","e4","e5","e27","e26","e25","e24","e23","e22","e21","46","47","48","49","50","51","52","53","54","55","56","57","58","59","60","61","62","63","64","65","66","67","68","69","70","71","72","73","74","75","76"))
tRNA <- trna[!grepl("mt",trna$chr),]

trna_site <- trna[trna$mean>0,] ## 13,20B, 27,28,32,38,39,40,e12,54,55    ###609
trna_site <- trna_site[trna_site$trna_pos=="13" | trna_site$trna_pos=="20B" | trna_site$trna_pos=="27" |trna_site$trna_pos=="28" |trna_site$trna_pos=="32"| trna_site$trna_pos=="38"|trna_site$trna_pos=="39"|trna_site$trna_pos=="40"|trna_site$trna_pos=="e12"|trna_site$trna_pos=="54"|trna_site$trna_pos=="55", ]
##nrow(trna_site) 569
p <- ggplot(trna_site,aes(x=trna_pos,y=mean*100)) + geom_jitter(size=2,shape=16, colour = "grey",position=position_jitter(0.2))+ geom_boxplot(outlier.shape = NA, color="#2166AC",fill="#006699",width=0.5,alpha = 0.2, lwd=1, fatten=1) +#geom_jitter(color="black", size=0.4) +##outlier.shape = NA, width=0.8 alpha=0.9
  theme_light() + xlab("pos") +  theme(axis.text = element_text(color = "black")) +# ylim(0,100)+
  ylab("level") +
  #scale_fill_manual(values=c("#67A9CF","#EF8A62"))+
  #scale_color_manual(values=c("#67A9CF","#EF8A62"))+
  theme(legend.position = "none") 
ggsave("plot/modification_level_trna_position.pdf",p,width=16, height = 4) ##, width=16, height = 4



#mtindex <- tRNA_index[grepl("mt", tRNA_index$chr),]
mttrna <- merge(tRNA_index[grepl("mt", tRNA_index$chr) & !grepl("e",tRNA_index$trna_pos),], mtsite, by.x=c("chr", "pos"), by.y=c("chr", "end"), all.x=TRUE)
mttrna[is.na(mttrna)] <- 0
mttrna$mean <- (mttrna$treat1_level+mttrna$treat2_level)/2
mttrna <- mttrna[,c("chr", "trna_pos", "mean")]
mttrna$trna_pos <- factor(mttrna$trna_pos, level=c("0","1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","17A","18","19","20","20A","20B","21","22","23","24","25","26","27","28","29","30","31","32","33","34","35","36","37","38","39","40","41","42","43","44","45","46","47","48","49","50","51","52","53","54","55","56","57","58","59","60","61","62","63","64","65","66","67","68","69","70","71","72","73","74","75","76"))

mttrna_site <- mttrna[mttrna$mean>0,] ##54
mttrna_site[mttrna_site$mean>1,]$mean <- 1
mttrna_site <- mttrna_site[mttrna_site$trna_pos=="27" |mttrna_site$trna_pos=="28"|mttrna_site$trna_pos=="39"|mttrna_site$trna_pos=="40"|mttrna_site$trna_pos=="55", ]
# mttrna_site %>% group_by(trna_pos) %>% summarize(n=n())
#`summarise()` ungrouping output (override with `.groups` argument)
# A tibble: 5 x 2
  trna_pos     n
  <fct>    <int>
1 27          10
2 28          11
3 39           9
4 40           5
5 55           6
mttrna_site$name <- gsub("27", "first", mttrna_site$trna_pos)
mttrna_site$name <- gsub("28", "first", mttrna_site$name)
mttrna_site$name <- gsub("39", "second", mttrna_site$name)
mttrna_site$name <- gsub("40", "second", mttrna_site$name)
mttrna_site$name <- gsub("first", "27_28", mttrna_site$name)
mttrna_site$name <- gsub("second", "39_40", mttrna_site$name)
p <- ggplot(mttrna_site,aes(x=name,y=mean*100)) +   geom_jitter(size=2,shape=16, colour = "grey",position=position_jitter(0.2)) + geom_boxplot(outlier.shape = NA, color="#EF8A62",fill="#F4A582",width=0.5,alpha = 0.2, lwd=1, fatten=1) +
  theme_light() + xlab("pos") +  theme(axis.text = element_text(color = "black")) +#ylim(0,102)+
  ylab("level") +
  #scale_fill_manual(values=c("#67A9CF","#EF8A62"))+
  #scale_color_manual(values=c("#67A9CF","#EF8A62"))+
  theme(legend.position = "none") 
ggsave("plot/modification_level_mttrna_position.pdf",p, width=4, height = 4)

mttrna$name <- c("mttRNA")
trna$name <- c("cytotRNA")
combinesite <- rbind(mttrna[mttrna$mean>0,], trna[trna$mean>0,]) ###663

p <- ggplot(combinesite,aes(x=name,y=mean*100,color=name)) +  geom_boxplot(width=0.8) + geom_point(aes(fill = name), size = 1, shape = 21, position = position_jitterdodge()) +#geom_jitter(color="black", size=0.4) +##outlier.shape = NA, width=0.8 alpha=0.9
  #stat_summary(fun = mean, geom = "point",color = "black") +
  #stat_summary(fun = mean, geom="text", aes(label = round(..y.., 2)), hjust = -0.2, vjust=-1.5)+
  theme_light() + xlab("pos27") +  theme(axis.text = element_text(color = "black")) + #ylim(0,100)+
  ylab("level") +
  scale_fill_manual(values=c("#67A9CF","#EF8A62"))+
  scale_color_manual(values=c("#67A9CF","#EF8A62"))+
  theme(legend.position = "none") 
ggsave("plot/modification_level_cy_vs_mt.pdf",p, width=4, height = 4)

##############draw overlapping venn plot between bacs and ms
###get merged mt-tRNA list from Hela, c666.1, raji, elijah
mt_merge <- read.table("site/bacs_mt-trna_allcelllines_merge_sites.txt", header=T, sep=" ", colClasses=c("ref"="character") )
nrow(mt_merge) ###65
mt_merge <- merge(tRNA_index[grepl("mt", tRNA_index$chr) & !grepl("e",tRNA_index$trna_pos),], mt_merge[,c("chr", "start", "end")], by.x=c("chr", "pos"), by.y=c("chr", "end"), all.y=TRUE) ###65
write.table(mt_merge, "site/mtsite_merge_list.txt", quote = FALSE, col.names = T, row.names = F)

library(VennDiagram)##########
pdf("plot/bacs_ms_venn_on_mtRNA.pdf")
draw.pairwise.venn(area1=65 , area2=52,cross.area=48, 
                   category=c("bacs","MS"),col=c("#FE767C","#666666"), cat.pos=c(-27, 27), cat.dist=c(0.03, 0.03), cat.col=c("#FE767C","#666666"), height = 300 , 
        width = 300)
dev.off()

pdf("plot/bacs_praise_venn_on_mtRNA.pdf")
draw.pairwise.venn(area1 = 65,    # Draw pairwise venn diagram
                   area2 = 34,
                   cross.area = 33,category=c("bacs","praise_seq"),col=c("#FE767C","#9ECE88"),fill=c("#FE767C","#9ECE88"), cat.cex=0,cat.pos=c(-27, 27), cat.dist=c(0.03, 0.03), cat.col=c("#FE767C","#666666"), height = 50 , 
        width = 50, lwd=0.5,ext.text=FALSE)
dev.off()

##plot modification level for trna_pos in different cell line
mt_merge <- read.table("site/bacs_mt-trna_allcelllines_merge_sites.txt", header=T, sep=" ", colClasses=c("ref"="character") )
mt_merge$hela_level <- (mt_merge$treat1_level+mt_merge$treat2_level)/2
mt_merge$c666_level <- (mt_merge$treat_level.x+mt_merge$treat_level.y)/2
mt_merge$raji_level <- mt_merge$treat_level.x.x
mt_merge$elijah_level <- mt_merge$treat_level.y.y
mt_merge <- merge(tRNA_index[grepl("mt", tRNA_index$chr) & !grepl("e",tRNA_index$trna_pos),], mt_merge[,c("chr", "start", "end","hela_level", "c666_level", "raji_level", "elijah_level")], by.x=c("chr", "pos"), by.y=c("chr", "end"), all.y=TRUE) ###65
# table(mt_merge$trna_pos)
#20 25 27 28 31 32 38 39 40 50 55 62 65 66 67 68  8 
#1  1 10 11  1  3  2  9  5  2  7  1  1  5  1  4  1

select_site <- mt_merge[mt_merge$trna_pos=="27" |mt_merge$trna_pos=="28"|mt_merge$trna_pos=="38"|mt_merge$trna_pos=="39"|mt_merge$trna_pos=="40"|mt_merge$trna_pos=="55" |mt_merge$trna_pos=="66"|mt_merge$trna_pos=="67"|mt_merge$trna_pos=="68", ]
select_site$name <- gsub("27", "first", select_site$trna_pos)
select_site$name <- gsub("28", "first", select_site$name)
select_site$name <- gsub("38", "second", select_site$name)
select_site$name <- gsub("39", "second", select_site$name)
select_site$name <- gsub("40", "second", select_site$name)
select_site$name <- gsub("66", "third", select_site$name)
select_site$name <- gsub("67", "third", select_site$name)
select_site$name <- gsub("68", "third", select_site$name)
select_site$name <- gsub("first", "27_28", select_site$name)
select_site$name <- gsub("second", "38_40", select_site$name)
select_site$name <- gsub("third", "66_68", select_site$name)

## table(select_site$name)
#27_28 38_40    55 66_68 
#   21    16     7    10 

melt_site <- melt(select_site[,c("name", "hela_level", "c666_level", "raji_level", "elijah_level")], id.vars="name")
melt_site[melt_site$value>1,]$value <- 1
melt_site[melt_site$value<0,]$value <- 0

pdf("plot/mt_site_pos_cell_line.pdf",width=8,height=3)
ggplot(melt_site, aes(x = name, y=value*100, fill = variable)) + geom_boxplot(outlier.shape=NA,width=0.6,position=position_dodge(width=0.7))+geom_point(size=2,pch = 21, position = position_jitterdodge(jitter.width =0.2, dodge.width=0.70))+scale_fill_manual(values=c("#77B661","#3B83C5","#F36D21", "#FFD504"))+
theme_bw()
dev.off()


#############mt comparision with ct-tRNA pseudouridine level
tRNA_site <- read.table("/users/ludwig/ebu571/ebu571/27Oct2023/site/snoRNA_tRNA_called_sites_fdrpassed_tRNA_filtered.txt",header=T, sep="\t",colClasses=c("ref"="character"))
tRNA_site <- tRNA_site[,c("chr", "treat1_level", "treat2_level")]
tRNA_site$mean_level <- (tRNA_site$treat1_level + tRNA_site$treat2_level)/2
tRNA_site$name <- c("ct_tRNA")
tRNA_site[tRNA_site$mean_level>=1,]$mean_level <- 1

site <- read.table("/users/ludwig/ebu571/ebu571/27Oct2023/site/snoRNA_tRNA_called_sites.txt",header=TRUE, sep=" ",colClasses=c("ref"="character")) ###1100
mtsite <- site[grepl("mt", site$chr),] ##54
mtsite <- mtsite[,c("chr", "treat1_level", "treat2_level")]
mtsite$mean_level <- (mtsite$treat1_level+mtsite$treat2_level)/2
mtsite$name <- c("mt_tRNA")

finaldat <- rbind(tRNA_site[,c("name", "mean_level")],mtsite[,c("name", "mean_level")])


p <- ggplot(finaldat,aes(x=name,y=mean_level*100, color=name)) + geom_point(aes(fill = name), size = 1, shape = 21, position = position_jitterdodge()) + geom_boxplot(width=0.6, alpha=0.2)+#geom_jitter(color="black", size=0.4) +##outlier.shape = NA, width=0.8 alpha=0.9
  theme_light() + xlab("") +  theme(axis.text = element_text(color = "black")) + #ylim(0,100) +
  ylab("level") +
  scale_fill_manual(values=c("#67A9CF","#EF8A62"))+
  scale_color_manual(values=c("#67A9CF","#EF8A62"))+
  theme(legend.position = "none") 
ggsave("plot/modification_level_cttrna_vsmttrna.pdf",p, width=3, height = 3) ##geom_boxplot(outlier.shape = NA, width=0.8)









################################################heatmap on all A of mt-tRNA and cyto-tRNA
input1 <- read.table("/users/ludwig/ebu571/ebu571/27Oct2023/align/H9_2023Oct26_S9.tRNA.filter.sort_A.mpile.txt",header=F, sep="\t",colClasses=c("V3"="character")) 
input2 <- read.table("/users/ludwig/ebu571/ebu571/27Oct2023/align/H10_2023Oct26_S10.tRNA.filter.sort_A.mpile.txt",header=F, sep="\t",colClasses=c("V3"="character")) 
treat1 <- read.table("/users/ludwig/ebu571/ebu571/27Oct2023/align/H7_2023Oct26_S7.tRNA.filter.sort_A.mpile.txt",header=F, sep="\t",colClasses=c("V3"="character")) 
treat2 <- read.table("/users/ludwig/ebu571/ebu571/27Oct2023/align/H8_2023Oct26_S8.tRNA.filter.sort_A.mpile.txt",header=F, sep="\t",colClasses=c("V3"="character")) 


mergelist <- list(input1[grepl("tRNA", input1$V1),c("V1", "V2", "V3", "V4", "V5")], input2[grepl("tRNA", input2$V1),c("V1", "V2", "V3", "V4", "V5")], treat1[grepl("tRNA", treat1$V1),c("V1", "V2", "V3", "V4", "V5")], treat2[grepl("tRNA", treat2$V1),c("V1", "V2", "V3", "V4", "V5")])
merge <- mergelist %>% reduce(full_join, by=c("V1", "V2", "V3")) %>% na.omit()
colnames(merge) <- c("chr", "pos", "ref", "input1_depth", "input1_mis", "input2_depth", "input2_mis", "treat1_depth", "treat1_mis", "treat2_depth", "treat2_mis")

merge$delta_mis1 <- merge$treat1_mis/merge$treat1_depth -(merge$input1_mis/merge$input1_depth)
merge$delta_mis2 <- merge$treat2_mis/merge$treat2_depth - (merge$input2_mis/merge$input2_depth)
merge$mean <- (merge$delta_mis1+merge$delta_mis2)/2
write.table(merge,"site/tRNA_m1A_site.txt",sep='\t',quote=F,row.names=F)

#merge <- read.table("site/tRNA_m1A_site.txt",header=TRUE)
tRNA_index <- read.table("resource/tRNA_position_index.txt", header=T, sep="\t",colClasses="character")
colnames(tRNA_index) <- c("chr", "trna_pos", "ref", "pos")

trna_A <- merge(tRNA_index, merge, by.x=c("chr", "pos"), by.y=c("chr", "pos"), all.x=T)
trna_A[is.na(trna_A)] <- 0

cyA <- trna_A[!grepl("mt", trna_A$chr),c("chr","trna_pos", "mean")]
cyA$trna_pos <- factor(cyA$trna_pos, level=c("0","1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","17A","18","19","20","20A","20B","21","22","23","24","25","26","27","28","29","30","31","32","33","34","35","36","37","38","39","40","41","42","43","44","45","e11","e12","e13","e14","e15","e16","e17","e1","e2","e3","e4","e5","e27","e26","e25","e24","e23","e22","e21","46","47","48","49","50","51","52","53","54","55","56","57","58","59","60","61","62","63","64","65","66","67","68","69","70","71","72","73","74","75","76"))

mtA <- trna_A[grepl("mt", trna_A$chr) & !grepl("e",trna_A$trna_pos),c("chr","trna_pos", "mean")]
mtA$trna_pos <- factor(mtA$trna_pos, level=c("0","1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","17A","18","19","20","20A","20B","21","22","23","24","25","26","27","28","29","30","31","32","33","34","35","36","37","38","39","40","41","42","43","44","45","46","47","48","49","50","51","52","53","54","55","56","57","58","59","60","61","62","63","64","65","66","67","68","69","70","71","72","73","74","75","76"))

###################only plot representatives
seltrna <- read.table("resource/representative_tRNA.txt", header=F, sep="\t",colClasses="character")
colnames(seltrna) <- c("chr")
selA <- merge(cyA, seltrna, by=c("chr"))


pdf("plot/tRNA_A_misincorporation_heatmap.pdf", width=80, height=80)
ggplot(cyA, aes(x = trna_pos, y = chr, fill = mean*100)) + geom_tile(color = "white", size=0.01)+ coord_fixed()+#+ ##to make Square tiles 
  #geom_text(aes(label = delta_conversion), color = "black", size = 6) +
  scale_fill_gradient2(low = "#075AFF",mid = "#FFFFCC",high = "#FF0000", guide = guide_colorbar(frame.colour = "black",ticks = TRUE,))+ theme(panel.border=element_rect(fill = NA, colour='black',size=1))+theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()

pdf("plot/tRNA_A_misincorporation_heatmap_representative.pdf", width=80, height=80)
ggplot(selA, aes(x = trna_pos, y = chr, fill = mean*100)) + geom_tile(color = "white", size=0.01)+ coord_fixed()+#+ ##to make Square tiles 
  #geom_text(aes(label = delta_conversion), color = "black", size = 6) +
  scale_fill_gradient2(low = "#075AFF",mid = "#FFFFCC",high = "#FF0000", guide = guide_colorbar(frame.colour = "black",ticks = TRUE,))+ theme(panel.border=element_rect(fill = NA, colour='black',size=1))+theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()

pdf("plot/mtRNA_A_misincorporation_heatmap.pdf", width=20, height=10)
ggplot(mtA, aes(x = trna_pos, y = chr, fill = mean*100)) + geom_tile(color = "white", size=0.01)+ coord_flip()+ #+ ##to make Square tiles 
  #geom_text(aes(label = delta_conversion), color = "black", size = 6) +
  scale_fill_gradient2(low = "blue",mid = "white",high = "red", limits=c(-60,60), guide = guide_colorbar(frame.colour = "black",ticks = TRUE,))+ theme(panel.border=element_rect(fill = NA, colour='black',size=1))+theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()



#####boxplot for position 55 for tRNA and pos 9 for mt-tRNA
trna_A <- merge(tRNA_index, merge, by.x=c("chr", "pos"), by.y=c("chr", "pos"))
pos58 <- trna_A[grepl("tRNA", trna_A$chr) & !grepl("mt", trna_A$chr) & trna_A$trna_pos=="58",] ###201  ##199
pos58$input_mis <- ((pos58$input1_mis/pos58$input1_depth) + (pos58$input2_mis/pos58$input2_depth))/2
pos58$treat_mis <- ((pos58$treat1_mis/pos58$treat1_depth)+(pos58$treat2_mis/pos58$treat2_depth))/2

filteredlist <- read.table("/users/ludwig/ebu571/ebu571/27Oct2023/featureCounts/tRNA_counts_filtered.txt", header=T, sep="\t")

filtered_pos58 <- merge(filteredlist[,c("chr", "anticodon")], pos58,by=c("chr")) ##192

pos58 <- filtered_pos58[,c("trna_pos", "input_mis", "treat_mis")]
pos9 <- trna_A[grepl("mt", trna_A$chr) & trna_A$trna_pos=="9",] ##22 ###14
pos9$input_mis <- ((pos9$input1_mis/pos9$input1_depth) + (pos9$input2_mis/pos9$input2_depth))/2
pos9$treat_mis <- ((pos9$treat1_mis/pos9$treat1_depth)+(pos9$treat2_mis/pos9$treat2_depth))/2

pos58 <- rbind(pos58[pos58$input_mis>=0.20,],pos9[,c("trna_pos", "input_mis", "treat_mis")]) ###223


pos58[!is.na(pos58$input_mis),] %>% nrow() ###213 ###184

 #table(pos58[!is.na(pos58$input_mis),]$trna_pos)
 # 58   9 
##170  14
melt58 <- melt(pos58[!is.na(pos58$input_mis),], id.vars=c("trna_pos"))
write.table(pos58,"site/tRNA_m1A_test.txt",sep='\t',quote=F,row.names=F)

p <- ggplot(melt58,aes(x=trna_pos,y=value*100, fill=variable)) +geom_boxplot(width=0.8, outlier.shape = NA)+# geom_point(aes(fill = variable), size = 1, shape = 21, position = position_jitterdodge()) #geom_jitter(color="black", size=0.4) +##outlier.shape = NA, width=0.8 alpha=0.9
  theme_light() + xlab("") +  theme(axis.text = element_text(color = "black")) + ylim(0,100) +
  ylab("level") +
  scale_fill_manual(values=c("#1b9e77","#7570b3"))
#theme(legend.position = "none") 
ggsave("plot/m1A_misincorporation_boxplot.pdf",p, width=4, height = 4) ##geom_boxplot(outlier.shape = NA, width=0.8)

pdf("plot/pos58.histagram.pdf")
ggplot(pos58[pos58$trna_pos==58,], aes(x=input_mis)) + geom_histogram()
dev.off()


###pie chart for mt
dat <- data.frame(c("A", "B", "C"),c(8,13,13))
colnames(dat) <- c("name", "n")
nb.cols <- 3
mycolors <- colorRampPalette(brewer.pal(3, "RdYlBu"))(nb.cols)

pdf("plot/mtRNA_pie.pdf")
ggplot(dat, aes(x="", y=n, fill=name))+
geom_bar(width = 1, stat = "identity")+coord_polar("y") + scale_fill_manual(values = mycolors) +  ###or RdYlBu
  geom_text(aes(label = n),position = position_stack(vjust = 0.5))+
  theme_void()
dev.off()




###orgsnize tRNA and mt-tRNA data
site <- read.table("/users/ludwig/ebu571/ebu571/27Oct2023/site/snoRNA_tRNA_called_sites.txt",header=TRUE, sep=" ",colClasses=c("ref"="character")) 
mtsite <- site[grepl("mt", site$chr),] ##54

tRNA_site <- read.table("/users/ludwig/ebu571/ebu571/27Oct2023/site/snoRNA_tRNA_called_sites_fdrpassed_tRNA_filtered.txt",header=T, sep="\t",colClasses=c("ref"="character")) ###609

tRNA_index <- read.table("resource/tRNA_position_index.txt", header=T, sep="\t",colClasses="character")
colnames(tRNA_index) <- c("chr", "trna_pos", "ref", "pos")

cytrna <- merge(tRNA_index, tRNA_site, by.x=c("chr", "pos"), by.y=c("chr", "end"), all.y=TRUE) ##609
mttrna <- merge(tRNA_index, mtsite, by.x=c("chr", "pos"), by.y=c("chr", "end"), all.y=TRUE) ##54
write.table(cytrna,"site/hela_cy_trna_site.txt",sep='\t',quote=F,row.names=F)
write.table(mttrna,"site/hela_mt_trna_site.txt",sep='\t',quote=F,row.names=F)