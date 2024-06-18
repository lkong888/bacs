#module load R/4.0.3-foss-2020b 
library(GenomicAlignments)
library(dplyr)
library(reshape2)
library(tidyr)
library(tidyverse)

###Site in Hela cell
allsite1 <- read.table("/gpfs3/well/ludwig/users/ebu571/27Oct2023/site/treat1_input1_smallrna_allsites.txt", header=T, sep=" ", colClasses=c("ref"="character") ) 
allsite2 <- read.table("/gpfs3/well/ludwig/users/ebu571/27Oct2023/site/treat2_input2_smallrna_allsites.txt", header=T, sep=" ", colClasses=c("ref"="character") ) 

merge1 <- list(allsite1, allsite2)
merge1 <- merge1 %>% reduce(full_join, by=c("chr", "start", "end", "ref", "motif")) %>% na.omit()
colnames(merge1) <- c("motif","chr", "start", "end", "ref","treat1_depth","treat1_T","treat1_C", "treat1_conversion","treat1_gap","treat1_gap_rate","ctrl1_depth","ctrl1_T","ctrl1_C","ctrl1_conversion","ctrl1_gap","ctrl1_gap_rate","conversion1", "treat1_fp", "treat1_level","treat1_p", "treat1_fdr","treat2_depth","treat2_T","treat2_C", "treat2_conversion","treat2_gap","treat2_gap_rate","ctrl2_depth","ctrl2_T","ctrl2_C","ctrl2_conversion","ctrl2_gap","ctrl2_gap_rate","conversion2", "treat2_fp", "treat2_level","treat2_p", "treat2_fdr")
nrow(merge1) #14967

high_confident_site <- merge1[(merge1$treat1_T+merge1$treat1_C)>=20 &(merge1$treat2_T+merge1$treat2_C)>=20 & merge1$treat1_level>=0.05 & merge1$treat2_level>=0.05 & merge1$treat1_fdr<=0.001 & merge1$treat2_fdr<=0.001, ]
high_confident_site <- high_confident_site[(high_confident_site$ctrl1_conversion<=0.01 |high_confident_site$ctrl1_C<=2) &(high_confident_site$ctrl2_conversion<=0.01 |high_confident_site$ctrl2_C<=2), ]
nrow(high_confident_site) ##1093

tRNA_index <- read.table("/users/ludwig/ebu571/ebu571/20May2023_final/resource/tRNA_position_index.txt", header=T, sep="\t",colClasses="character")
colnames(tRNA_index) <- c("chr", "trna_pos", "ref", "end")
merge1 <- merge(tRNA_index[,c("chr", "trna_pos", "end")], merge1, by.x=c("chr", "end"), by.y=c("chr", "end"), all.y=TRUE)
merge1[!grepl("tRNA", merge1$chr),]$trna_pos <- merge1[!grepl("tRNA", merge1$chr),]$end



####read KO data for three TRUB1 KO cell lines
##/users/ludwig/ebu571/ebu571/20Nov2023/ko_sample/treat1_input1_rrna_smallrna_allT_motif.bed for g2_9
##/users/ludwig/ebu571/ebu571/20Nov2023/ko_sample/treat2_input2_rrna_smallrna_allT_motif.bed for g2_11
ko1 <- read.table("/users/ludwig/ebu571/ebu571/20Nov2023/ko_sample/treat1_input1_rrna_smallrna_allsites.txt", header=T, sep=" ", colClasses=c("ref"="character") ) 
ko2 <- read.table("/users/ludwig/ebu571/ebu571/20Nov2023/ko_sample/treat2_input2_rrna_smallrna_allsites.txt", header=T, sep=" ", colClasses=c("ref"="character") ) 
ko3 <- read.table("/users/ludwig/ebu571/ebu571/14Oct2023/ko_sample/treat2_input2_rrna_smallrna_allsites.txt", header=T, sep=" ", colClasses=c("ref"="character") ) 

ko_trub1 <- list(ko1, ko2, ko3)
ko_trub1 <- ko_trub1 %>% reduce(full_join, by=c("chr", "start", "end", "ref", "motif")) %>% na.omit()
colnames(ko_trub1) <- c("motif","chr", "start", "end", "ref","treat1_depth","treat1_T","treat1_C", "treat1_conversion","treat1_gap","treat1_gap_rate","ctrl1_depth","ctrl1_T","ctrl1_C","ctrl1_conversion","ctrl1_gap","ctrl1_gap_rate","conversion1", "treat1_fp", "treat1_level","treat1_p", "treat1_fdr","treat2_depth","treat2_T","treat2_C", "treat2_conversion","treat2_gap","treat2_gap_rate","ctrl2_depth","ctrl2_T","ctrl2_C","ctrl2_conversion","ctrl2_gap","ctrl2_gap_rate","conversion2", "treat2_fp", "treat2_level","treat2_p", "treat2_fdr","treat3_depth","treat3_T","treat3_C", "treat3_conversion","treat3_gap","treat3_gap_rate","ctrl3_depth","ctrl3_T","ctrl3_C","ctrl3_conversion","ctrl3_gap","ctrl3_gap_rate","conversion1", "treat3_fp", "treat3_level","treat3_p", "treat3_fdr")
nrow(ko_trub1) #14027

wt_ko_trub <- merge(merge1[,c("chr", "start", "end","trna_pos", "ref","motif", "treat1_level", "treat2_level")], ko_trub1[,c("chr", "start", "end", "treat1_level", "treat2_level","treat3_level")], by=c("chr", "start", "end"))
colnames(wt_ko_trub) <- c("chr", "start", "end", "trna_pos","ref","motif", "wt1_level", "wt2_level", "ko1_level", "ko2_level","ko3_level")###11658
wt_ko_trub$wt_mean <- (wt_ko_trub$wt1_level+wt_ko_trub$wt2_level)/2
wt_ko_trub$ko_mean <- (wt_ko_trub$ko1_level+wt_ko_trub$ko2_level+wt_ko_trub$ko3_level)/3

wt_ko_trub[grepl("mt_tRNA", wt_ko_trub$chr),c("group")] <- "mt_tRNA"
wt_ko_trub[grepl("^tRNA", wt_ko_trub$chr),c("group")] <- "tRNA"
wt_ko_trub[grepl("SNO", wt_ko_trub$chr) | grepl("SCA", wt_ko_trub$chr) ,c("group")] <- "snoRNA"
write.table(wt_ko_trub,"plot/wt_trub1ko_table.txt",sep='\t',quote=F,row.names=F)
for (i in 1:nrow(wt_ko_trub)){
    if(abs(wt_ko_trub[i,]$wt_mean-wt_ko_trub[i,]$ko_mean)>=0.05 & wt_ko_trub[i,]$ko_mean/wt_ko_trub[i,]$wt_mean>=2){
        wt_ko_trub$type[i] <- "up"
    }else if(abs(wt_ko_trub[i,]$wt_mean-wt_ko_trub[i,]$ko_mean)>=0.05 & wt_ko_trub[i,]$ko_mean/wt_ko_trub[i,]$wt_mean<=1/2){
        wt_ko_trub$type[i] <- "down"
    }else{
        wt_ko_trub$type[i] <- "nochanged"
    }
}



wt_ko_trub[wt_ko_trub$wt_mean>=1,]$wt_mean <- 1
wt_ko_trub[wt_ko_trub$wt_mean<0,]$wt_mean <- 0
wt_ko_trub[wt_ko_trub$ko_mean>=1,]$ko_mean <- 1
wt_ko_trub[wt_ko_trub$ko_mean<0,]$ko_mean <- 0


###all site in mt
###all site in cytosolic tRNA
mt_site1 <- wt_ko_trub[wt_ko_trub$group=="mt_tRNA",]
mt_site1$label <- paste(mt_site1$chr, mt_site1$trna_pos, sep="_")
mt_site1 <- merge(mt_site1,high_confident_site[grepl("mt_tRNA", high_confident_site$chr),c("chr", "start", "end")], by=c("chr", "start", "end")) ###54
###all sites in cy-tRNA
trna_site <- read.table("/users/ludwig/ebu571/ebu571/27Oct2023/site/snoRNA_tRNA_called_sites_fdrpassed_tRNA_filtered.txt", header=T, sep="\t", colClasses=c("ref"="character") )
ct_site1 <- wt_ko_trub[wt_ko_trub$group=="tRNA" & !is.na(wt_ko_trub$trna_pos),] ##2741
ct_site1$label <-  paste(ct_site1$chr, ct_site1$trna_pos, sep="_")
ct_site1 <- merge(ct_site1,trna_site[grepl("^tRNA", trna_site$chr),c("chr", "start", "end")], by=c("chr", "start", "end")) ###541




##############################################################pus7

ko4 <- read.table("/users/ludwig/ebu571/ebu571/07Nov2023/site/treat8_input8_rrna_smallrna_allsites.txt", header=T, sep=" ", colClasses=c("ref"="character") ) 
ko5 <- read.table("/users/ludwig/ebu571/ebu571/07Nov2023/site/treat9_input9_rrna_smallrna_allsites.txt", header=T, sep=" ", colClasses=c("ref"="character") ) 

ko_pus7 <- list(ko4, ko5)
ko_pus7 <- ko_pus7 %>% reduce(full_join, by=c("chr", "start", "end", "ref", "motif")) %>% na.omit()
colnames(ko_pus7) <- c("motif","chr", "start", "end", "ref","treat1_depth","treat1_T","treat1_C", "treat1_conversion","treat1_gap","treat1_gap_rate","ctrl1_depth","ctrl1_T","ctrl1_C","ctrl1_conversion","ctrl1_gap","ctrl1_gap_rate","conversion1", "treat1_fp", "treat1_level","treat1_p", "treat1_fdr","treat2_depth","treat2_T","treat2_C", "treat2_conversion","treat2_gap","treat2_gap_rate","ctrl2_depth","ctrl2_T","ctrl2_C","ctrl2_conversion","ctrl2_gap","ctrl2_gap_rate","conversion2", "treat2_fp", "treat2_level","treat2_p", "treat2_fdr")
nrow(ko_pus7) #15961

wt_ko_pus7 <- merge(merge1[,c("chr", "start", "end","trna_pos", "ref","motif", "treat1_level", "treat2_level")], ko_pus7[,c("chr", "start", "end", "treat1_level", "treat2_level")], by=c("chr", "start", "end"))
colnames(wt_ko_pus7) <- c("chr", "start", "end", "trna_pos","ref","motif", "wt1_level", "wt2_level", "ko1_level", "ko2_level")###11658
wt_ko_pus7$wt_mean <- (wt_ko_pus7$wt1_level+wt_ko_pus7$wt2_level)/2
wt_ko_pus7$ko_mean <- (wt_ko_pus7$ko1_level+wt_ko_pus7$ko2_level)/2

for (i in 1:nrow(wt_ko_pus7)){
    if(abs(wt_ko_pus7[i,]$wt_mean-wt_ko_pus7[i,]$ko_mean)>=0.05 & wt_ko_pus7[i,]$ko_mean/wt_ko_pus7[i,]$wt_mean>=2){
        wt_ko_pus7$type[i] <- "up"
    }else if(abs(wt_ko_pus7[i,]$wt_mean-wt_ko_pus7[i,]$ko_mean)>=0.05 & wt_ko_pus7[i,]$ko_mean/wt_ko_pus7[i,]$wt_mean<=1/2){
        wt_ko_pus7$type[i] <- "down"
    }else{
        wt_ko_pus7$type[i] <- "nochanged"
    }
}

wt_ko_pus7[grepl("mt_tRNA", wt_ko_pus7$chr),c("group")] <- "mt_tRNA"
wt_ko_pus7[grepl("^tRNA", wt_ko_pus7$chr),c("group")] <- "tRNA"
wt_ko_pus7[grepl("SNO", wt_ko_pus7$chr) | grepl("SCA", wt_ko_pus7$chr) ,c("group")] <- "snoRNA"
write.table(wt_ko_pus7,"plot/wt_pus7ko_table.txt",sep='\t',quote=F,row.names=F)

wt_ko_pus7[wt_ko_pus7$wt_mean>=1,]$wt_mean <- 1
wt_ko_pus7[wt_ko_pus7$wt_mean<0,]$wt_mean <- 0
wt_ko_pus7[wt_ko_pus7$ko_mean>=1,]$ko_mean <- 1
wt_ko_pus7[wt_ko_pus7$ko_mean<0,]$ko_mean <- 0


###plot heatmap
###all site in mt
###all site in cytosolic tRNA
mt_site2 <- wt_ko_pus7[wt_ko_pus7$group=="mt_tRNA",]
mt_site2$label <-  paste(mt_site2$chr, mt_site2$trna_pos, sep="_")
mt_site2 <- merge(mt_site2,high_confident_site[grepl("mt_tRNA", high_confident_site$chr),c("chr", "start", "end")], by=c("chr", "start", "end")) ###54
mt_site2$type <- factor(mt_site2$type, level=c("down", "up", "nochanged"))


ct_site2 <- wt_ko_pus7[wt_ko_pus7$group=="tRNA" & !is.na(wt_ko_pus7$trna_pos),] ##2741
ct_site2$label <-  paste(ct_site2$chr, ct_site2$trna_pos, sep="_")
ct_site2 <- merge(ct_site2,trna_site[grepl("^tRNA", trna_site$chr),c("chr", "start", "end")], by=c("chr", "start", "end")) ###586
ct_site2$type <- factor(ct_site2$type, level=c("down", "up", "nochanged"))


###############pus1
ko6 <- read.table("/users/ludwig/ebu571/ebu571/07Nov2023/site/treat10_input10_rrna_smallrna_allsites.txt", header=T, sep=" ", colClasses=c("ref"="character") ) 
ko7 <- read.table("/users/ludwig/ebu571/ebu571/07Nov2023/site/treat11_input11_rrna_smallrna_allsites.txt", header=T, sep=" ", colClasses=c("ref"="character") ) 

ko_pus1 <- list(ko6, ko7)
ko_pus1 <- ko_pus1 %>% reduce(full_join, by=c("chr", "start", "end", "ref", "motif")) %>% na.omit()
colnames(ko_pus1) <- c("motif","chr", "start", "end", "ref","treat1_depth","treat1_T","treat1_C", "treat1_conversion","treat1_gap","treat1_gap_rate","ctrl1_depth","ctrl1_T","ctrl1_C","ctrl1_conversion","ctrl1_gap","ctrl1_gap_rate","conversion1", "treat1_fp", "treat1_level","treat1_p", "treat1_fdr","treat2_depth","treat2_T","treat2_C", "treat2_conversion","treat2_gap","treat2_gap_rate","ctrl2_depth","ctrl2_T","ctrl2_C","ctrl2_conversion","ctrl2_gap","ctrl2_gap_rate","conversion2", "treat2_fp", "treat2_level","treat2_p", "treat2_fdr")
nrow(ko_pus1) #15690

wt_ko_pus1 <- merge(merge1[,c("chr", "start", "end","trna_pos", "ref","motif", "treat1_level", "treat2_level")], ko_pus1[,c("chr", "start", "end", "treat1_level", "treat2_level")], by=c("chr", "start", "end"))
colnames(wt_ko_pus1) <- c("chr", "start", "end", "trna_pos","ref","motif", "wt1_level", "wt2_level", "ko1_level", "ko2_level")###11658
wt_ko_pus1$wt_mean <- (wt_ko_pus1$wt1_level+wt_ko_pus1$wt2_level)/2
wt_ko_pus1$ko_mean <- (wt_ko_pus1$ko1_level+wt_ko_pus1$ko2_level)/2

for (i in 1:nrow(wt_ko_pus1)){
    if(abs(wt_ko_pus1[i,]$wt_mean-wt_ko_pus1[i,]$ko_mean)>=0.05 & wt_ko_pus1[i,]$ko_mean/wt_ko_pus1[i,]$wt_mean>=2){
        wt_ko_pus1$type[i] <- "up"
    }else if(abs(wt_ko_pus1[i,]$wt_mean-wt_ko_pus1[i,]$ko_mean)>=0.05 & wt_ko_pus1[i,]$ko_mean/wt_ko_pus1[i,]$wt_mean<=1/2){
        wt_ko_pus1$type[i] <- "down"
    }else{
        wt_ko_pus1$type[i] <- "nochanged"
    }
}

wt_ko_pus1[grepl("mt_tRNA", wt_ko_pus1$chr),c("group")] <- "mt_tRNA"
wt_ko_pus1[grepl("^tRNA", wt_ko_pus1$chr),c("group")] <- "tRNA"
wt_ko_pus1[grepl("SNO", wt_ko_pus1$chr) | grepl("SCA", wt_ko_pus1$chr) ,c("group")] <- "snoRNA"
write.table(wt_ko_pus1,"plot/wt_pus1ko_table.txt",sep='\t',quote=F,row.names=F)

wt_ko_pus1[wt_ko_pus1$wt_mean>=1,]$wt_mean <- 1
wt_ko_pus1[wt_ko_pus1$wt_mean<0,]$wt_mean <- 0
wt_ko_pus1[wt_ko_pus1$ko_mean>=1,]$ko_mean <- 1
wt_ko_pus1[wt_ko_pus1$ko_mean<0,]$ko_mean <- 0

###all site in mt
###all site in cytosolic tRNA
mt_site3 <- wt_ko_pus1[wt_ko_pus1$group=="mt_tRNA",]
mt_site3$label <-  paste(mt_site3$chr, mt_site3$trna_pos, sep="_")
mt_site3 <- merge(mt_site3,high_confident_site[grepl("mt_tRNA", high_confident_site$chr),c("chr", "start", "end")], by=c("chr", "start", "end")) ###54
mt_site3$type <- factor(mt_site3$type, level=c("down", "up", "nochanged"))

ct_site3 <- wt_ko_pus1[wt_ko_pus1$group=="tRNA" & !is.na(wt_ko_pus1$trna_pos),] ##2741
ct_site3$label <-  paste(ct_site3$chr, ct_site3$trna_pos, sep="_")
ct_site3 <- merge(ct_site3,trna_site[grepl("^tRNA", trna_site$chr),c("chr", "start", "end")], by=c("chr", "start", "end")) ###611
ct_site3$type <- factor(ct_site3$type, level=c("down", "up", "nochanged"))

####plot for mt-tRNA 55, 27/28, 66/67/68
###merge all mt-tRNA site
####


mt_merge <- list(mt_site1, mt_site2,mt_site3)
mt_merge <- mt_merge %>% reduce(full_join, by=c("chr", "start", "trna_pos","end", "ref", "motif", "group", "label")) %>% na.omit() %>% dplyr::select(chr,start,end,trna_pos,ref,motif,group,label,wt1_level.x,wt2_level.x,wt_mean.x,ko1_level.x,ko2_level.x,ko3_level,ko_mean.x,type.x,ko1_level.y,ko2_level.y,ko_mean.y,type.y,ko1_level,ko2_level,ko_mean,type)
colnames(mt_merge) <- c("chr", "start", "end", "trna_pos","ref","motif","group","label","wt1_level", "wt2_level","wt_mean","trub1_ko1","trub1_ko2","trub1_ko3","trub1_ko_mean","trub1_type","pus7_ko1","pus7_ko2","pus7_ko_mean", "pus7_ko_type", "pus1_ko1","pus1_ko2", "pus1_ko_mean","pus1_ko_type")
nrow(mt_merge) #54
write.table(mt_merge,"plot/mt_site.txt",sep='\t',quote=F,row.names=F)


library(ggpubr)
seldat <- melt(mt_merge[mt_merge$trna_pos==55 |mt_merge$trna_pos==27 | mt_merge$trna_pos==28 | mt_merge$trna_pos==66 | mt_merge$trna_pos==67 | mt_merge$trna_pos==68 ,c("trna_pos","wt_mean", "trub1_ko_mean", "pus7_ko_mean", "pus1_ko_mean")], id.vars=c("trna_pos"))
seldat[seldat$trna_pos==55,c("xlabel")] <- "pos55"
seldat[seldat$trna_pos==27 | seldat$trna_pos==28,c("xlabel")] <- "pos27/28"
seldat[seldat$trna_pos==66 | seldat$trna_pos==67 | seldat$trna_pos==68,c("xlabel")] <- "pos66/67/68"


ct_merge <- list(ct_site1, ct_site2,ct_site3)
ct_merge <- ct_merge %>% reduce(full_join, by=c("chr", "start", "trna_pos","end", "ref", "motif", "group", "label"))
write.table(ct_merge,"plot/ct_site.txt",sep='\t',quote=F,row.names=F)


######plot 55 position
mt_55 <-seldat[seldat$xlabel=="pos55" & (seldat$variable=="wt_mean" | seldat$variable=="trub1_ko_mean"),]
mt_55$xlabel <- "mt_55"
ct_55 <- melt(ct_site1[ct_site1$trna_pos=="55",c("trna_pos", "wt_mean", "ko_mean")], id.vars=c("trna_pos"))
ct_55$variable <- gsub("ko_mean", "trub1_ko_mean",ct_55$variable)
ct_55$xlabel <- "cy_pos55"
pos55 <- rbind(mt_55, ct_55)
pos55$variable <- factor(pos55$variable, level=c("wt_mean", "trub1_ko_mean"))

pdf("plot/pos55_trub1_wt_ko_boxplot.pdf",width=4,height=4)
ggplot(pos55, aes(x = xlabel, y=value*100, color=variable,fill = variable)) + geom_boxplot(outlier.shape=NA,width = 0.5,position=position_dodge(width=0.7))+geom_point(pch = 21, position = position_jitterdodge(jitter.width =0.3, dodge.width=0.7))+scale_color_manual(values=c('#404040','#006699'))+scale_fill_manual(values=c('#DBDBDB','#E5F4FF'))+
#stat_summary(fun = mean, geom="text", aes(label = round(..y.., 3)),vjust=-12)+
theme_bw()+theme(legend.position="none")#+stat_compare_means(aes(group = variable), label = "p.signif",method = "t.test")
dev.off()


####make bar graph
agg_55 <- aggregate(value ~ xlabel+variable,pos55, mean)
pdf("plot/pos55_trub1_wt_ko_barplot.pdf",width=4,height=4)
ggplot(agg_55, aes(x = xlabel, y=value*100, fill = variable))+geom_point(data=pos55,aes(x=xlabel, y=value*100),pch = 21, position = position_jitterdodge())+scale_fill_manual(values=c("grey40","#006699"))+
#stat_summary(fun = mean, geom="text", aes(label = round(..y.., 3)),vjust=-12)+ + geom_bar(stat="identity", position=position_dodge())
theme_bw()+theme(legend.position="none")#+stat_compare_means(aes(group = variable), label = "p.signif",method = "t.test")
dev.off()



###plot ct_27/28 mt_27/28, 
ct_pos1_pos <- melt(ct_site3[(ct_site3$trna_pos=="27" | ct_site3$trna_pos=="28"),c("trna_pos", "wt_mean", "ko_mean")],id.vars=c("trna_pos"))
ct_pos1_pos$xlabel <- "ct_27/28"
ct_pos1_pos$variable <- gsub("ko_mean", "pus1_ko_mean",ct_pos1_pos$variable)
mt_pos1_pos <-seldat[(seldat$xlabel=="pos27/28" | seldat$xlabel=="pos66/67/68") & (seldat$variable=="wt_mean" | seldat$variable=="pus1_ko_mean"),]
pus1_pos <- rbind(mt_pos1_pos, ct_pos1_pos)
pus1_pos$variable <- factor(pus1_pos$variable, level=c("wt_mean", "pus1_ko_mean"))

pdf("plot/pus1_wt_ko_boxplot.pdf",width=5,height=4)
ggplot(pus1_pos, aes(x = xlabel, y=value*100,color=variable, fill = variable)) + geom_boxplot(outlier.shape=NA, width = 0.5, position=position_dodge(width=0.7))+geom_point(pch = 21, position = position_jitterdodge(jitter.width =0.25, dodge.width=0.7))+scale_color_manual(values=c('#404040','#228B22'))+scale_fill_manual(values=c("#DBDBDB","#CBFFCB"))+
#stat_summary(fun = mean, geom="text", aes(label = round(..y.., 3)),vjust=-12)+
theme_bw()+theme(legend.position="none")#+stat_compare_means(aes(group = variable), label = "p.signif",method = "t.test")
dev.off()


#####plot site for PUS7 
pos_pus7 <-melt(ct_site2[(ct_site2$trna_pos=="13" | ct_site2$trna_pos=="20B" | ct_site2$trna_pos=="35" | ct_site2$trna_pos=="36" | ct_site2$trna_pos=="50"),c("trna_pos", "wt_mean", "ko_mean")],id.vars=c("trna_pos"))
pos_pus7$xlabel <- pos_pus7$trna_pos


pus7_dat <-ct_site2[ct_site2$trna_pos=="50",]
mt_melt1 <- melt(pus7_dat[,c("label", "wt1_level", "wt2_level", "ko1_level","ko2_level")],id.var=c("label"))
mt_melt1[mt_melt1$value>=1,]$value <- 1
mt_melt1[mt_melt1$value<0,]$value <- 0
mdf1 <- spread(mt_melt1, key = variable, value = value) %>%  as.data.frame()
final_matrix1 <- data.frame(mdf1[,-1], row.names=as.character(mdf1[,1])) %>% as.matrix(rownames = TRUE)


pdf("plot/pus7.wt_ko.heatmap.pdf", width=5, height=2.5)
pheatmap(final_matrix1,cluster_cols = FALSE, cluster_rows = FALSE)
dev.off()

pdf("plot/pus7_wt_ko_boxplot.pdf",,width=6,height=4)
ggplot(pos_pus7, aes(x = xlabel, y=value*100, color=variable,fill = variable)) + geom_boxplot(outlier.shape=NA, width = 0.5, position=position_dodge(width=0.7))+geom_point(pch = 21, position = position_jitterdodge(jitter.width =0.25, dodge.width=0.7))+scale_color_manual(values=c('#404040','#DC143C'))+scale_fill_manual(values=c("#DBDBDB","#FFBFCA"))+ylim(0,100)+
#stat_summary(fun = mean, geom="text", aes(label = round(..y.., 3)),vjust=-12)+
theme_bw()+theme(legend.position="none")#+stat_compare_means(aes(group = variable), label = "p.signif",method = "t.test")
dev.off()


p1 <- ggplot(pos55, aes(x = xlabel, y=value*100, color=variable,fill = variable)) + geom_boxplot(outlier.shape=NA,width = 0.5,position=position_dodge(width=0.7))+geom_point(pch = 21, position = position_jitterdodge(jitter.width =0.3, dodge.width=0.7))+scale_color_manual(values=c('#404040','#006699'))+scale_fill_manual(values=c('#DBDBDB','#E5F4FF'))+
#stat_summary(fun = mean, geom="text", aes(label = round(..y.., 3)),vjust=-12)+
theme_bw()
p2 <- ggplot(pus1_pos, aes(x = xlabel, y=value*100,color=variable, fill = variable)) + geom_boxplot(outlier.shape=NA, width = 0.5, position=position_dodge(width=0.7))+geom_point(pch = 21, position = position_jitterdodge(jitter.width =0.25, dodge.width=0.7))+scale_color_manual(values=c('#404040','#228B22'))+scale_fill_manual(values=c("#DBDBDB","#CBFFCB"))+
#stat_summary(fun = mean, geom="text", aes(label = round(..y.., 3)),vjust=-12)+
theme_bw()

p3 <- ggplot(pos_pus7, aes(x = xlabel, y=value*100, color=variable,fill = variable)) + geom_boxplot(outlier.shape=NA, width = 0.5, position=position_dodge(width=0.7))+geom_point(pch = 21, position = position_jitterdodge(jitter.width =0.25, dodge.width=0.7))+scale_color_manual(values=c('#404040','#DC143C'))+scale_fill_manual(values=c("#DBDBDB","#FFBFCA"))+ylim(0,100)+
#stat_summary(fun = mean, geom="text", aes(label = round(..y.., 3)),vjust=-12)+
theme_bw()

library("gridExtra")
pdf("plot/wt_ko_boxplot.pdf",width=8, height=4)
grid.arrange(p1,p2, p3,
        ncol=2, nrow=2)
dev.off()


######plot combining sites in mt_tRNA, tRNA and mRNA
trub1_mrna <- read.table("site/wt_trub1.site.txt", header=T, sep="\t",colClasses=c("ref.x"="character")) ##986
pus7_mrna <- read.table("site/wt_pus7.site.txt", header=T, sep="\t", colClasses=c("ref.x"="character") ) #656
pus1_mrna <- read.table("site/wt_pus1.site.txt", header=T, sep="\t", colClasses=c("ref.x"="character") ) ##590



###for trub1
###all snoRNA, tRNA, mt-tRNA, mRNA site
snosite <- read.table("/users/ludwig/ebu571/ebu571/27Oct2023/site/snoRNA_tRNA_called_sites.txt", header=TRUE, sep=" ",colClasses=c("ref"="character")) 
snosite <- snosite[!grepl("tRNA", snosite$chr),] ###304

sno_trub <- merge(snosite[,c("chr", "start", "end", "ref","motif", "treat1_level", "treat2_level")], ko_trub1[,c("chr", "start", "end", "treat1_level", "treat2_level","treat3_level")], by=c("chr", "start", "end"))
colnames(sno_trub) <- c("chr", "start", "end","ref","motif", "wt1_level", "wt2_level", "ko1_level", "ko2_level","ko3_level")###221
sno_trub$ko_mean <- (sno_trub$ko1_level+sno_trub$ko2_level+sno_trub$ko3_level)/3
sno_trub$wt_mean <- (sno_trub$wt1_level+sno_trub$wt2_level)/2

sno_trub[sno_trub$ko_mean<0.05,] ##3
for (i in 1:nrow(sno_trub)){
    if(sno_trub[i,]$ko1_level<=0.01 & sno_trub[i,]$ko2_level<=0.01& sno_trub[i,]$ko3_level<=0.01){
        sno_trub$change[i] <- "depleted"
    }else if(sno_trub[i,]$wt_mean-sno_trub[i,]$ko_mean>=0.20){
        sno_trub$change[i] <- "down"
    }else{
        sno_trub$change[i] <- "non"
    }
}
sno_trub_dat <- sno_trub[,c("chr", "end","wt_mean", "ko_mean", "change")]
colnames(sno_trub_dat) <-c("chr", "pos","wt_mean", "ko_mean","change")
sno_trub_dat$type <- "snoRNA"


##> table(sno_trub_dat$change)

non 
221 


##trub1_mRNA_site
##only retain site treat_T.x+treat_C.x>=20, ctrl_T.x+ctrl_C.x>=20, treat_level2.x>0.05, treat_C.x>=10
mrna_trub_dat <- trub1_mrna[(trub1_mrna$treat_T.x+trub1_mrna$treat_C.x)>=20 & (trub1_mrna$ctrl_T.x+trub1_mrna$ctrl_C.x)>=20 & trub1_mrna$treat_level2.x>=0.05 & trub1_mrna$treat_C.x>=10 & (trub1_mrna$treat_T.y+trub1_mrna$treat_C.y)>=20,]
nrow(mrna_trub_dat) ##309

for (i in 1:nrow(mrna_trub_dat)){
    if(mrna_trub_dat[i,]$treat_level2.y<=0.01){
        mrna_trub_dat$change[i] <- "depleted"
    }else if((mrna_trub_dat[i,]$treat_level2.x-mrna_trub_dat[i,]$treat_level2.y)>=0.20){
        mrna_trub_dat$change[i] <- "down"
    }else{
        mrna_trub_dat$change[i] <- "non"
    }
}

write.table(mrna_trub_dat,"plot/trub1.mrna.site.txt",sep='\t',quote=F,row.names=F)

mrna_trub_dat <- mrna_trub_dat[,c("chr","end","treat_level2.x", "treat_level2.y", "change")]
colnames(mrna_trub_dat) <- c("chr", "pos","wt_mean", "ko_mean", "change")
mrna_trub_dat$type <- "mrna"
###trna_site for trub1 :ct_site1 mt_site1 
nrow(ct_site1) ###554
ct_trub1_dat <- ct_site1

for (i in 1:nrow(ct_trub1_dat)){
    if(ct_trub1_dat[i,]$ko1_level<=0.01 & ct_trub1_dat[i,]$ko2_level<=0.01& ct_trub1_dat[i,]$ko3_level<=0.01){
        ct_trub1_dat$change[i] <- "depleted"
    }else if(ct_trub1_dat[i,]$wt_mean-ct_trub1_dat[i,]$ko_mean>=0.20){
        ct_trub1_dat$change[i] <- "down"
    }else{
        ct_trub1_dat$change[i] <- "non"
    }
}
ct_trub1_dat <- ct_trub1_dat[,c("chr", "trna_pos","wt_mean", "ko_mean", "change")]
colnames(ct_trub1_dat) <- c("chr", "pos","wt_mean", "ko_mean", "change")
ct_trub1_dat$type <- "cy_trna"



###mt_site: 54
nrow(mt_site1) ##54
mt_trub1_dat <- mt_site1

mt_trub1_dat[mt_trub1_dat$trna_pos=="55",c("change")] <- "depleted"
mt_trub1_dat[mt_trub1_dat$trna_pos!="55",c("change")] <- "non"

mt_trub1_dat <- mt_trub1_dat[,c("chr", "trna_pos","wt_mean", "ko_mean", "change")]
colnames(mt_trub1_dat) <- c("chr", "pos","wt_mean", "ko_mean", "change")
mt_trub1_dat$type <- "mt_trna"

trub1_dat <- rbind(sno_trub_dat, mrna_trub_dat,ct_trub1_dat,mt_trub1_dat)
#> table(trub1_dat$type, trub1_dat$change)
         
          depleted down non
  cy_trna        0   88 466
  mrna          41    4 264
  mt_trna        6    0  48
  snoRNA         0    0 221

trub1_dat[trub1_dat$wt_mean>=1,]$wt_mean <- 1
trub1_dat[trub1_dat$wt_mean<0,]$wt_mean <- 0
trub1_dat[trub1_dat$ko_mean>=1,]$ko_mean <- 1
trub1_dat[trub1_dat$ko_mean<0,]$ko_mean <- 0
trub1_dat$group <- paste(trub1_dat$type,trub1_dat$change, sep="_")
trub1_dat$group <- factor(trub1_dat$group, level=c("snoRNA_non","mrna_non","cy_trna_non","mt_trna_non","mrna_down","cy_trna_down","mrna_depleted", "mt_trna_depleted"))
library(ggrepel)
pdf("plot/trub1_ko_scatter_figure.pdf",width=10,height=10)
ggplot(trub1_dat, aes(x=round(wt_mean*100,2), y=round(as.numeric(ko_mean)*100,2),color=group))+
geom_point(size = 3)+
scale_color_manual(values=c("grey","grey","grey","grey",'#1B9E77','#1B9E77','red', "red"))+
#geom_text_repel(data=wt_ko_trub[wt_ko_trub$type=="up"| wt_ko_trub$type=="down",],aes(x=round(wt_mean*100,2), y=round(as.numeric(ko_mean)*100,2),label=trna_pos),size = 5)+
ylab("ko_mean level(%)")+xlab("wt mean_level level(%)")+scale_x_continuous(limits = c(-1,101), expand = c(0, 0)) +
scale_y_continuous(limits = c(-1,101), expand = c(0, 0))+theme(axis.text = element_text(color = "black"))+theme(legend.position="none") +geom_abline(intercept = 0,slope = 1,linetype=2)+
theme_bw()
dev.off()

pdf("plot/trub1_ko_scatter_figure_nolegend.pdf",width=8,height=8)
ggplot(trub1_dat, aes(x=round(wt_mean*100,2), y=round(as.numeric(ko_mean)*100,2),color=group))+
geom_point(size = 5)+
scale_color_manual(values=c("grey","grey","grey","grey",'#1B9E77','#1B9E77','red', "red"))+
#geom_text_repel(data=wt_ko_trub[wt_ko_trub$type=="up"| wt_ko_trub$type=="down",],aes(x=round(wt_mean*100,2), y=round(as.numeric(ko_mean)*100,2),label=trna_pos),size = 5)+
ylab("ko_mean level(%)")+xlab("wt mean_level level(%)")+scale_x_continuous(limits = c(-1.5,101.5), expand = c(0, 0)) +
scale_y_continuous(limits = c(-1.5,101.5), expand = c(0, 0))+theme(axis.text = element_text(color = "black"))+theme(legend.position="none") +geom_abline(intercept = 0,slope = 1,linetype=2)+
theme_bw()+theme(legend.position="none")
dev.off()



########################for pus7
sno_pus7 <- merge(snosite[,c("chr", "start", "end", "ref","motif", "treat1_level", "treat2_level")], ko_pus7[,c("chr", "start", "end", "treat1_level", "treat2_level")], by=c("chr", "start", "end"))
colnames(sno_pus7) <- c("chr", "start", "end","ref","motif", "wt1_level", "wt2_level", "ko1_level", "ko2_level")###268
sno_pus7$ko_mean <- (sno_pus7$ko1_level+sno_pus7$ko2_level)/2
sno_pus7$wt_mean <- (sno_pus7$wt1_level+sno_pus7$wt2_level)/2

sno_pus7[sno_pus7$ko_mean<0.05,] ##3
for (i in 1:nrow(sno_pus7)){
    if(sno_pus7[i,]$ko1_level<=0.01 & sno_pus7[i,]$ko2_level<=0.01){
        sno_pus7$change[i] <- "depleted"
    }else{
        sno_pus7$change[i] <- "non"
    }
}
write.table(sno_pus7,"plot/pus7.snorna.site.txt",sep='\t',quote=F,row.names=F)
sno_pus7_dat <- sno_pus7[,c("chr", "end","wt_mean", "ko_mean", "change")]
colnames(sno_pus7_dat) <-c("chr", "pos","wt_mean", "ko_mean","change")
sno_pus7_dat$type <- "snoRNA"


mrna_pus7_dat <- pus7_mrna[(pus7_mrna$treat_T.x+pus7_mrna$treat_C.x)>=20 & (pus7_mrna$ctrl_T.x+pus7_mrna$ctrl_C.x)>=20 & pus7_mrna$treat_level2.x>=0.05 & pus7_mrna$treat_C.x>=10 & (pus7_mrna$treat_T.y+pus7_mrna$treat_C.y)>=20,]
nrow(mrna_pus7_dat) ##263

for (i in 1:nrow(mrna_pus7_dat)){
    if(mrna_pus7_dat[i,]$treat_level2.y<=0.01){
        mrna_pus7_dat$change[i] <- "depleted"
    }else{
        mrna_pus7_dat$change[i] <- "non"
    }
}

write.table(mrna_pus7_dat,"plot/pus7.mrna.site.txt",sep='\t',quote=F,row.names=F)
mrna_pus7_dat <- mrna_pus7_dat[,c("chr","end","treat_level2.x", "treat_level2.y", "change")]
colnames(mrna_pus7_dat) <- c("chr", "pos","wt_mean", "ko_mean", "change")
mrna_pus7_dat$type <- "mrna"
###trna_site for pus7 :ct_site1 mt_site1 
nrow(ct_site2) ###602
ct_pus7_dat <- ct_site2
ct_pus7_dat[ct_pus7_dat$trna_pos=="13"| ct_pus7_dat$trna_pos=="20B" |ct_pus7_dat$trna_pos=="35" | ct_pus7_dat$trna_pos=="36" | ct_pus7_dat$trna_pos=="50",c("change")] <- "depleted"
ct_pus7_dat[ct_pus7_dat$trna_pos!="13"& ct_pus7_dat$trna_pos!="20B" & ct_pus7_dat$trna_pos!="35" & ct_pus7_dat$trna_pos!="36" & ct_pus7_dat$trna_pos!="50",c("change")] <- "non"

ct_pus7_dat <- ct_pus7_dat[,c("chr", "trna_pos","wt_mean", "ko_mean", "change")]
colnames(ct_pus7_dat) <- c("chr", "pos","wt_mean", "ko_mean", "change")
ct_pus7_dat$type <- "cy_trna"



###mt_site: 54
nrow(mt_site2) ##54
mt_pus7_dat <- mt_site2
mt_pus7_dat[mt_pus7_dat$trna_pos=="50",c("change")] <- "depleted"
mt_pus7_dat[mt_pus7_dat$trna_pos!="50",c("change")] <- "non"

mt_pus7_dat <- mt_pus7_dat[,c("chr", "trna_pos","wt_mean", "ko_mean", "change")]
colnames(mt_pus7_dat) <- c("chr", "pos","wt_mean", "ko_mean", "change")
mt_pus7_dat$type <- "mt_trna"

#write.table(mrna_trub_dat,"plot/trub1.mrna.site.txt",sep='\t',quote=F,row.names=F)
pus7_dat <- rbind(sno_pus7_dat, mrna_pus7_dat,ct_pus7_dat,mt_pus7_dat)
#> table(pus7_dat$type, pus7_dat$change)
         
          depleted non
  cy_trna       69 533
  mrna          22 241
  mt_trna        1  53
  snoRNA         1 267

pus7_dat[pus7_dat$wt_mean>=1,]$wt_mean <- 1
pus7_dat[pus7_dat$wt_mean<0,]$wt_mean <- 0
pus7_dat[pus7_dat$ko_mean>=1,]$ko_mean <- 1
pus7_dat[pus7_dat$ko_mean<0,]$ko_mean <- 0
pus7_dat$group <- paste(pus7_dat$type,pus7_dat$change, sep="_")
pus7_dat$group <- factor(pus7_dat$group, level=c("snoRNA_non","mrna_non","cy_trna_non","mt_trna_non","snoRNA_depleted","mrna_depleted","cy_trna_depleted", "mt_trna_depleted"))
library(ggrepel)
pdf("plot/pus7_ko_scatter_figure_legend.pdf",width=15,height=15)
ggplot(pus7_dat, aes(x=round(wt_mean*100,2), y=round(as.numeric(ko_mean)*100,2),color=group))+
geom_point(size = 3)+
scale_color_manual(values=c("grey","grey","grey","grey","red",'red','red', "red"))+
#geom_text_repel(data=wt_ko_pus7[wt_ko_pus7$type=="up"| wt_ko_pus7$type=="down",],aes(x=round(wt_mean*100,2), y=round(as.numeric(ko_mean)*100,2),label=trna_pos),size = 5)+
ylab("ko_mean level(%)")+xlab("wt mean_level level(%)")+
scale_y_continuous(limits = c(-1,100), expand = c(0, 0))+theme(axis.text = element_text(color = "black"))+theme(legend.position="none") +geom_abline(intercept = 0,slope = 1,linetype=2)+
theme_bw()
dev.off()

pdf("plot/pus7_ko_scatter_figure_nolegend.pdf",width=8,height=8)
ggplot(pus7_dat, aes(x=round(wt_mean*100,2), y=round(as.numeric(ko_mean)*100,2),color=group))+
geom_point(size = 5)+
scale_color_manual(values=c("grey","grey","grey","grey","red",'red','red', "red"))+
#geom_text_repel(data=wt_ko_pus7[wt_ko_pus7$type=="up"| wt_ko_pus7$type=="down",],aes(x=round(wt_mean*100,2), y=round(as.numeric(ko_mean)*100,2),label=trna_pos),size = 5)+
ylab("ko_mean level(%)")+xlab("wt mean_level level(%)")+scale_x_continuous(limits = c(-1.5,101.5), expand = c(0, 0))+
scale_y_continuous(limits = c(-1.5,101.5), expand = c(0, 0))+theme(axis.text = element_text(color = "black"))+theme(legend.position="none") +geom_abline(intercept = 0,slope = 1,linetype=2)+
theme_bw()+theme(legend.position="none")
dev.off()



####

########################for pus1
sno_pus1 <- merge(snosite[,c("chr", "start", "end", "ref","motif", "treat1_level", "treat2_level")], ko_pus1[,c("chr", "start", "end", "treat1_level", "treat2_level")], by=c("chr", "start", "end"))
colnames(sno_pus1) <- c("chr", "start", "end","ref","motif", "wt1_level", "wt2_level", "ko1_level", "ko2_level")###251
sno_pus1$ko_mean <- (sno_pus1$ko1_level+sno_pus1$ko2_level)/2
sno_pus1$wt_mean <- (sno_pus1$wt1_level+sno_pus1$wt2_level)/2

sno_pus1[sno_pus1$ko_mean<0.05,] ##3
for (i in 1:nrow(sno_pus1)){
    if(sno_pus1[i,]$ko1_level<=0.01 & sno_pus1[i,]$ko2_level<=0.01){
        sno_pus1$change[i] <- "depleted"
    }else{
        sno_pus1$change[i] <- "non"
    }
}
sno_pus1_dat <- sno_pus1[,c("chr", "end","wt_mean", "ko_mean", "change")]
colnames(sno_pus1_dat) <-c("chr", "pos","wt_mean", "ko_mean","change")
sno_pus1_dat$type <- "snoRNA"


mrna_pus1_dat <- pus1_mrna[(pus1_mrna$treat_T.x+pus1_mrna$treat_C.x)>=20 & (pus1_mrna$ctrl_T.x+pus1_mrna$ctrl_C.x)>=20 & pus1_mrna$treat_level2.x>=0.05 & pus1_mrna$treat_C.x>=10 & (pus1_mrna$treat_T.y+pus1_mrna$treat_C.y)>=20,]
nrow(mrna_pus1_dat) ##259

for (i in 1:nrow(mrna_pus1_dat)){
    if(mrna_pus1_dat[i,]$treat_level2.y<=0.01){
        mrna_pus1_dat$change[i] <- "depleted"
    }else{
        mrna_pus1_dat$change[i] <- "non"
    }
}
write.table(mrna_pus1_dat,"plot/pus1.mrna.site.txt",sep='\t',quote=F,row.names=F)

mrna_pus1_dat <- mrna_pus1_dat[,c("chr","end","treat_level2.x", "treat_level2.y", "change")]
colnames(mrna_pus1_dat) <- c("chr", "pos","wt_mean", "ko_mean", "change")
mrna_pus1_dat$type <- "mrna"
###trna_site for pus1 :ct_site1 mt_site1 
nrow(ct_site3) ###600
ct_pus1_dat <- ct_site3
ct_pus1_dat[ct_pus1_dat$trna_pos=="27"| ct_pus1_dat$trna_pos=="28", c("change")] <- "depleted"
ct_pus1_dat[ct_pus1_dat$trna_pos!="27"& ct_pus1_dat$trna_pos!="28",c("change")] <- "non"

ct_pus1_dat <- ct_pus1_dat[,c("chr", "trna_pos","wt_mean", "ko_mean", "change")]
colnames(ct_pus1_dat) <- c("chr", "pos","wt_mean", "ko_mean", "change")
ct_pus1_dat$type <- "cy_trna"



###mt_site: 54
nrow(mt_site3) ##54
mt_pus1_dat <- mt_site3
mt_pus1_dat[mt_pus1_dat$trna_pos=="27" | mt_pus1_dat$trna_pos=="28" | mt_pus1_dat$trna_pos=="66" | mt_pus1_dat$trna_pos=="67" | mt_pus1_dat$trna_pos=="68" | mt_pus1_dat$trna_pos=="20" | mt_pus1_dat$trna_pos=="25",c("change")] <- "depleted"
mt_pus1_dat[mt_pus1_dat$trna_pos!="27" & mt_pus1_dat$trna_pos!="28" & mt_pus1_dat$trna_pos!="66" & mt_pus1_dat$trna_pos!="67" & mt_pus1_dat$trna_pos!="68" & mt_pus1_dat$trna_pos!="20" & mt_pus1_dat$trna_pos!="25",c("change")] <- "non"

mt_pus1_dat <- mt_pus1_dat[,c("chr", "trna_pos","wt_mean", "ko_mean", "change")]
colnames(mt_pus1_dat) <- c("chr", "pos","wt_mean", "ko_mean", "change")
mt_pus1_dat$type <- "mt_trna"

pus1_dat <- rbind(sno_pus1_dat, mrna_pus1_dat,ct_pus1_dat,mt_pus1_dat)
#> table(pus1_dat$type, pus1_dat$change)
         
          depleted non
  cy_trna      130 470
  mrna           8 251
  mt_trna       27  27
  snoRNA         0 251

pus1_dat[pus1_dat$wt_mean>=1,]$wt_mean <- 1
pus1_dat[pus1_dat$wt_mean<0,]$wt_mean <- 0
pus1_dat[pus1_dat$ko_mean>=1,]$ko_mean <- 1
pus1_dat[pus1_dat$ko_mean<0,]$ko_mean <- 0
pus1_dat$group <- paste(pus1_dat$type,pus1_dat$change, sep="_")
pus1_dat$group <- factor(pus1_dat$group, level=c("snoRNA_non","mrna_non","cy_trna_non","mt_trna_non","snoRNA_depleted","mrna_depleted","cy_trna_depleted", "mt_trna_depleted"))
library(ggrepel)
pdf("plot/pus1_ko_scatter_figure_legend.pdf",width=15,height=15)
ggplot(pus1_dat, aes(x=round(wt_mean*100,2), y=round(as.numeric(ko_mean)*100,2),color=group))+
geom_point(size = 3)+
scale_color_manual(values=c("grey","grey","grey","grey","red",'red','red', "red"))+
#geom_text_repel(data=wt_ko_pus1[wt_ko_pus1$type=="up"| wt_ko_pus1$type=="down",],aes(x=round(wt_mean*100,2), y=round(as.numeric(ko_mean)*100,2),label=trna_pos),size = 5)+
ylab("ko_mean level(%)")+xlab("wt mean_level level(%)")+
scale_y_continuous(limits = c(-1,100), expand = c(0, 0))+theme(axis.text = element_text(color = "black"))+theme(legend.position="none") +geom_abline(intercept = 0,slope = 1,linetype=2)+
theme_bw()
dev.off()

pdf("plot/pus1_ko_scatter_figure_nolegend.pdf",width=8,height=8)
ggplot(pus1_dat, aes(x=round(wt_mean*100,2), y=round(as.numeric(ko_mean)*100,2),color=group))+
geom_point(size = 5)+
scale_color_manual(values=c("grey","grey","grey","grey","red",'red','red', "red"))+
#geom_text_repel(data=wt_ko_pus1[wt_ko_pus1$type=="up"| wt_ko_pus1$type=="down",],aes(x=round(wt_mean*100,2), y=round(as.numeric(ko_mean)*100,2),label=trna_pos),size = 5)+
ylab("ko_mean level(%)")+xlab("wt mean_level level(%)")+scale_x_continuous(limits = c(-1.5,101.5), expand = c(0, 0)) +
scale_y_continuous(limits = c(-1.5,101.5), expand = c(0, 0))+theme(axis.text = element_text(color = "black"))+theme(legend.position="none") +geom_abline(intercept = 0,slope = 1,linetype=2)+
theme_bw()+theme(legend.position="none")
dev.off()


write.table(ct_site1,"plot/trub1.txt",sep='\t',quote=F,row.names=F)
write.table(ct_site2,"plot/pus7.txt",sep='\t',quote=F,row.names=F)
write.table(ct_site3,"plot/pus1.txt",sep='\t',quote=F,row.names=F)



##############plot PUS7 mt_50, pus1, mt25, and mt-20
wt_ko_pus7 <- read.table("plot/wt_pus7ko_table.txt", header=T, sep="\t", colClasses=c("ref"="character") )
wt_ko_pus7_mt50 <- wt_ko_pus7[wt_ko_pus7$chr=="mt_tRNA_Met", ] ###16
wt_ko_pus7_mt50[wt_ko_pus7_mt50$wt_mean>=1,]$wt_mean <- 1
wt_ko_pus7_mt50[wt_ko_pus7_mt50$wt_mean<0,]$wt_mean <- 0
wt_ko_pus7_mt50[wt_ko_pus7_mt50$ko_mean>=1,]$ko_mean <- 1
wt_ko_pus7_mt50[wt_ko_pus7_mt50$ko_mean<0,]$ko_mean <- 0
melt_pus7 <- melt(wt_ko_pus7_mt50[,c("trna_pos", "wt_mean", "ko_mean")], id.vars=c("trna_pos"))
melt_pus7$trna_pos <- as.numeric(melt_pus7$trna_pos)
pdf("plot/PUS7_mt_50.pdf", width=12,height = 4)
ggplot(melt_pus7, aes(group=variable, col=variable))+geom_segment(aes(x=trna_pos, y=round(value,4)*100,xend=trna_pos,yend=0), size=0.65)+geom_point(aes(x=trna_pos, y=round(value,4)*100),size=4,alpha = 0.8) +ylab("level(%)")+xlab("trna_pos")+
theme_bw() +scale_x_continuous(limits = c(0,75), expand = c(0, 0))+ scale_y_continuous(limits = c(-2,102), expand = c(0, 0))+scale_color_manual(values=c("grey40","#DC143C"))
dev.off()

wt_ko_pus7_mt50$delta <- wt_ko_pus7_mt50$ko_mean-wt_ko_pus7_mt50$wt_mean
wt_ko_pus7_mt50$trna_pos <- as.numeric(wt_ko_pus7_mt50$trna_pos)
#wt_ko_pus7_mt50[wt_ko_pus7_mt50$delta<0,]$delta <- 0
pdf("plot/PUS7_mt_delta_level.pdf", width=12,height = 4)
ggplot(wt_ko_pus7_mt50, aes(group=trna_pos))+geom_segment(aes(x=trna_pos, y=round(delta,4)*100,xend=trna_pos,yend=0), size=0.65)+geom_point(aes(x=trna_pos, y=round(delta,4)*100),size=4, color="red") +ylab("delta(%)")+xlab("trna_pos")+
theme_bw() +scale_x_continuous(limits = c(0,75), expand = c(0, 0))+ scale_y_continuous(limits = c(-101,5), expand = c(0, 0))#+scale_color_manual(values=c("grey40","#DC143C"))
dev.off()




wt_ko_pus1 <- read.table("plot/wt_pus1ko_table.txt", header=T, sep="\t", colClasses=c("ref"="character") )
wt_ko_pus1_mt20 <- wt_ko_pus1[wt_ko_pus1$chr=="mt_tRNA_Leu_UUR", ] ###19
wt_ko_pus1_mt20[wt_ko_pus1_mt20$wt_mean>=1,]$wt_mean <- 1
wt_ko_pus1_mt20[wt_ko_pus1_mt20$wt_mean<0,]$wt_mean <- 0
wt_ko_pus1_mt20[wt_ko_pus1_mt20$ko_mean>=1,]$ko_mean <- 1
wt_ko_pus1_mt20[wt_ko_pus1_mt20$ko_mean<0,]$ko_mean <- 0
melt_pus1_20 <- melt(wt_ko_pus1_mt20[,c("trna_pos", "wt_mean", "ko_mean")], id.vars=c("trna_pos"))
melt_pus1_20$trna_pos <- as.numeric(melt_pus1_20$trna_pos)
pdf("plot/PUS1_mt_20.pdf", width=12,height = 4)
ggplot(melt_pus1_20, aes(group=variable, col=variable))+geom_segment(aes(x=trna_pos, y=round(value,4)*100,xend=trna_pos,yend=0), size=0.65)+geom_point(aes(x=trna_pos, y=round(value,4)*100),size=4,alpha = 0.8) +ylab("level(%)")+xlab("trna_pos")+
theme_bw() +scale_x_continuous(limits = c(0,75),breaks = seq(0,75,25), expand = c(0, 0))+ scale_y_continuous(limits = c(-2,102), expand = c(0, 0))+scale_color_manual(values=c("grey40","#228B22"))
dev.off()

wt_ko_pus1_mt20$delta <- wt_ko_pus1_mt20$ko_mean-wt_ko_pus1_mt20$wt_mean
wt_ko_pus1_mt20$trna_pos <- as.numeric(wt_ko_pus1_mt20$trna_pos)
#wt_ko_pus1_mt20[wt_ko_pus1_mt20$delta<0,]$delta <- 0
pdf("plot/PUS1_mt_delta_level.pdf", width=12,height = 4)
ggplot(wt_ko_pus1_mt20, aes(group=trna_pos))+geom_segment(aes(x=trna_pos, y=round(delta,4)*100,xend=trna_pos,yend=0), size=0.65)+geom_point(aes(x=trna_pos, y=round(delta,4)*100),size=4, color="red") +ylab("delta(%)")+xlab("trna_pos")+
theme_bw() +scale_x_continuous(limits = c(0,75), expand = c(0, 0))+ scale_y_continuous(limits = c(-101,5), expand = c(0, 0))#+scale_color_manual(values=c("grey40","#DC143C"))
dev.off()


wt_ko_pus1_mt25 <- wt_ko_pus1[wt_ko_pus1$chr=="mt_tRNA_Asn", ] ###26
wt_ko_pus1_mt25[wt_ko_pus1_mt25$wt_mean>=1,]$wt_mean <- 1
wt_ko_pus1_mt25[wt_ko_pus1_mt25$wt_mean<0,]$wt_mean <- 0
wt_ko_pus1_mt25[wt_ko_pus1_mt25$ko_mean>=1,]$ko_mean <- 1
wt_ko_pus1_mt25[wt_ko_pus1_mt25$ko_mean<0,]$ko_mean <- 0
melt_pus1_25 <- melt(wt_ko_pus1_mt25[,c("trna_pos", "wt_mean", "ko_mean")], id.vars=c("trna_pos"))

melt_pus1_25$trna_pos <- as.numeric(melt_pus1_25$trna_pos)
pdf("plot/PUS1_mt_25.pdf", width=12,height = 4)
ggplot(melt_pus1_25, aes(group=variable, col=variable))+geom_segment(aes(x=trna_pos, y=round(value,4)*100,xend=trna_pos,yend=0), size=0.65)+geom_point(aes(x=trna_pos, y=round(value,4)*100),size=4,alpha = 0.8) +ylab("level(%)")+xlab("trna_pos")+
theme_bw() +scale_x_continuous(limits = c(0,75),breaks = seq(0,75,25), expand = c(0, 0))+ scale_y_continuous(limits = c(-2,102), expand = c(0, 0))+scale_color_manual(values=c("grey40","#228B22"))
dev.off()

wt_ko_pus1_mt25$delta <- wt_ko_pus1_mt25$ko_mean-wt_ko_pus1_mt25$wt_mean
wt_ko_pus1_mt25$trna_pos <- as.numeric(wt_ko_pus1_mt25$trna_pos)
#wt_ko_pus1_mt25[wt_ko_pus1_mt25$delta<0,]$delta <- 0
pdf("plot/pus1_mt25_delta_level.pdf", width=12,height = 4)
ggplot(wt_ko_pus1_mt25, aes(group=trna_pos))+geom_segment(aes(x=trna_pos, y=round(delta,4)*100,xend=trna_pos,yend=0), size=0.65)+geom_point(aes(x=trna_pos, y=round(delta,4)*100),size=4, color="red") +ylab("delta(%)")+xlab("trna_pos")+
theme_bw() +scale_x_continuous(limits = c(0,75), expand = c(0, 0))+ scale_y_continuous(limits = c(-101,5), expand = c(0, 0))#+scale_color_manual(values=c("grey40","#DC143C"))
dev.off()



####plot mt-mrna site
chrm <- read.table("site/bacs_chrM.txt", header=F, sep="\t",colClasses=c("V5"="character"))
chrm_site <- chrm[chrm$V41=="protein_coding" & chrm$V55!="improper_strand",] ###7
chrm_site <- chrm_site[,c("V1", "V2", "V3", "V4", "V5", "V6", "V26")]
colnames(chrm_site) <- c("chr", "start", "end", "motif", "strand", "ref", "wt_level")

ribo_chrm <- read.table("site/ribo_merge_treat_input_chrM_ivt_pvalue1", header=T, sep=" ", colClasses=c("ref"="character") ) 
pus1_chrm <- read.table("site/pus1_treat_input_chrM_ivt_pvalue1", header=T, sep=" ", colClasses=c("ref"="character") ) 
pus7_chrm <- read.table("site/pus7_treat_input_chrM_ivt_pvalue1", header=T, sep=" ", colClasses=c("ref"="character") ) 
trub1_chrm <- read.table("site/trub1_treat_input_chrM_ivt_pvalue1", header=T, sep=" ", colClasses=c("ref"="character") ) 


mergechrm <- list(chrm_site, ribo_chrm[,c("chr", "start", "end", "treat_level2")],pus1_chrm[,c("chr", "start", "end", "treat_level2")], pus7_chrm[,c("chr", "start", "end", "treat_level2")], trub1_chrm[,c("chr", "start", "end", "treat_level2")])
#mergechrm <- list(chrm_site, ribo_chrm,pus1_chrm, pus7_chrm, trub1_chrm)

mergechrm <- mergechrm %>% reduce(full_join, by=c("chr", "start", "end")) %>% na.omit()
nrow(mergechrm) #7
write.table(mergechrm,"plot/mt_mrna_wt_vs_ko.site.txt",sep='\t',quote=F,row.names=F)
> mergechrm[,7:11]
    wt_level treat_level2.x treat_level2.y treat_level2.x.x treat_level2.y.y
1 0.17627017     0.19336994   0.0013460521       0.19095093       0.28663036
2 0.09802562     0.10592030   0.0011240496       0.11658249       0.16103449
3 0.08911599     0.07062650  -0.0013260441       0.07597526       0.11363132
4 0.05513533     0.05010173  -0.0048283565       0.04841059       0.09350979
5 0.25696798     0.24032701  -0.0046729865       0.23767527       0.33254941
6 0.06565096     0.07015237   0.0025247714       0.05783962       0.10380120
7 0.16294760     0.15426722  -0.0008895296       0.15444488       0.23541034

mergechrm$xlabel <- "mt_mrna"
melt_chrm <- melt(mergechrm[,c("xlabel", "treat_level2.x", "treat_level2.y", "treat_level2.x.x", "treat_level2.y.y")], id.vars="xlabel")

melt_chrm[melt_chrm$value<0,]$value <- 0

pdf("plot/wt_ko_chrM_mrna_boxplot.pdf",width=4,height=4)
ggplot(melt_chrm, aes(x = xlabel, y=value*100, fill = variable)) + geom_boxplot(outlier.shape=NA, width = 0.5, position=position_dodge(width=0.7))+geom_point(pch = 21, position = position_jitterdodge(jitter.width =0.25, dodge.width=0.7))+scale_fill_manual(values=c("grey40","#228B22","#DC143C","#006699"))+ylim(0,50)+
theme_bw()+theme(legend.position="none")
dev.off()

melt_chrm1 <- melt(mergechrm[,c("xlabel", "treat_level2.x", "treat_level2.y")], id.vars="xlabel")

melt_chrm1[melt_chrm1$value<0,]$value <- 0

pdf("plot/pus1_wt_ko_chrM_boxplot.pdf",width=2,height=3)
ggplot(melt_chrm1, aes(x = xlabel, y=value*100, fill = variable)) + geom_boxplot(outlier.shape=NA, width = 0.5, position=position_dodge(width=0.7))+geom_point(pch = 21, position = position_jitterdodge(jitter.width =0.25, dodge.width=0.7))+scale_fill_manual(values=c("grey40","#228B22"))+ylim(0,50)+
theme_bw()+theme(legend.position="none")
dev.off()

mergechrm[mergechrm$treat_level2.x<0,]$treat_level2.x <- 0
mergechrm[mergechrm$treat_level2.y<0,]$treat_level2.y <- 0
pdf("plot/pus1_ko_chrM_scatter_figure_nolegend.pdf",width=8,height=8)
ggplot(mergechrm, aes(x=round(treat_level2.x*100,2), y=round(as.numeric(treat_level2.y )*100,2),color=xlabel))+
geom_point(size = 5)+
scale_color_manual(values=c("red"))+
#geom_text_repel(data=wt_ko_pus1[wt_ko_pus1$type=="up"| wt_ko_pus1$type=="down",],aes(x=round(wt_mean*100,2), y=round(as.numeric(ko_mean)*100,2),label=trna_pos),size = 5)+
ylab("ko_mean level(%)")+xlab("wt mean_level level(%)")+scale_x_continuous(limits = c(-0.3,25.3), expand = c(0, 0)) +
scale_y_continuous(limits = c(-0.3,25.3), expand = c(0, 0))+theme(axis.text = element_text(color = "black"))+theme(legend.position="none") +geom_abline(intercept = 0,slope = 1,linetype=2)+
theme_bw()+theme(legend.position="none")
dev.off()






##############################################compare pseudou sites on tRNA and rRNA
###trub1 dependent mt and mrna sites
trub1_mrna1 <- read.table("plot/trub1.mrna.site.txt", header=T, sep="\t",) ##309
pus7_mrna1 <- read.table("plot/pus7.mrna.site.txt", header=T, sep="\t") #263
pus1_mrna1 <- read.table("plot/pus1.mrna.site.txt", header=T, sep="\t") ##259

trub1_mrna1 <- trub1_mrna1[trub1_mrna1$change=="depleted",] ##41
trub1_mrna1$type <- "mrna"
trub1_mrna1 <- trub1_mrna1[,c("type","treat_level2.x")]
colnames(trub1_mrna1) <- c("type", "value")
trub1_mrna1$group <- "trub1"
trub1_trna <- seldat[seldat$xlabel=="pos55" & (seldat$variable=="wt_mean" | seldat$variable=="trub1_ko_mean"),]
trub1_trna$type <- "mt_trna"
trub1_trna <- trub1_trna[trub1_trna$variable=="wt_mean",c("type", "value")]
trub1_trna$group <- "trub1"
##PUS7_cytrna and pus7_mrna1
pus7_mrna1 <- pus7_mrna1[pus7_mrna1$change=="depleted",] ##22
pus7_mrna1$type <- "mrna"
pus7_mrna1 <- pus7_mrna1[,c("type","treat_level2.x")]
colnames(pus7_mrna1) <- c("type", "value")
pus7_mrna1$group <- "pus7"
pus7_trna <- pos_pus7
pus7_trna$type <- "trna"
pus7_trna <- pus7_trna[pus7_trna$variable=="wt_mean",c("type", "value")]
pus7_trna$group <- "pus7"
###pus1_mrna1
pus1_mrna1 <- pus1_mrna1[pus1_mrna1$change=="depleted",] ##8
pus1_mrna1$type <- "mrna"
pus1_mrna1 <- pus1_mrna1[,c("type","treat_level2.x")]
colnames(pus1_mrna1) <- c("type", "value")
pus1_mrna1$group <- "pus1"
pus1_trna <- pus1_pos ##310
pus1_trna[pus1_trna$xlabel=="ct_27/28","type"] <- "trna"
pus1_trna[pus1_trna$xlabel!="ct_27/28","type"] <- "mt_rna"
pus1_trna <- pus1_trna[pus1_trna$variable=="wt_mean",c("type", "value")]
pus1_trna$group <- "pus1"

pus1_mt <- mt_pus1_dat[mt_pus1_dat$change=="depleted",] ##27
pus1_mt <- pus1_mt[,c("type", "wt_mean")]
colnames(pus1_mt) <- c("type", "value")
pus1_mt$group <- "pus1"

combine <- rbind(trub1_mrna1, trub1_trna, pus7_mrna1, pus7_trna,pus1_mrna1, pus1_trna[pus1_trna$type=="trna",], pus1_mt)
          pus1 pus7 trub1
  mrna       8   22    41
  mt_trna   27    0     6
  trna     130   69     0

combine[combine$value<0,]$value <- 0
combine[combine$value>1,]$value <- 1
combine$group <- factor(combine$group, level=c("trub1", "pus7", "pus1"))
combine$type <- factor(combine$type, level=c("trna", "mt_trna", "mrna"))
p <- ggplot(combine,aes(x=group,y=value*100, fill=type)) +geom_boxplot(outlier.shape=NA,width=0.9,position = position_dodge2(preserve = "single",padding = 0.3))+# geom_point(aes(fill = variable), size = 1, shape = 21, position = position_jitterdodge()) #geom_jitter(color="black", size=0.4) +##outlier.shape = NA, width=0.8 alpha=0.9
  theme_light() + xlab("") +  theme(axis.text = element_text(color = "black")) + #scale_y_continuous(limits = c(-1,100), expand = c(0, 0))+
  ylab("level") +
  scale_fill_manual(values=c("#EE7621","#006699","#7FC97F"))
#theme(legend.position = "none") 
ggsave("plot/mRNA_tRNA_compare_boxplot.pdf",p, width=6, height = 4) 
write.table(combine,"plot/mRNA_tRNA_compare_boxplot.txt",sep='\t',quote=F,row.names=F)