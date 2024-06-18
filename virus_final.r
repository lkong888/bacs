#module load R/4.0.3-foss-2020b 
library(dplyr)
library(tidyr)
library(reshape2)

####plot for c666.1, raji, elijah
c666_rep1 <- read.table("/users/ludwig/ebu571/ebu571/07Nov2023/virus_sample/treat1_input1_virus_allsites.txt", header=T, sep=" ", colClasses=c("ref"="character") ) 
c666_rep2 <- read.table("/users/ludwig/ebu571/ebu571/14Oct2023/virus_sample/treat1_input1_virus_allsites.txt", header=T, sep=" ", colClasses=c("ref"="character") ) 
raji <- read.table("/users/ludwig/ebu571/ebu571/20Nov2023/virus_sample/treat10_input10_virus_allsites.txt", header=T, sep=" ", colClasses=c("ref"="character") ) 
eligah <- read.table("/users/ludwig/ebu571/ebu571/20Nov2023/virus_sample/treat11_input11_virus_allsites.txt", header=T, sep=" ", colClasses=c("ref"="character") ) 

###check cooverage
mean(c666_rep1[(c666_rep1$end>=6629 &c666_rep1$end<=6795) |(c666_rep1$end>=6956 &c666_rep1$end<=7128),]$treat_T+c666_rep1[(c666_rep1$end>=6629 &c666_rep1$end<=6795) |(c666_rep1$end>=6956 &c666_rep1$end<=7128),]$treat_C) ##66080.82
mean(c666_rep1[(c666_rep1$end>=6629 &c666_rep1$end<=6795) |(c666_rep1$end>=6956 &c666_rep1$end<=7128),]$ctrl_T+c666_rep1[(c666_rep1$end>=6629 &c666_rep1$end<=6795) |(c666_rep1$end>=6956 &c666_rep1$end<=7128),]$ctrl_C) ##85132.04
mean(c666_rep2[(c666_rep2$end>=6629 &c666_rep2$end<=6795) |(c666_rep2$end>=6956 &c666_rep2$end<=7128),]$treat_T+c666_rep2[(c666_rep2$end>=6629 &c666_rep2$end<=6795) |(c666_rep2$end>=6956 &c666_rep2$end<=7128),]$treat_C) ##32214.82
mean(c666_rep2[(c666_rep2$end>=6629 &c666_rep2$end<=6795) |(c666_rep2$end>=6956 &c666_rep2$end<=7128),]$ctrl_T+c666_rep2[(c666_rep2$end>=6629 &c666_rep2$end<=6795) |(c666_rep2$end>=6956 &c666_rep2$end<=7128),]$ctrl_C) ##34844.64
mean(raji[(raji$end>=6629 &raji$end<=6795) |(raji$end>=6956 &raji$end<=7128),]$treat_T+raji[(raji$end>=6629 &raji$end<=6795) |(raji$end>=6956 &raji$end<=7128),]$treat_C) ##87397.75
mean(raji[(raji$end>=6629 &raji$end<=6795) |(raji$end>=6956 &raji$end<=7128),]$ctrl_T+raji[(raji$end>=6629 &raji$end<=6795) |(raji$end>=6956 &raji$end<=7128),]$ctrl_C) ##52287.33
mean(eligah[(eligah$end>=6629 &eligah$end<=6795) |(eligah$end>=6956 &eligah$end<=7128),]$treat_T+eligah[(eligah$end>=6629 &eligah$end<=6795) |(eligah$end>=6956 &eligah$end<=7128),]$treat_C) ##10161.67
mean(eligah[(eligah$end>=6629 &eligah$end<=6795) |(eligah$end>=6956 &eligah$end<=7128),]$ctrl_T+eligah[(eligah$end>=6629 &eligah$end<=6795) |(eligah$end>=6956 &eligah$end<=7128),]$ctrl_C) ##6167.965


c666 <- merge(c666_rep1[,c("chr", "start", "end", "treat_level","ctrl_conversion")], c666_rep2[,c("chr", "start", "end", "treat_level","ctrl_conversion")], by=c("chr", "start", "end"))
c666$treat_level <- (c666$treat_level.x+c666$treat_level.y)/2
c666$ctrl_conversion <- (c666$ctrl_conversion.x+c666$ctrl_conversion.y)/2

c666_melt <- melt(c666[(c666$end>=6629 &c666$end<=6795) |(c666$end>=6956 &c666$end<=7128) ,c("end", "treat_level", "ctrl_conversion")], ,id.var="end")
###6999 is a T to C SNP
c666_melt[c666_melt$end==6999 & c666_melt$variable=="treat_level",]$value <- 0
c666_melt[c666_melt$value< 0,]$value <- 0
#c666_melt[c666_melt$value> 1,]$value <- 1
raji_melt <- melt(raji[(raji$end>=6629 &raji$end<=6795) |(raji$end>=6956 &raji$end<=7128) ,c("end", "treat_level", "ctrl_conversion")], ,id.var="end")
raji_melt[raji_melt$value< 0,]$value <- 0
#raji_melt[raji_melt$value> 1,]$value <- 1
eligah_melt <- melt(eligah[(eligah$end>=6629 &eligah$end<=6795) |(eligah$end>=6956 &eligah$end<=7128) ,c("end", "treat_level", "ctrl_conversion")], ,id.var="end")
eligah_melt[eligah_melt$value< 0,]$value <- 0
#eligah_melt[eligah_melt$value> 1,]$value <- 1
p1 <- ggplot(c666_melt, aes(group=variable, col=variable))+geom_segment(aes(x=end, y=round(value,4)*100,xend=end,yend=0), size=0.75)+geom_point(aes(x=end, y=round(value,4)*100),size=4) +ylab("ratio(%)")+xlab("")+ geom_hline(yintercept=5, linetype="dashed", color = "black")+
theme_bw()+scale_x_continuous() + scale_y_continuous(limits = c(0,100), expand = c(0, 0))+scale_color_manual(values=c('#FE767C','#ABABAB'))+ggtitle("c666.1_EBV")


p2 <- ggplot(raji_melt, aes(group=variable, col=variable))+geom_segment(aes(x=end, y=round(value,4)*100,xend=end,yend=0), size=0.75)+geom_point(aes(x=end, y=round(value,4)*100),size=4) +ylab("ratio(%)")+xlab("")+ geom_hline(yintercept=5, linetype="dashed", color = "black")+
theme_bw()+scale_x_continuous() + scale_y_continuous(limits = c(0,100), expand = c(0, 0))+scale_color_manual(values=c('#FE767C','#ABABAB'))+ggtitle("raji_EBV")

p3 <- ggplot(eligah_melt, aes(group=variable, col=variable))+geom_segment(aes(x=end, y=round(value,4)*100,xend=end,yend=0), size=0.75)+geom_point(aes(x=end, y=round(value,4)*100),size=4) +ylab("ratio(%)")+xlab("")+ geom_hline(yintercept=5, linetype="dashed", color = "black")+
theme_bw()+scale_x_continuous() + scale_y_continuous(limits = c(0,100), expand = c(0, 0))+scale_color_manual(values=c('#FE767C','#ABABAB'))+ggtitle("eligah_EBV")

library("ggpubr")
pdf("plot/ebv.conversion.pdf",width=10, height=10)
ggarrange(p1, p2, p3,
          ncol = 1, nrow = 3)
dev.off()

p1 <- ggplot(c666_melt[c666_melt$variable=="treat_level",], aes(group=variable, col=variable))+geom_segment(aes(x=end, y=round(value,4)*100,xend=end,yend=0), size=0.75)+geom_point(aes(x=end, y=round(value,4)*100),size=4) +ylab("ratio(%)")+xlab("")+ geom_hline(yintercept=5, linetype="dashed", color = "black")+
theme_bw()+scale_x_continuous() + scale_y_continuous(limits = c(0,100), expand = c(0, 0))+scale_color_manual(values=c('#77A9D7'))+ggtitle("c666.1_EBV")


p2 <- ggplot(raji_melt[raji_melt$variable=="treat_level",], aes(group=variable, col=variable))+geom_segment(aes(x=end, y=round(value,4)*100,xend=end,yend=0), size=0.75)+geom_point(aes(x=end, y=round(value,4)*100),size=4) +ylab("ratio(%)")+xlab("")+ geom_hline(yintercept=5, linetype="dashed", color = "black")+
theme_bw()+scale_x_continuous() + scale_y_continuous(limits = c(0,100), expand = c(0, 0))+scale_color_manual(values=c('#77A9D7'))+ggtitle("raji_EBV")

p3 <- ggplot(eligah_melt[eligah_melt$variable=="treat_level",], aes(group=variable, col=variable))+geom_segment(aes(x=end, y=round(value,4)*100,xend=end,yend=0), size=0.75)+geom_point(aes(x=end, y=round(value,4)*100),size=4) +ylab("ratio(%)")+xlab("")+ geom_hline(yintercept=5, linetype="dashed", color = "black")+
theme_bw()+scale_x_continuous() + scale_y_continuous(limits = c(0,100), expand = c(0, 0))+scale_color_manual(values=c('#77A9D7'))+ggtitle("eligah_EBV")

library("ggpubr")
pdf("plot/ebv.level.pdf",width=10, height=10)
ggarrange(p1, p2, p3,
          ncol = 1, nrow = 3)
dev.off()

#################virus plot for other RNA virus
covid <- read.table("/users/ludwig/ebu571/ebu571/20Nov2023/virus_sample/treat1_input1_virus_allsites.txt", header=T, sep=" ", colClasses=c("ref"="character") ) 
covid$ctrl_depth <- covid$ctrl_T+covid$ctrl_C
mean(covid[covid$chr=="NC_045512.2" & (covid$treat_T+covid$treat_C)>50 & (covid$ctrl_T+covid$ctrl_C)>50 & covid$ctrl_conversion<0.01,]$treat_T+covid[covid$chr=="NC_045512.2" & (covid$treat_T+covid$treat_C)>50 & (covid$ctrl_T+covid$ctrl_C)>50 & covid$ctrl_conversion<0.01,]$treat_C) ##457.0872
mean(covid[covid$chr=="NC_045512.2" & (covid$treat_T+covid$treat_C)>50 & (covid$ctrl_T+covid$ctrl_C)>50 & covid$ctrl_conversion<0.01,]$ctrl_T+covid[covid$chr=="NC_045512.2" & (covid$treat_T+covid$treat_C)>50 & (covid$ctrl_T+covid$ctrl_C)>50 & covid$ctrl_conversion<0.01,]$ctrl_C) ###272.9478
###read coverage
covid_ctrl <- read.table("/users/ludwig/ebu571/ebu571/20Nov2023/align/H9_2023Nov19_S9.virus.clean.filtered.mpile.txt", header=T, sep="\t",colClasses=c("ref_base"="character"))
covid_ctrl$coverage <- covid_ctrl$A+covid_ctrl$C+covid_ctrl$G+covid_ctrl$T+covid_ctrl$Gap+covid_ctrl$a+covid_ctrl$c+covid_ctrl$g+covid_ctrl$t+covid_ctrl$gap

pdf("plot/covid_coverage.pdf",width=6,height = 2)
p1 <- ggplot(covid_ctrl[covid_ctrl$chr=="NC_045512.2",], aes(x=pos,y=coverage)) +
  geom_line(size=0.5)+
  theme(legend.position = "bottom")+
  #ylim(0,100) +
  ylab("ctrl_coverage") +
  xlab("covid_pos") + scale_y_continuous(limits = c(0,4500), expand = c(0, 0))+
  #scale_color_manual(values=c("#23b177","#446db4"),name="") +
  theme(
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    panel.background = element_blank(), 
    axis.line = element_line(colour = "black"))
dev.off()

covid[covid$treat_level<0,]$treat_level <- 0
covid$type <- "covid"
pdf("plot/covid_treat_ctrl_point.pdf", width=6,height = 2)
p2 <- ggplot(covid[covid$chr=="NC_045512.2" & (covid$treat_T+covid$treat_C)>50 & (covid$ctrl_T+covid$ctrl_C)>50 & covid$ctrl_conversion<0.01,],aes(colour=type))+geom_segment(aes(x=end, y=round(treat_level,4)*100,xend=end,yend=0), size=0.65)+geom_point(aes(x=end, y=round(treat_level,4)*100),size=2,alpha = 0.8) +ylab("level(%)")+xlab("pos")+
theme_bw()+scale_y_continuous(limits = c(-2,102), expand = c(0, 0))+geom_hline(yintercept=10, linetype="dashed", color = "black")+theme(legend.position="none")+scale_color_manual(values=c('#77A9D7'))
dev.off()

library(ggpubr)
pdf("plot/covid_merge.pdf", width=6,height = 4)
ggarrange(p1,p2,ncol = 1, nrow = 2)
dev.off()

pdf("plot/covid_utoc_figure.pdf",width=5,height=5)
ggplot(covid[covid$chr=="NC_045512.2" & (covid$treat_T+covid$treat_C)>50 & (covid$ctrl_T+covid$ctrl_C)>50,], aes(x=round(ctrl_conversion*100,2), y=round(as.numeric(treat_conversion)*100,2)))+
geom_point(size = 2, color="#ABABAB")+scale_x_continuous(limits = c(-1,101), expand = c(0, 0)) +
  scale_y_continuous(limits = c(-1,101), expand = c(0, 0))+geom_abline(intercept = 0,slope = 1,linetype=2)+#+theme(axis.text = element_text(color = "#ABABAB"))
xlab("ctrl_utoc_ratio(%)")+ylab("treat_utoc_ratio(%)")+theme_bw()
dev.off()
##nrow(covid[covid$chr=="NC_045512.2" & (covid$treat_T+covid$treat_C)>200 & (covid$ctrl_T+covid$ctrl_C)>200,]) 1769

hcv <- read.table("/users/ludwig/ebu571/ebu571/20Nov2023/virus_sample/treat3_input3_virus_allsites.txt", header=T, sep=" ", colClasses=c("ref"="character") ) 
hcv$ctrl_depth <- hcv$ctrl_T+hcv$ctrl_C
mean(hcv[hcv$chr=="JF343782.1",]$treat_T+hcv[hcv$chr=="JF343782.1",]$treat_C) ##1740.219
mean(hcv[hcv$chr=="JF343782.1",]$ctrl_T+hcv[hcv$chr=="JF343782.1",]$ctrl_C) ##637.2646

hcv_ctrl <- read.table("/users/ludwig/ebu571/ebu571/20Nov2023/align/H3_2023Nov19_S3.virus.clean.filtered.mpile.txt", header=T, sep="\t",colClasses=c("ref_base"="character"))
hcv_ctrl$coverage <- hcv_ctrl$A+hcv_ctrl$C+hcv_ctrl$G+hcv_ctrl$T+hcv_ctrl$Gap+hcv_ctrl$a+hcv_ctrl$c+hcv_ctrl$g+hcv_ctrl$t+hcv_ctrl$gap

p1 <- ggplot(hcv_ctrl[hcv_ctrl$chr=="JF343782.1",], aes(x=pos,y=coverage)) +
  geom_line(size=0.5)+
  theme(legend.position = "bottom")+
  #ylim(0,100) +
  ylab("ctrl_coverage") +
  xlab("hcv_pos") + scale_y_continuous(limits = c(0,4000), expand = c(0, 0))+
  #scale_color_manual(values=c("#23b177","#446db4"),name="") +
  theme(
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    panel.background = element_blank(), 
    axis.line = element_line(colour = "black"))

hcv[hcv$treat_level<0,]$treat_level <- 0
hcv$type <- "hcv"

p2 <- ggplot(hcv[hcv$chr=="JF343782.1" & (hcv$treat_T+hcv$treat_C)>50 & (hcv$ctrl_T+hcv$ctrl_C)>50 & hcv$ctrl_conversion<0.01,],aes(colour=type))+geom_segment(aes(x=end, y=round(treat_level,4)*100,xend=end,yend=0), size=0.65)+geom_point(aes(x=end, y=round(treat_level,4)*100),size=2,alpha = 0.8) +ylab("level(%)")+xlab("pos")+
theme_bw()+scale_y_continuous(limits = c(-2,102), expand = c(0, 0))+geom_hline(yintercept=10, linetype="dashed", color = "black")+theme(legend.position="none")+scale_color_manual(values=c('#77A9D7'))

library(ggpubr)
pdf("plot/hcv_merge.pdf", width=6,height = 4)
ggarrange(p1,p2,ncol = 1, nrow = 2)
dev.off()

pdf("plot/hcv_utoc_figure.pdf",width=5,height=5)
ggplot(hcv[hcv$chr=="JF343782.1" & (hcv$treat_T+hcv$treat_C)>50 & (hcv$ctrl_T+hcv$ctrl_C)>50,], aes(x=round(ctrl_conversion*100,2), y=round(as.numeric(treat_conversion)*100,2)))+
geom_point(size = 2, color="#ABABAB")+scale_x_continuous(limits = c(-1,101), expand = c(0, 0)) +
  scale_y_continuous(limits = c(-1,101), expand = c(0, 0))+theme(axis.text = element_text(color = "black"))+geom_abline(intercept = 0,slope = 1,linetype=2)+
xlab("ctrl_utoc_ratio(%)")+ylab("treat_utoc_ratio(%)")+theme_bw()
dev.off()
##nrow(hcv[hcv$chr=="JF343782.1",]) ##2048

zika <- read.table("/users/ludwig/ebu571/ebu571/20Nov2023/virus_sample/treat4_input4_virus_allsites.txt", header=T, sep=" ", colClasses=c("ref"="character") ) 
zika$ctrl_depth <- zika$ctrl_T+zika$ctrl_C
mean(zika[zika$chr=="NC_035889.1",]$treat_T+zika[zika$chr=="NC_035889.1",]$treat_C) ##11918.64
mean(zika[zika$chr=="NC_035889.1",]$ctrl_T+zika[zika$chr=="NC_035889.1",]$ctrl_C) ##7226.122

zika_ctrl <- read.table("/users/ludwig/ebu571/ebu571/20Nov2023/align/H4_2023Nov19_S4.virus.clean.filtered.mpile.txt", header=T, sep="\t",colClasses=c("ref_base"="character"))
zika_ctrl$coverage <- zika_ctrl$A+zika_ctrl$C+zika_ctrl$G+zika_ctrl$T+zika_ctrl$Gap+zika_ctrl$a+zika_ctrl$c+zika_ctrl$g+zika_ctrl$t+zika_ctrl$gap


p1 <- ggplot(zika_ctrl[zika_ctrl$chr=="NC_035889.1",], aes(x=pos,y=coverage)) +
  geom_line(size=0.5)+
  theme(legend.position = "bottom")+
  #ylim(0,100) +
  ylab("ctrl_coverage") +
  xlab("zika_pos") + scale_y_continuous(limits = c(0,30000), expand = c(0, 0))+
  #scale_color_manual(values=c("#23b177","#446db4"),name="") +
  theme(
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    panel.background = element_blank(), 
    axis.line = element_line(colour = "black"))

zika[zika$treat_level<0,]$treat_level <- 0
zika$type <- "zika"

p2 <- ggplot(zika[zika$chr=="NC_035889.1" & (zika$treat_T+zika$treat_C)>50 & (zika$ctrl_T+zika$ctrl_C)>50 & zika$ctrl_conversion<0.01,],aes(colour=type))+geom_segment(aes(x=end, y=round(treat_level,4)*100,xend=end,yend=0), size=0.65)+geom_point(aes(x=end, y=round(treat_level,4)*100),size=2,alpha = 0.8) +ylab("level(%)")+xlab("pos")+
theme_bw()+scale_y_continuous(limits = c(-2,102), expand = c(0, 0))+geom_hline(yintercept=10, linetype="dashed", color = "black")+theme(legend.position="none")+scale_color_manual(values=c('#77A9D7'))

library(ggpubr)
pdf("plot/zika_merge.pdf", width=6,height = 4)
ggarrange(p1,p2,ncol = 1, nrow = 2)
dev.off()
pdf("plot/zika_utoc_figure.pdf",width=5,height=5)
ggplot(zika[zika$chr=="NC_035889.1" & (zika$treat_T+zika$treat_C)>50 & (zika$ctrl_T+zika$ctrl_C)>50,], aes(x=round(ctrl_conversion*100,2), y=round(as.numeric(treat_conversion)*100,2)))+
geom_point(size = 2, color="#ABABAB")+scale_x_continuous(limits = c(-1,101), expand = c(0, 0)) +
  scale_y_continuous(limits = c(-1,101), expand = c(0, 0))+theme(axis.text = element_text(color = "black"))+geom_abline(intercept = 0,slope = 1,linetype=2)+
xlab("ctrl_utoc_ratio(%)")+ylab("treat_utoc_ratio(%)")+theme_bw()
dev.off()
##nrow(zika[zika$chr=="NC_035889.1",]) ##2309
##########hdv
hdv <- read.table("/users/ludwig/ebu571/ebu571/20Nov2023/virus_sample/treat6_input6_virus_allsites.txt", header=T, sep=" ", colClasses=c("ref"="character") ) 
hdv$ctrl_depth <- hdv$ctrl_T+hdv$ctrl_C
mean(hdv[hdv$chr=="AJ000558.1",]$treat_T+hdv[hdv$chr=="AJ000558.1",]$treat_C) ##6639.882
mean(hdv[hdv$chr=="AJ000558.1",]$ctrl_T+hdv[hdv$chr=="AJ000558.1",]$ctrl_C) ##3683.87

hdv_ctrl <- read.table("/users/ludwig/ebu571/ebu571/20Nov2023/align/H6_2023Nov19_S6.virus.clean.filtered.mpile.txt", header=T, sep="\t",colClasses=c("ref_base"="character"))
hdv_ctrl$coverage <- hdv_ctrl$A+hdv_ctrl$C+hdv_ctrl$G+hdv_ctrl$T+hdv_ctrl$Gap+hdv_ctrl$a+hdv_ctrl$c+hdv_ctrl$g+hdv_ctrl$t+hdv_ctrl$gap


p1 <- ggplot(hdv_ctrl[hdv_ctrl$chr=="AJ000558.1",], aes(x=pos,y=coverage)) +
  geom_line(size=0.5)+
  theme(legend.position = "bottom")+
  #ylim(0,100) +
  ylab("ctrl_coverage") +
  xlab("hdv_pos") + scale_y_continuous(limits = c(0,15000), expand = c(0, 0))+
  #scale_color_manual(values=c("#23b177","#446db4"),name="") +
  theme(
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    panel.background = element_blank(), 
    axis.line = element_line(colour = "black"))

hdv[hdv$treat_level<0,]$treat_level <- 0
hdv$type <- "hdv"

p2 <- ggplot(hdv[hdv$chr=="AJ000558.1" & (hdv$treat_T+hdv$treat_C)>50 & (hdv$ctrl_T+hdv$ctrl_C)>50 & hdv$ctrl_conversion<0.01,],aes(colour=type))+geom_segment(aes(x=end, y=round(treat_level,4)*100,xend=end,yend=0), size=0.65)+geom_point(aes(x=end, y=round(treat_level,4)*100),size=2,alpha = 0.8) +ylab("level(%)")+xlab("pos")+
theme_bw()+scale_y_continuous(limits = c(-2,102), expand = c(0, 0))+geom_hline(yintercept=10, linetype="dashed", color = "black")+theme(legend.position="none")+scale_color_manual(values=c('#77A9D7'))

library(ggpubr)
pdf("plot/hdv_merge.pdf", width=6,height = 4)
ggarrange(p1,p2,ncol = 1, nrow = 2)
dev.off()

pdf("plot/hdv_utoc_figure.pdf",width=5,height=5)
ggplot(hdv[hdv$chr=="AJ000558.1" & (hdv$treat_T+hdv$treat_C)>50 & (hdv$ctrl_T+hdv$ctrl_C)>50,], aes(x=round(ctrl_conversion*100,2), y=round(as.numeric(treat_conversion)*100,2)))+
geom_point(size = 2, color="#ABABAB")+scale_x_continuous(limits = c(-1,101), expand = c(0, 0)) +
  scale_y_continuous(limits = c(-1,101), expand = c(0, 0))+theme(axis.text = element_text(color = "black"))+geom_abline(intercept = 0,slope = 1,linetype=2)+
xlab("ctrl_utoc_ratio(%)")+ylab("treat_utoc_ratio(%)")+theme_bw()
dev.off()
##nrow(hdv[hdv$chr=="AJ000558.1",]) ##339


sinv <- read.table("/users/ludwig/ebu571/ebu571/14Oct2023/virus_sample/treat1_input1_sinv_allsites.txt", header=T, sep=" ", colClasses=c("ref"="character") ) 
nrow(sinv[sinv$chr=="NC_001547.1" & sinv$ctrl_conversion<0.01 & sinv$treat_level>0.05,]) ##2 but they are snp
mean(sinv[sinv$chr=="NC_001547.1",]$treat_T+sinv[sinv$chr=="NC_001547.1",]$treat_C) ##16182.95
mean(sinv[sinv$chr=="NC_001547.1",]$ctrl_T+sinv[sinv$chr=="NC_001547.1",]$ctrl_C) ##11966.25

sinv_ctrl <- read.table("/users/ludwig/ebu571/ebu571/14Oct2023/align/H9_2023Oct13_S9.virus.clean.filtered.mpile.txt", header=T, sep="\t",colClasses=c("ref_base"="character"))
sinv_ctrl$coverage <- sinv_ctrl$A+sinv_ctrl$C+sinv_ctrl$G+sinv_ctrl$T+sinv_ctrl$Gap+sinv_ctrl$a+sinv_ctrl$c+sinv_ctrl$g+sinv_ctrl$t+sinv_ctrl$gap


p1 <- ggplot(sinv_ctrl[sinv_ctrl$chr=="NC_001547.1",], aes(x=pos,y=coverage)) +
  geom_line(size=0.5)+
  theme(legend.position = "bottom")+
  #ylim(0,100) +
  ylab("ctrl_coverage") +
  xlab("sinv_pos") + scale_y_continuous(limits = c(0,60000), expand = c(0, 0))+
  #scale_color_manual(values=c("#23b177","#446db4"),name="") +
  theme(
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    panel.background = element_blank(), 
    axis.line = element_line(colour = "black"))

sinv[sinv$treat_level<0,]$treat_level <- 0
sinv$type <- "sinv"

p2 <- ggplot(sinv[sinv$chr=="NC_001547.1" & (sinv$treat_T+sinv$treat_C)>50 & (sinv$ctrl_T+sinv$ctrl_C)>50 & sinv$ctrl_conversion<0.01,],aes(colour=type))+geom_segment(aes(x=end, y=round(treat_level,4)*100,xend=end,yend=0), size=0.65)+geom_point(aes(x=end, y=round(treat_level,4)*100),size=2,alpha = 0.8) +ylab("level(%)")+xlab("pos")+
theme_bw()+scale_y_continuous(limits = c(-2,102), expand = c(0, 0))+geom_hline(yintercept=10, linetype="dashed", color = "black")+theme(legend.position="none")+scale_color_manual(values=c('#77A9D7'))

library(ggpubr)
pdf("plot/sinv_merge.pdf", width=6,height = 4)
ggarrange(p1,p2,ncol = 1, nrow = 2)
dev.off()
pdf("plot/sinv_utoc_figure.pdf",width=5,height=5)
###one snp site was filtered out: (sinv$treat_T+sinv$treat_C)/sinv$treat_depth>0.90
ggplot(sinv[sinv$chr=="NC_001547.1" & (sinv$treat_T+sinv$treat_C)>50 & (sinv$ctrl_T+sinv$ctrl_C)>50 &(sinv$treat_T+sinv$treat_C)/sinv$treat_depth>0.90 ,], aes(x=round(ctrl_conversion*100,2), y=round(as.numeric(treat_conversion)*100,2)))+
geom_point(size = 2, color="#ABABAB")+scale_x_continuous(limits = c(-1,101), expand = c(0, 0)) +
  scale_y_continuous(limits = c(-1,101), expand = c(0, 0))+theme(axis.text = element_text(color = "black"))+geom_abline(intercept = 0,slope = 1,linetype=2)+
xlab("ctrl_utoc_ratio(%)")+ylab("treat_utoc_ratio(%)")+theme_bw()
dev.off()
##nrow(sinv[sinv$chr=="NC_001547.1",]) ##2048

####only plot E6: 37-474; E7:446-712 L1: 5333-6841
hpv1 <- read.table("/users/ludwig/ebu571/ebu571/27Oct2023/virus_sample/treat3_input3_virus_allsites.txt", header=T, sep=" ", colClasses=c("ref"="character") ) 
mean(hpv1[(hpv1$end>=37&hpv1$end<=712) | (hpv1$end>=5333 & hpv1$end<=6841),]$treat_T+hpv1[(hpv1$end>=37&hpv1$end<=712) | (hpv1$end>=5333 & hpv1$end<=6841),]$treat_C) ##1104.23
mean(hpv1[(hpv1$end>=37&hpv1$end<=712) | (hpv1$end>=5333 & hpv1$end<=6841),]$ctrl_T+hpv1[(hpv1$end>=37&hpv1$end<=712) | (hpv1$end>=5333 & hpv1$end<=6841),]$ctrl_C) ##740.6384
hpv2 <- read.table("/users/ludwig/ebu571/ebu571/27Oct2023/virus_sample/treat4_input4_virus_allsites.txt", header=T, sep=" ", colClasses=c("ref"="character") ) 
mean(hpv2[(hpv2$end>=37&hpv2$end<=712) | (hpv2$end>=5333 & hpv2$end<=6841),]$treat_T+hpv2[(hpv2$end>=37&hpv2$end<=712) | (hpv2$end>=5333 & hpv2$end<=6841),]$treat_C) ##1162.516
mean(hpv2[(hpv2$end>=37&hpv2$end<=712) | (hpv2$end>=5333 & hpv2$end<=6841),]$ctrl_T+hpv2[(hpv2$end>=37&hpv2$end<=712) | (hpv2$end>=5333 & hpv2$end<=6841),]$ctrl_C) #654.9641

hpv <- merge(hpv1[(hpv1$treat_T+hpv1$treat_C)>=50 &(hpv1$ctrl_T+hpv1$ctrl_C)>=50 ,c("chr", "start", "end", "treat_conversion", "ctrl_conversion","treat_level")],hpv2[(hpv2$treat_T+hpv2$treat_C)>=50 &(hpv2$ctrl_T+hpv2$ctrl_C)>=50,c("chr", "start", "end", "treat_conversion", "ctrl_conversion","treat_level")], by=c("chr", "start", "end"))
hpv$treat_conversion <- (hpv$treat_conversion.x+hpv$treat_conversion.y)/2
hpv$ctrl_conversion <- (hpv$ctrl_conversion.x+hpv$ctrl_conversion.y)/2
hpv$treat_level <- (hpv$treat_level.x+hpv$treat_level.y)/2
hpv <- hpv[(hpv$end>=37&hpv$end<=712) | (hpv$end>=5333 & hpv$end<=6841),]

hpv1_ctrl <- read.table("/users/ludwig/ebu571/ebu571/27Oct2023/align/H1_2023Oct26_S1.virus.clean.filtered.mpile.txt", header=T, sep="\t",colClasses=c("ref_base"="character"))
hpv1_ctrl$coverage <- hpv1_ctrl$A+hpv1_ctrl$C+hpv1_ctrl$G+hpv1_ctrl$T+hpv1_ctrl$Gap+hpv1_ctrl$a+hpv1_ctrl$c+hpv1_ctrl$g+hpv1_ctrl$t+hpv1_ctrl$gap

hpv2_ctrl <- read.table("/users/ludwig/ebu571/ebu571/27Oct2023/align/H2_2023Oct26_S2.virus.clean.filtered.mpile.txt", header=T, sep="\t",colClasses=c("ref_base"="character"))
hpv2_ctrl$coverage <- hpv2_ctrl$A+hpv2_ctrl$C+hpv2_ctrl$G+hpv2_ctrl$T+hpv2_ctrl$Gap+hpv2_ctrl$a+hpv2_ctrl$c+hpv2_ctrl$g+hpv2_ctrl$t+hpv2_ctrl$gap

hpv_ctrl <- merge(hpv1_ctrl,hpv2_ctrl, by=c("chr", "pos"))
hpv_ctrl$coverage <- (hpv_ctrl$coverage.x+hpv_ctrl$coverage.y)/2
hpv_ctrl <- hpv_ctrl[((hpv_ctrl$pos>=37&hpv_ctrl$pos<=712) | (hpv_ctrl$pos>=5333 & hpv_ctrl$pos<=6841)),]

p1 <- ggplot(hpv_ctrl[hpv_ctrl$chr=="HPV18REF",], aes(x=pos,y=coverage)) +
  geom_line(size=0.5)+
  theme(legend.position = "bottom")+
  #ylim(0,100) +
  ylab("ctrl_coverage") +
  xlab("hpv_pos") + #scale_y_continuous(limits = c(0,15000), expand = c(0, 0))+
  #scale_color_manual(values=c("#23b177","#446db4"),name="") +
  theme(
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    panel.background = element_blank(), 
    axis.line = element_line(colour = "black"))

hpv[hpv$treat_level<0,]$treat_level <- 0
hpv$type <- "hpv"

p2 <- ggplot(hpv[hpv$chr=="HPV18REF" & hpv$ctrl_conversion<0.01,],aes(colour=type))+geom_segment(aes(x=end, y=round(treat_level,4)*100,xend=end,yend=0), size=0.65)+geom_point(aes(x=end, y=round(treat_level,4)*100),size=2,alpha = 0.8) +ylab("level(%)")+xlab("pos")+
theme_bw()+scale_y_continuous(limits = c(-2,102), expand = c(0, 0))+geom_hline(yintercept=10, linetype="dashed", color = "black")+theme(legend.position="none")+scale_color_manual(values=c('#77A9D7'))

library(ggpubr)
pdf("plot/hpv_merge.pdf", width=6,height = 4)
ggarrange(p1,p2,ncol = 1, nrow = 2)
dev.off()

nrow(hpv[hpv$chr=="HPV18REF" & hpv$ctrl_conversion<0.01 & hpv$treat_level>0.05,]) ##0
nrow(hpv[hpv$chr=="HPV18REF" ,]) ##612
pdf("plot/hpv_utoc_figure.pdf",width=5,height=5)
ggplot(hpv[hpv$chr=="HPV18REF" ,], aes(x=round(ctrl_conversion*100,2), y=round(as.numeric(treat_conversion)*100,2)))+
geom_point(size = 2, color="#ABABAB")+scale_x_continuous(limits = c(-1,101), expand = c(0, 0)) +
  scale_y_continuous(limits = c(-1,101), expand = c(0, 0))+theme(axis.text = element_text(color = "black"))+geom_abline(intercept = 0,slope = 1,linetype=2)+
xlab("ctrl_utoc_ratio(%)")+ylab("treat_utoc_ratio(%)")+theme_bw()
dev.off()



