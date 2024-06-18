##module load R/4.0.3-foss-2020b
bid_only_in_bacs <- read.table("others_data/bid_uniq_site_in_bacs.txt") ##276

nrow(bid_only_in_bacs[bid_only_in_bacs$V41 <0.05, ]) ##245
min(bid_only_in_bacs[bid_only_in_bacs$V41>=0.05, ]$V43 ) ##0.01028826

praise_only <- read.table("others_data/praise_uniq_site_in_bacs_coverage20.txt") ##1718
unique(praise_only[,c("V1","V2","V3")]) %>% nrow() ##909
 
 praise_agg <- aggregate(V29 ~ V1+V2+V3, praise_only, mean) ##909
 nrow(praise_agg[praise_agg$V29<0.05,]) ##760

 praise_dat1 <- merge(praise_only, praise_agg[praise_agg$V29>=0.05,], by=c("V1","V2","V3")) ##187
  praise_dat1[,c("V1","V2","V3")] %>% unique() %>% nrow() ##149
 nrow(praise_dat1)
praise_only[praise_only$V29<0.05,c("V1","V2","V3")] %>% unique() %>% nrow() ##796
praise_only[praise_only$V29>=0.05,c("V1","V2","V3")] %>% unique() %>% nrow() ##256
praise_only[praise_only$V29>=0.05,]$V31 %>% min() ##0.01005673


library(VennDiagram)
##########################Final plot for comparision with others
pdf(file="plot/mRNA_bacs_bid_venn.pdf",width=4, height=4)
draw.pairwise.venn(area1 = 1335,    # Draw pairwise venn diagram
                   area2 = 575,
                   cross.area = 230,category=c("bacs","bid_seq"),col=c("#FE767C","#53D2EE"),fill=c("#FE767C","#53D2EE"), cat.cex=0,cat.pos=c(-27, 27), cat.dist=c(0.03, 0.03), height = 50 , 
        width = 50, lwd=0.5,ext.text=FALSE)
dev.off()

pdf(file="plot/mRNA_bacs_praise_venn.pdf",width=4, height=4)
draw.pairwise.venn(area1 = 1335,    # Draw pairwise venn diagram
                   area2 = 1995,
                   cross.area = 651,category=c("bacs","praise_seq"),col=c("#FE767C","#9ECE88"),fill=c("#FE767C","#9ECE88"), cat.cex=0,cat.pos=c(-27, 27), cat.dist=c(0.03, 0.03), height = 50 ,inverted=TRUE, 
        width = 50, lwd=0.5,ext.text=FALSE)
dev.off()


pdf(file="plot/mRNA_bacs_grhighest_venn.pdf",width=4, height=4)
draw.pairwise.venn(area1 = 1335,    # Draw pairwise venn diagram
                   area2 = 70,
                   cross.area = 61,category=c("bacs","highest"),col=c("#FE767C","#666666"),fill=c("#FE767C","#666666"), cat.cex=0,cat.pos=c(-27, 27), cat.dist=c(0.03, 0.03), height = 50 , 
        width = 50, lwd=0.5,ext.text=FALSE)
dev.off()

pdf(file="plot/mRNA_bacs_grhigh_venn.pdf",width=4, height=4)
draw.pairwise.venn(area1 = 1335,    # Draw pairwise venn diagram
                   area2 = 813,
                   cross.area = 196,category=c("bacs","high"),col=c("#FE767C","#666666"),fill=c("#FE767C","#666666"), cat.cex=0,cat.pos=c(-27, 27), cat.dist=c(0.03, 0.03), height = 50 , 
        width = 50, lwd=0.5,ext.text=FALSE)
dev.off()


pdf(file="plot/mRNA_bacs_grhigh_covered9smp_venn.pdf",width=4, height=4)
draw.pairwise.venn(area1 = 1335,    # Draw pairwise venn diagram
                   area2 = 320,
                   cross.area = 177,category=c("bacs","high_9smp"),col=c("#FE767C","#666666"),fill=c("#FE767C","#666666"), cat.cex=0,cat.pos=c(-27, 27), cat.dist=c(0.03, 0.03), height = 50 , 
        width = 50, lwd=0.5,ext.text=FALSE)
dev.off()


bid_site <- data.frame(group=c("low_coverage","low_level","low_p"),value=c(69,245,31))
bid_site$group <- factor(bid_site$group,level=c("low_coverage","low_level","low_p"))
praise_site <- data.frame(group=c("low_coverage","low_level","low_p"),value=c(435,760,149))
praise_site$group <- factor(praise_site$group,level=c("low_coverage","low_level","low_p"))

p1 <- ggplot(bid_site, aes(x="", y=as.numeric(value), fill=group))+
geom_bar(width = 1, stat = "identity")+coord_polar("y") + scale_fill_brewer(palette="Set2") +  ###or RdYlBu
  geom_text(aes(label = as.numeric(value)),position = position_stack(vjust = 0.5))+
  theme_void()

p2 <- ggplot(praise_site, aes(x="", y=as.numeric(value), fill=group))+
geom_bar(width = 1, stat = "identity")+coord_polar("y") + scale_fill_brewer(palette="Set2") +  ###or RdYlBu
  geom_text(aes(label = as.numeric(value)),position = position_stack(vjust = 0.5))+
  theme_void()
library(ggpubr)
pdf("plot/bis_praise.pdf")
ggarrange(p1,p2,ncol = 2, nrow = 1)
dev.off()
