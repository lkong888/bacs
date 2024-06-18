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

###################70mer spikein 
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
dat$smp <- gsub("*/align/","",dat$sample)
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



###write intergrated output for calibration curve 
h1 <- read.table("/users/ludwig/ebu571/ebu571/07Nov2023/align/H1_2023Nov06_S1.clean.nnunn.sort_NNUNN1.contex.table", header=T, sep=" ")
h2 <- read.table("/users/ludwig/ebu571/ebu571/07Nov2023/align/H2_2023Nov06_S2.clean.nnunn.sort_NNUNN1.contex.table", header=T, sep=" ")
h3 <- read.table("/users/ludwig/ebu571/ebu571/07Nov2023/align/H3_2023Nov06_S3.clean.nnunn.sort_NNUNN1.contex.table", header=T, sep=" ")
h4 <- read.table("/users/ludwig/ebu571/ebu571/07Nov2023/align/H4_2023Nov06_S4.clean.nnunn.sort_NNUNN1.contex.table", header=T, sep=" ")
h5 <- read.table("/users/ludwig/ebu571/ebu571/07Nov2023/align/H5_2023Nov06_S5.clean.nnunn.sort_NNUNN1.contex.table", header=T, sep=" ")
h6 <- read.table("/users/ludwig/ebu571/ebu571/07Nov2023/align/H6_2023Nov06_S6.clean.nnunn.sort_NNUNN1.contex.table", header=T, sep=" ")
merge1 <- list(h1[,c("first", "last", "ratio")], h2[,c("first", "last", "ratio")],h3[,c("first", "last", "ratio")],h4[,c("first", "last", "ratio")],h5[,c("first", "last", "ratio")],h6[,c("first", "last", "ratio")])
merge1 <- merge1 %>% reduce(full_join, by=c("first", "last"))
colnames(merge1) <- c("first", "last", "h1_ratio", "h2_ratio", "h3_ratio", "h4_ratio", "h5_ratio", "h6_ratio")
write.table(merge1, "site/calibration_nnunn1.context_all.table", quote = FALSE, col.names = T, row.names = F)

####plot conversion rate for NNUNN1
h7 <- read.table("/users/ludwig/ebu571/ebu571/27Oct2023/align/H7_2023Oct26_S7.clean.nnunn.sort_NNUNN1.contex.table", header=T, sep=" ")
h8 <- read.table("/users/ludwig/ebu571/ebu571/27Oct2023/align/H8_2023Oct26_S8.clean.nnunn.sort_NNUNN1.contex.table", header=T, sep=" ")

h3 <- read.table("/users/ludwig/ebu571/ebu571/27Oct2023/align/H3_2023Oct26_S3.clean.nnunn.sort_NNUNN1.contex.table", header=T, sep=" ")
h4 <- read.table("/users/ludwig/ebu571/ebu571/27Oct2023/align/H4_2023Oct26_S4.clean.nnunn.sort_NNUNN1.contex.table", header=T, sep=" ")

h5 <- read.table("/users/ludwig/ebu571/ebu571/07Nov2023/align/H5_2023Nov06_S5.clean.nnunn.sort_NNUNN1.contex.table", header=T, sep=" ")
h6 <- read.table("/users/ludwig/ebu571/ebu571/07Nov2023/align/H6_2023Nov06_S6.clean.nnunn.sort_NNUNN1.contex.table", header=T, sep=" ")

mergemdf <- list(h3[,c("first", "last", "ratio")], h4[,c("first", "last", "ratio")],h7[,c("first", "last", "ratio")], h8[,c("first", "last", "ratio")],h5[,c("first", "last", "ratio")], h6[,c("first", "last", "ratio")])
mergemdf <- mergemdf %>% reduce(full_join, by=c("first", "last"))
colnames(mergemdf) <- c("first", "last", "h3_ratio", "h4_ratio","h7_ratio", "h8_ratio","h5_ratio", "h6_ratio")

mergemdf$ratio <- round((mergemdf$h3_ratio+mergemdf$h4_ratio+mergemdf$h7_ratio+mergemdf$h8_ratio+mergemdf$h5_ratio+mergemdf$h6_ratio)/6,2)*100
mdf1 <- spread(mergemdf[,c("first", "last", "ratio")], key = last, value = ratio) %>%  as.data.frame()
final_matrix <- data.frame(mdf1[,-1], row.names=as.character(mdf1[,1])) %>% as.matrix(rownames = TRUE)

write.table(mergemdf, "plot/modified_nnunn.context.table", quote = FALSE, col.names = T, row.names = F)

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




###for unmodified nnunn (nnunn2)
h7 <- read.table("/users/ludwig/ebu571/ebu571/27Oct2023/align/H7_2023Oct26_S7.clean.nnunn.sort_NNUNN2.contex.table", header=T, sep=" ")
h8 <- read.table("/users/ludwig/ebu571/ebu571/27Oct2023/align/H8_2023Oct26_S8.clean.nnunn.sort_NNUNN2.contex.table", header=T, sep=" ")

h3 <- read.table("/users/ludwig/ebu571/ebu571/27Oct2023/align/H3_2023Oct26_S3.clean.nnunn.sort_NNUNN2.contex.table", header=T, sep=" ")
h4 <- read.table("/users/ludwig/ebu571/ebu571/27Oct2023/align/H4_2023Oct26_S4.clean.nnunn.sort_NNUNN2.contex.table", header=T, sep=" ")

h5 <- read.table("/users/ludwig/ebu571/ebu571/07Nov2023/align/H5_2023Nov06_S5.clean.nnunn.sort_NNUNN2.contex.table", header=T, sep=" ")
h6 <- read.table("/users/ludwig/ebu571/ebu571/07Nov2023/align/H6_2023Nov06_S6.clean.nnunn.sort_NNUNN2.contex.table", header=T, sep=" ")

merge1 <- list(h3[,c("first", "last", "ratio")], h4[,c("first", "last", "ratio")],h7[,c("first", "last", "ratio")], h8[,c("first", "last", "ratio")],h5[,c("first", "last", "ratio")], h6[,c("first", "last", "ratio")])
merge1 <- merge1 %>% reduce(full_join, by=c("first", "last"))
colnames(merge1) <- c("first", "last", "h3_ratio", "h4_ratio","h7_ratio", "h8_ratio","h5_ratio", "h6_ratio")

merge1$ratio <- round((merge1$h3_ratio+merge1$h4_ratio+merge1$h7_ratio+merge1$h8_ratio+merge1$h5_ratio+merge1$h6_ratio)/6,2)*100
mdf2 <- spread(merge1[,c("first", "last", "ratio")], key = last, value = ratio) %>%  as.data.frame()
final_matrix1 <- data.frame(mdf2[,-1], row.names=as.character(mdf2[,1])) %>% as.matrix(rownames = TRUE)

write.table(merge1, "plot/unmodified_nnunn2.context.table", quote = FALSE, col.names = T, row.names = F)


library(pheatmap)
pdf("plot/unmodified_nnunn2.ratio.pdf")
pheatmap(final_matrix1,cluster_cols = FALSE, cluster_rows = FALSE, display_numbers = TRUE, number_format="%.1f", fontsize_number = 12, color=colorRampPalette(c("white", "red"))(length(breaksList)), breaks = breaksList, legend_breaks = seq(0, 100, 50), legend_labels = c("0", "50", "100"), border_color = "grey", number_color="grey0")
dev.off()

pdf("plot/unmodified_nnunn2.ratio.pdf")
ggplot(merge1[,c("first", "last", "ratio")], aes(x = first, y = last, fill = ratio)) + geom_tile(color = "gray", size=0.1)+ coord_fixed()+ ##to make Square tiles
  geom_text(aes(label = ratio), color = "black", size = 6) +
  scale_fill_gradient2(low = "white", high = "red", limits=c(0,100))
dev.off()

