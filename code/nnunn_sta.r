#!/users/ludwig/cfo155/miniconda2/bin/Rscript --vanilla
#module load R/4.0.3-foss-2020b 
library(GenomicAlignments)
library(dplyr)
library(reshape2)
library(tidyr)
options <- commandArgs(trailingOnly = TRUE)
fn0 <- options[1]
fn <- gsub(".bam$","",fn0)
#write.table(stack1, "H3_2023Mar27_S3.trimmed_R2.nnunn.sort.txt", quote = FALSE, col.names = F, row.names = F)

motif <- function(input, ref){
    ref <- paste0(ref)
    data <- data.frame(input)
    colnames(data) <- c("stack1")
    data$first <- substring(data$stack1, 1, 2)
    data$mid <- substring(data$stack1, 3, 3)
    data$last <- substring(data$stack1, 4, 5)
    dat <- data %>% group_by(first,mid,last) %>% summarize(n = n()) ## sum(dat$n) 
    #dat2 <- dat[grep("\\+",dat$first),] ##
    #dat2 <- dat[grep("\\+",dat$first, invert = TRUE),]   or dat2 <- dat[!grepl("\\+",dat$first),] 
    dat1 <- dat[!grepl("\\+",dat$first)& !grepl("\\+",dat$mid) &!grepl("\\+",dat$last)& !grepl("\\-",dat$first)& !grepl("\\-",dat$mid) &!grepl("\\-",dat$last),] ###
    discard <- dat[grepl("\\+",dat$first) | grepl("\\+",dat$mid) | grepl("\\+",dat$last)  | grepl("\\-",dat$first) | grepl("\\-",dat$mid) | grepl("\\-",dat$last),] ###
    write.table(discard, paste0(fn, "_", ref,".discarded.table"), quote = FALSE, col.names = T, row.names = F)

    redat <- spread(dat1, key = mid, value = n)
    redat <- redat[!grepl("N",redat$first) & !grepl("N",redat$last),]
    redat[is.na(redat)] <- 0                                # Applying spread function
    redat$ratio <- redat$C/(redat$C+redat$T)
    redat1 <- redat %>% dplyr::select(first, last, ratio)
    write.table(redat, paste0(fn, "_", ref, ".contex.table"), quote = FALSE, col.names = T, row.names = F)

    final <- spread(redat1, key = last, value = ratio) %>%  as.data.frame()
    write.table(final, paste0(fn, "_", ref, ".motif.table"), quote = FALSE, col.names = T, row.names = F)
    final_matrix <- data.frame(final[,-1], row.names=as.character(final[,1])) %>% as.matrix(rownames = TRUE)
    library(pheatmap)
    pdf(paste0(fn, "_", ref, ".ratio.pdf"))
    pheatmap(final_matrix,cluster_cols = FALSE, cluster_rows = FALSE, display_numbers = TRUE, fontsize_number = 12)
    dev.off()

    redat$depth <- redat$C+redat$T
    redat2 <- redat %>% dplyr::select(first, last, depth)
    final2 <- spread(redat2, key = last, value = depth) %>%  as.data.frame()
    final_matrix2 <- data.frame(final2[,-1], row.names=as.character(final2[,1])) %>% as.matrix(rownames = TRUE)
    pdf(paste0(fn, "_", ref, ".depth.pdf"))
    pheatmap(final_matrix2,cluster_cols = FALSE, cluster_rows = FALSE, display_numbers = TRUE, number_format="%.0i",fontsize_number = 12)
    dev.off()

    redat$left1 <- gsub("AA|CA|GA|TA|AG|CG|GG|TG|AT|CT|GT|TT","ND",redat$first)
    redat$left2 <- gsub("AC|GC|TC","DC",redat$left1)
    redat$left <- gsub("CC","CC",redat$left2)


    redat$right1 <- gsub("AA|AC|AG|AT|GA|GC|GG|GT|TA|TC|TG|TT","DN",redat$last)
    redat$right2 <- gsub("CA|CG|CT","CD",redat$right1)
    redat$right <- gsub("CC","CC",redat$right2)

    redat_depth<- redat %>% as.data.frame() %>% dplyr::select(left, right, depth)
    redat_depth1 <- aggregate(depth ~ left+right, redat_depth, sum)
    redat_depth <- redat_depth1[order(redat_depth1$left,redat_depth1$right),]

    redat_C<- redat %>% as.data.frame() %>% dplyr::select(left, right, C)
    redat_C1 <- aggregate(C ~ left+right, redat_C, sum)
    redat_C<- redat_C1[order(redat_C1$left,redat_C1$right),]

    redat_C$ratio <- redat_C$C/redat_depth$depth
    redat_C <- redat_C %>% dplyr::select(left,right,ratio)
    final3 <- spread(redat_C, key = right, value = ratio) %>%  as.data.frame()
    final_matrix <- data.frame(final3[,-1], row.names=as.character(final3[,1])) %>% as.matrix(rownames = TRUE)
    library(pheatmap)
    pdf(paste0(fn, "_", ref, ".ratio1.pdf"))
    pheatmap(final_matrix,cluster_cols = FALSE, cluster_rows = FALSE, display_numbers = TRUE, fontsize_number = 12)
    dev.off()

    final4 <- spread(redat_depth, key = right, value = depth) %>%  as.data.frame()
    final_matrix <- data.frame(final4[,-1], row.names=as.character(final4[,1])) %>% as.matrix(rownames = TRUE)
    library(pheatmap)
    pdf(paste0(fn, "_", ref, ".depth1.pdf"))
    pheatmap(final_matrix,cluster_cols = FALSE, cluster_rows = FALSE, display_numbers = TRUE,number_format="%.0i", fontsize_number = 12)
    dev.off()
}

stack1 <- stackStringsFromBam(paste0(fn, ".bam"), param=GRanges("NNUNN1:13-17")) 
if(length(stack1)>10){
    motif(stack1, "NNUNN1")
}

stack2 <- stackStringsFromBam(paste0(fn, ".bam"), param=GRanges("NNUNN2:13-17")) 
if(length(stack2)>10){
    motif(stack2, "NNUNN2")
}
