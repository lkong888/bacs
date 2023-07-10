#!/users/ludwig/cfo155/miniconda2/bin/Rscript --vanilla
options <- commandArgs(trailingOnly = TRUE)

.libPaths("/well/ludwig/users/ebu571/R/4.0/skylake")
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(parallel))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(RColorBrewer))
library(readr)
options(bitmapType='cairo-png')
library(matrixStats)

plot_mutation <- function(input){

    dat <- read.table(input, header=TRUE)
    dat$total_insert <- dat$Insert+dat$insert
    dat$total_del <- dat$Delete+dat$delete
    dat$total_A <- dat$A+dat$a
    dat$total_C <- dat$C+dat$c
    dat$total_G <- dat$G+dat$g
    dat$total_T <- dat$T+dat$t
    dat$depth <- dat$total_del+dat$total_A+dat$total_C+dat$total_G+dat$total_T
    dat1 <- dat[,c("chr", "pos", "ref_base", "depth", "total_insert","total_del", "total_A", "total_C", "total_G", "total_T")]
    colnames(dat1) <- c("reference", "pos", "ref", "depth", "insert","del", "A", "C", "G", "T")

    datA <- dat1[dat1$reference=="2kb_spikein" & dat1$ref=="A",]
    datC <- dat1[dat1$reference=="2kb_spikein" & dat1$ref=="C",]
    datG <- dat1[dat1$reference=="2kb_spikein" & dat1$ref=="G",]
    datT <- dat1[dat1$reference=="2kb_spikein" & dat1$ref=="T",]
    ###for mean and median statistics
    sta <-c()

    #### ref: A
    datA$Adel_ratio <- datA$del/datA$depth
    datA$Avariant <- (datA$del + datA$C + datA$G +datA$T)
    datA$Avariant_ratio <- datA$Avariant/datA$depth
    datA$Amutation_ratio <- 1- (datA$A/(datA$A+datA$C+datA$G+datA$T))
    datA$AtoA <- datA$A/(datA$A+datA$C+datA$G+datA$T)
    datA$AtoC <- datA$C/(datA$A+datA$C+datA$G+datA$T)
    datA$AtoG <- datA$G/(datA$A+datA$C+datA$G+datA$T)
    datA$AtoT <- datA$T/(datA$A+datA$C+datA$G+datA$T)
    datA1 <- datA[,c("reference", "pos", "Adel_ratio", "Avariant_ratio", "Amutation_ratio", "AtoC", "AtoG", "AtoT")] %>% melt(id.vars=c("reference", "pos"))
    sta$meanA <- colMeans(datA[,c("Adel_ratio", "Avariant_ratio", "Amutation_ratio", "AtoA", "AtoG", "AtoC", "AtoT")])
    sta$medianA <- apply(datA[,c("Adel_ratio", "Avariant_ratio", "Amutation_ratio", "AtoA", "AtoG", "AtoC", "AtoT")],2,median)

    #### ref: G
    datG$Gdel_ratio <- datG$del/datG$depth
    datG$Gvariant <- datG$del + datG$A + datG$C +datG$T
    datG$Gvariant_ratio <- datG$Gvariant/datG$depth
    datG$Gmutation_ratio <- 1- (datG$G/(datG$A+datG$C+datG$G+datG$T))
    datG$GtoA <- datG$A/(datG$A+datG$C+datG$G+datG$T)
    datG$GtoC <- datG$C/(datG$A+datG$C+datG$G+datG$T)
    datG$GtoG <- datG$G/(datG$A+datG$C+datG$G+datG$T)
    datG$GtoT <- datG$T/(datG$A+datG$C+datG$G+datG$T)
    datG1 <- datG[,c("reference", "pos", "Gdel_ratio", "Gvariant_ratio", "Gmutation_ratio", "GtoA", "GtoC", "GtoT")] %>% melt(id.vars=c("reference", "pos"))
    sta$meanG <- colMeans(datG[,c("Gdel_ratio", "Gvariant_ratio", "Gmutation_ratio", "GtoA", "GtoC", "GtoG", "GtoT")])
    sta$medianG <- apply(datG[,c("Gdel_ratio", "Gvariant_ratio", "Gmutation_ratio", "GtoA", "GtoC", "GtoG", "GtoT")],2,median)

    #### ref: C
    datC$Cdel_ratio <- datC$del/datC$depth
    datC$Cvariant <- datC$del + datC$A + datC$G +datC$T
    datC$Cvariant_ratio <- datC$Cvariant/datC$depth
    datC$Cmutation_ratio <- 1- (datC$C/(datC$A+datC$C+datC$G+datC$T))
    datC$CtoA <- datC$A/(datC$A+datC$C+datC$G+datC$T)
    datC$CtoC <- datC$C/(datC$A+datC$C+datC$G+datC$T)
    datC$CtoG <- datC$G/(datC$A+datC$C+datC$G+datC$T)
    datC$CtoT <- datC$T/(datC$A+datC$C+datC$G+datC$T)
    datC1 <- datC[,c("reference", "pos", "Cdel_ratio", "Cvariant_ratio", "Cmutation_ratio", "CtoA", "CtoG", "CtoT")] %>% melt(id.vars=c("reference", "pos"))
    sta$meanC <- colMeans(datC[,c("Cdel_ratio", "Cvariant_ratio", "Cmutation_ratio", "CtoA", "CtoC", "CtoG", "CtoT")])
    sta$medianC <- apply(datC[,c("Cdel_ratio", "Cvariant_ratio", "Cmutation_ratio", "CtoA", "CtoC", "CtoG", "CtoT")],2,median)

    #### ref: T
    datT$Tdel_ratio <- datT$del/datT$depth
    datT$Tvariant <- datT$del + datT$A + datT$G +datT$C
    datT$Tvariant_ratio <- datT$Tvariant/datT$depth
    datT$Tmutation_ratio <- 1- (datT$T/(datT$A+datT$C+datT$G+datT$T))
    datT$TtoA <- datT$A/(datT$A+datT$C+datT$G+datT$T)
    datT$TtoC <- datT$C/(datT$A+datT$C+datT$G+datT$T)
    datT$TtoG <- datT$G/(datT$A+datT$C+datT$G+datT$T)
    datT$TtoT <- datT$T/(datT$A+datT$C+datT$G+datT$T)
    datT1 <- datT[,c("reference", "pos", "Tdel_ratio", "Tvariant_ratio", "Tmutation_ratio", "TtoA","TtoC", "TtoG")] %>% melt(id.vars=c("reference", "pos"))
    sta$meanT <- colMeans(datT[,c("Tdel_ratio", "Tvariant_ratio", "Tmutation_ratio", "TtoA", "TtoC", "TtoG","TtoT")])
    sta$medianT <- apply(datT[,c("Tdel_ratio", "Tvariant_ratio", "Tmutation_ratio", "TtoA", "TtoC", "TtoG","TtoT")],2,median)

     write.table(sta, file=paste0(gsub(".txt","",input),".sta.txt"), quote = FALSE, col.names = T, row.names = F)

    colnames(datA) <- c("reference", "pos", "ref", "depth", "insert","del", "A", "C", "G", "T", "del_ratio", "variant", "variant_ratio", "mutation_ratio", "toA", "toC", "toG", "toT")
    colnames(datG) <- c("reference", "pos", "ref", "depth", "insert","del", "A", "C", "G", "T", "del_ratio", "variant", "variant_ratio", "mutation_ratio", "toA", "toC", "toG", "toT")
    colnames(datC) <- c("reference", "pos", "ref", "depth", "insert","del", "A", "C", "G", "T", "del_ratio", "variant", "variant_ratio", "mutation_ratio", "toA", "toC", "toG", "toT")
    colnames(datT) <- c("reference", "pos", "ref", "depth", "insert","del", "A", "C", "G", "T", "del_ratio", "variant", "variant_ratio", "mutation_ratio", "toA", "toC", "toG", "toT")
    output <- rbind(datA, datG, datC, datT)
    write.table(output, file=paste0(gsub(".txt","",input),".2kb_mutation_sta.txt"), quote = FALSE, col.names = T, row.names = F)

    #output <- rbind(datA,setNames(datC, names(datA)))
    #output <- rbind(output,setNames(datG, names(output)))
    #output <- rbind(output,setNames(datT, names(output)))
    #output <- rbind(output,setNames(dat_pseu, names(output)))
    #write.table(output, "H5_14Dec2022_S5.fastp.filter.sort.mutation_sta1.txt", quote = FALSE, col.names = F, row.names = F)

    alldat <- rbind(datA1,datC1,datG1,datT1)
    alldat$mut_type <- factor(alldat$variable, levels=c("Adel_ratio", "Avariant_ratio", "Amutation_ratio", "AtoG", "AtoC", "AtoT","Gdel_ratio","Gvariant_ratio", "Gmutation_ratio", "GtoA", "GtoC", "GtoT","Cdel_ratio", "Cvariant_ratio", "Cmutation_ratio", "CtoA", "CtoG", "CtoT","Tdel_ratio", "Tvariant_ratio", "Tmutation_ratio", "TtoA", "TtoG", "TtoC"))
    write.table(alldat, file=paste0(gsub(".txt","",input),".2kb_mutation_sta_plot.txt"), quote = FALSE, col.names = T, row.names = F)

    p <- ggplot(alldat, aes(x=mut_type,y=value, fill=mut_type))+ geom_boxplot()
    ggsave(file=paste0(gsub(".txt","",input),".2kb_mutation_sta_box.pdf"),p, width=30,height = 5)
}

treatment <- options[1]
plot_mutation(treatment)




####Combine Plot
#H7 <- read.table("align/H7_2023May19_S7.rRNA.filter.dedup.mpile.2kb_mutation_sta_plot.txt", header=TRUE)
#H9 <- read.table("align/H9_2023May19_S9.rRNA.filter.dedup.mpile.2kb_mutation_sta_plot.txt", header=TRUE)
#mergedat <- merge(H7,H9, by=c("reference", "pos", "mut_type")) 
#colnames(mergedat) <- c("reference", "pos", "mut_type", "treat_variable","treat_value", "ctrl_variable","ctrl_value")
#mergedat$mut_type <- factor(mergedat$mut_type, levels=c("Adel_ratio", "Avariant_ratio", "Amutation_ratio", "AtoG", "AtoC", "AtoT","Gdel_ratio","Gvariant_ratio", "Gmutation_ratio", "GtoA", "GtoC", "GtoT","Cdel_ratio", "Cvariant_ratio", "Cmutation_ratio", "CtoA", "CtoG", "CtoT","Tdel_ratio", "Tvariant_ratio", "Tmutation_ratio", "TtoA", "TtoG", "TtoC", "pseu_del_ratio", "pseu_variant_ratio", "pseu_mutation_ratio", "PSEUtoA", "PSEUtoG", "PSEUtoC"))
#seldat <- mergedat %>% dplyr::select("mut_type", "treat_value", "ctrl_value") %>% melt(id.vars=c("mut_type"))
#p <- ggplot(seldat, aes(x=mut_type, y=value*100, fill=variable))+ geom_boxplot()+ylim(0,50)+theme_bw()
#ggsave("plot/2kb_mutation_sta_box.pdf",p, width=30,height = 5)

