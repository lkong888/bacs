##module load R/4.0.3-foss-2020b
library(dplyr)
library(reshape2)
library(tidyr)
library(ggplot2)
library(cowplot)
library(tidyverse)
library(parallel)
setwd("/users/ludwig/ebu571/ebu571/20May2023_final")

#####get the tRNA indexing
trnaalign <- readLines("/users/ludwig/ebu571/ebu571/test/hg38-trnaalign.stk")
output_file <- file("resource/hg38-trnaalign.position.txt", "w")

for (i in 1:length(trnaalign)) {
    line <- trnaalign[[i]]
    if (grepl("^tRNA", line)){
    fields <- strsplit(line, "\\s+")[[1]]
    chr <- fields[1]
    sequence <- sub("-", "",fields[2])

    for (n in 1:nchar(sequence)){
        position <- n
        ref <- substring(sequence, n, n)
        output <- paste(chr, position, toupper(ref), sep=" ")
        cat(output, "\n", file = output_file, append=TRUE)
    }
    }
}

close(output_file)

###cat /users/ludwig/ebu571/ebu571/20May2023_final/resource/tRNA_high_confident_tmp.fa | paste - - | sed 's/>//g' | sed 's/ /\t/g' > resource/tRNA_high_confident.txt
trnaalign <- readLines("resource/tRNA_high_confident.txt")
output_file <- file("resource/trna.position.txt", "w")

for (i in 1:length(trnaalign)) {
    line <- trnaalign[[i]]
    if (grepl("^tRNA", line)){
    fields <- strsplit(line, "\\s+")[[1]]
    chr <- fields[1]
    sequence <- sub("-", "",fields[2])

    for (n in 1:nchar(sequence)){
        position <- n
        ref <- substring(sequence, n, n)
        output <- paste(chr, position, toupper(ref), sep=" ")
        cat(output, "\n", file = output_file, append=TRUE)
    }
    }
}

close(output_file)



#####join -1 2 -2 1 <(cat hg38-trnaalign.position.txt |sed 's/ /\t/g' | awk '$3=="A" || $3=="C" || $3=="G" || $3=="U"' | awk '{OFS="\t"}{print $1, $2,$3}'| sort -k2) <(cat /users/ludwig/ebu571/ebu571/test/hg38-alignnum.txt | grep -v ^0|cut -f1,4| sort -k1) | awk '{OFS="\t"}{print $2, $1, $3, $4}' | sort -k 1,1 -k2,2n > tRNA.index.txt

###tRNA_Ser_AGA keep 1 and 4
library(Biostrings)

# Read sequences from the FASTA file
sequences <- readDNAStringSet("resource/tRNA_high_confident.fa")
output_file <- file("resource/tRNA_high_confident_pairwise.txt", "w")

for (i in 1:(length(sequences)-1)) {
    seq1 <- sequences[[i]]
    names1 <- names(sequences)[[i]]

  for (j in (i+1):length(sequences)) {
    seq2 <- sequences[[j]]
    names2 <- names(sequences)[[j]]

    alignment <- pairwiseAlignment(seq1, seq2)
    nedit <- nedit(alignment)
    if (nedit<5){
        output <- paste(names1, seq1, names2, seq2, nedit,sep=" ")
        cat(output, "\n", file = output_file, append=TRUE)
    }
}
}


###do filtering
aligndat <- read.table("resource/tRNA_high_confident_pairwise.txt", stringsAsFactors = FALSE)
colnames(aligndat) <- c("seq1_name", "seq1", "seq2_name", "seq2", "nedit")

dat1 <- aligndat[aligndat$nedit==1 & sub("tRNA_(\\w+_\\w+)_\\w+", "\\1", aligndat$seq1_name) == sub("tRNA_(\\w+_\\w+)_\\w+", "\\1", aligndat$seq2_name),]##90
anticodon <- unique(sub("tRNA_(\\w+_\\w+)_\\w+", "\\1", dat1$seq1_name))
for (i in anticodon){
    dat <- dat1[sub("tRNA_(\\w+_\\w+)_\\w+", "\\1", dat1$seq1_name)==i,]
    name <- c(dat$seq1_name, dat$seq2_name)

    max_frequency <- max(table(name))
    if (max_frequency <2){
        if(length(dat$seq1)>length(dat$seq2)){
            filtered=

        }
    }


}



###########################count TPM for tRNA (Mim-tRNA-seq output)
a <- read.table("featureCounts/input_feature.txt",header=F, sep="\t")
colnames(a) <- c("parent", "isodecoder", "h9", "h10", "length")
tRNA <- a[!grepl("mito", a$isodecoder) & !grepl("Escherichia", a$isodecoder),]

tRNA[,3] = tRNA[,3]/tRNA$length *10**3/sum(tRNA[,3]/tRNA$length *10**3) *10**6
tRNA[,4] = tRNA[,4]/tRNA$length *10**3/sum(tRNA[,4]/tRNA$length *10**3) *10**6
write.table(tRNA,"featureCounts/input_feature_tpm.txt",sep='\t',quote=F,row.names=F)

tRNA1 <- tRNA[!grepl("tRX",tRNA$isodecoder) & tRNA$h9>0 & tRNA$h10>0,] ##291

summary(tRNA$h9)
#
   Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
     0.00      0.81     23.65   2577.32   1020.07 125849.68

summary(tRNA$h10)
     Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
     0.00      0.63     23.96   2577.32    945.89 134816.26 



> summary(tRNA1$h9)
     Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
     0.22      6.81    146.04   3428.11   2036.09 125849.68

> summary(tRNA1$h10)
     Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
     0.21      6.34    132.30   3427.89   1856.70 134816.26 



pdf("processed_plot/h9_input_TPM_distribution.pdf")
hist(tRNA$h9)
dev.off()


tRNA1$h9_range <- cut(tRNA1$h9, breaks=c(-1,50,100,200,400,800,1000,150000), labels=c("level1", "level2", "level3","level4", "level5", "level6","level7"))
dat1 <- tRNA1 %>% group_by(h9_range) %>% summarize(n = n())

nrow(tRNA1[tRNA1$h9>=146.04,]) ##145
nrow(tRNA1[tRNA1$h10>=132.30,]) ##145

nrow(merge(tRNA1[tRNA1$h9>=146.04,], tRNA1[tRNA1$h10>=132.30,], by=c("isodecoder"))) ##144


#################################get the high-expressed isodecoder list
######raw read counts > 1% of total reads in the same anticodon catogary
###cat featureCounts/input1.counts | grep -v ^# |cut -f2,7,8 | grep tRNA | grep -v mt >featureCounts/tRNA.counts.txt
counts <- read.table("featureCounts/tRNA.counts.txt",header=T, sep="\t")
colnames(counts) <- c("chr", "input1", "input2")
counts$tmp <- sub("tRNA_(\\w+_\\w+)_(\\w+_\\w+)", "\\1", counts$chr)
counts$anticodon <- sub("tRNA_(\\w+_\\w+)_\\w+", "\\1", counts$tmp)


counts_agg1 <- aggregate(input1~anticodon, counts, sum)
colnames(counts_agg1) <- c("anticodon", "input1_sum")
counts_agg2 <- aggregate(input2~anticodon, counts, sum)
colnames(counts_agg2) <- c("anticodon", "input2_sum")
counts_merge <- merge(counts, counts_agg1, by=c("anticodon"))
counts_merge <- merge(counts_merge, counts_agg2, by=c("anticodon"))

excludelist <- counts_merge[counts_merge$input1 <= (counts_merge$input1_sum*0.02) &counts_merge$input2 <= (counts_merge$input2_sum*0.02), ] ##60
filteredlist <- counts_merge[counts_merge$input1 > (counts_merge$input1_sum*0.02) | counts_merge$input2 > (counts_merge$input2_sum*0.02), ] ##200


write.table(excludelist,"featureCounts/tRNA_counts_lessthan1per.txt",sep='\t',quote=F,row.names=F)
write.table(filteredlist,"featureCounts/tRNA_counts_filtered.txt",sep='\t',quote=F,row.names=F)


#############################tRNA filtered####################################
###fdr filtered and counts filtered 
###10% level in at least two replicates
site <- read.table("site/snoRNA_tRNA_called_sites_fdrpassed.txt",header=T, sep="\t",colClasses=c("ref"="character")) ###1100
tRNA_site <- merge(site, filteredlist[,c("chr", "input1","input2", "input1_sum", "input2_sum")], by=c("chr")) ##650
tRNA_site1 <- tRNA_site[(tRNA_site$treat1_level >= 0.10 & tRNA_site$treat2_level >= 0.10 ) | (tRNA_site$treat1_level >= 0.10 & tRNA_site$treat3_level >= 0.10) | (tRNA_site$treat2_level >= 0.10 & tRNA_site$treat3_level >= 0.10), ]
##tRNA_site1 <- tRNA_site[(tRNA_site$treat1_level >= 0.10 & tRNA_site$treat2_level >= 0.10 &tRNA_site$treat3_level < 0.10) | (tRNA_site$treat1_level >= 0.10 & tRNA_site$treat3_level >= 0.10 & tRNA_site$treat2_level < 0.10) | (tRNA_site$treat2_level >= 0.10 & tRNA_site$treat3_level >= 0.10 & tRNA_site$treat1_level < 0.10) |  (tRNA_site$treat2_level >= 0.10 & tRNA_site$treat3_level >= 0.10 & tRNA_site$treat1_level >= 0.10), ]
nrow(tRNA_site1) ##625
exclude_site <- tRNA_site[(tRNA_site$treat1_level < 0.10 & tRNA_site$treat2_level < 0.10 & tRNA_site$treat3_level < 0.10) | (tRNA_site$treat1_level >= 0.10 & tRNA_site$treat2_level < 0.10 & tRNA_site$treat3_level < 0.10) | (tRNA_site$treat1_level < 0.10 & tRNA_site$treat2_level >= 0.10 & tRNA_site$treat3_level < 0.10) | (tRNA_site$treat1_level < 0.10 & tRNA_site$treat2_level < 0.10 & tRNA_site$treat3_level >= 0.10), ]

nrow(exclude_site) ##25
write.table(tRNA_site1,"site/snoRNA_tRNA_called_sites_fdrpassed_tRNA_filtered.txt",sep='\t',quote=F,row.names=F)
write.table(exclude_site,"site/snoRNA_tRNA_called_sites_fdrpassed_tRNA_10%excludesite.txt",sep='\t',quote=F,row.names=F)



##################get tRNA table 
##cat resource/tRNA_table.txt | sed 's/ /_/g' | sed 's/?/U/g' | sed 's/?1/U/g' | cut -f1,6-$NF | sed 's/G\tG\tG\tC\tU\t_\tU\tG\tA\tA/G\tG\tG\tC\tU\tN\tU\tG\tA\tA/g' >resource/tRNA_table1.txt
tRNA <- read.table("resource/tRNA_table1.txt", header=T, sep="\t",colClasses="character")
melttrna <- melt(tRNA, id.vars="Full_name")
isodecoder <- unique(melttrna$Full_name)

final <- data.frame(matrix(ncol = 4, nrow = 0))
colnames(final) <- c("Full_name", "variable", "value", "V4")
for (i in isodecoder){
    dat <- melttrna[melttrna$Full_name==i,]
    x <- 1
    for (j in 1:nrow(dat)){
        if (dat[j,]$value=="A" | dat[j,]$value=="C" | dat[j,]$value=="G" | dat[j,]$value=="U" | dat[j,]$value=="Ul" | dat[j,]$value=="N"){
        dat[j,4] <- x
        x <- x+1
        }
    }
    final <- rbind(final,dat)
}
nrow(final) ###25308
nrow(melttrna) ###25308

write.table(final,"resource/tRNA_position_table.txt",sep='\t',quote=F,row.names=F)

##cat resource/tRNA_position_table.txt | sort -k 1 | grep -v X_ | sed 's/X//g' | sed 's/\.1/0/g' | sed 's/Ul/U/g' | sed 's/Ul/U/g' > resource/tRNA_position_index.txt
