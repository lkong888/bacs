#!/users/ludwig/cfo155/miniconda2/bin/Rscript --vanilla

# Specify a job name
#$ -N fdr

# Project name and target queue
#$ -P ludwig.prjc
#$ -q short.qc

# Run the job in the current working directory
#$ -cwd -j y

suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(parallel))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(RColorBrewer))

options <- commandArgs(trailingOnly = TRUE)
fn <- options[1]
outputname <- options[2]

# hmC site ----
site <- fread(fn,  header = FALSE, nThread = 10, stringsAsFactors = FALSE, data.table = FALSE, sep = "\t", verbose = T)
colnames(site) <- c("chr", "start", "end", "strand", "ref", "motif", "treat_depth", "treat_T", "treat_C", "treat_conversion", "treat_gap", "treat_gap_rate", "control_ref","ctrl_depth", "ctrl_T", "ctrl_C", "ctrl_conversion", "ctrl_gap", "ctrl_gap_rate","conversion", "fp", "treat_level", "p", "fdr")

site$fdr <- p.adjust(site$p, method = "BH")
print(paste0("nsite", nrow(site)))
write.table(site, outputname , quote=FALSE, sep="\t", row.names = FALSE, col.names = TRUE)