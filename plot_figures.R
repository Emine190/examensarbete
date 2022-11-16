
install.packages("BiocManager")
install.packages("tidyverse")
#BiocManager::install("rtracklayer")
setwd("/media/god/jellyfish/Emil/R/")
BiocManager::install(version = "3.16")
library(BiocManager)
library(cowplot)
library(ggplot2)
library(data.table)
library(tidyverse)
library(rstatix)
library(EnvStats)
library(ggbeeswarm)
library(ggpubr)
library(karyoploteR)
library(rtracklayer)

#temp = list.files(pattern="*.vcf")
#myfiles = lapply(temp, read.delim)
#imports teh data of the three non-mosaic women
UPIC <- fread("het_UPIC.bed")
PLJ <- fread("het_13PLJ.bed")
ZZPU <- fread("het_ZZPU.bed")

#check which SNPs are in shared
kek1 <- intersect(UPIC$V2, PLJ$V2)
kek2 <- intersect(kek1, ZZPU$V2)

#Make a file with all of them together
bur <- rbind(UPIC, PLJ, ZZPU)

bur_filterOG <- as.data.frame(bur[bur$V2 %in% kek2,])

colnames(bur_filterOG) <- c("chr", "start", "ref", "alt", "depth", "sample")

bur_filterOG$end <- bur_filterOG$start + 1
bur_filterOG$width <- 1

shortz <- dplyr::select(bur_filterOG, chr, start, end, sample)

gr_bur_filterOG <- toGRanges(shortz)

kp <- plotKaryotype(plot.type=3, main="plot.type=3", genome = "hg38", chromosomes = "chrX")
kpPoints(karyoplot = kp , data=gr_bur_filterOG[gr_bur_filterOG$sample == "GTEX-UPIC-0004-SM-5SOEF",], chr =  seqnames(gr_bur_filterOG[gr_bur_filterOG$sample == "GTEX-UPIC-0004-SM-5SOEF",]), x = start(gr_bur_filterOG[gr_bur_filterOG$sample == "GTEX-UPIC-0004-SM-5SOEF",]), col = "blue", y=0.05)
kpPoints(karyoplot = kp , data=gr_bur_filterOG[gr_bur_filterOG$sample == "GTEX-13PLJ-0003-SM-6WSSCN",], chr =  seqnames(gr_bur_filterOG[gr_bur_filterOG$sample == "GTEX-13PLJ-0003-SM-6WSSCN",]), x = start(gr_bur_filterOG[gr_bur_filterOG$sample == "GTEX-13PLJ-0003-SM-6WSSCN",]), col = "red", y=0.15)
kpPoints(karyoplot = kp , data=gr_bur_filterOG[gr_bur_filterOG$sample == "GTEX-ZZPU-0003-SM-6WBUC",], chr =  seqnames(gr_bur_filterOG[gr_bur_filterOG$sample == "GTEX-ZZPU-0003-SM-6WBUC",]), x = start(gr_bur_filterOG[gr_bur_filterOG$sample == "GTEX-ZZPU-0003-SM-6WBUC",]), col = "grey", y=0.25)
kpAxis()


#Import the Controls
Control <- fread("het_Control.bed")
#apply the control to solo out the one unique for the 3 women.
kek3 <- setdiff(kek2, Control$V2)


Unique_snp <- as.data.frame(bur[bur$V2 %in% kek3,])


colnames(Unique_snp) <- c("chr", "start", "ref", "alt", "depth", "sample")

Unique_snp$end <- Unique_snp$start + 1
Unique_snp$width <- 1

thingforgr <- dplyr::select(Unique_snp, chr, start, end, sample)

gr_Unique_snp <- toGRanges(thingforgr)

#buren <- data.frame(bur[bur$V2 %in% kek3,])
#mcols(buren) <- cbind(mcols(buren), newMcols)

kp <- plotKaryotype(plot.type=3, main="plot.type=3", genome = "hg38", chromosomes = "chrX")
  kpPoints(karyoplot = kp , data=gr_bur_filterOG, chr =  seqnames(gr_bur_filterOG), x = start(gr_bur_filterOG), col = "blue", y=0.05)
  kpPoints(karyoplot = kp , data=gr_Unique_snp, chr =  seqnames(gr_Unique_snp), x = start(gr_Unique_snp), col = "red", y=0.15)
  kpAxis()
  kpPlotMarkers(karyoplot = kp, data=gr_Unique_snp, chr = seqnames(gr_Unique_snp), y=0.3, labels=)

