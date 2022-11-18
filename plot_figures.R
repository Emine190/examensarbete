
install.packages("BiocManager")
install.packages("tidyverse")
#BiocManager::install("rtracklayer")
setwd("/media/MY/files/R/")
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
Woman1 <- fread("het_Woman1.bed")
Woman2 <- fread("het_Woman2.bed")
Woman3 <- fread("het_Woamn3.bed")

#check which SNPs are in shared between all three women
kek1 <- intersect(Woman1$V2, Woman2$V2)
kek2 <- intersect(kek1, Woman3$V2)

#Make a file with all of them together
bur <- rbind(Woman1, Woman2, Woman3)

bur_filterOG <- as.data.frame(bur[bur$V2 %in% kek2,])

colnames(bur_filterOG) <- c("chr", "start", "ref", "alt", "depth", "sample")

bur_filterOG$end <- bur_filterOG$start + 1
bur_filterOG$width <- 1

shortz <- dplyr::select(bur_filterOG, chr, start, end, sample)

gr_bur_filterOG <- toGRanges(shortz)

kp <- plotKaryotype(plot.type=3, main="plot.type=3", genome = "hg38", chromosomes = "chrX")
kpPoints(karyoplot = kp , data=gr_bur_filterOG[gr_bur_filterOG$sample == "Woman1",], chr =  seqnames(gr_bur_filterOG[gr_bur_filterOG$sample == "Woman1",]), x = start(gr_bur_filterOG[gr_bur_filterOG$sample == "GTEX-UPIC-0004-SM-5SOEF",]), col = "blue", y=0.05)
kpPoints(karyoplot = kp , data=gr_bur_filterOG[gr_bur_filterOG$sample == "Woman2",], chr =  seqnames(gr_bur_filterOG[gr_bur_filterOG$sample == "Woman2",]), x = start(gr_bur_filterOG[gr_bur_filterOG$sample == "GTEX-13PLJ-0003-SM-6WSSCN",]), col = "red", y=0.15)
kpPoints(karyoplot = kp , data=gr_bur_filterOG[gr_bur_filterOG$sample == "Woman3",], chr =  seqnames(gr_bur_filterOG[gr_bur_filterOG$sample == "Woman3",]), x = start(gr_bur_filterOG[gr_bur_filterOG$sample == "GTEX-ZZPU-0003-SM-6WBUC",]), col = "grey", y=0.25)




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

kp <- plotKaryotype(plot.type=4, main="plot.type=4", genome = "hg38", chromosomes = "chrX")
kpPlotRainfall(kp, data=gr_Unique_snp[gr_Unique_snp$sample == "Woman1",] , r0 = 0.1, r1=0.2)


kp <- plotKaryotype(plot.type=3, main="plot.type=3", genome = "hg38", chromosomes = "chrX")
  kpPoints(karyoplot = kp , data=gr_bur_filterOG, chr =  seqnames(gr_bur_filterOG), x = start(gr_bur_filterOG), col = "blue", y=0.05)
  kpPoints(karyoplot = kp , data=gr_Unique_snp, chr =  seqnames(gr_Unique_snp), x = start(gr_Unique_snp), col = "red", y=0.15)

 # try and plot the ACGT snp to have a show of the normal split of it. 
  
  
  
#Import the bed files from the COV files 
# first set a new working directory so the importing of the files gets easier.
  setwd("/media/MY/files/R/sv/cov/")

  Woman1_cov <- fread("Woman1_cov..bed")
  Woman1_cov_log <- Woman1_cov
  Woman1_cov_log$coverage <- log10(Woman1_cov$coverage +1)
  Woman1_COV_gr <- toGRanges(Woman1_cov)
  Woman1_cov_loggr <- toGRanges(Woman1_cov_log)
  
  Woman2_COV <-fread("Woman2_cov..bed")
  Woman2_COV_GR <- toGRanges(Woman2_COV)
  Woman2_cov_log <- Woman2_COV
  Woman2_cov_log$coverage <- log10(Woman2_COV$coverage)
  Woman2_cov_loggr <- toGRanges(Woman2_cov_log)
  
 Woman3_COV <- fread("Woman3_cov..bed")
  Woman3_COV_GR <- toGRanges(Woman3_COV)
  Woman3_cov_log <- Woman3_COV
  Woman3_cov_log$coverage <- log10(Woman3_COV$coverage)
  Woman3_cov_loggr <- toGRanges(Woman3_cov_log)
  
  #Creates a list for importing the COV files 

  #ListforCOV = list.files(pattern="*.bed")
#covfiles = lapply(ListforCOV, read.delim)


#Goal is a plot to show the coverage of the whole genome and then plot it against the Control


kp_covUWoman1 <- plotKaryotype(plot.type = 2, main="plot of COV for Woman1", genome = "hg38", chromosomes =c("chr17"))
#kpPlotCoverage(kp_covWoman1, data=Woman1_cov_loggr, show.0.cov = TRUE, data.panel=1, r0=0, r1=0.2, col="red", border = "blue")
kpPoints(karyoplot = kp_covWoman1 , data=Woman1_cov_loggr, chr =  seqnames(Woman1_cov_loggr), x = start(Woman1_cov_loggr), col = "black", y=Woman1_cov_loggr$coverage)

kp_covWoman2 <- plotKaryotype(plot.type = 2, main="plot of COV for Woman2", genome = "hg38", chromosomes =c("chr17","chr16", "chrX"))
kpPlotCoverage(kp_covWoman2, data=Woman2_cov_loggr, show.0.cov = TRUE, data.panel=1, r0=0, r1=0.2, col="red", border="blue")

kp_covWoman3 <- plotKaryotype(plot.type = 2, main="plot of COV for Woman3", genome = "hg38", chromosomes =c("chr17","chr16", "chrX"))
kpPlotCoverage(kp_covWoman3, data=Woman3_cov_loggr, show.0.cov = TRUE, data.panel=1, r0=0, r1=0.2, col="red", border="blue")

kp_Region <- plotKaryotype("hg38", plot.type=2, chromosomes = c("chrX"))
kpPlotRegions(kp_Region, data = Woman1_cov_loggr,r0=0, r1=0.1, col="red")






#the SV files to be read in.
setwd("/media/MY/files/R")

Woman1_cov5k <- fread("Woman1_cov_5K.bed")
Woman1_cov5k_log <- Woman1_cov5k
Woman1_cov5k_log$coverage <- log10(Woman1_cov5k$coverage +1)
Woman1_COV5k_gr <- toGRanges(Woman1_cov5k)
Woman1_cov5k_loggr <- toGRanges(Woman1_cov5k_log)

Woman2_COV5k <-fread("Woman2_cov_5K.bed")
Woman2_COV5k_GR <- toGRanges(Woman2_COV5k)
Woman2_cov5k_log <- Woman2_COV5k
Woman2_cov5k_log$coverage <- log10(Woman2_COV5k$coverage +1)
Woman2_cov5k_loggr <- toGRanges(Woman2_cov5k_log)

ZZPU_COV5k <- fread("GTEX-ZZPU_cov_5K.bed")
ZZPU_COV5k_GR <- toGRanges(ZZPU_COV5k)
ZZPU_cov5k_log <- ZZPU_COV5k
ZZPU_cov5k_log$coverage <- log10(ZZPU_COV5k$coverage +1)
ZZPU_cov5k_loggr <- toGRanges(ZZPU_cov_log)


kp_cov <- plotKaryotype(plot.type = 4, main="plot of COV for Woman1", genome = "hg38", chromosomes =c("chr17"))
#kpPlotCoverage(kp_cov, data=UPIC_cov5k_loggr, show.0.cov = TRUE, data.panel=1, r0=0, r1=0.2, col="red", border = "blue")
kpPoints(karyoplot = kp_cov , data=Woman1_cov5k_loggr, chr =  seqnames(Woman1_cov5k_loggr), x = start(Woman1_cov5k_loggr), col = "black", y=Woman1_cov5k_loggr$coverage)




