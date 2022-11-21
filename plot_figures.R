
install.packages("BiocManager")
install.packages("tidyverse")
install.packages("SNPfiltR")
install.packages("vcfR")
#BiocManager::install("rtracklayer")
setwd("/media/god/jellyfish/Emil/R/")
BiocManager::install(version = "3.16")
library(BiocManager)
library(SNPfiltR)
library(vcfR)
library(cowplot)
library(ggplot2)
library(data.table)
library(tidyr)
library(stringr)
library(ggbeeswarm)
library(ggpubr)
library(karyoploteR)
library(rtracklayer)
library(dplyr)

#temp = list.files(pattern="*.vcf")
#myfiles = lapply(temp, read.delim)
#imports the data of the three non-mosaic women
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

kp <- plotKaryotype(plot.type=3, main="Blue= Woman1,RED ", genome = "hg38", chromosomes = "chrX")
kpPoints(karyoplot = kp , data=gr_bur_filterOG[gr_bur_filterOG$sample == "GTEX-UPIC-0004-SM-5SOEF",], chr =  seqnames(gr_bur_filterOG[gr_bur_filterOG$sample == "GTEX-UPIC-0004-SM-5SOEF",]), x = start(gr_bur_filterOG[gr_bur_filterOG$sample == "GTEX-UPIC-0004-SM-5SOEF",]), col = "blue", y=0.05)
kpPoints(karyoplot = kp , data=gr_bur_filterOG[gr_bur_filterOG$sample == "GTEX-13PLJ-0003-SM-6WSSCN",], chr =  seqnames(gr_bur_filterOG[gr_bur_filterOG$sample == "GTEX-13PLJ-0003-SM-6WSSCN",]), x = start(gr_bur_filterOG[gr_bur_filterOG$sample == "GTEX-13PLJ-0003-SM-6WSSCN",]), col = "red", y=0.15)
kpPoints(karyoplot = kp , data=gr_bur_filterOG[gr_bur_filterOG$sample == "GTEX-ZZPU-0003-SM-6WBUC",], chr =  seqnames(gr_bur_filterOG[gr_bur_filterOG$sample == "GTEX-ZZPU-0003-SM-6WBUC",]), x = start(gr_bur_filterOG[gr_bur_filterOG$sample == "GTEX-ZZPU-0003-SM-6WBUC",]), col = "grey", y=0.25)


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

# Redoing and splitting the depth so all zeros are removed.
bur_split <- Unique_snp %>% separate(depth, c('sampleX', 'refX'))
BurUNiqueNoZero <- filter(bur_split, sampleX > 0, refX > 0)

NoZerogr <- dplyr::select(BurUNiqueNoZero, chr, start, end, sample)


gr_NoZero_snp <- toGRanges(NoZerogr)

#PLot of all the different SNPs with the first being no zeros  

kp <- plotKaryotype(plot.type=3, main="SNPs", genome = "hg38", chromosomes = "chrX")
kpPlotRainfall(kp, data=gr_NoZero_snp, r0=0.25, r1=0.5, col="blue")
kpPoints(karyoplot = kp , data=gr_NoZero_snp, chr =  seqnames(gr_NoZero_snp), x = start(gr_NoZero_snp), col = "blue", y=-0.15)

kp <- plotKaryotype(plot.type=3, main="Unique SNPs", genome = "hg38", chromosomes = "chrX")
kpPlotRainfall(kp, data=gr_Unique_snp , r0 = 0.1, r1=0.3, col = "red")
kpPoints(karyoplot = kp , data=gr_Unique_snp, chr =  seqnames(gr_Unique_snp), x = start(gr_Unique_snp), col = "red", y=-0.15)

kp <- plotKaryotype(plot.type=3, main="SNPs for non-mosaic", genome = "hg38", chromosomes = "chrX")
kpPlotRainfall(karyoplot = kp , data=gr_bur_filterOG, r0=0.1, r1=0.3, col= "black")
kpPoints(karyoplot = kp , data=gr_bur_filterOG, chr =  seqnames(gr_bur_filterOG), x = start(gr_bur_filterOG), col = "black", y=-0.15)


 # try and plot the ACGT snp to have a show of the normal split of it. 
  
  
  
#Import the bed files from the COV files 
# first set a new working directory so the importing of the files gets easier.
  setwd("/media/god/jellyfish/Emil/R/sv/cov/")

  UPIC_cov <- fread("GTEX-UPIC_cov..bed")
  UPIC_cov_log <- UPIC_cov
  UPIC_cov_log$coverage <- log10(UPIC_cov$coverage +1)
  UPIC_COV_gr <- toGRanges(UPIC_cov)
  UPIC_cov_loggr <- toGRanges(UPIC_cov_log)
  
  PLJ_COV <-fread("GTEX-13PLJ_cov..bed")
  PLJ_COV_GR <- toGRanges(PLJ_COV)
  PLJ_cov_log <- PLJ_COV
  PLJ_cov_log$coverage <- log10(PLJ_COV$coverage)
  PLJ_cov_loggr <- toGRanges(PLJ_cov_log)
  
  ZZPU_COV <- fread("GTEX-ZZPU_cov..bed")
  ZZPU_COV_GR <- toGRanges(ZZPU_COV)
  ZZPU_cov_log <- ZZPU_COV
  ZZPU_cov_log$coverage <- log10(ZZPU_COV$coverage)
  ZZPU_cov_loggr <- toGRanges(ZZPU_cov_log)
  
  
  UPIC_NOZero <- filter(UPIC_cov, coverage > 0)
  UPIC_NOZero_log <- UPIC_NOZero
  UPIC_NOZero_log$coverage <- log10(UPIC_NOZero$coverage)
  UPIC_NOzero_gr <- toGRanges(UPIC_NOZero)
 
  
  UPIC_filt_500 <- subset(UPIC_NOZero, start < 55000001 | start > 64995001)
  UPIC_filt500_GR <- toGRanges(UPIC_filt_500)

  PLJ_NOZero <- filter(PLJ_COV, coverage > 0)
  PLJ_filt_500 <- subset(PLJ_NOZero, start < 55000001 | start > 64995001)
  PLJ_filt500_GR <- toGRanges(PLJ_filt_500)
  
  ZZPU_NOZero <- filter(ZZPU_COV, coverage > 0)
  ZZPU_filt_500 <- subset(ZZPU_NOZero, start < 55000001 | start > 64995001)
  ZZPU_filt500_GR <- toGRanges(ZZPU_filt_500)
  
  kp_cov <- plotKaryotype(plot.type = 4, main="plot of COV for Woman1_500", genome = "hg38", chromosomes =c("chrX"))
  kpPoints(karyoplot = kp_cov , data=UPIC_filt500_GR, chr =  seqnames(UPIC_filt500_GR), x = start(UPIC_filt500_GR), col = "black", y=UPIC_filt500_GR$coverage)
  
  kp_cov <- plotKaryotype(plot.type = 4, main="plot of COV for Woman2_500", genome = "hg38", chromosomes =c("chrX"))
  kpPoints(karyoplot = kp_cov , data=PLJ_filt500_GR, chr =  seqnames(PLJ_filt500_GR), x = start(PLJ_filt500_GR), col = "black", y=PLJ_filt500_GR$coverage)
  
  kp_cov <- plotKaryotype(plot.type = 4, main="plot of COV for Woman3_500", genome = "hg38", chromosomes =c("chrX"))
  kpPoints(karyoplot = kp_cov , data=ZZPU_filt500_GR, chr =  seqnames(ZZPU_filt500_GR), x = start(ZZPU_filt500_GR), col = "black", y=ZZPU_filt500_GR$coverage)
  #Goal is a plot to show the coverage of the whole genome and then plot it against the Control


kp_covUPIC <- plotKaryotype(plot.type = 4, main="plot of COV for Woman1", genome = "hg38", chromosomes =c("chrX"))
#kpPlotCoverage(kp_covUPIC, data=UPIC_cov_loggr, show.0.cov = TRUE, data.panel=1, r0=0, r1=0.2, col="red", border = "blue")
kpPoints(karyoplot = kp_covUPIC , data=UPIC_cov_loggr, chr =  seqnames(UPIC_cov_loggr), x = start(UPIC_cov_loggr), col = "black", y=UPIC_cov_loggr$coverage)

kp_covPLJ <- plotKaryotype(plot.type = 2, main="plot of COV for Woman2", genome = "hg38", chromosomes =c("chr17","chr16", "chrX"))
kpPlotCoverage(kp_covPLJ, data=PLJ_cov_loggr, show.0.cov = TRUE, data.panel=1, r0=0, r1=0.2, col="red", border="blue")

kp_covZZPU <- plotKaryotype(plot.type = 2, main="plot of COV for Woman3", genome = "hg38", chromosomes =c("chr17","chr16", "chrX"))
kpPlotCoverage(kp_covZZPU, data=ZZPU_cov_loggr, show.0.cov = TRUE, data.panel=1, r0=0, r1=0.2, col="red", border="blue")

#kp_Region <- plotKaryotype("hg38", plot.type=2, chromosomes = c("chrX"))
#kpPlotRegions(kp_Region, data = UPIC_cov_loggr,r0=0, r1=0.1, col="red")




#the SV files to be read in.
setwd("/media/god/jellyfish/Emil/R")

UPIC_cov5k <- subset(fread("GTEX-UPIC_cov_5K.bed"), `#chromosome` == "chrX")
UPIC_COV5k_gr <- toGRanges(UPIC_cov5k)

#UPIC_cov5k_log <- UPIC_cov5k
#UPIC_cov5k_log$coverage <- log10(UPIC_cov5k$coverage +1)
#UPIC_cov5k_loggr <- toGRanges(UPIC_cov5k_log)

PLJ_COV5k <-subset(fread("GTEX-13PLJ_cov_5K.bed"),  `#chromosome` == "chrX")
PLJ_COV5k_GR <- toGRanges(PLJ_COV5k)

#PLJ_cov5k_log <- PLJ_COV5k
#PLJ_cov5k_log$coverage <- log10(PLJ_COV5k$coverage +1)
#PLJ_cov5k_loggr <- toGRanges(PLJ_cov5k_log)

ZZPU_COV5k <- fread("GTEX-ZZPU_cov_5K.bed")
ZZPU_COV5k_GR <- toGRanges(ZZPU_COV5k)
ZZPU_cov5k_log <- ZZPU_COV5k
ZZPU_cov5k_log$coverage <- log10(ZZPU_COV5k$coverage +1)
ZZPU_cov5k_loggr <- toGRanges(ZZPU_cov_log)
#filtering out the centromere data. for all three samples. 
#Woman 1 
UPIC_filt <- UPIC_cov5k
UPIC_filt <- subset(UPIC_cov5k, start < 55000001 | start > 64995001)
UPIC_filt_GR <- toGRanges(UPIC_filt)

#Woman 2
PLJ_filter <- PLJ_COV5k
PLJ_filter <-  subset(PLJ_COV5k, start < 55000001 | start > 64995001)
PLJ_filter_GR <- toGRanges(PLJ_filter)

#woman 3
ZZPU_filt <- subset(ZZPU_COV5k, start < 55000001 | start > 64995001)
ZZPU_filter_GR <- toGRanges(ZZPU_filt)

#kprect plots
kp_cov <- plotKaryotype(plot.type = 4, main="plot of COV for Woman1_5k", genome = "hg38", chromosomes =c("chrX"))
kpRect(kp_cov, data=UPIC_filt_GR, chr =  seqnames(UPIC_filt_GR), y0=UPIC_filt_GR$coverage, y1=UPIC_filt_GR$coverage)

kp_cov <- plotKaryotype(plot.type = 4, main="plot of COV for Woman2_5k", genome = "hg38", chromosomes =c("chrX"))
kpRect(kp_cov, data=PLJ_filter_GR, chr =  seqnames(PLJ_filter_GR), y0=PLJ_filter_GR$coverage, y1=PLJ_filter_GR$coverage)

kp_cov <- plotKaryotype(plot.type = 4, main="plot of COV for Woman3_5k", genome = "hg38", chromosomes =c("chrX"))
kpRect(kp_cov, data=ZZPU_filter_GR, chr =  seqnames(ZZPU_filter_GR), y0=ZZPU_filter_GR$coverage, y1=ZZPU_filter_GR$coverage)

# KP points plots
kp_cov <- plotKaryotype(plot.type = 4, main="plot of COV for Woman3_5k", genome = "hg38", chromosomes =c("chrX"))
kpPoints(karyoplot = kp_cov , data=ZZPU_filter_GR, chr =  seqnames(ZZPU_filter_GR), x = start(ZZPU_filter_GR), col = "black", y=ZZPU_filter_GR$coverage)

kp_cov <- plotKaryotype(plot.type = 4, main="plot of COV for Woman2_5k", genome = "hg38", chromosomes =c("chrX"))
kpPoints(karyoplot = kp_cov , data=PLJ_filter_GR, chr =  seqnames(PLJ_filter_GR), x = start(PLJ_filter_GR), col = "black", y=PLJ_filter_GR$coverage)

kp_cov <- plotKaryotype(plot.type = 4, main="plot of COV for Woman1_5k", genome = "hg38", chromosomes=c("chrX"))
kpPoints(karyoplot = kp_cov, data =UPIC_filt_GR, chr = seqnames(UPIC_filt_GR), x = start(UPIC_filt_GR), col="black", y=UPIC_filt_GR$coverage )


UPIC_NOZero <- filter(UPIC_cov5k, coverage > 0.3)
UPIC_NOZero_log <- UPIC_NOZero
UPIC_NOZero_log$coverage <- log10(UPIC_NOZero$coverage)
UPIC_NOzero_gr <- toGRanges(UPIC_NOZero)

kp_cov <- plotKaryotype(plot.type = 7, main="     Woman1_5k", genome = "hg38", chromosomes =c("chrX"))
kpPoints(karyoplot = kp_cov , data=UPIC_NOzero_gr, chr =  seqnames(UPIC_NOzero_gr), x = start(UPIC_NOzero_gr), col = "blue", y=UPIC_NOzero_gr$coverage)


UPICVCFR <- read.vcfR( 'GTEX-UPIC..vcf', verbose= FALSE )
assess_missing_data_tsne(UPICVCFR, popmap=, thresholds = 3, perplexity = 5, iterations = 500, initial_dims = 5, clustering=TRUE)


