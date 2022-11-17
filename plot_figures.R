#Install all the packages that are needed.
install.packages("BiocManager")
install.packages("tidyverse")
#BiocManager::install("rtracklayer")

#Set the working directory to the correct directory.
setwd("/media/god/jellyfish/Emil/R/")

#Check the version of biocManager and set it to the latest. and set all the packages that are needed to the library.
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

#an example of a creating a list with all the vcf files and and import function for them.
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

#Filter out so only the SNPs that are in all samples are in the files. 
bur_filterOG <- as.data.frame(bur[bur$V2 %in% kek2,])

#Give the column names for the further processing in the change to a GRanges object. 
colnames(bur_filterOG) <- c("chr", "start", "ref", "alt", "depth", "sample")

#Sets a new column and then assigns it the value of start of snp +1.
bur_filterOG$end <- bur_filterOG$start + 1
#Sets the width to 1.
bur_filterOG$width <- 1
#sets and renames the selected and makes it easier for thetransformation to a GRanges object.
shortz <- dplyr::select(bur_filterOG, chr, start, end, sample)

#Creates the GRanges for the data of the SNPS that exists in all of the women.
gr_bur_filterOG <- toGRanges(shortz)

#plots the data and separetes the samples to 3 different lines, acts as a form of control 
kp <- plotKaryotype(plot.type=3, main="plot.type=3", genome = "hg38", chromosomes = "chrX")
kpPoints(karyoplot = kp , data=gr_bur_filterOG[gr_bur_filterOG$sample == "GTEX-UPIC-0004-SM-5SOEF",], chr =  seqnames(gr_bur_filterOG[gr_bur_filterOG$sample == "GTEX-UPIC-0004-SM-5SOEF",]), x = start(gr_bur_filterOG[gr_bur_filterOG$sample == "GTEX-UPIC-0004-SM-5SOEF",]), col = "blue", y=0.05)
kpPoints(karyoplot = kp , data=gr_bur_filterOG[gr_bur_filterOG$sample == "GTEX-13PLJ-0003-SM-6WSSCN",], chr =  seqnames(gr_bur_filterOG[gr_bur_filterOG$sample == "GTEX-13PLJ-0003-SM-6WSSCN",]), x = start(gr_bur_filterOG[gr_bur_filterOG$sample == "GTEX-13PLJ-0003-SM-6WSSCN",]), col = "red", y=0.15)
kpPoints(karyoplot = kp , data=gr_bur_filterOG[gr_bur_filterOG$sample == "GTEX-ZZPU-0003-SM-6WBUC",], chr =  seqnames(gr_bur_filterOG[gr_bur_filterOG$sample == "GTEX-ZZPU-0003-SM-6WBUC",]), x = start(gr_bur_filterOG[gr_bur_filterOG$sample == "GTEX-ZZPU-0003-SM-6WBUC",]), col = "grey", y=0.25)
kpAxis()


#Import the Controls
Control <- fread("het_Control.bed")
#apply the control to solo out the one unique for the 3 women.
kek3 <- setdiff(kek2, Control$V2)

#Isolates the samples data with the SNPS that only exists in the 3 women by checking wich SNPs that are in the control and subtracting thoose.
Unique_snp <- as.data.frame(bur[bur$V2 %in% kek3,])

#Give the column names for the further processing in the change to a GRanges object. 
colnames(Unique_snp) <- c("chr", "start", "ref", "alt", "depth", "sample")

#Sets a new column and then assigns it the value of start of snp +1.
Unique_snp$end <- Unique_snp$start + 1
#Sets the width to 1.
Unique_snp$width <- 1

#sets and renames the selected and makes it easier for thetransformation to a GRanges object.
thingforgr <- dplyr::select(Unique_snp, chr, start, end, sample)

#Creates the GRanges for the data of the SNPS that are unique to the non-mosaic women.
gr_Unique_snp <- toGRanges(thingforgr)

#buren <- data.frame(bur[bur$V2 %in% kek3,])
#mcols(buren) <- cbind(mcols(buren), newMcols)

#plots the data and separetes the snps from the 3 women that is not subtracted with the control SNPs. And the second line is SNPs with the COntrol SNPs subtracted.
kp <- plotKaryotype(plot.type=3, main="plot.type=3", genome = "hg38", chromosomes = "chrX")
  kpPoints(karyoplot = kp , data=gr_bur_filterOG, chr =  seqnames(gr_bur_filterOG), x = start(gr_bur_filterOG), col = "blue", y=0.05)
  kpPoints(karyoplot = kp , data=gr_Unique_snp, chr =  seqnames(gr_Unique_snp), x = start(gr_Unique_snp), col = "red", y=0.15)
  kpAxis()
  kpPlotMarkers(karyoplot = kp, data=gr_Unique_snp, chr = seqnames(gr_Unique_snp), y=0.3, labels=)



#change the working directory
setwd("/media/god/jellyfish/Emil/R/SV/cov/")

cov_list = list.files(pattern="*.bed")
cov_files = lapply(temp, read.delim)

#Fix the samples to Granges objects. ask Björn about this if i am supposed to do the same as above or change it.
UPIC_COV <- as.data.frame(cov_files$sample == "GTEX-UPIC-0004-SM-5SOEF")
PLJ_COV <- as.data.frame(cov_files$sample == "GTEX-13PLJ-0003-SM-6WSSCN")
ZZPU_COV <- as.data.frame(cov_files$sample == "GTEX-ZZPU-0003-SM-6WBUC")

# Fix the combined file for the 11 controls.



#plot the coverage of the non-mosaic and compare to the control 



