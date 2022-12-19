
install.packages("BiocManager")
install.packages("tidyverse")

#BiocManager::install("rtracklayer")
setwd("/path/to/dir/R/")
#BiocManager::install(version = "3.16")
BiocManager::install("AnnotationDbi")
library(BiocManager) #Version = 1.30.10
library(cowplot) #Version= 1.1.1
library(ggplot2) # Version = 3.3.3
library(data.table) # Version = 1.13.6
library(tidyr)# Version = 1.1.2
library(stringr)# Version =1.4.0
library(ggbeeswarm) # Version = 0.6.0
library(ggpubr) # Version = 0.4.0
library(karyoploteR)# Version =1.16.0
library(rtracklayer) # Version = 1.50.0
library(dplyr) # Version = 1.0.4


#imports the data of the three non-mosaic women

Woman1 <- fread("het_Woman1.bed")
Woman2 <- fread("het_Woman2.bed")
Woman3 <- fread("het_Woman3.bed")

#check which SNPs are in shared
Cmpare1 <- intersect(Woman1$V2, Woman2$V2)
Compared <- intersect(Compare1, Woman3$V2)

#Make a file with all of them together
Combined <- rbind(Woman1, Woman2, Woman3)

SNPs_inAll <- as.data.frame(Combined[Combined$V2 %in% Compared,])

colnames(SNPs_inAll) <- c("chr", "start", "ref", "alt", "depth", "sample")

SNPs_inAll$end <- SNPs_inAll$start + 1
SNPs_inAll$width <- 1

shortz <- dplyr::select(SNPs_inAll, chr, start, end, sample)

gr_SNPs_inAll <- toGRanges(shortz)

#Plots the SNPs of the three women. This is a control to check so the SNPs for them are the same or if eny sticks out to redo the part above.
kp <- plotKaryotype(plot.type=3, main="Blue= Woman1,RED ", genome = "hg38", chromosomes = "chrX")
kpPoints(karyoplot = kp , data=gr_SNPs_inAll[gr_SNPs_inAll$sample == "Woman1",], chr =  seqnames(gr_SNPs_inAll[gr_SNPs_inAll$sample == "Woman1",]), x = start(gr_SNPs_inAll[gr_SNPs_inAll$sample == "Woman1",]), col = "blue", y=0.05)
kpPoints(karyoplot = kp , data=gr_SNPs_inAll[gr_SNPs_inAll$sample == "Woman2",], chr =  seqnames(gr_SNPs_inAll[gr_SNPs_inAll$sample == "Woman2",]), x = start(gr_SNPs_inAll[gr_SNPs_inAll$sample == "Woman2",]), col = "red", y=0.15)
kpPoints(karyoplot = kp , data=gr_SNPs_inAll[gr_SNPs_inAll$sample == "Woman3",], chr =  seqnames(gr_SNPs_inAll[gr_SNPs_inAll$sample == "Woman3",]), x = start(gr_SNPs_inAll[gr_SNPs_inAll$sample == "Woman3",]), col = "grey", y=0.25)


#Import the Controls
Control <- fread("het_Control.bed")
#apply the control to solo out the one unique for the 3 women.
Unique <- setdiff(Compared, Control$V2)

#Creates a set with the SNPs that are not present in the control.
Unique_snp <- as.data.frame(>Combined[Combined$V2 %in% Unique,])

colnames(Unique_snp) <- c("chr", "start", "ref", "alt", "depth", "sample")

#Creates a end to the SNPs so it can be plotted.
Unique_snp$end <- Unique_snp$start + 1
Unique_snp$width <- 1

thingforgr <- dplyr::select(Unique_snp, chr, start, end, sample)

#turn the data table to a GRanges object.
gr_Unique_snp <- toGRanges(thingforgr)

# Redoing and splitting the depth so all zeros are removed.
bur_split <- Unique_snp %>% separate(depth, c('sampleX', 'refX'))
BurUNiqueNoZero <- filter(bur_split, sampleX > 0, refX > 0)

NoZerogr <- dplyr::select(BurUNiqueNoZero, chr, start, end, sample)
gr_NoZero_snp <- toGRanges(NoZerogr)

#PLot of all the different SNPs with the first being no zeros  

#PLots the Snps with nonZero values that are unique for the three Women 
kp <- plotKaryotype(plot.type=3, main="SNPs", genome = "hg38", chromosomes = "chrX")
kpPlotRainfall(kp, data=gr_NoZero_snp, r0=0.25, r1=0.5, col="blue")
kpPoints(karyoplot = kp , data=gr_NoZero_snp, chr =  seqnames(gr_NoZero_snp), x = start(gr_NoZero_snp), col = "blue", y=-0.15)

#Plots all the Unique SNPs from the three Women that have gone thru filtering away the SNPs that are present in the controls 
kp <- plotKaryotype(plot.type=3, main="Unique SNPs", genome = "hg38", chromosomes = "chrX")
kpPlotRainfall(kp, data=gr_Unique_snp , r0 = 0.1, r1=0.3, col = "red")
kpPoints(karyoplot = kp , data=gr_Unique_snp, chr =  seqnames(gr_Unique_snp), x = start(gr_Unique_snp), col = "red", y=-0.15)

#plot the SNPs from the three women, these are NOT filtered.
kp <- plotKaryotype(plot.type=3, main="SNPs for non-mosaic", genome = "hg38", chromosomes = "chrX")
kpPlotRainfall(karyoplot = kp , data=gr_SNPs_inAll, r0=0.1, r1=0.3, col= "black")
kpPoints(karyoplot = kp , data=gr_SNPs_inAll, chr =  seqnames(gr_SNPs_inAll), x = start(gr_SNPs_inAll), col = "black", y=-0.15)

# check if you will get any value from plotting the SNPs from just the XIC region and how that can further impact the other figures.


#Loads in data if it is wanted to check on the gene names. IMportant that you then create a subset where you only have the genes you are intrested in.
gennames <- loadDb("txdb.hg38.knownGene.sqlite")
genes.data <- makeGenesDataFromTxDb(gennames)
genes.data <- addGeneNames(genes.data)
genes.data <- mergeTranscripts(genes.data)


#Create the Zoom file for the XIC regions
zoom.region <- toGRanges(data.frame("chrX", 72e6, 74.3e6))


Just_XIC72k <-subset(NoZerogr, NoZerogr$start > 73792000)
Just_XICSNP74k <- subset(Just_XIC72k, Just_XIC72k$start < 74295000)
XICSNPS_GR <- toGRanges(Just_XICSNP74k)

# PLots a figure with the l,ocation on the selected region.
kp <- plotKaryotype(plot.type=3, main="SNPs for XIC-region 71mb to 74mb", genome = "hg38", chromosomes = "chrX", zoom = zoom.region)
kpPlotRainfall(karyoplot = kp , data=XICSNPS_GR, r0=0.1, r1=0.3, col= "black")
#Zoom in further to the SNPs that are present.
zoom.region <- toGRanges(data.frame("chrX", 74e6, 74.3e6))


kp <- plotKaryotype(plot.type=3, main="SNPs for XIC-region 74mb to 74.3mb", genome = "hg38", chromosomes = "chrX", zoom = zoom.region)
kpPlotRainfall(karyoplot = kp , data=XICSNPS_GR, r0=0.1, r1=0.3, col= "black")

#Xact region plot 
#Create the Zoom file for the Xact regions
zoom.XACT <- toGRanges(data.frame("chrX", 113e6, 114.5e6))
Just_XICSNP113k <- subset(NoZerogr, NoZerogr$start > 113112000)
Just_XICSNP114k <- subset(Just_XICSNP113k, Just_XICSNP113k$start < 114512000)
XACTSNPs <- toGRanges(Just_XICSNP114k)

#plots the SNps over the Xact gene 
kp <- plotKaryotype(plot.type=3, main="SNPs for XACT 113mb to 114mb", genome = "hg38", chromosomes = "chrX", zoom = zoom.XACT)
kpPlotRainfall(karyoplot = kp , data=XACTSNPs, r0=0.1, r1=0.3, col= "black")

#kpPlotGenes(karyoplot = kp , data= gennames, add.transcript.names = TRUE, plot.transcripts = T, gene.names = NULL, r1=0.6, data.panel = 2)


#Import the bed files from the COV files 
# first set a new working directory so the importing of the files gets easier.
  setwd("/path/to/directory/R/")

  Woman1_cov <- fread("Woman1_cov..bed")
  Woman1_COV_gr <- toGRanges(Woman1_cov)

#If to high Cov values use the log outputs instead
  Woman1_cov_log <- Woman1_cov
  Woman1_cov_log$coverage <- log10(Woman1_cov$coverage +1)
  Woman1_cov_loggr <- toGRanges(Woman1_cov_log)
  
  Woman2_COV <-fread("Woman2_cov..bed")
  Woman2_COV_GR <- toGRanges(Woman2_COV)
  
#If to high Cov values use the log outputs instead
Woman2_cov_log <- Woman2_COV
  Woman2_cov_log$coverage <- log10(Woman2_COV$coverage)
  Woman2_cov_loggr <- toGRanges(Woman2_cov_log)
  
  Woman3_COV <- fread("Woman3_cov..bed")
  Woman3_COV_GR <- toGRanges(Woman3_COV)


#If to high Cov values use the log outputs instead
Woman3_cov_log <- Woman3_COV
 Woman3_cov_log$coverage <- log10(Woman3_COV$coverage)
  Woman3_cov_loggr <- toGRanges(Woman3_cov_log)
  
  #removes all teh tiles with 0 coverage
  Woman1_NOZero <- filter(Woman1_cov, coverage > 0)
  Woman1_NOZero_log <- Woman1_NOZero
 #If the coverage values are high then apply the 10log over the coverage column to get more managable values
   Woman1_NOZero_log$coverage <- log10(Woman1_NOZero$coverage)
  
   Woman1_NOzero_gr <- toGRanges(Woman1_NOZero)
 
  # Remove the reads over the Centromere. 
  Woman1_NOZero <- filter(Woman1_cov, coverage > 0)
  Woman1_filt_500 <- subset(Woman1_NOZero, start < 55000001 | start > 64995001)
  Woman1_filt500_GR <- toGRanges(Woman1_filt_500)

  Woman2_NOZero <- filter(Woman2_COV, coverage > 0)
  Woman2_filt_500 <- subset(Woman2_NOZero, start < 55000001 | start > 64995001)
  Woman2_filt500_GR <- toGRanges(Woman2_filt_500)
  
  Woman3_NOZero <- filter(Woman3_COV, coverage > 0)
  Woman3_filt_500 <- subset(ZZPU_NOZero, start < 55000001 | start > 64995001)
 Woman3_filt500_GR <- toGRanges(Woman3_filt_500)
  
  #PLots the coverage of all the non-mosaic women.
  kp_cov <- plotKaryotype(plot.type = 4, main="plot of COV for Woman1_500", genome = "hg38", chromosomes =c("chrX"))
  kpPoints(karyoplot = kp_cov , data=Woman1_filt500_GR, chr =  seqnames(Woman1_filt500_GR), x = start(Woman1_filt500_GR), col = "black", y=Woman1_filt500_GR$coverage)
  
  kp_cov <- plotKaryotype(plot.type = 4, main="plot of COV for Woman2_500", genome = "hg38", chromosomes =c("chrX"))
  kpPoints(karyoplot = kp_cov , data=Woman2_filt500_GR, chr =  seqnames(Woman2_filt500_GR), x = start(Woman2_filt500_GR), col = "black", y=Woman2_filt500_GR$coverage)
  
  kp_cov <- plotKaryotype(plot.type = 4, main="plot of COV for Woman3_500", genome = "hg38", chromosomes =c("chrX"))
  kpPoints(karyoplot = kp_cov , data=Woman3_filt500_GR, chr =  seqnames(Woman3_filt500_GR), x = start(Woman3_filt500_GR), col = "black", y=Woman3_filt500_GR$coverage)
  
  
#here plot the specific Cov and Sv from the XIC region from the 500 tile files.
kp_cov <- plotKaryotype(plot.type = 4, main="plot of COV for Woman1_500", genome = "hg38", chromosomes =c("chrX"), zoom=zoom.region)
kpPoints(karyoplot = kp_cov , data=Woman1_filt500_GR, chr =  seqnames(Woman1_filt500_GR), x = start(Woman1_filt500_GR), col = "black", y=Woman1_filt500_GR$coverage)

kp_cov <- plotKaryotype(plot.type = 4, main="plot of COV for Woman2_500", genome = "hg38", chromosomes =c("chrX"))
kpPoints(karyoplot = kp_cov , data=Woman2_filt500_GR, chr =  seqnames(Woman2_filt500_GR), x = start(Woman2_filt500_GR), col = "black", y=Woman2_filt500_GR$coverage)

kp_cov <- plotKaryotype(plot.type = 4, main="plot of COV for Woman3_500", genome = "hg38", chromosomes =c("chrX"))
kpPoints(karyoplot = kp_cov , data=Woman3_filt500_GR, chr =  seqnames(Woman3_filt500_GR), x = start(Woman3_filt500_GR), col = "black", y=Woman3_filt500_GR$coverage)


#PLots the 10log of the three Women if needed due to high values from before. 
kp_covWoman1 <- plotKaryotype(plot.type = 4, main="plot of COV for Woman1", genome = "hg38", chromosomes =c("chrX"))
#kpPlotCoverage(kp_covUPIC, data=Woman1_cov_loggr, show.0.cov = TRUE, data.panel=1, r0=0, r1=0.2, col="red", border = "blue")
kpPoints(karyoplot = kp_covWoman1 , data=Woman1_cov_loggr, chr =  seqnames(Woman1_cov_loggr), x = start(Woman1_cov_loggr), col = "black", y=Woman1_cov_loggr$coverage)

kp_covWoman2 <- plotKaryotype(plot.type = 2, main="plot of COV for Woman2", genome = "hg38", chromosomes =c("chr17","chr16", "chrX"))
kpPlotCoverage(kp_covWoman2, data=Woman2_cov_loggr, show.0.cov = TRUE, data.panel=1, r0=0, r1=0.2, col="red", border="blue")

kp_covWoman3 <- plotKaryotype(plot.type = 2, main="plot of COV for Woman3", genome = "hg38", chromosomes =c("chr17","chr16", "chrX"))
kpPlotCoverage(kp_covWoman3, data=Woman3_cov_loggr, show.0.cov = TRUE, data.panel=1, r0=0, r1=0.2, col="red", border="blue")






#the COV files with 5k tiles to be read in.
setwd("/path/to/directory/R")
#This reads in only the chrX from the file.
#Woman1 
Woman1_cov5k <- subset(fread("Woman1_cov_5K.bed"), `#chromosome` == "chrX")
Woman1_COV5k_gr <- toGRanges(Woman1_cov5k)
#filtering out the centromere data.
Woman1_filt <- Woman1_cov5k
Woman1_filt <- subset(Woman1_cov5k, start < 55000001 | start > 64995001)
Woman1_filt_GR <- toGRanges(Woman1_filt)

#Woman2
Woman2_COV5k <-subset(fread("Woman2_cov_5K.bed"),  `#chromosome` == "chrX")
Woman2_COV5k_GR <- toGRanges(Woman2_COV5k)
#filtering out the centromere data.
Woman2_filter <- Woman2_COV5k
Woman2_filter <-  subset(Woman2_COV5k, start < 55000001 | start > 64995001)
Woman2_filter_GR <- toGRanges(Woman2_filter)
#Woman3
Woman3_COV5k <- subset(fread("Woman3_cov_5K.bed"), `#chromosome` == "chrX")
Woman3_COV5k_GR <- toGRanges(Woman3_COV5k)
#filtering out the centromere data.
Woman3_filter <- Woman3_COV5k
Woman3_filt <- subset(Woman3_COV5k, start < 55000001 | start > 64995001)
Woman3_filter_GR <- toGRanges(Woman3_filt)









#kpRect plots for the 5k tiles so it can be obeserved in the region and on other chromosomes if needed.
kp_cov <- plotKaryotype(plot.type = 4, main="plot of COV for Woman1_5k", genome = "hg38", chromosomes =c("chrX"), zoom=zoom.region)
kpRect(kp_cov, data=Woman1_filt_GR, chr =  seqnames(Woman1_filt_GR), y0=Woman1_filt_GR$coverage, y1=Woman1_filt_GR$coverage)

kp_cov <- plotKaryotype(plot.type = 4, main="plot of COV for Woman2_5k", genome = "hg38", chromosomes =c("chrX"), zoom = zoom.region)
kpRect(kp_cov, data=Woman2_filter_GR, chr =  seqnames(Woman2_filter_GR), y0=Woman2_filter_GR$coverage, y1=Woman2_filter_GR$coverage)

kp_cov <- plotKaryotype(plot.type = 4, main="plot of COV for Woman3_5k", genome = "hg38", chromosomes =c("chrX"), zoom = zoom.region)
kpRect(kp_cov, data=Woman3_filter_GR, chr =  seqnames(Woman3_filter_GR), y0=Woman3_filter_GR$coverage, y1=Woman3_filter_GR$coverage)

# KP points plots
kp_cov <- plotKaryotype(plot.type = 4, main="plot of COV for Woman3_5k", genome = "hg38", chromosomes =c("chrX"))
kpPoints(karyoplot = kp_cov , data=Woman3_filter_GR, chr =  seqnames(Woman3_filter_GR), x = start(Woman3_filter_GR), col = "black", y=Woman3_filter_GR$coverage)

kp_cov <- plotKaryotype(plot.type = 4, main="plot of COV for Woman2_5k", genome = "hg38", chromosomes =c("chrX"))
kpPoints(karyoplot = kp_cov , data=Woman2_filter_GR, chr =  seqnames(Woman2_filter_GR), x = start(Woman2_filter_GR), col = "black", y=Woman2_filter_GR$coverage)

kp_cov <- plotKaryotype(plot.type = 4, main="plot of COV for Woman1_5k", genome = "hg38", chromosomes=c("chrX"))
kpPoints(karyoplot = kp_cov, data =Woman1_filt_GR, chr = seqnames(Woman1_filt_GR), x = start(Woman1_filt_GR), col="black", y=Woman1_filt_GR$coverage )

#changing the filter for the coverage to a slightly higher value. 
Woman1_NOZero <- filter(Woman1_cov5k, coverage > 0.3)
Woman1_NOZero_log <- Woman1_NOZero
Woman1_NOZero_log$coverage <- log10(Woman1_NOZero$coverage)
Woman1_NOzero_gr <- toGRanges(Woman1_NOZero)

kp_cov <- plotKaryotype(plot.type = 7, main="     Woman1_5k", genome = "hg38", chromosomes =c("chrX"))
kpPoints(karyoplot = kp_cov , data=Woman1_NOzero_gr, chr =  seqnames(Woman1_NOzero_gr), x = start(Woman1_NOzero_gr), col = "blue", y=Woman1_NOzero_gr$coverage)


Woman1_ploidies <- read.table("/path/to/directory/R/Woman1..ploidies.tab", header = T)
Woman2_ploidies <- read.table("/path/to/directory/R/Woman2..ploidies.tab", header = T)
Woman3_ploidies <- read.table("/path/to/directory/R/Woman3..ploidies.tab", header = T)

Woman1_ploidies <- Woman1_ploidies[1:23,]
Woman2_ploidies <- Woman2_ploidies[1:23,]
Woman3_ploidies <- Woman3_ploidies[1:23,]

k="Woman1"
p="Woman2"
m="Woman3"

Woman1_ploidies$sample <- k
Woman2_ploidies$sample <- p
Woman3_ploidies$sample <- m
Women <- rbind(Woman1_ploidies,Woman2_ploidies,Woman3_ploidies)


standard_contigs <- c(paste0(rep("chr", 22), seq(1,22,1)), "chrX")

C1 <- read.table("/path/to/directory/R/C1..ploidies.tab", header = T)
C2 <- read.table("/path/to/directory/R/C2..ploidies.tab", header = T)
C3 <- read.table("/path/to/directory/R/C3..ploidies.tab", header = T)
C4 <- read.table("/path/to/directory/R/C4..ploidies.tab", header = T)
C5 <- read.table("/path/to/directory/R/C5..ploidies.tab", header = T)
C6 <- read.table("/path/to/directory/R/C6..ploidies.tab", header = T)
C7 <- read.table("/path/to/directory/R/C7..ploidies.tab", header = T)
C8 <- read.table("/path/to/directory/R/C8..ploidies.tab", header = T)
C9 <- read.table("/path/to/directory/R/C9..ploidies.tab", header = T)
C10 <- read.table("/path/to/directory/R/C10..ploidies.tab", header = T)
C11 <- read.table("/path/to/directory/R/C11..ploidies.tab", header = T)
C12 <- read.table("/path/to/directory/R/C12..ploidies.tab", header = T)

C1 <-C1[1:23,]
C2 <-C2[1:23,]
C3 <-C3[1:23,]
C4 <-C4[1:23,]
C5 <-C5[1:23,]
C6 <-C6[1:23,]
C7 <-C7[1:23,]
C8 <-C8[1:23,]
C9 <-C9[1:23,]
C10 <-C10[1:23,]
C11 <-C11[1:23,]
C12 <-C12[1:23,]
C="Control"
Control_Cov <- rbind(C1,C2,C3,C4,C5,C6,C7,C8,C9,C10,C11,C12)
Control_Cov$sample <- C

Women <- rbind(Women,Control_Cov)
source("//path/to/directory/R/AL_ggplot_themes.R")
library(rstatix)

women_stats <- Women[Women$Chromosome %in% standard_contigs,] %>% dplyr::group_by(Chromosome) %>% rstatix::get_summary_stats(Ploidy, type = "common")

Control_stats <- Control_Cov[Control_Cov$Chromosome %in% standard_contigs,] %>% dplyr::group_by(Chromosome) %>% rstatix::get_summary_stats(Ploidy, type = "common")

ggplot(Women[Women$Chromosome %in% standard_contigs,], aes(x= factor(Chromosome, levels = standard_contigs),  y= Ploidy ))+
  stat_summary(geom="line", fun="mean", aes(group=1))+
  geom_pointrange(data = women_stats, aes(y=mean, ymin=mean-sd, ymax=mean+sd, group=1), lty=2)+
  #geom_line(data = women_stats, aes(y=mean+sd, group=1), lty=2)+
  geom_point(aes(color=sample))+
  theme_AL_box_rotX()+
  theme(
    legend.position = "right", plot.title = element_text(size=12))+
  ggtitle("PLoidy plot")+
    labs(x="")+
  geom_hline(yintercept = 2, lty=2)



Women[Women$Chromosome %in% standard_contigs,] %>% dplyr::group_by(Chromosome) %>% kruskal_test(Ploidy ~ Chromosome)

ggplot(Women[Women$Chromosome %in% standard_contigs,], aes(x=factor(Chromosome, levels = standard_contigs), y=Ploidy)) + geom_boxplot() + ggbeeswarm::geom_quasirandom(aes(col= sample))+
  #geom_line(data = Control_stats, aes(y=mean+sd, group=1), lty=2)+
  geom_point(data = Control_stats,aes(x=Chromosome, y=median))+
  theme_AL_box_rotX()
