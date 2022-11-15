install.packages("BiocManager")

BiocManager::install("karyoploteR")
BiocManager::install("rtracklayer")


library(chromomap)
library(cowplot)
library(ggplot2)
library(rstatix)
library(EnvStats)
library(ggbeeswarm)
library(ggpubr)
library(tidyverse)
library(karyoploteR)

temp = list.files(pattern="*.bed")
myfiles = lapply(temp, read.delim)

kp <- plotKaryotype(plot.type=6, main="plot.type=1")
kpDataBackground(kp, data.panel = 1)
kpText(kp, chr="chr1", x=60e6, y=0.5, labels="data.panel=1", data.panel = 1)
kpDataBackground(kp, data.panel = 2)
kpText(kp, chr="chr1", x=60e6, y=0.5, labels="data.panel=2", data.panel = 2)


regions <- createRandomRegions(nregions=400, length.mean = 3e6, mask=NA, non.overlapping = FALSE)
kpPlotRegions(kp, data=regions,chromosomes="chrX")



kp4 <- plotKaryotype("hg19", plot.type=1, chromosomes=c("chr17", "chrX"))

regs <- createRandomRegions(nregions = 1000, length.mean = 10000000, length.sd = 1000000,
                            non.overlapping = FALSE, genome = "hg19", mask=NA)
kpPlotRegions(kp4, regs, r0 = 0, r1 = 0.8, col="#AAAAAA")

kpPlotCoverage(kp4, regs, ymax = 20, r0=0.8,  r1=1, col="#CCCCFF")
kpAxis(kp, ymin = 0, ymax= 20, numticks = 2, r0 = 0.8, r1=1)


### T Test
stat_test_t.test <- df_raw %>% 
  t_test(formula = reformulate(colnames(df_raw)[2], (colnames(df_raw)[1])), paired = F, alternative = "two.sided")

### Wilcoxon
stat_test_Wilcoxon <- df_raw %>% 
  wilcox_test(formula = reformulate(colnames(df_raw)[2], (colnames(df_raw)[1])), paired = F, alternative = "two.sided")

### One-way ANOVA with post-hoc analysis
stat_test_anova_tukey <- df_raw %>% 
  aov(formula = reformulate(colnames(df_raw)[2], (colnames(df_raw)[1]))) %>% 
  tukey_hsd(formula = reformulate(colnames(df_raw)[2], (colnames(df_raw)[1])))

### Kruskal-Wallis with post-hoc analysis
stat_test_KS <- df_raw %>% 
  kruskal_test(formula = reformulate(colnames(df_raw)[2], (colnames(df_raw)[1])))
stat.test.dunn <- df_raw %>% dunn_test(reformulate(colnames(df_raw)[2], (colnames(df_raw)[1])))


base_plot <- ggplot(data = datasets, mapping = aes(x = expression))

base_plot + geom_boxplot()

base_plot + geom_point()

base_plot + geom_histogram()

kp <- plotKaryotype(chromosomes=(c("autosomal"))
                   

kp <- plotKaryotype(chromosomes=c("chrX")

 kp <- plotKaryotype(chromosomes = c("chr1", "chr2"), plot.type = 2)
 kpDataBackground(kp, data.panel = 1)
 kpDataBackground(kp, data.panel = 2)

