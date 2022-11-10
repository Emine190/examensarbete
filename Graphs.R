install.packages("chromoMap")

library(chromomap)
library(cowplot)
library(ggplot2)
library(rstatix)
library(EnvStats)
library(ggbeeswarm)
library(ggpubr)
library(tidyverse)
library(karyoploteR)



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

