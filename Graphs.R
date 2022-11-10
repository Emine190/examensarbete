library(cowplot)
library(ggplot2)
library(rstatix)
library(EnvStats)
library(ggbeeswarm)
library(ggpubr)
library(tidyverse)


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
