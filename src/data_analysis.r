#!/usr/bin/env Rscript

# required 
library(readxl)
library(ggplot2)
library(ggfortify) # required for autoplot
library(tidyr)
library(dplyr)
library(visdat) # required for vis_miss
library(tibble)
library(limma)
library(RColorBrewer)
library(pheatmap)
library(PMCMRplus) # posthoc test

source('src/helper_functions.R')

# load raw data from excel
prot_exp_raw <- read_excel("data_052020/Protein_Exp__FeatureLevel_artificial.xlsx", 
                       sheet = "Sheet 1",
                       col_names = TRUE,
                       na = "NA")

# split feature pattern and sum by protein, ungroup after operation
prot_exp <- prot_exp_raw %>%
  separate(Feature, c('Protein', 'Peptide', 'Charge'), sep = '_', remove = TRUE) %>%
  group_by(Protein) %>%
  summarise_at(vars(Pool1_126:Pool2_131N), sum, na.rm = TRUE) %>%
  ungroup()

# visualization of missing values
vis_miss(prot_exp[,2:20]) # no missing values after aggregation
# vis_miss(prot_exp_raw[,2:20]) # a batch pattern visible here, too

# box plot to explore data distribution
pr_expr_box_plot(prot_exp)

########### prep for PCA
# read annotation
sample_annotation <- read_excel("data_052020/ExpDesign.xlsx",
                                col_names = TRUE)
sample_annotation <- sample_annotation %>%
  unite(sample_name, Pool, Channel, remove = FALSE)


# transpose (Protein is the name of the first column) and add annotation
# tr_prot_exp_annotated <- prot_exp %>% 
#   gather(sample_name, values, 2:ncol(prot_exp)) %>%
#   spread(Protein, values) %>%
#   add_column(sample_annotation, .after = "sample_name")

# transpose (Protein is the name of the first column) and add annotation
tr_prot_exp_annotated <- prot_exp %>% 
  gather(sample_name, values, 2:ncol(prot_exp)) %>%
  spread(Protein, values) %>%
  right_join(sample_annotation, by = 'sample_name') %>%
  select(Sample, sample_name, Pool, Channel, BioRep, Condition, everything())

# compute and plot PCA
my_scaled_pca <- prcomp(tr_prot_exp_annotated[7:ncol(tr_prot_exp_annotated)], scale. = TRUE)
autoplot(my_scaled_pca, data=tr_prot_exp_annotated, colour = 'Pool', shape = 'Condition', scale = 0, size = 5)


# remove batch effect
prot_exp_matrix <- as.matrix(prot_exp[,2:20])
prot_exp_matrix_corrected <- removeBatchEffect(prot_exp_matrix, sample_annotation$Pool)

# copy back to tibble
prot_exp_cor <- prot_exp
prot_exp_cor[,2:20] <- prot_exp_matrix_corrected


# Control batch effect correction
pr_expr_box_plot(prot_exp_cor)

# add additional vector for replicates 1-3 and 4-6
RepBinary <- c("Rep1-3", "Rep1-3", "Rep1-3", "Rep1-3", "Rep4-6", "Rep1-3", "Rep1-3", "Rep4-6", "Rep4-6", "Rep4-6", "Rep4-6", "Rep1-3", "Rep1-3", "Rep1-3", "Rep1-3", "Rep4-6", "Rep1-3", "Rep1-3", "Rep4-6")

# transpose (Feature is the name of the first column) and add annotation
tr_prot_exp_cor_annotated <- prot_exp_cor %>% 
  gather(sample_name, values, 2:ncol(prot_exp_cor)) %>%
  spread(Protein, values) %>%
  right_join(sample_annotation, by = "sample_name") %>%
  select(Sample, sample_name, Pool, Channel, BioRep, Condition, everything()) %>%
  add_column(RepBinary, .after = "BioRep")


# compute and plot PCA
my_scaled_pca <- prcomp(tr_prot_exp_cor_annotated[8:ncol(tr_prot_exp_cor_annotated)], scale. = TRUE)
autoplot(my_scaled_pca, data=tr_prot_exp_cor_annotated, colour = 'Condition', shape = 'BioRep', scale = 0, size = 5)
# possible mislabeling of sample 17 (Rep3/Cond3)
autoplot(my_scaled_pca, data=tr_prot_exp_cor_annotated, colour = 'Condition', shape = 'RepBinary', scale = 0, size = 5)


# compute variance and keep top 2% only
prot_exp_cor_top2perc <- prot_exp_cor %>%
  rowwise() %>%
  mutate(variance = var(c_across(starts_with("Pool")))) %>%
  ungroup() %>%
  slice_max(variance, prop = 0.02)

# check mean / variance distribution to check variance increase in low value
prot_exp_cor_mnvar <- prot_exp_cor %>%
  rowwise() %>%
  mutate(variance = var(c_across(starts_with("Pool"))),
         mean = mean(c_across(starts_with("Pool"))),
         median = median(c_across(starts_with("Pool"))))

ggplot(prot_exp_cor_mnvar, aes(x=median, y= variance)) +
  geom_point() +
  scale_x_continuous(trans = 'log10') +
  scale_y_continuous(trans = 'log10') +
  geom_smooth(method=lm)


# heatmap plotting
prot_exp_matrix_cor_top2perc <- as.matrix(prot_exp_cor_top2perc[2:20])
rownames(prot_exp_matrix_cor_top2perc) <- prot_exp_cor_top2perc$Protein
colnames(prot_exp_matrix_cor_top2perc) <- paste(sample_annotation$Condition, '/', sample_annotation$BioRep, sep='')
# color_vec <- colorRampPalette(brewer.pal(9, "YlOrRd"))(255)
color_vec <- colorRampPalette(brewer.pal(9, "Blues"))(255)
# as expected, not very expressive
pheatmap(prot_exp_matrix_cor_top2perc,col=color_vec)

# two values are less than 0, replace with a values slightly lower than the next lowest
prot_exp_matrix_cor_top2perc[prot_exp_matrix_cor_top2perc<0] <- 6
pheatmap(log10(prot_exp_matrix_cor_top2perc),col=color_vec)

# not sure what this was for
ggplot(prot_exp_cor, aes(x=Pool1_126, y=Pool2_130N)) +
  geom_point() +
  scale_x_continuous(trans = 'log10') +
  scale_y_continuous(trans = 'log10')


## stats
prot_exp_cor_gathered <- prot_exp_cor %>%
  gather(sample_name, values, 2:ncol(prot_exp_cor)) %>%
  right_join(sample_annotation, by = 'sample_name') %>%
  select(Sample, Pool, Channel, BioRep, Condition, Protein, values) %>%
  mutate_at(vars(Sample, Pool, Channel, BioRep, Condition, Protein), as.factor)

# Control vs cond1
prot_exp_cor_gathered_VEH_cond1 <- prot_exp_cor_gathered %>%
  filter(Condition == 'VEH' | Condition == 'Cond1') %>%
  unite(CndPrt, Condition, Protein, remove = FALSE) %>%
  mutate_at(vars(CndPrt), as.factor)

kruskal_result_1 <- kruskal.test(values ~ CndPrt, data = prot_exp_cor_gathered_VEH_cond1)
dunn_result_1 <- kwAllPairsDunnTest(values ~ CndPrt, data = prot_exp_cor_gathered_VEH_cond1, p.adjust.methods = 'BH')

# Control vs cond2
prot_exp_cor_gathered_VEH_cond2 <- prot_exp_cor_gathered %>%
  filter(Condition == 'VEH' | Condition == 'Cond2') %>%
  unite(CndPrt, Condition, Protein, remove = FALSE) %>%
  mutate_at(vars(CndPrt), as.factor)

kruskal_result_2 <- kruskal.test(values ~ CndPrt, data = prot_exp_cor_gathered_VEH_cond2)
dunn_result_2 <- kwAllPairsDunnTest(values ~ CndPrt, data = prot_exp_cor_gathered_VEH_cond2, p.adjust.methods = 'BH')

# Control vs cond3
prot_exp_cor_gathered_VEH_cond3 <- prot_exp_cor_gathered %>%
  filter(Condition == 'VEH' | Condition == 'Cond3') %>%
  unite(CndPrt, Condition, Protein, remove = FALSE) %>%
  mutate_at(vars(CndPrt), as.factor)

kruskal_result_3 <- kruskal.test(values ~ CndPrt, data = prot_exp_cor_gathered_VEH_cond3)
dunn_result_3 <- kwAllPairsDunnTest(values ~ CndPrt, data = prot_exp_cor_gathered_VEH_cond3, p.adjust.methods = 'BH')
