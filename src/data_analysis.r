# data analysis
library(readxl)
library(ggplot2)
library(ggfortify) # required for autoplot
library(tidyr)
library(dplyr)
library(visdat)
library(tibble)
prot_exp_raw <- read_excel("data_052020/Protein_Exp__FeatureLevel_artificial.xlsx", 
                       sheet = "Sheet 1",
                       col_names = TRUE,
                       na = "NA")

prot_exp <- prot_exp_raw %>%
  separate(Feature, c('Protein', 'Peptide', 'Charge'), sep = '_', remove = TRUE) %>%
    group_by(Protein) %>%
    summarise_at(vars(Pool1_126:Pool2_131N), sum, na.rm = TRUE)

# heatmap for visualization of missing values
vis_miss(prot_exp[,2:20])

# box plot
gthrd_prot_exp <- prot_exp %>% 
  gather(sample_name, values, 2:ncol(prot_exp))

ggplot(gthrd_prot_exp, aes(x=sample_name, y=values), la) + 
  geom_boxplot(notch = TRUE) +
  labs( x = "Sample", y = "Protein expression (Sum of features)",
      title ="Cumulative protein expression per sample")

#tr_prot_exp <- as_tibble(cbind(nms = names(prot_exp), t(prot_exp)))
# testData = prot_exp[1:20, 1:5]
# testData %>% 
#   gather(var, values, 2:ncol(testData)) %>%
#   spread(Feature, values)

# add annotation
sample_annotation <- read_excel("data_052020/ExpDesign.xlsx",
                                col_names = TRUE)

# transpose (Feature is the name of the first column) and add annotation
tr_prot_exp_annotated <- prot_exp %>% 
  gather(sample_name, values, 2:ncol(prot_exp)) %>%
  spread(Protein, values) %>%
  add_column(sample_annotation, .after = "sample_name")


my_scaled_pca <- prcomp(tr_prot_exp_annotated[7:ncol(tr_prot_exp_annotated)], scale. = TRUE)
autoplot(my_scaled_pca, data=tr_prot_exp_annotated, colour = 'Pool', shape = 'Condition', scale = 0)
# autoplot(my_scaled_pca, data=tr_prot_exp_annotated, colour = 'Condition', scale = 0)
# autoplot(my_scaled_pca, data=tr_prot_exp_annotated, colour = 'BioRep', scale = 0)
# autoplot(my_scaled_pca, data=tr_prot_exp_annotated, colour = 'Channel', scale = 0)


# NA problem solved after summarizing
tr_prot_exp_no_NA <- tr_prot_exp_annotated %>%
  gather(protein, expression_values, 7:ncol(tr_prot_exp_annotated)) %>%
  mutate(expression_values = replace_na(expression_values, 0)) %>%
  spread(protein, expression_values)

tr_prot_exp_track_NA <- tr_prot_exp_annotated %>%
  gather(protein, expression_values, 7:ncol(tr_prot_exp_annotated)) %>%
  mutate(expression_values, was_NA = is.na(expression_values), .before = expression_values) %>%
  select(sample_name:was_NA) %>%
  spread(protein, wamy_scaled_pca_2s_NA)


my_scaled_pca_2 <- prcomp(tr_prot_exp_no_NA[7:ncol(tr_prot_exp_no_NA)], scale. = T)
autoplot(my_scaled_pca_2, data=tr_prot_exp_no_NA, colour = 'Pool', pch = '4', scale = 0)
autoplot(my_scaled_pca_2, data=tr_prot_exp_no_NA, colour = 'Condition', scale = 0)
autoplot(my_scaled_pca_2, data=tr_prot_exp_no_NA, colour = 'BioRep', scale = 0)
autoplot(my_scaled_pca_2, data=tr_prot_exp_no_NA, colour = 'Channel', scale = 0)

# depends on vqv-ggbiplot
# ggbiplot(my_scaled_pca_2, ellipse=TRUE,  labels=tr_prot_exp_no_NA$Sample, groups=tr_prot_exp_no_NA$Pool)
