#!/usr/bin/env Rscript

pr_expr_box_plot <- function(my_data) {
  gthrd_data <- my_data %>% 
    gather(sample_name, values, 2:ncol(my_data))
  
  # zoom to visualize mild batch effect
  ggplot(gthrd_data, aes(x=sample_name, y=values)) + 
    geom_boxplot(notch = FALSE) +
    coord_cartesian(ylim = c(5,40)) +
    labs( x = "Sample", y = "Protein expression (Sum of features)",
          title ="Cumulative protein expression per sample")
  
}

verify_NA_effect <- function(my_data){
  ###### NA handling
  # NA problem solved by summarizing by proteins
  # leaving here for educational purpose
  
  # replace NA with zero
  data_NA_to_zero <- my_data %>%
    gather(protein, expression_values, 7:ncol(my_data)) %>%
    mutate(expression_values = replace_na(expression_values, 0)) %>%
    spread(protein, expression_values)
  
  # track NA positions
  tr_prot_exp_track_NA <- my_data %>%
    gather(protein, expression_values, 7:ncol(my_data)) %>%
    mutate(expression_values, was_NA = is.na(expression_values), .before = expression_values) %>%
    select(sample_name:was_NA) %>%
    spread(protein, was_NA)
  
  # PCA and plot
  pcaObj <- prcomp(data_NA_to_zero[7:ncol(data_NA_to_zero)], scale. = T)
  autoplot(pcaObj, data=data_NA_to_zero, colour = 'Pool', shape = 'Condition', scale = 0, size = 5)
}