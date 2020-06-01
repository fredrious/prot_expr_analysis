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