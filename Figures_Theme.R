#!/usr/bin/env Rscript
#
# Custom ggplot2 theme

library(ggplot2)

# plotting parameters
text_col <- "black"
file_theme.size.minor <- 7
file_theme.size.major <- 8

# ggplot theme for plots intended for DIN A4 sized pages / Word docs
file_theme <- theme_minimal() +
  theme(legend.text = element_text(size = file_theme.size.minor, color = text_col),
        legend.title = element_text(size = file_theme.size.minor, color = text_col),
        legend.position = "bottom",
        legend.margin=margin(t=-5),
        axis.text = element_text(size = file_theme.size.minor, color = text_col),
        axis.title = element_text(size = file_theme.size.major, color = text_col),
        plot.title = element_text(size = file_theme.size.major, hjust = 0.5, 
                                  color = text_col, face = "bold"),
        strip.text = element_text(size = file_theme.size.major, 
                                  color = text_col, face = "bold",
                                  margin = margin( b = 1.5, t = 1.5,
                                                   l = 1.5, r = 1.5)))
