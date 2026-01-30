### Figure 1 heritability and PGS prediction ###

library(tidyverse)
library(data.table)
library(patchwork)

source("Figures_Theme.R")

#### Figure 1: Heritability ####

fig1_data <- fread("data/heritability.tsv", data.table = FALSE) %>% 
  pivot_longer(cols = c(h2_twin, h2_SNP, PGS_R2), names_to = "Estimate") %>% 
  mutate(Estimate = factor(Estimate, levels = c("h2_twin", "h2_SNP", "PGS_R2"),
                           labels = c("Twin Heritability",
                                      "SNP Heritability",
                                      "PGS Prediction")))

fig1 <- fig1_data %>% 
    mutate(Disorder = fct_reorder(Disorder, value, .fun = max, .desc = TRUE)) %>% 
    ggplot(aes(x = Disorder,  y = value, fill = Estimate)) +
    geom_bar(stat = "identity", width = 0.9, position = position_dodge()) +
    labs(x = NULL, y = "% Variance Explained") +
    scale_fill_viridis_d(option = "D", end = 0.85) +
    ylim(0, 1) +
    file_theme +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          legend.position = c(0.8, 0.85),
          legend.key.size = unit(0.25, "cm"))

ggsave("Figure_1.png", path = "output", plot = fig1, bg = "white",
       width = 8.9, height = 10, units = "cm")
ggsave("Figure_1.pdf", path = "output", plot = fig1, bg = "white",
       width = 8.9, height = 10, units = "cm")


#### Figure 2: Relative and absolute risk ####

# calculations under liability threshold model for given prevalence k, 
# heritability h2, and genetic correlation r
ltm <- function(k, h2, r) {
  tidyr::crossing(prevalence = k, h2 = h2, gen_cor = r) %>% rowwise() %>% 
    mutate(
      # liability threshold for a binary disorder (disease when L>T)
      threshold = qnorm(1 - prevalence),
      
      # mean liability of affected individuals
      phi_T = dnorm(threshold),
      mean_affected = phi_T / prevalence,
      
      # mean liability shift in a first-degree relative (FDR)
      mu_fdr = gen_cor * h2 * mean_affected,
      
      # absolute risk for a first-degree relative
      K_fdr = 1 - pnorm(threshold - mu_fdr),
      
      # relative risk for a first-degree relative compared to general population
      RR = K_fdr / prevalence
    ) %>% 
    ungroup()
}

# illustration for given prevalence, heritability and genetic correlation
plot_illustration <- function(k, h2, r) {
  params <- ltm(k, h2, r)
  
  ggplot(data.frame(x = c(-4, 4)), aes(x = x)) + 
    stat_function(fun = dnorm, args = list(mean = params$mu_fdr, sd = 1),
                  color = "#440154FF") +
    stat_function(fun = dnorm, args = list(mean = params$mu_fdr, sd = 1),
                  xlim = c(params$threshold, 4), geom = "area", 
                  fill = "#440154FF", color = "#440154FF", alpha = 1) +
    stat_function(fun = dnorm, args = list(mean = 0, sd = 1), color = "#84CA72") +
    stat_function(fun = dnorm, args = list(mean = 0, sd = 1),
                  xlim = c(params$threshold, 4), geom = "area", 
                  fill = "#84CA72", color = "#84CA72", alpha = 1) +
    labs(title = paste0("h2: ", h2), x = "Liability", y = "Density") +
    file_theme 
}

# set of heritabilites, prevalences and genetic correlations
h2 <- c(0.2, 0.4, 0.8)
prevalence <- seq(0, 0.2, 0.02)
gen_cor <- c(0, 0.25, 0.5)
prev_highlight <-  c(0.02, 0.2)

# calculations for different combinations of prevalence, heritability, relatedness
fig2_data <- ltm(prevalence, h2, gen_cor) %>% 
  mutate(highlight = prevalence %in% prev_highlight,
         Relatedness = case_when(gen_cor == 0 ~ "unrelated",
                                 gen_cor == 0.25 ~ "2nd degree",
                                 gen_cor == 0.5 ~ "1st degree",
                                 TRUE ~ NA_character_) %>% factor()) %>% na.omit()

# panel a: illustration of liability threshold model, with prevalence = 2% (~BD)
illustration <- lapply(h2, plot_illustration, k = 0.02, r = 0.5) %>% 
    wrap_plots(nrow = 1)

# panel b: risk of individuals with first- or second-degree affected relative in
# relation to the general population (across prevalence range 2% - 20%)
relrisk <- fig2_data  %>% 
  filter(gen_cor != 0) %>% 
  ggplot(aes(x = prevalence, y = RR, fill = Relatedness,
             alpha = highlight)) +
  geom_bar(stat = "identity", position = position_dodge2(width = 0.5)) +
  scale_fill_manual(values = c("2nd degree" = "#2A788EFF",
                               "1st degree" = "#440154FF"), guide = "none") +
  scale_alpha_manual(values = c(0.25, 1), guide = "none") +
  facet_wrap(~ h2, labeller = label_both, nrow = 1) +
  geom_hline(yintercept = 1, col = "#7AD151FF") +
  labs(x = "Prevalence", y = "Relative risk") +
  scale_x_continuous(breaks = seq(0.04, 0.2, 0.04)) +
  coord_cartesian(xlim=c(0.02, max(prevalence))) +
  file_theme

# panel b: absolute risk across prevalence range 2% - 20%
absrisk <- fig2_data %>% 
  ggplot(aes(x = prevalence, y = K_fdr, color = Relatedness)) +
  geom_line(linewidth = 0.5) +
  geom_point(data = function(df) filter(df, highlight),
             aes(shape = Relatedness), size = 2) +
  facet_wrap(~ h2, labeller = label_both, nrow = 1) +
  labs(x = "Prevalence", y = "Absolute risk", fill = "Prevalence") +
  scale_colour_manual(values = c(unrelated = "#7AD151FF",
                                 "2nd degree" = "#2A788EFF",
                                 "1st degree" = "#440154FF")) +
  scale_x_continuous(breaks = seq(0.04, 0.2, 0.04)) +
  coord_cartesian(xlim = c(0.02, max(prevalence))) +
  file_theme

fig2 <- illustration / relrisk / absrisk

ggsave("Figure_2.png", path = "output", plot = fig2, bg = "white",
       width = 16, height = 16, units = "cm", dpi = 600)
ggsave("Figure_2.pdf", path = "output", plot = fig2, bg = "white",
       width = 16, height = 16, units = "cm")

