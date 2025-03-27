########## FIGURE 1 - PHENOTYPIC DATA ##########

##### BASICS #####

# load packages
library(readxl)
library(tidyverse)
library(viridis)
library(magick)
library(grid)
library(ggpubr)
library(psych)
library(rstatix)
library(emmeans)
library(ggbeeswarm)
library(patchwork)

# check and save session info
sessionInfo() %>%
  capture.output(file = "results/Figure1/Figure1_session_info.txt")


##### DATA LOADING, CLEANING & SCALING #####

# laod data
pheno_data <- read_csv("data/phenotypic_data/pheno_data.csv")

# exclude contaminated and dead plants
pheno_data_final <- subset(pheno_data, Exclude == "no")

# create numeric dataframe excluding rows columns only "0" entries and character columns
pheno_data_num <- pheno_data_final %>%
  dplyr::select(-c("N_Other", "N_No")) %>%
  dplyr::select(where(is.numeric))


##### DEFINE BASIC PLOTTING PARAMETERS #####

# define treatment labels
treatment_labels <- c(
  "Ceno" = expression(italic("Cenococcum geophilum")),
  "Pilo" = expression(italic("Piloderma croceum")),
  "Co_Inoc" = "Co-Inoculation",
  "Control" = "Control"
)

# define legend order
treatment_order <- c("Control", "Pilo", "Ceno", "Co_Inoc")

# define custom colors
custom_colors <- c(
  "Control" = "#55C667FF", # Light green
  "Ceno" = "#440154FF", # Purple
  "Pilo" = "#FDE725FF", # Yellow
  "Co_Inoc" = "#39568CFF" # Light blue
)


##### PLOT DRY WEIGHT COMPARISONS #####

### do statistical comparison between treatments for durations of each stage ###

# do wilcoxon test with bonferroni p adjustment for plant dry weight
stat_test_dw_plant <- pheno_data_final %>%
  wilcox_test(DW_Plant ~ Treatment) %>%
  adjust_pvalue(method = "bonferroni") %>%
  add_significance()

# save_results
dw_plant_results <- subset(stat_test_dw_plant, select = -c(.y.))
write.csv(
  dw_plant_results,
  "results/Figure1/dw_plant_wilcoxon.csv",
  row.names = FALSE
)


### plot results ###

# set plot order
pheno_data_final$Treatment <- factor(
  pheno_data_final$Treatment,
  levels = treatment_order
)

# create plot
dw_plant_plot <- ggplot(
  pheno_data_final,
  aes(x = Treatment, y = DW_Plant, fill = Treatment)
) +
  geom_boxplot(
    outlier.shape = 16,
    outlier.size = 2,
    alpha = 0.2,
    lwd = 0.5,
    aes(color = Treatment)
  ) +
  scale_fill_manual(
    values = custom_colors,
    labels = treatment_labels,
    breaks = treatment_order
  ) +
  scale_color_manual(
    values = custom_colors,
    labels = treatment_labels,
    breaks = treatment_order
  ) +
  scale_x_discrete(labels = treatment_labels) + # format x-axis
  theme_classic() +
  geom_signif(
    y_position = c(1.66, 1.8),
    xmin = c(2, 2),
    xmax = c(3, 4),
    annotation = c("< 0.001", "0.004")
  ) +
  labs(
    title = "B    Total Biomass Produced",
    subtitle = "        Only significant p-values are shown (Wilcoxon test with Bonferroni correction) ",
    x = "Treatment",
    y = "Dry Weight (g)"
  ) +
  theme_pubr() +
  theme(
    legend.position = "none",
    axis.text.x = element_text(size = 10),
    axis.text.y = element_text(size = 10),
    axis.title = element_text(size = 12),
    plot.title = element_text(face = "bold", size = 14),
    plot.subtitle = element_text(size = 12, face = "italic"),
    panel.grid.major.y = element_line(color = "gray90"),
    panel.grid.minor.y = element_line(color = "gray95")
  )

# show plot
dw_plant_plot


##### PLOT DRY WEIGHT ABOVE AND BELOW GROUND #####

### do statistical comparison between treatments for durations of each stage ###

# do wilcoxon test with bonferroni p adjustment for plant dry weight
stat_test_above <- pheno_data_final %>%
  wilcox_test(DW_Shoot ~ Treatment) %>%
  adjust_pvalue(method = "bonferroni") %>%
  add_significance()
stat_test_above$Part <- "above ground"
stat_test_below <- pheno_data_final %>%
  wilcox_test(DW_Roots_Total ~ Treatment) %>%
  adjust_pvalue(method = "bonferroni") %>%
  add_significance()
stat_test_below$Part <- "below ground"

# combine results
wilcoxon_results_above_below <- rbind(stat_test_above, stat_test_below)
print(wilcoxon_results_above_below)

# save_results
wilxcoxon_reuslts_above_below <- subset(
  wilcoxon_results_above_below,
  select = -c(.y.)
)
write.csv(
  wilcoxon_results_above_below,
  "results/Figure1/above_below_ground_wilcoxon.csv",
  row.names = FALSE
)


### plot results ###

# create long format
pheno_data_long_above_below <- pivot_longer(
  pheno_data_final,
  cols = c("DW_Roots_Total", "DW_Shoot"),
  names_to = "Part",
  values_to = "Part_DW"
)

# set treatment labels
treatment_labels_above_below <- c("Shoot", "Root")

# set plot order
pheno_data_long_above_below$Part <- factor(
  pheno_data_long_above_below$Part,
  levels = c("DW_Shoot", "DW_Roots_Total")
)
treatment_order_above_below <- c("Shoot", "Root")

# create plot
dw_above_below_plot <- ggplot(
  pheno_data_long_above_below,
  aes(x = Part, y = Part_DW, fill = Treatment)
) +
  geom_boxplot(
    outlier.shape = 16,
    outlier.size = 2,
    alpha = 0.2,
    lwd = 0.5,
    aes(color = Treatment)
  ) +
  scale_fill_manual(
    values = custom_colors,
    labels = treatment_labels,
    breaks = treatment_order
  ) +
  scale_color_manual(
    values = custom_colors,
    labels = treatment_labels,
    breaks = treatment_order
  ) +
  scale_x_discrete(labels = treatment_labels_above_below) +
  geom_signif(
    y_position = c(0.98, 1.05, 0.86, 0.93),
    xmin = c(0.906, 0.906, 1.906, 1.906),
    xmax = c(1.094, 1.282, 2.0935, 2.282),
    annotation = c("< 0.001", "0.012", "< 0.001", "< 0.001")
  ) +
  labs(
    title = "C    Shoot vs. Root Biomass Produced",
    subtitle = "        Only significant p-values are shown (Wilcoxon test with Bonferroni correction) ",
    x = "Plant Part",
    y = "Dry Weight (g)"
  ) +
  theme_pubr() +
  theme(
    legend.position = "right",
    axis.text.x = element_text(size = 10),
    axis.text.y = element_text(size = 10),
    axis.title = element_text(size = 12),
    plot.title = element_text(face = "bold", size = 14),
    plot.subtitle = element_text(size = 12, face = "italic"),
    panel.grid.major.y = element_line(color = "gray90"),
    panel.grid.minor.y = element_line(color = "gray95")
  )

# show plot
dw_above_below_plot


##### BIOMASS ALLOCATION PATTERNS ######

### do statistical test ###

# create long format
pheno_data_long_allocation <- pivot_longer(
  pheno_data_final,
  cols = c("DW_Shoot", "DW_Roots_Total"),
  names_to = "Part",
  values_to = "Part_DW"
)

# create data subsets
subset_DW_Shoot <- subset(pheno_data_long_allocation, Part == "DW_Shoot")
subset_DW_Roots <- subset(pheno_data_long_allocation, Part == "DW_Roots_Total")

# calculate GLMs (sum)
model_DW_Shoot <- glm(
  Part_DW ~ Treatment + DW_Plant,
  data = subset_DW_Shoot,
  family = gaussian()
)
model_DW_Roots <- glm(
  Part_DW ~ Treatment + DW_Plant,
  data = subset_DW_Roots,
  family = gaussian()
)

# check distributions of residuals
hist(residuals(model_DW_Shoot))
qqnorm(residuals(model_DW_Shoot))
qqline(residuals(model_DW_Shoot))
hist(residuals(model_DW_Roots))
qqnorm(residuals(model_DW_Roots))
qqline(residuals(model_DW_Roots))

# get summary statistics
summary(model_DW_Shoot)
summary(model_DW_Roots)

# calculate and plot estimated means
emm_DW_Shoot <- emmeans(
  model_DW_Shoot,
  ~ Treatment | DW_Plant,
  type = "response",
  adjust = "sidak"
)
plot(emm_DW_Shoot, comparisons = TRUE)
emm_DW_Shoot_results <- as.data.frame(emm_DW_Shoot)
emm_DW_Shoot_results$Part <- "Shoot"
emm_DW_Roots <- emmeans(
  model_DW_Roots,
  ~ Treatment | DW_Plant,
  type = "response",
  adjust = "sidak"
)
plot(emm_DW_Roots, comparisons = TRUE)
emm_DW_Roots_results <- as.data.frame(emm_DW_Roots)
emm_DW_Roots_results$Part <- "Roots"

# combine results
emm_results <- rbind(emm_DW_Shoot_results, emm_DW_Roots_results)

# save_results
emm_results <- subset(emm_results, select = -c(DW_Plant))
write.csv(
  emm_results,
  "results/Figure1/emm_resource_allocation.csv",
  row.names = FALSE
)

# calculate pairwise comparison statistics with sidak correction
pairs_Shoot <- as.data.frame(pairs(emm_DW_Shoot, adjust = "sidak"))
pairs_Shoot$Part <- "Shoot"
pairs_Roots <- as.data.frame(pairs(emm_DW_Roots, adjust = "sidak"))
pairs_Roots$Part <- "Roots"

# combine results
emm_pairwise <- rbind(pairs_Shoot, pairs_Roots)

# save_results
emm_pairwise <- subset(emm_pairwise, select = -c(DW_Plant))
write.csv(
  emm_pairwise,
  "results/Figure1/pairwise_comparison_resource_allocation.csv",
  row.names = FALSE
)


### plot results ###

# define plot order
plot_order_allocation <- c("Shoot", "Roots")

# set treatment order
pheno_data_long_allocation$Part <- factor(pheno_data_long_allocation$Part)
levels(pheno_data_long_allocation$Part) <- c("Roots", "Shoot")

# set part order
emm_results$Part <- factor(emm_results$Part, levels = plot_order_allocation)

# set treatment labels
treatment_labels_allocation <- c("Shoot", "Root")

# create plot
emm_plot <- ggplot() +
  geom_point(
    data = emm_results,
    aes(x = Part, y = emmean, color = Treatment),
    position = position_dodge(width = 0.75),
    size = 3
  ) +
  geom_linerange(
    data = emm_results,
    aes(
      x = Part,
      y = emmean,
      color = Treatment,
      ymin = lower.CL,
      ymax = upper.CL
    ),
    position = position_dodge(width = 0.75),
    size = 1.2
  ) +
  geom_beeswarm(
    data = pheno_data_long_allocation,
    aes(x = Part, y = Part_DW, color = Treatment),
    alpha = 0.25,
    dodge.width = 0.75,
    shape = 16
  ) +
  scale_color_manual(
    values = custom_colors,
    labels = treatment_labels,
    breaks = treatment_order
  ) +
  scale_x_discrete(labels = c("Shoot", "Root")) +
  labs(
    title = "D    Treatment Effects on Plant Biomass Allocation",
    subtitle = "        Plant dry weight is used as a model co-factor; only significant p-values are shown (pairwise t-statistic test with Sidak correction)",
    x = "Plant Part",
    y = "Dry Weight Estimated Marginal Mean (g)",
    color = "Treatment"
  ) +
  annotate("segment", x = 0.905, xend = 1.095, y = 0.95, yend = 0.95) +
  annotate("segment", x = 0.905, xend = 0.905, y = 0.922, yend = 0.95) +
  annotate("segment", x = 1.095, xend = 1.095, y = 0.922, yend = 0.95) +
  annotate("text", x = 1.0, y = 0.98, label = "0.024") +
  annotate("segment", x = 0.725, xend = 1.095, y = 1.04, yend = 1.04) +
  annotate("segment", x = 0.725, xend = 0.725, y = 1.012, yend = 1.04) +
  annotate("segment", x = 1.095, xend = 1.095, y = 1.012, yend = 1.04) +
  annotate("text", x = 0.91, y = 1.07, label = "0.049") +
  theme_pubr() +
  theme(
    axis.text.x = element_text(size = 10),
    axis.text.y = element_text(size = 10),
    axis.title = element_text(size = 12),
    legend.position = "right",
    plot.title = element_text(face = "bold", size = 14),
    plot.subtitle = element_text(size = 12, face = "italic"),
    panel.grid.major.y = element_line(color = "gray90"),
    panel.grid.minor.y = element_line(color = "gray95")
  )

# show plot
emm_plot


##### FINAL VISUALIZATION #####

# create combined plot
combined_plot_pheno <- dw_plant_plot /
  dw_above_below_plot /
  emm_plot +
  plot_layout(widths = c(5), heights = (10)) +
  plot_annotation(
    title = "Phenotypic Variation Analysis",
    theme = theme(
      plot.title = element_text(size = 16, face = "bold", hjust = 0.5)
    )
  )

# show plot
combined_plot_pheno

# save combined plot (has to be combined with experimental setup plot)
ggsave(
  "results/Figure1/Figure1_R.pdf",
  plot = combined_plot_pheno,
  width = 12.5,
  height = 15
)
dev.off()
