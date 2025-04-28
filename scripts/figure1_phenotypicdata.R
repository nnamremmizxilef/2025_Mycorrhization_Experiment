########## FIGURE 1 - PHENOTYPIC DATA ##########

##### BASICS #####

# load packages
library(tidyverse)
library(viridis)
library(magick)
library(grid)
library(ggpubr)
library(psych)
library(rstatix)
library(MASS)
library(emmeans)
library(ggbeeswarm)
library(patchwork)

# check and save session info (only run in new R session)
#sessionInfo() %>% capture.output(file = "results/Figure1/Figure1_session_info.txt")

##### DATA LOADING, CLEANING & SCALING #####

# laod data
pheno_data <- read_csv("data/phenotypic_data/pheno_data.csv")

# exclude contaminated and dead plants
pheno_data_final <- subset(pheno_data, Exclude == "no")


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
  "Control" = "#4B8C6F",
  "Ceno" = "#8C5E8C",
  "Pilo" = "#D6C9A0",
  "Co_Inoc" = "#5E96A6"
)


##### PLOT DRY WEIGHT COMPARISONS #####

### do statistical comparison between treatments for durations of each stage ###

# do wilcoxon test with Benjamini-Hochberg p adjustment for plant dry weight
stat_test_dw_plant <- pheno_data_final %>%
  wilcox_test(DW_Plant ~ Treatment) %>%
  adjust_pvalue(method = "BH") %>%
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
    alpha = 0.5,
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
    y_position = c(1.8, 1.94, 1.66),
    xmin = c(2, 2, 1),
    xmax = c(3, 4, 3),
    annotation = c("< 0.001", "0.002", "0.036"),
    textsize = 4
  ) +
  labs(
    title = "b",
    x = "Treatment",
    y = "Dry Weight (g)"
  ) +
  theme_pubr() +
  theme(
    legend.position = "none",
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12),
    axis.title = element_text(size = 14),
    plot.title = element_text(
      face = "bold",
      size = 20,
      vjust = 2,
      hjust = -0.075
    ),
    panel.grid.major.y = element_line(color = "gray90"),
    panel.grid.minor.y = element_line(color = "gray95")
  )

# show plot
dw_plant_plot


##### PLOT DRY WEIGHT ABOVE AND BELOW GROUND #####

### do statistical comparison between treatments for durations of each stage ###

# do wilcoxon test with Benjamini-Hochberg p adjustment for plant dry weight
stat_test_above <- pheno_data_final %>%
  wilcox_test(DW_Shoot ~ Treatment) %>%
  adjust_pvalue(method = "BH") %>%
  add_significance()
stat_test_above$Part <- "Shoot"
stat_test_below <- pheno_data_final %>%
  wilcox_test(DW_Roots_Total ~ Treatment) %>%
  adjust_pvalue(method = "BH") %>%
  add_significance()
stat_test_below$Part <- "Root"

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
    alpha = 0.5,
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
    y_position = c(1.05, 1.12, 0.98, 0.93, 1, 0.86),
    xmin = c(0.906, 0.906, 0.716, 1.906, 1.906, 1.72),
    xmax = c(1.094, 1.282, 1.094, 2.0935, 2.282, 1.906),
    annotation = c("< 0.001", "< 0.001", "0.030", "< 0.001", "0.006", "0.048"),
    textsize = 4
  ) +
  labs(
    title = "c",
    x = "Plant Part",
    y = "Dry Weight (g)"
  ) +
  theme_pubr() +
  theme(
    legend.position = "right",
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12),
    axis.title = element_text(size = 14),
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 12),
    plot.title = element_text(
      face = "bold",
      size = 20,
      hjust = -0.075,
      vjust = 2
    ),
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
model_DW_Shoot_no <- glm(
  (Part_DW) ~ Treatment + DW_Plant,
  data = subset_DW_Shoot,
  family = gaussian()
)
model_DW_Roots_no <- glm(
  (Part_DW) ~ Treatment + DW_Plant,
  data = subset_DW_Roots,
  family = gaussian()
)
bc_shoot <- boxcox(Part_DW ~ Treatment + DW_Plant, data = subset_DW_Shoot)
lambda_shoot <- bc_shoot$x[which.max(bc_shoot$y)]
model_DW_Shoot_bc <- glm(
  ((Part_DW^lambda_shoot - 1) / lambda_shoot) ~ Treatment + DW_Plant,
  data = subset_DW_Shoot,
  family = gaussian()
)
bc_roots <- boxcox(Part_DW ~ Treatment + DW_Plant, data = subset_DW_Roots)
lambda_roots <- bc_roots$x[which.max(bc_roots$y)]
model_DW_Roots_bc <- glm(
  ((Part_DW^lambda_roots - 1) / lambda_roots) ~ Treatment + DW_Plant,
  data = subset_DW_Roots,
  family = gaussian()
)
model_DW_Shoot_log <- glm(
  log(Part_DW) ~ Treatment + DW_Plant,
  data = subset_DW_Shoot,
  family = gaussian()
)
model_DW_Roots_log <- glm(
  log(Part_DW) ~ Treatment + DW_Plant,
  data = subset_DW_Roots,
  family = gaussian()
)

# check AIC and BIC
AIC(model_DW_Shoot_bc, model_DW_Shoot_log, model_DW_Shoot_no)
AIC(model_DW_Roots_bc, model_DW_Roots_log, model_DW_Roots_no)
BIC(model_DW_Shoot_bc, model_DW_Shoot_log, model_DW_Shoot_no)
BIC(model_DW_Roots_bc, model_DW_Roots_log, model_DW_Roots_no)
# use no transformation

# check distributions of residuals
hist(residuals(model_DW_Shoot_no))
qqnorm(residuals(model_DW_Shoot_no))
qqline(residuals(model_DW_Shoot_no))
hist(residuals(model_DW_Roots_no))
qqnorm(residuals(model_DW_Roots_no))
qqline(residuals(model_DW_Roots_no))

# get summary statistics
summary(model_DW_Shoot_no)
summary(model_DW_Roots_no)

# calculate and plot estimated means
emm_DW_Shoot <- emmeans(
  model_DW_Shoot_no,
  ~ Treatment | DW_Plant,
  type = "response",
  adjust = "BH"
)
plot(emm_DW_Shoot, comparisons = TRUE)
emm_DW_Shoot_results <- as.data.frame(emm_DW_Shoot)
emm_DW_Shoot_results$Part <- "Shoot"
emm_DW_Roots <- emmeans(
  model_DW_Roots_no,
  ~ Treatment | DW_Plant,
  type = "response",
  adjust = "BH"
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

# calculate pairwise comparison statistics with Benjamini-Hochberg correction
pairs_Shoot <- as.data.frame(pairs(
  emm_DW_Shoot,
  adjust = "BH",
  type = "response"
))
pairs_Shoot$Part <- "Shoot"
pairs_Roots <- as.data.frame(pairs(
  emm_DW_Roots,
  adjust = "BH",
  type = "response"
))
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
    aes(x = Part, y = response, color = Treatment),
    position = position_dodge(width = 0.75),
    size = 3
  ) +
  geom_linerange(
    data = emm_results,
    aes(
      x = Part,
      y = response,
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
    alpha = 0.5,
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
    title = "d",
    x = "Plant Part",
    y = "Dry Weight Estimated Marginal Mean (g)",
    color = "Treatment"
  ) +
  annotate("segment", x = 0.905, xend = 1.095, y = 0.95, yend = 0.95) +
  annotate("segment", x = 0.905, xend = 0.905, y = 0.922, yend = 0.95) +
  annotate("segment", x = 1.095, xend = 1.095, y = 0.922, yend = 0.95) +
  annotate("text", x = 1.0, y = 0.98, label = "0.025", size = 4) +
  annotate("segment", x = 0.725, xend = 1.095, y = 1.04, yend = 1.04) +
  annotate("segment", x = 0.725, xend = 0.725, y = 1.012, yend = 1.04) +
  annotate("segment", x = 1.095, xend = 1.095, y = 1.012, yend = 1.04) +
  annotate("text", x = 0.91, y = 1.07, label = "0.025", size = 4) +
  theme_pubr() +
  theme(
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12),
    axis.title = element_text(size = 14),
    legend.position = "right",
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 12),
    plot.title = element_text(
      face = "bold",
      size = 20,
      hjust = -0.075,
      vjust = 2
    ),
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
  plot_layout(widths = c(5), heights = (10))

# show plot
combined_plot_pheno

# save combined plot (has to be combined with experimental setup plot)
ggsave(
  "results/Figure1/Figure1.pdf",
  plot = combined_plot_pheno,
  width = 12.5,
  height = 15
)
dev.off()
