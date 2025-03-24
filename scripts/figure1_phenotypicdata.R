########## FIGURE 1 - PHENOTYPIC DATA ##########


##### BASICS #####


# load packages
library(readxl)
library(tidyverse)
library(viridis)
library(pdftools)
library(magick)
library(grid)
library(ggpubr)
library(psych)
library(rstatix)
library(patchwork)

# check package versions
sessionInfo()



##### DATA LOADING, CLEANING & SCALING #####


# laod data (dead plants are excluded in time_seriess.xlsx file)
pheno_data <- read_excel("data/phenotypic_data/pheno_data.xlsx")

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
  "Control" = "#55C667FF",      # Light green
  "Ceno" = "#440154FF",         # Purple
  "Pilo" = "#FDE725FF",         # Yellow
  "Co_Inoc" = "#39568CFF"       # Light blue
)



#### PLOT DRY WEIGHT COMPARISONS ####


### do statistical comparison between treatments for durations of each stage ###


# do wilcoxon test with bonferroni p adjustment for plant dry weight
stat_test_dw_plant <- pheno_data_final %>%
  wilcox_test(DW_Plant ~ Treatment) %>%
  adjust_pvalue(method = "bonferroni") %>%
  add_significance()

# save_results
dw_plant_results <- subset(stat_test_dw_plant, select = -c(.y.))
write.csv(dw_plant_results, "results/Figure1/dw_plant_wilcoxon.csv", row.names = FALSE)




### plot results ###


# set treatment order
pheno_data_final$Treatment <- factor(pheno_data_final$Treatment, levels = treatment_order)

# create plot
dw_plant_plot <- ggplot(pheno_data_final, aes(x = Treatment, y = DW_Plant, fill = Treatment)) +
  geom_boxplot(outlier.shape = 16, outlier.size = 2, alpha = 0.2, lwd = 0.5, aes(color = Treatment)) +
  scale_fill_manual(values = custom_colors, labels = treatment_labels, breaks = treatment_order) +
  scale_color_manual(values = custom_colors, labels = treatment_labels, breaks = treatment_order) +
  scale_x_discrete(labels = treatment_labels) + # format x-axis
  theme_classic() +
  geom_signif(
    y_position = c(1.66, 1.8), xmin = c(2, 2), xmax = c(3, 4),
    annotation = c("< 0.001", "0.004")
  ) +
  labs(title = "B    Total Biomass Produced",
       subtitle = "        Only significant p-values are shown (Wilcoxon test with Bonferroni correction) ",
       x = "Treatment",
       y = "Plant Dry Weight (g)"
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
dw_plant_plot



#### PLOT DRY WEIGHT ABOVE AND BELOW GROUND ####


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
wilxcoxon_reuslts_above_below <- subset(wilcoxon_results_above_below, select = -c(.y.))
write.csv(wilcoxon_results_above_below, "results/Figure1/above_below_ground_wilcoxon.csv", row.names = FALSE)



### plot results ###


# create long format
pheno_data_long <- pivot_longer(pheno_data_final, cols = c("DW_Roots_Total", "DW_Shoot"), names_to = "Part", values_to = "Part_DW")

# set treatment labels
treatment_labels_above_below <- c("Below Ground", "Above Ground")

# set treatment labels
treatment_order_above_below <- c("Below Ground", "Above Ground")

# create plot
dw_above_below_plot <- ggplot(pheno_data_long, aes(x = Part, y = Part_DW, fill = Treatment)) +
  geom_boxplot(outlier.shape = 16, outlier.size = 2, alpha = 0.2, lwd = 0.5, aes(color = Treatment)) +
  scale_fill_manual(values = custom_colors, labels = treatment_labels_above_below, breaks = treatment_order_above_below) +
  scale_color_manual(values = custom_colors, labels = treatment_labels_above_below, breaks = treatment_order_above_below) +
  scale_x_discrete(labels = treatment_labels_above_below) + # format x-axis
  theme_classic() +
  geom_signif(
    y_position = c(0.86, 0.93, 0.98, 1.05), xmin = c(0.906, 0.906, 1.906, 1.906), xmax = c(1.094, 1.282, 2.0935, 2.282),
    annotation = c("< 0.001", "0.012", "< 0.001", "< 0.001")
  ) +
  labs(title = "C    Above- vs. Below-Ground Biomass Produced",
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



##### MYCORRHIZATION EFFECTS ######


# create long format
pheno_data_weight_long <- pivot_longer(pheno_data_final, cols = c("DW_Shoot", "DW_Roots_Total", "DW_LR", "DW_PR"), names_to = "Part", values_to = "Part_DW")

# define histogram order
part_order <- c("DW_Shoot", "DW_Roots_Total", "DW_PR", "DW_LR")

# set treatment order
pheno_data_weight_long$Treatment <- factor(pheno_data_weight_long$Treatment, levels = treatment_order)

# set treatment labels
treatment_labels_weight <- c("DW Shoot", "DW Total Roots", "DW Principal Roots", "DW Lateral Roots")

# set treatment labels
treatment_order_weight <- c("DW Shoot", "DW Total Roots", "DW Principal Roots", "DW Lateral Roots")

# create plot
dw_comparison <- ggplot(pheno_data_weight_long, aes(x = Part, y = Part_DW, fill = Treatment)) +
  geom_boxplot(outlier.shape = 16, outlier.size = 2, alpha = 0.2, lwd = 0.5, aes(color = Treatment)) +
  scale_fill_manual(values = custom_colors, labels = treatment_labels_weight, breaks = treatment_order_weight) +
  scale_color_manual(values = custom_colors, labels = treatment_labels_weight, breaks = treatment_order_weight) +
  scale_x_discrete(limits = part_order, labels = treatment_labels_weight) + # format x-axis
  theme_classic() +
  geom_signif(
    y_position = c(0.86, 0.93, 0.98, 1.05), xmin = c(0.906, 0.906, 1.906, 1.906), xmax = c(1.094, 1.282, 2.0935, 2.282),
    annotation = c("< 0.001", "0.012", "< 0.001", "< 0.001")
  ) +
  labs(title = "C    Above- vs. Below-Ground Biomass Produced",
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
dw_comparison







# calculate mycorrhization rates
pheno_data_final$rate_ceno <- pheno_data_final$N_Ceno / (pheno_data_final$N_Ceno + pheno_data_final$N_Pilo + pheno_data_final$N_No)
pheno_data_final$water_cont_total <- (pheno_data_final$FW_Plant - pheno_data_final$DW_Plant) / pheno_data_final$FW_Plant
pheno_data_final$water_cont_roots <- (pheno_data_final$FW_Roots_Total - pheno_data_final$DW_Roots_Total) / pheno_data_final$FW_Roots_Total
pheno_data_ceno <- subset(pheno_data_final, Treatment == c("Ceno", "Co_Inoc"))



model <- lm(water_cont_roots ~ rate_ceno, data = pheno_data_ceno)
r2 <- summary(model)$r.squared
pval <- summary(model)$coefficients[2,4]
r2
pval






##### FINAL VISUALIZATION #####

# create combined plot
combined_plot_pheno <- dw_plant_plot /
  dw_above_below_plot +
  plot_layout(widths = c(5), heights = (10)) +
  plot_annotation(
    title = "Phenotypic Variation Analysis",
    theme = theme(
      plot.title = element_text(size = 16, face = "bold", hjust = 0.5)
    )
  )

# show plot
combined_plot_pheno

# save combined plot
ggsave("results/Figure2/Figure2.pdf", plot = combined_plot_pheno, width = 12.5, height = 20)
dev.off()


