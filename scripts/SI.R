########## SI FIGURES ##########

##### BASICS #####

# load packages
library(tidyverse)
library(viridis)
library(psych)
library(corrplot)
library(ggpubr)
library(patchwork)
library(forcats)

# check and save session info (only run in new R session)
#sessionInfo() %>% capture.output(file = "results/SI/SI_session_info.txt")

##### DATA LOADING, CLEANING & SCALING #####

# laod data (dead plants are excluded in time_seriess.xlsx file)
pheno_data <- read_csv("data/phenotypic_data/pheno_data.csv")
randomization_data <- read_csv("data/experimental_design/randomization.csv")

# exclude contaminated and dead plants
pheno_data_final <- subset(pheno_data, Exclude == "no")

# create numeric dataframe excluding rows columns only "0" entries and character columns
pheno_data_num <- pheno_data_final %>%
  dplyr::select(-c("N_Other", "N_No", "Max_Stages_Reached")) %>%
  dplyr::select(where(is.numeric))

# load estimated marginal means data for resource allocation
resource_allocation_data <- read_csv(
  "results/Figure1/emm_resource_allocation.csv"
)


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


##### PLANT RANDOMIZATION PLOT #####

### RANDOMIZATION BY STAGE ###

# set plot order
randomization_data$Treatment <- factor(
  randomization_data$Treatment,
  levels = treatment_order
)

# calculate stage counts
stage_counts <- randomization_data %>%
  count(Treatment, Stage) %>%
  pivot_wider(
    names_from = Stage,
    values_from = n,
    values_fill = 0
  )

# convert to long format
stage_counts_long <- stage_counts %>%
  pivot_longer(
    cols = -Treatment,
    names_to = "Stage",
    values_to = "N"
  )

# plot stage counts
stage_plot <- ggplot(stage_counts_long) +
  geom_bar(aes(x = Treatment, y = N, fill = Stage), stat = "identity") +
  theme_pubr() +
  scale_fill_viridis_d(option = "C") +
  scale_y_continuous(limits = c(0, 24), breaks = seq(0, 24, 1)) +
  scale_x_discrete(labels = treatment_labels) +
  labs(
    title = "a",
    x = "Treatment",
    y = "N Plants"
  ) +
  theme(
    legend.position = "right",
    axis.text.x = element_text(size = 10),
    axis.text.y = element_text(size = 10),
    axis.title = element_text(size = 12),
    plot.title = element_text(face = "bold", size = 14),
    plot.subtitle = element_text(size = 12, face = "italic")
  )

# show plot
stage_plot


### RANDOMIZATION BY SIZE ###

# plot initial stem length
stem_length_plot <- ggplot(
  randomization_data %>%
    mutate(Treatment = factor(Treatment, levels = treatment_order)),
  aes(x = Treatment, y = Length_Init, fill = Treatment, color = Treatment)
) +
  geom_boxplot(outlier.alpha = 0.2, alpha = 0.2) +
  theme_pubr() +
  scale_fill_manual(values = custom_colors) +
  scale_color_manual(values = custom_colors) +
  scale_x_discrete(labels = treatment_labels) + # Add this line for custom x-axis labels
  theme(axis.title.x = element_blank()) +
  labs(
    title = "b",
    x = "Treatment",
    y = "Stem Length (cm)"
  ) +
  theme(
    legend.position = "none",
    axis.text.x = element_text(size = 10),
    axis.text.y = element_text(size = 10),
    axis.title = element_text(size = 12),
    plot.title = element_text(face = "bold", size = 14),
    plot.subtitle = element_text(size = 12, face = "italic")
  )

# show plot
stem_length_plot


### RANDOMIZATION BY ROOTING DATE ###

# plot initial stem length
date_plot <- ggplot(
  randomization_data %>%
    mutate(Treatment = factor(Treatment, levels = treatment_order)),
  aes(
    x = Treatment,
    y = Rooting_Date_Divergence,
    fill = Treatment,
    color = Treatment
  )
) +
  geom_boxplot(outlier.alpha = 0.2, alpha = 0.2) +
  theme_pubr() +
  scale_fill_manual(values = custom_colors) +
  scale_color_manual(values = custom_colors) +
  scale_x_discrete(labels = treatment_labels) + # Add this line for custom x-axis labels
  theme(axis.title.x = element_blank()) +
  labs(
    title = "c",
    x = "Treatment",
    y = "N Days Divergence from Mean Rooting Date"
  ) +
  theme(
    legend.position = "none",
    axis.text.x = element_text(size = 10),
    axis.text.y = element_text(size = 10),
    axis.title = element_text(size = 12),
    plot.title = element_text(face = "bold", size = 14),
    plot.subtitle = element_text(size = 12, face = "italic")
  )

# show plot
date_plot


#### FINAL VISUALIZATION ####

# create combined plot
combined_plot_randomization <- stage_plot /
  stem_length_plot /
  date_plot +
  plot_layout(widths = c(5), heights = (10)) +
  plot_annotation(
    title = "",
    theme = theme(
      plot.title = element_text(size = 16, face = "bold", hjust = 0.5)
    )
  )

# show plot
combined_plot_randomization

# save combined plot (has to be combined with experimental setup plot)
ggsave(
  "results/SI/Randomization_Plot.pdf",
  plot = combined_plot_randomization,
  width = 12.5,
  height = 15
)
dev.off()


##### CORRELATION MATRIX HEATMAP #####

# calculate correlation matrix and p-values
corr_results <- corr.test(na.omit(pheno_data_num), method = "spearman")
corr_matrix <- as.matrix(corr_results$r)
p_matrix <- as.matrix(corr_results$p)

# create pdf to save plot
pdf("results/SI/Correlation_Plot.pdf", width = 10, height = 10)

# plot all correlations in upper triangle
corrplot(
  corr_matrix,
  type = "full",
  method = "circle",
  tl.col = "black",
  tl.srt = 90,
  tl.cex = 1,
  cl.cex = 1,
  diag = FALSE,
  title = "",
  mar = c(0, 0, 0, 0)
)

# plot only significant correlations in lower triangle
corrplot(
  corr_matrix,
  add = TRUE,
  type = "lower",
  method = "circle",
  diag = FALSE,
  tl.pos = "n",
  cl.pos = "n",
  p.mat = p_matrix,
  sig.level = 0.05,
  insig = "blank"
)

# save pdf
dev.off()

# convert correlation matrices to dataframe
corr_df <- as.data.frame(corr_matrix)
corr_df <- corr_df %>%
  rownames_to_column(var = "var")
p_df <- as.data.frame(p_matrix)
p_df <- p_df %>%
  rownames_to_column(var = "var")

# save results
write_csv(corr_df, "results/SI/Spearman_R.csv")
write_csv(p_df, "results/SI/Spearman_p.csv")


##### RESOURCE ALLOCATION PROPORTIONS #####

# calculate proportions
resource_allocation_data <- resource_allocation_data %>%
  group_by(Treatment) %>%
  mutate(proportion = response / sum(response)) %>%
  ungroup()

# convert proportions to percent
resource_allocation_data <- resource_allocation_data %>%
  mutate(percentage = proportion * 100)


# create reversed factor for plotting shoot over root
resource_allocation_data <- resource_allocation_data %>%
  mutate(Part = fct_rev(Part))

# create stacked barplot
resource_proportions_plot <- ggplot(
  resource_allocation_data,
  aes(
    x = factor(Treatment, levels = treatment_order),
    y = percentage,
    fill = Part
  )
) +
  geom_bar(stat = "identity", position = "stack", width = 0.7) +
  geom_hline(
    yintercept = 50,
    linetype = "dashed",
    color = "black",
    size = 0.5
  ) +
  scale_y_continuous(
    name = "Proportion (%)",
    limits = c(0, 100),
    breaks = seq(0, 100, by = 25)
  ) +
  scale_fill_viridis_d(option = "A") +
  scale_x_discrete(labels = treatment_labels) +
  labs(title = "", x = "Treatment") +
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
resource_proportions_plot

# save plot

ggsave(
  "results/SI/Resource_Allocation_Proportions.pdf",
  plot = resource_proportions_plot,
  width = 12.5,
  height = 5
)
dev.off()
