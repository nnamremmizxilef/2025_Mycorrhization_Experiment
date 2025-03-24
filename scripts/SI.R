########## SI FIGURES ##########


##### BASICS #####


# load packages
library(readxl)
library(tidyverse)
library(viridis)
library(psych)
library(corrplot)
library(BMS)

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
  "Co_Inok" = "Co-Inoculation",
  "Control" = "Control"
)

# define legend order
treatment_order <- c("Control", "Pilo", "Ceno", "Co_Inok")

# define custom colors
custom_colors <- c(
  "Control" = "#55C667FF",      # Light green
  "Ceno" = "#440154FF",         # Purple
  "Pilo" = "#FDE725FF",         # Yellow
  "Co_Inok" = "#39568CFF"       # Light blue
)


##### CREATE CORRELATION MATRIX HEATMAP #####


# calculate correlation matrix and p-values
corr_results <- corr.test(na.omit(pheno_data_num), method = "spearman")
corr_matrix <- as.matrix(corr_results$r)
p_matrix <- as.matrix(corr_results$p)

# create pdf to save plot
pdf("results/SI/correlation_plot.pdf", width = 10, height = 10)

# plot all correlations in upper triangle
corrplot(corr_matrix, 
         type = "full",
         method = "circle", 
         tl.col = "black",
         tl.srt = 90,
         tl.cex = 1,
         cl.cex = 1,
         diag = FALSE,
         title = "Spearman Correlation Matrix",
         mar = c(0, 0, 4, 0))

# plot only significant correlations in lower triangle
corrplot(corr_matrix, 
         add = TRUE,
         type = "lower",
         method = "circle", 
         diag = FALSE,
         tl.pos = "n",
         cl.pos = "n",
         p.mat = p_matrix,
         sig.level = 0.05,
         insig = "blank")

# add subtitle
mtext("The upper triangle shows all Spearman correlation R values;
      the lower triangle shows only significant correlations.", 
      side = 3,           # 3 = top
      line = -0.3,         # Position below the main title
      cex = 0.8,          # Smaller text size than title
      col = "black",
      font = 3)

# save pdf
dev.off()

# convert correlation matrices to dataframe
corr_df <- as.data.frame(corr_matrix)
corr_df <- corr_df %>%
  rownames_to_column(var="var")
p_df <- as.data.frame(p_matrix)
p_df <- p_df %>%
  rownames_to_column(var="var")

# save results
write_csv(corr_df, "results/SI/Spearman_R.csv")
write_csv(p_df, "results/SI/Spearman_p.csv")

