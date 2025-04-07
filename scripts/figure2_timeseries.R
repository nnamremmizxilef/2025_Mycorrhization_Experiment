########## FIGURE 2 - TIME SERIES ANALYSIS ##########

##### BASICS #####

# load packages
library(tidyverse)
library(viridis)
library(ggpubr)
library(brms)
library(bayesnec)
library(bayestestR)
library(tidybayes)
library(ggtext)
library(rstatix)
library(ggsignif)
library(grid)
library(patchwork)

# check and save session info (only run in new R session)
#sessionInfo() %>% capture.output(file = "results/Figure2/Figure2_session_info.txt")

##### DATA LOADING, CLEANING & SCALING #####

# laod data (dead plants are excluded in time_seriess.csv file)
time_data <- read_csv("data/time_series_data/time_series.csv")

# create stage_reference for calculation of stages reached
stage_reference <- data.frame(
  Stage = c(paste0(rep(1:5, each = 4), rep(c("A", "B", "C", "D"), 5))),
  stage_index = 1:20
)

# calculate the number of stages reached for each date
time_data <- time_data %>%
  left_join(stage_reference, by = "Stage") %>%
  group_by(ID) %>%
  arrange(ID, Days) %>%
  mutate(
    initial_index = first(stage_index),
    Stages_Reached = ifelse(Days == 0, 0, stage_index - initial_index)
  ) %>%
  dplyr::select(-c(stage_index, initial_index))

# create a column where the total number of stages reached is added for each ID
time_data <- time_data %>%
  group_by(ID) %>%
  mutate(Max_Stages_Reached = max(Stages_Reached)) %>%
  ungroup()

# exclude contaminated plants
time_data_final <- subset(time_data, Contaminated == "no")

# scale time variable, remove T0 where 0 stages are reached (cumulative distribution family requires positive integers), define plant ID & Treatment as factor
time_data_scaled <- time_data_final
time_data_scaled$days_scaled <- scale(time_data_final$Days)
time_data_scaled$ID <- factor(time_data_scaled$ID)
time_data_scaled$Treatment <- factor(time_data_scaled$Treatment)
time_data_scaled <- subset(time_data_scaled, Stages_Reached > 0)


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


##### CREATE BASIC GROWTH STAGE DEVELOPMENT PLOT #####

# remove week 10 (only 1 set reached week 10)
time_data_basic_plot <- subset(time_data_final, Week < 9)

# create basic plot for growth stage development (mean and 95% confidence interval)
basic_growth_plot <- ggplot(
  time_data_basic_plot,
  aes(x = Week, y = Stages_Reached, color = Treatment)
) +
  stat_summary(fun = mean, geom = "line", linewidth = 1) +
  stat_summary(
    fun.data = mean_se,
    geom = "ribbon",
    alpha = 0.2,
    aes(fill = Treatment)
  ) +
  scale_color_manual(
    values = custom_colors,
    labels = treatment_labels,
    breaks = treatment_order
  ) +
  scale_fill_manual(
    values = custom_colors,
    labels = treatment_labels,
    breaks = treatment_order
  ) +
  scale_x_continuous(
    expand = c(0, 0),
    limits = c(0, 8),
    breaks = (c(0, 1, 2, 3, 4, 5, 6, 7, 8, 9))
  ) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 7.25)) +
  labs(
    title = "b",
    y = "N Stages Reached",
    x = "N Weeks",
    color = "Treatment",
    fill = "Treatment"
  ) +
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
basic_growth_plot


##### GROWTH CURVE ANALYSIS - BRM (Cumulative Distribution Family) #####

### FIT AND EVALUATE MODELS ###

# set control treatment as reference
time_data_scaled$Treatment <- relevel(
  time_data_scaled$Treatment,
  ref = "Control"
)

# fit growth curve model with cumulative distribution family using bayesian regression (Treatment * days_scaled + random)
bay_model <- brm(
  Stages_Reached ~ Treatment * days_scaled + (1 | ID), # plant ID as random factor to account for non-independent measurements (1 wo random effect, 1 sum, 1 treatment only, )
  family = cumulative(), # cumulative distribution family to account for cumulative nature of growth-stage succession
  data = time_data_scaled, # use scaled numeric explanatory variable in model
  chains = 4, # use 4 independent fitting procedures
  iter = 10000, # use 10000 iterations
  warmup = 5000, # use 50% of iterations as warm-up
  seed = 2025 # use seed
)

# fit growth curve model with cumulative distribution family using bayesian regression (Treatment * days_scaled)
bay_model_no_random <- brm(
  Stages_Reached ~ Treatment * days_scaled,
  family = cumulative(), # cumulative distribution family to account for cumulative nature of growth-stage succession
  data = time_data_scaled, # use scaled numeric explanatory variable in model
  chains = 4, # use 4 independent fitting procedures
  iter = 10000, # use 10000 iterations
  warmup = 5000, # use 50% of iterations as warm-up
  seed = 2025 # use seed
)

# fit growth curve model with cumulative distribution family using bayesian regression (Treatment + days_scaled + random)
bay_model_sum <- brm(
  Stages_Reached ~ Treatment + days_scaled + (1 | ID), # plant ID as random factor to account for non-independent measurements
  family = cumulative(), # cumulative distribution family to account for cumulative nature of growth-stage succession
  data = time_data_scaled, # use scaled numeric explanatory variable in model
  chains = 4, # use 4 independent fitting procedures
  iter = 10000, # use 10000 iterations
  warmup = 5000, # use 50% of iterations as warm-up
  seed = 2025 # use seed
)

# fit growth curve model with cumulative distribution family using bayesian regression (Treatment + days_scaled)
bay_model_sum_no_random <- brm(
  Stages_Reached ~ Treatment + days_scaled,
  family = cumulative(), # cumulative distribution family to account for cumulative nature of growth-stage succession
  data = time_data_scaled, # use scaled numeric explanatory variable in model
  chains = 4, # use 4 independent fitting procedures
  iter = 10000, # use 10000 iterations
  warmup = 5000, # use 50% of iterations as warm-up
  seed = 2025 # use seed
)

# compare all models
loo(bay_model, bay_model_no_random, bay_model_sum, bay_model_sum_no_random) # lower looic is desired, higher elpd is desired
waic(bay_model, bay_model_no_random, bay_model_sum, bay_model_sum_no_random) # lower waic is desired, higher elpd is desired

# check model summary statistics and plot posterior distribution (Treatment * days_scaled + random)
summary(bay_model)
fixef(bay_model) # fixed effects
ranef(bay_model)
bayes_R2(bay_model)
plot(bay_model)
plot(conditional_effects(bay_model))

# check model summary statistics and plot posterior distribution (Treatment * days_scaled + random)
summary(bay_model_no_random, waic = TRUE)
fixef(bay_model_no_random) # fixed effects
bayes_R2(bay_model_no_random)
plot(bay_model_no_random)
plot(conditional_effects(bay_model_no_random))

# check model summary statistics and plot posterior distribution (Treatment + days_scaled + random)
summary(bay_model_sum)
fixef(bay_model_sum) # fixed effects
ranef(bay_model_sum)
bayes_R2(bay_model_sum)
plot(bay_model_sum)
plot(conditional_effects(bay_model_sum))

# check model summary statistics and plot posterior distribution (Treatment + days_scaled)
summary(bay_model_sum_no_random)
fixef(bay_model_sum_no_random) # fixed effects
bayes_R2(bay_model_sum_no_random)
plot(bay_model_sum_no_random)
plot(conditional_effects(bay_model_sum_no_random))


### plot predictions ###

# extract conditional effects data and convert to dataframe
ce_data <- conditional_effects(bay_model, effects = "days_scaled:Treatment")
ce_df <- as.data.frame(ce_data$`days_scaled:Treatment`)

# calculate scaling parameters for scaled days and scale back to original
days_mean <- mean(time_data_final$Days)
days_sd <- sd(time_data_final$Days)
ce_df$days_original <- ce_df$days_scaled * days_sd + days_mean

# reorder factor levels
ce_df$Treatment <- factor(
  ce_df$Treatment,
  levels = c("Control", "Ceno", "Pilo", "Co_Inoc")
)

# create growth model plot
model_plot <- ggplot(
  ce_df,
  aes(x = days_original, y = estimate__, color = Treatment)
) +
  geom_line(size = 1) +
  geom_ribbon(
    aes(ymin = lower__, ymax = upper__, fill = Treatment),
    alpha = 0.2,
    color = NA
  ) +
  geom_line(aes(y = lower__), linetype = "solid", size = 0.5) +
  geom_line(aes(y = upper__), linetype = "solid", size = 0.5) +
  scale_color_manual(
    values = custom_colors,
    labels = treatment_labels,
    name = "Treatment"
  ) + # apply custom color scheme
  scale_fill_manual(
    values = custom_colors,
    labels = treatment_labels,
    name = "Treatment"
  ) +
  scale_x_continuous(expand = c(0, 0), limits = c(0, 56)) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 7.25)) +
  labs(
    title = "c",
    y = "N Stages Reached",
    x = "N Days"
  ) +
  theme_pubr() +
  theme(
    axis.text.x = element_text(size = 10),
    axis.text.y = element_text(size = 10),
    axis.title = element_text(size = 12),
    legend.position = "none",
    plot.title = element_text(face = "bold", size = 14),
    plot.subtitle = element_text(size = 12, face = "italic"),
    panel.grid.major.y = element_line(color = "gray90"),
    panel.grid.minor.y = element_line(color = "gray95")
  )

# show plot
model_plot

# save results data
model_predictions <- subset(
  ce_df,
  select = -c(Stages_Reached, ID, cond__, effect1__, effect2__, days_scaled)
)
model_predictions <- model_predictions %>% rename(estimate = estimate__)
model_predictions <- model_predictions %>% rename(upper = upper__)
model_predictions <- model_predictions %>% rename(lower = lower__)
model_predictions <- model_predictions %>% rename(se = se__)
model_predictions <- model_predictions %>% rename(days = days_original)
model_predictions <- model_predictions %>% rename(treatment = Treatment)
write.csv(
  model_predictions,
  "results/Figure2/model_predictions_timeseries.csv",
  row.names = FALSE
)


### INVESTIGATE SIGNIFICANCE OF TREATMENT EFFECTS ###

# extract treatment effects
treatment_effects <- conditional_effects(bay_model, effects = "Treatment")
treatment_effect_data <- as.data.frame(treatment_effects$Treatment)

# reorder factor levels
treatment_effect_data$Treatment <- factor(
  treatment_effect_data$Treatment,
  levels = c("Control", "Ceno", "Pilo", "Co_Inoc")
)

# calculate grand mean
avg_data <- aggregate(
  estimate__ ~ days_scaled,
  data = treatment_effect_data,
  FUN = mean
)
names(avg_data)[2] <- "avg_estimate"

# merge average estimate with main data
treatment_effect_data <- merge(
  treatment_effect_data,
  avg_data,
  by = "days_scaled"
)

# calculate estimate, upper and lower relative to grand mean
treatment_effect_data$estimate__rel <- treatment_effect_data$estimate__ -
  treatment_effect_data$avg_estimate
treatment_effect_data$lower__rel <- treatment_effect_data$lower__ -
  treatment_effect_data$avg_estimate
treatment_effect_data$upper__rel <- treatment_effect_data$upper__ -
  treatment_effect_data$avg_estimate

# remove temporary column
treatment_effect_data$avg_estimate <- NULL

# create treatment order
treatment_effect_data$Treatment <- factor(
  treatment_effect_data$Treatment,
  levels = treatment_order
)

# get significance
treatment_effect_data$Significant <- (treatment_effect_data$lower__rel > 0 |
  treatment_effect_data$upper__rel < 0)

# create plot for treatment effects on stage development (mean and 95% credible interval)
overall_treatment_plot <- ggplot(
  treatment_effect_data,
  aes(x = Treatment, y = estimate__rel, color = Treatment)
) +
  geom_point(aes(shape = Significant), size = 3) +
  geom_linerange(
    aes(color = Treatment, ymin = lower__rel, ymax = upper__rel),
    size = 1.2
  ) +
  scale_color_manual(
    values = custom_colors,
    labels = treatment_labels,
    guide = "none" # Remove the color legend
  ) +
  scale_x_discrete(labels = treatment_labels) +
  geom_hline(yintercept = 0, size = 0.2, col = "black") +
  labs(
    title = "d",
    y = "Relative Response",
    x = "Treatment"
  ) +
  theme_pubr() +
  theme(
    axis.text.x = element_text(angle = 25, hjust = 1, size = 10),
    axis.text.y = element_text(size = 10),
    axis.title = element_text(size = 12),
    legend.position = "right",
    plot.title = element_text(face = "bold", size = 14),
    plot.subtitle = element_text(size = 12, face = "italic"),
    panel.grid.major.y = element_line(color = "gray90"),
    panel.grid.minor.y = element_line(color = "gray95")
  )

# show plot
overall_treatment_plot

# save results data
overall_effects <- subset(
  treatment_effect_data,
  select = -c(Stages_Reached, ID, cond__, effect1__, days_scaled)
)
overall_effects <- overall_effects %>% rename(se = se__)
overall_effects <- overall_effects %>% rename(estimate_absolute = estimate__)
overall_effects <- overall_effects %>% rename(lower_absolute = lower__)
overall_effects <- overall_effects %>% rename(upper_absolute = upper__)
overall_effects <- overall_effects %>% rename(estimate_relative = estimate__rel)
overall_effects <- overall_effects %>% rename(lower_relative = lower__rel)
overall_effects <- overall_effects %>% rename(upper_relative = upper__rel)
write.csv(
  overall_effects,
  "results/Figure2/overall_treatment_effects.csv",
  row.names = FALSE
)


### IDENTIFY DAYS OF SIGNIFICANT DIFFERENCE ###

# extract unique treatments from original data
treatments <- unique(bay_model$data$Treatment)

# identify treatment pairs
treatment_pairs <- combn(treatments, 2, simplify = FALSE)

# create data frame for significantly different days
significance_days <- data.frame(
  Treatment1 = character(),
  Treatment2 = character(),
  First_Significant_Day = numeric(),
  stringsAsFactors = FALSE
)

# extract data
for (pair in treatment_pairs) {
  t1 <- pair[1]
  t2 <- pair[2]
  first_sig_day <- NA
  # check each day
  for (day in 1:56) {
    # revert scaling of days
    if (is.null(days_mean) || is.null(days_sd)) {
      day_scaled <- day
    } else {
      day_scaled <- (day - days_mean) / days_sd
    }
    # create prediction data for treatment pairs
    pred_data <- data.frame(
      Treatment = c(t1, t2),
      days_scaled = rep(day_scaled, 2),
      ID = rep("new_subject", 2)
    )
    # get posterior predictions
    preds <- posterior_epred(
      bay_model,
      newdata = pred_data,
      allow_new_levels = TRUE
    )
    # compute expected stages for each treatment
    expected_stages <- array(0, dim = c(dim(preds)[1], 2))
    for (i in 1:dim(preds)[3]) {
      # for each outcome category/stage
      expected_stages[, 1] <- expected_stages[, 1] + i * preds[, 1, i]
      expected_stages[, 2] <- expected_stages[, 2] + i * preds[, 2, i]
    }
    # compute differences and check significance
    stage_diffs <- expected_stages[, 1] - expected_stages[, 2]
    # check if 95% CI excludes zero
    ci_lower <- quantile(stage_diffs, 0.025)
    ci_upper <- quantile(stage_diffs, 0.975)
    is_significant <- (ci_lower > 0) | (ci_upper < 0)
    # record first significant day
    if (is_significant && is.na(first_sig_day)) {
      first_sig_day <- day
      break
    }
  }
  # add results to data frame
  significance_days <- rbind(
    significance_days,
    data.frame(
      Treatment1 = t1,
      Treatment2 = t2,
      First_Significant_Day = first_sig_day
    )
  )
}

# check results
print(significance_days)

# create data frame to store results
diff_over_time <- data.frame(
  Treatment1 = character(),
  Treatment2 = character(),
  Day = numeric(),
  Mean_Diff = numeric(),
  Lower_CI = numeric(),
  Upper_CI = numeric(),
  Significant = logical(),
  stringsAsFactors = FALSE
)

# for each pair, calculate differences over time
for (pair in treatment_pairs) {
  t1 <- pair[1]
  t2 <- pair[2]
  for (day in 1:56) {
    # scale day parameter
    if (is.null(days_mean) || is.null(days_sd)) {
      day_scaled <- day
    } else {
      day_scaled <- (day - days_mean) / days_sd
    }
    # create prediction data for treatment paris
    pred_data <- data.frame(
      Treatment = c(t1, t2),
      days_scaled = rep(day_scaled, 2),
      ID = rep("new_subject", 2)
    )
    # get posterior predictions
    preds <- posterior_epred(
      bay_model,
      newdata = pred_data,
      allow_new_levels = TRUE
    )
    # compute expected stages for each treatment
    expected_stages <- array(0, dim = c(dim(preds)[1], 2))
    for (i in 1:dim(preds)[3]) {
      # for each outcome category/stage
      expected_stages[, 1] <- expected_stages[, 1] + i * preds[, 1, i]
      expected_stages[, 2] <- expected_stages[, 2] + i * preds[, 2, i]
    }
    # compute differences and check significance
    stage_diffs <- expected_stages[, 1] - expected_stages[, 2]
    # calculate mean and CI
    mean_diff <- mean(stage_diffs)
    ci_lower <- quantile(stage_diffs, 0.025)
    ci_upper <- quantile(stage_diffs, 0.975)
    is_significant <- (ci_lower > 0) | (ci_upper < 0)
    # add to results data frame
    diff_over_time <- rbind(
      diff_over_time,
      data.frame(
        Treatment1 = t1,
        Treatment2 = t2,
        Day = day,
        Mean_Diff = mean_diff,
        Lower_CI = ci_lower,
        Upper_CI = ci_upper,
        Significant = is_significant
      )
    )
  }
}

# define significant treatment pairs for final plot
treatment_pairs_final <- list(
  factor(
    c(2, 4),
    levels = 1:4,
    labels = c("Control", "Ceno", "Co_Inoc", "Pilo")
  ),
  factor(
    c(2, 3),
    levels = 1:4,
    labels = c("Control", "Ceno", "Co_Inoc", "Pilo")
  ),
  factor(
    c(2, 1),
    levels = 1:4,
    labels = c("Control", "Ceno", "Co_Inoc", "Pilo")
  ),
  factor(
    c(3, 1),
    levels = 1:4,
    labels = c("Control", "Ceno", "Co_Inoc", "Pilo")
  ),
  factor(
    c(3, 4),
    levels = 1:4,
    labels = c("Control", "Ceno", "Co_Inoc", "Pilo")
  )
)

# create data frame to store results for final plot
diff_over_time_final <- data.frame(
  Treatment1 = character(),
  Treatment2 = character(),
  Day = numeric(),
  Mean_Diff = numeric(),
  Lower_CI = numeric(),
  Upper_CI = numeric(),
  Significant = logical(),
  stringsAsFactors = FALSE
)

# for each final pair, calculate differences over time
for (pair in treatment_pairs_final) {
  t1 <- pair[1]
  t2 <- pair[2]
  for (day in 1:56) {
    # scale day parameter
    if (is.null(days_mean) || is.null(days_sd)) {
      day_scaled <- day
    } else {
      day_scaled <- (day - days_mean) / days_sd
    }
    # create prediction data for treatment paris
    pred_data <- data.frame(
      Treatment = c(t1, t2),
      days_scaled = rep(day_scaled, 2),
      ID = rep("new_subject", 2)
    )
    # get posterior predictions
    preds <- posterior_epred(
      bay_model,
      newdata = pred_data,
      allow_new_levels = TRUE
    )
    # compute expected stages for each treatment
    expected_stages <- array(0, dim = c(dim(preds)[1], 2))
    for (i in 1:dim(preds)[3]) {
      # for each outcome category/stage
      expected_stages[, 1] <- expected_stages[, 1] + i * preds[, 1, i]
      expected_stages[, 2] <- expected_stages[, 2] + i * preds[, 2, i]
    }
    # compute differences and check significance
    stage_diffs <- expected_stages[, 1] - expected_stages[, 2]
    # calculate mean and CI
    mean_diff <- mean(stage_diffs)
    ci_lower <- quantile(stage_diffs, 0.025)
    ci_upper <- quantile(stage_diffs, 0.975)
    is_significant <- (ci_lower > 0) | (ci_upper < 0)
    # add to results data frame
    diff_over_time_final <- rbind(
      diff_over_time_final,
      data.frame(
        Treatment1 = t1,
        Treatment2 = t2,
        Day = day,
        Mean_Diff = mean_diff,
        Lower_CI = ci_lower,
        Upper_CI = ci_upper,
        Significant = is_significant
      )
    )
  }
}

# save pairwise timeseries credible intervals
diff_over_time <- diff_over_time %>%
  dplyr::select(
    Treatment1,
    Treatment2,
    Day,
    Lower_CI,
    Upper_CI,
    Significant
  ) %>%
  arrange(Treatment1, Treatment2, Day)
write.csv(
  diff_over_time,
  "results/Figure2/pairwise_treatment_differences_over_time.csv",
  row.names = FALSE
)

# create final pairwise differences plot
pairwise_plot_final <- ggplot(
  diff_over_time_final,
  aes(
    x = Day,
    y = Mean_Diff,
    ymin = Lower_CI,
    ymax = Upper_CI,
    color = Significant
  )
) +
  geom_ribbon(aes(fill = Significant), alpha = 0.5, color = NA) +
  geom_line() +
  geom_hline(yintercept = 0, linetype = "dashed") +
  facet_wrap(
    ~ paste(Treatment1, "vs", Treatment2),
    labeller = labeller(
      .default = label_value,
      `paste(Treatment1, "vs", Treatment2)` = function(x) {
        x <- gsub(
          "Ceno vs Co_Inoc",
          "*Cenococcum geophilum* vs.<br>Co-Inoculation",
          x
        )
        x <- gsub(
          "Ceno vs Pilo",
          "*Cenococcum geophilum* vs.<br>*Piloderma croceum*",
          x
        )
        x <- gsub(
          "Ceno vs Control",
          "*Cenococcum geophilum* vs. <br>Control",
          x
        )
        x <- gsub("Co_Inoc vs Control", "Co-Inoculation vs.<br>Control", x)
        x <- gsub(
          "Co_Inoc vs Pilo",
          "Co-Inoculation vs.<br>*Piloderma croceum*",
          x
        )
      }
    )
  ) +
  scale_color_manual(values = c("gray60", "gray10")) +
  scale_fill_manual(values = c("gray60", "gray10")) +
  theme_minimal() +
  labs(
    title = "e",
    x = "N Days",
    y = "Expected Stage (Treatment 1 - Treatment 2)"
  ) +
  theme_pubr() +
  theme(
    axis.text.x = element_text(size = 10),
    axis.text.y = element_text(size = 10),
    axis.title = element_text(size = 12),
    legend.position = "right",
    plot.title = element_text(face = "bold", size = 14),
    plot.subtitle = element_text(size = 12, face = "italic"),
    panel.grid.major.y = element_line(color = "gray90"),
    panel.grid.minor.y = element_line(color = "gray95"),
    strip.background = element_rect(fill = "white"),
    strip.text = element_markdown(size = 10)
  )

# show final plot
pairwise_plot_final

# save significant day results
write.csv(
  significance_days,
  "results/Figure2/significance_days.csv",
  row.names = FALSE
)


###### CHECK FOR DIFFERENCES IN STAGE DURATIONS FOR EACH STAGE BETWEEN TREATMENTS #####

### CALCULATE STAGE DURATIONS ###

# write function for stage duration calculation
calculate_stage_duration <- function(df) {
  # convert date data
  df$Date <- as.Date(as.character(df$Date), format = "%Y%m%d")
  # sort data by ID and Date
  df <- df[order(df$ID, df$Date), ]
  # create results dataframe
  result <- data.frame(
    ID = integer(),
    Treatment = character(),
    Stage = character(),
    StartDate = as.Date(character()),
    EndDate = as.Date(character()),
    DaysInStage = integer(),
    stringsAsFactors = FALSE
  )
  # analyze each ID separately
  for (id in unique(df$ID)) {
    id_data <- df[df$ID == id, ]
    if (nrow(id_data) <= 1) next
    # track stages
    current_stage <- NULL
    stage_start_date <- NULL
    for (i in 1:nrow(id_data)) {
      # check if new stage was reached
      if (is.null(current_stage) || id_data$Stage[i] != current_stage) {
        # record duration for previous stage
        if (!is.null(current_stage)) {
          result <- rbind(
            result,
            data.frame(
              ID = id,
              Treatment = id_data$Treatment[i - 1],
              Stage = current_stage,
              StartDate = stage_start_date,
              EndDate = id_data$Date[i],
              DaysInStage = as.integer(id_data$Date[i] - stage_start_date),
              stringsAsFactors = FALSE
            )
          )
        }
        # start tracking new stage
        current_stage <- id_data$Stage[i]
        stage_start_date <- id_data$Date[i]
      }
    }
  }
  return(result)
}

# calculate stage durations
stage_durations <- calculate_stage_duration(time_data_final)
print(stage_durations)

# add stage letter for plotting
stage_durations <- stage_durations %>%
  mutate(StageLetter = substr(Stage, nchar(Stage), nchar(Stage)))


### COMPARISONS BETWEEN TREATMENTS FOR EACH STAGE ###

# subset stage durations for specific stages
stage_durations_A <- subset(stage_durations, grepl("A$", Stage))
stage_durations_B <- subset(stage_durations, grepl("B$", Stage))
stage_durations_C <- subset(stage_durations, grepl("C$", Stage))
stage_durations_D <- subset(stage_durations, grepl("D$", Stage))

# add stage letters for stat test
stage_durations_A$stage_letter <- "A"
stage_durations_B$stage_letter <- "B"
stage_durations_C$stage_letter <- "C"
stage_durations_D$stage_letter <- "D"

# do wilcoxon test with bonferroni p adjustment for each stage
stat_test_A <- stage_durations_A %>%
  wilcox_test(DaysInStage ~ Treatment) %>%
  adjust_pvalue(method = "bonferroni") %>%
  add_significance()
stat_test_A$stage <- "A"
stat_test_B <- stage_durations_B %>%
  wilcox_test(DaysInStage ~ Treatment) %>%
  adjust_pvalue(method = "bonferroni") %>%
  add_significance()
stat_test_B$stage <- "B"
stat_test_C <- stage_durations_C %>%
  wilcox_test(DaysInStage ~ Treatment) %>%
  adjust_pvalue(method = "bonferroni") %>%
  add_significance()
stat_test_C$stage <- "C"
stat_test_D <- stage_durations_D %>%
  wilcox_test(DaysInStage ~ Treatment) %>%
  adjust_pvalue(method = "bonferroni") %>%
  add_significance()
stat_test_D$stage <- "D"

# combine results
wilcoxon_results_durations <- rbind(
  stat_test_A,
  stat_test_B,
  stat_test_C,
  stat_test_D
)
print(wilcoxon_results_durations)

# save_results
stage_duration_results <- subset(wilcoxon_results_durations, select = -c(.y.))
write.csv(
  stage_duration_results,
  "results/Figure2/stage_duration_wilcoxon.csv",
  row.names = FALSE
)


### PLOT RESULTS ###

# set treatment order
stage_durations$Treatment <- factor(
  stage_durations$Treatment,
  levels = treatment_order
)

# create plot
stage_plot <- ggplot(
  stage_durations,
  aes(x = StageLetter, y = DaysInStage, fill = Treatment)
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
  labs(y = "Days in Stage", x = "Stage") +
  theme_classic() +
  geom_signif(
    y_position = c(64, 59, 24, 30),
    xmin = c(0.72, 0.912, 2.907, 3.1),
    xmax = c(1.095, 1.095, 3.1, 3.278),
    annotation = c("0.024", "0.001", "0.048", "0.012")
  ) +
  labs(
    title = "f    Treatment-Specifc Stage Durations",
    subtitle = "        Only significant p-values are shown (Wilcoxon test with Bonferroni correction) ",
    x = "Stage",
    y = "N Days in Stage"
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
stage_plot


##### FINAL VISUALIZATION #####

# create combined plot
combined_plot <- (plot_spacer() | basic_growth_plot) /
  (model_plot | overall_treatment_plot) /
  pairwise_plot_final /
  stage_plot +
  plot_layout() +
  plot_annotation(
    theme = theme(
      plot.title = element_text(size = 16, face = "bold", hjust = 0.5)
    )
  )

# show plot
combined_plot

# save combined plot
ggsave(
  "results/Figure2/Figure2.pdf",
  plot = combined_plot,
  width = 12.5,
  height = 20
)
dev.off()
