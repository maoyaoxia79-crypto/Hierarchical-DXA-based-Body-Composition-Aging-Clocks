# ==============================================================================
# Demo and Simulation Script: DXA-based Body Composition Aging Clocks
# ==============================================================================
# This script generates synthetic (mock) data to simulate the UK Biobank cohort
# characteristics and runs the biological clock building pipeline.
#
# It demonstrates that the clock construction algorithms work out-of-the-box.
#
# Requirements:
#   R libraries: survival, glmnet, caret, flexsurv, ggplot2, Metrics, ggExtra
# ==============================================================================

# 1. Load the core functions
source("clock_construction_functions.R")

cat("--- Preparing Directories ---\n")
if (!dir.exists("output")) dir.create("output")

# 2. Generate Synthetic (Mock) Data
# In actual analysis, this corresponds to UK Biobank derived clinical and DXA dataset
set.seed(42)
n_samples <- 1500

cat("--- Generating Synthetic Cohort (n =", n_samples, ") ---\n")

# Define clock prefixes (representing muscle, fat, and skeletal domains)
clock_prefixes <- c(
  "Central.Fat.Age", "Peripheral.Fat.Age", "Total.Fat.Age",
  "Axial.Muscle.Age", "Peripheral.Muscle.Age", "Total.Muscle.Age",
  "Axial.Skeletal.Age", "Peripheral.Skeletal.Age", "Total.Skeletal.Age",
  "Body.Composition.Age"
)

# For each prefix, we simulate 5 standardized biomarkers (features)
all_exposure_vars <- unlist(lapply(clock_prefixes, function(p) paste0(p, "_", 1:5, "_std")))

# Base clinical cohort data
mock_data <- data.frame(
  id = 1:n_samples,
  sex = sample(0:1, n_samples, replace = TRUE), # 0 = Female, 1 = Male
  age_at_ins2_std = rnorm(n_samples, 0, 1), # Standardized chronological age
  death = sample(0:1, n_samples, replace = TRUE, prob = c(0.88, 0.12)), # Event status
  death.ins2py = runif(n_samples, 1, 15), # Follow-up time (years)
  split_group = sample(c("Training", "Validation", "Disease.Test"), 
                       n_samples, replace = TRUE, prob = c(0.6, 0.2, 0.2))
)

# Add mock standardized DXA biomarkers (drawn from normal distribution)
for (var in all_exposure_vars) {
  mock_data[[var]] <- rnorm(n_samples, 0, 1)
}

# Simulate correlation/collinearity between fat markers and muscle markers
mock_data$Peripheral.Muscle.Age_2_std <- mock_data$Peripheral.Muscle.Age_1_std * 0.95 + rnorm(n_samples, 0, 0.05)
mock_data$Central.Fat.Age_2_std <- mock_data$Central.Fat.Age_1_std * 0.92 + rnorm(n_samples, 0, 0.08)

# 3. Perform Data Splitting
train_data_full <- mock_data[mock_data$split_group == "Training", ]
test_data_full <- mock_data[mock_data$split_group == "Validation", ]
disease_test_data_full <- mock_data[mock_data$split_group == "Disease.Test", ]

train_data_male <- train_data_full[train_data_full$sex == 1, ]
train_data_female <- train_data_full[train_data_full$sex == 0, ]
test_data_male <- test_data_full[test_data_full$sex == 1, ]
test_data_female <- test_data_full[test_data_full$sex == 0, ]
disease_test_data_male <- disease_test_data_full[disease_test_data_full$sex == 1, ]
disease_test_data_female <- disease_test_data_full[disease_test_data_full$sex == 0, ]

cat("Cohort Split Summary:\n")
cat(" - Male:   Train =", nrow(train_data_male), ", Val =", nrow(test_data_male), ", Disease =", nrow(disease_test_data_male), "\n")
cat(" - Female: Train =", nrow(train_data_female), ", Val =", nrow(test_data_female), ", Disease =", nrow(disease_test_data_female), "\n")

# 4. Define Clocks and Pipeline Configurations to Run
clocks_to_build <- list(
  "Central.Fat" = "Central.Fat.Age",
  "Peripheral.Fat" = "Peripheral.Fat.Age",
  "Total.Fat" = "Total.Fat.Age",
  "Axial.Muscle" = "Axial.Muscle.Age",
  "Peripheral.Muscle" = "Peripheral.Muscle.Age",
  "Total.Muscle" = "Total.Muscle.Age",
  "Axial.Skeletal" = "Axial.Skeletal.Age",
  "Peripheral.Skeletal" = "Peripheral.Skeletal.Age",
  "Total.Skeletal" = "Total.Skeletal.Age",
  "Body.Composition" = "Body.Composition.Age"
)

# Run 2 key flows for demonstration:
# Flow 1: No LASSO, Filter Collinearity (Collinear variables removed)
# Flow 2: Staged LASSO (LASSO variables selection, Age is NOT penalized)
analysis_flows <- list(
  "No_LASSO_Filtered" = list(params = list(use_lasso = FALSE, filter_collinearity = TRUE), suffix = "_no_LASSO_filtered"),
  "Staged_LASSO" = list(params = list(use_lasso = TRUE, age_in_lasso = FALSE, age_in_gompertz = TRUE), suffix = "_staged_LASSO")
)

all_performance_results <- list()

# 5. Run the Clock Building Loop
for (flow_name in names(analysis_flows)) {
  flow <- analysis_flows[[flow_name]]
  cat(paste("\n=============================================\n"))
  cat(paste(" RUNNING FLOW:", flow_name, "\n"))
  cat(paste("=============================================\n"))
  
  for (clock_short_name in names(clocks_to_build)) {
    prefix <- clocks_to_build[[clock_short_name]]
    dir_name <- file.path("output", paste0(gsub("\\.", "_", clock_short_name), "_Age_Results", flow$suffix))
    
    # Run Male model
    male_args <- c(list(train_data_male, test_data_male, disease_test_data_male, "Male", prefix, dir_name), flow$params)
    male_results <- do.call(build_aging_clock_main, male_args)
    
    # Run Female model
    female_args <- c(list(train_data_female, test_data_female, disease_test_data_female, "Female", prefix, dir_name), flow$params)
    female_results <- do.call(build_aging_clock_main, female_args)
    
    # Save model results locally
    save(male_results, female_results, file = file.path(dir_name, paste0(clock_short_name, "_Age_Final_Data.RData")))
    
    # Gather performance metrics
    for (dataset_type in names(male_results$performance_metrics)) {
      all_performance_results[[length(all_performance_results) + 1]] <- 
        data.frame(Flow = flow_name, Clock = clock_short_name, Sex = "Male", Dataset = dataset_type, 
                   as.data.frame(male_results$performance_metrics[[dataset_type]]))
      
      all_performance_results[[length(all_performance_results) + 1]] <- 
        data.frame(Flow = flow_name, Clock = clock_short_name, Sex = "Female", Dataset = dataset_type, 
                   as.data.frame(female_results$performance_metrics[[dataset_type]]))
    }
  }
}

# 6. Summarize & Save Performance Metrics
cat("\n=============================================\n")
cat(" SUMMARIZING RESULTS\n")
cat("=============================================\n")

performance_summary_df <- do.call(rbind, all_performance_results)
write.csv(performance_summary_df, file.path("output", "All_Models_Performance_Summary.csv"), row.names = FALSE)
cat("All model performance metrics saved to output/All_Models_Performance_Summary.csv\n\n")

# Display a subset of the results (Combined correlation)
cat("Sample Results (Validation Set Correlation):\n")
print(head(performance_summary_df[performance_summary_df$Dataset == "Validation", c("Flow", "Clock", "Sex", "r", "mae")], 10))

# 7. Generate Performance Visualization Plots
cat("\nGenerating Performance Plots...\n")
plot_data <- performance_summary_df %>%
  mutate(Sex = factor(Sex, levels = c("Female", "Male")),
         Dataset = factor(Dataset, levels = c("Training", "Validation", "Disease.Test")))

for (flow_name in names(analysis_flows)) {
  flow_plot_data <- filter(plot_data, Flow == flow_name)
  
  # Correlation (r) plot
  p_corr <- ggplot(flow_plot_data, aes(x = Clock, y = r, fill = Sex)) +
    geom_bar(stat = "identity", position = position_dodge(width = 0.9), width = 0.8) +
    facet_wrap(~Dataset, ncol = 3) +
    scale_fill_manual(values = c("Female" = "#F28080", "Male" = "#6B7EB9")) +
    labs(
      title = paste("Model Performance (Correlation) - Flow:", flow_name),
      x = "",
      y = "Correlation (r)",
      fill = "Sex"
    ) +
    theme_bw(base_size = 12) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 9),
      panel.grid.major.x = element_blank(),
      panel.grid.minor.y = element_blank(),
      legend.position = "top"
    )
  
  # MAE plot
  p_mae <- ggplot(flow_plot_data, aes(x = Clock, y = mae, fill = Sex)) +
    geom_bar(stat = "identity", position = position_dodge(width = 0.9), width = 0.8) +
    facet_wrap(~Dataset, ncol = 3) +
    scale_fill_manual(values = c("Female" = "#F28080", "Male" = "#6B7EB9")) +
    labs(
      title = paste("Model Performance (MAE) - Flow:", flow_name),
      x = "Aging Clock",
      y = "MAE (years)",
      fill = "Sex"
    ) +
    theme_bw(base_size = 12) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 9),
      panel.grid.major.x = element_blank(),
      panel.grid.minor.y = element_blank(),
      legend.position = "top"
    )
  
  # Save plots as PDF
  pdf_path <- file.path("output", paste0("Summary_Performance_Plots_", flow_name, ".pdf"))
  pdf(pdf_path, width = 12, height = 8)
  print(p_corr)
  print(p_mae)
  dev.off()
  cat(paste("Summary plots saved to", pdf_path, "\n"))
}

cat("\n=============================================\n")
cat(" DEMO RUN SUCCESSFULLY COMPLETED!\n")
cat("=============================================\n")
