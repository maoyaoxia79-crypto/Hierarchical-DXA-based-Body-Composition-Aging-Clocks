#' Modular Functions for DXA-derived Biological Age Clock Construction
#'
#' This script contains the primary utility and model building functions to construct
#' sex-specific Gompertz-Cox based biological age clocks, filter collinear features,
#' perform LASSO variable selection, evaluate clock performance, and generate plots.
#' 
#' Author: maoyaoxia79-crypto
#' Year: 2026
#' License: MIT

library(survival)
library(glmnet)
library(caret)
library(ggplot2)
library(Metrics)
library(flexsurv)
library(ggExtra)
library(tidyr)
library(dplyr)

#' Calculate Model Performance Metrics
#'
#' @param actual Numeric vector of chronological age.
#' @param predicted Numeric vector of predicted biological age.
#' @return A list containing Pearson correlation (r), R-squared (r_squared), 
#' Mean Absolute Error (mae), Root Mean Squared Error (rmse), and R2 p-value.
calculate_metrics <- function(actual, predicted) {
  valid_indices <- !is.na(actual) & !is.na(predicted)
  actual <- actual[valid_indices]
  predicted <- predicted[valid_indices]
  
  if (length(actual) < 2) {
    return(list(r = NA, r_squared = NA, mae = NA, rmse = NA, p_value_R2 = NA))
  }
  
  model <- lm(actual ~ predicted)
  model_summary <- summary(model)
  r_squared <- model_summary$r.squared
  p_value_R2 <- if (nrow(model_summary$coefficients) > 1) model_summary$coefficients[2, 4] else NA
  r <- cor(actual, predicted)
  mae <- mean(abs(actual - predicted))
  rmse <- sqrt(mean((actual - predicted)^2))
  
  return(list(r = r, r_squared = r_squared, mae = mae, rmse = rmse, p_value_R2 = p_value_R2))
}

#' Plot Chronological Age vs. Biological Age and Save as PDF
#'
#' @param data Data frame containing chronological age (`age`) and biological age (`Bio_Age`).
#' @param dataset_type Character string indicating the split name (e.g. "Training", "Validation", "Disease.Test").
#' @param clock_name_prefix Character string indicating the name of the clock.
#' @param output_dir Character string indicating the output directory path.
#' @param sex_label Character string for the cohort sex ("Male" or "Female").
#' @param plot_color HEX string for plot coloring.
#' @return A list containing calculated metrics.
plot_age_correlation <- function(data, dataset_type, clock_name_prefix, output_dir, sex_label, plot_color) {
  metrics <- calculate_metrics(data$age, data$Bio_Age)
  file_path <- file.path(output_dir, paste0(sex_label, "_", gsub("\\.", "_", dataset_type), "_Plot.pdf"))
  
  p <- ggplot(data, aes(x = age, y = Bio_Age)) +
    geom_point(color = plot_color, size = 2, alpha = 0.6) +
    stat_smooth(method = "lm", color = "gray30", linewidth = 1.2) +
    labs(title = paste(sex_label, "-", dataset_type), 
         x = "Chronological Age (years)", 
         y = paste0(clock_name_prefix, " (years)")) +
    theme_bw(base_size = 14) +
    theme(aspect.ratio = 1, 
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(), 
          panel.border = element_rect(linewidth = 1.5, color = "black", fill = NA), 
          plot.title = element_text(hjust = 0.5, face = "bold")) +
    annotate("text", x = min(data$age, na.rm = TRUE) + 1, y = max(data$Bio_Age, na.rm = TRUE) - 1, 
             label = paste0("r=", round(metrics$r, 3), 
                            "\nR²=", round(metrics$r_squared, 3), 
                            "\nMAE=", round(metrics$mae, 3), 
                            "\nRMSE=", round(metrics$rmse, 3)),
             hjust = 0, vjust = 1, size = 5, fontface = "bold", color = "gray30") +
    coord_fixed(ratio = 1, xlim = range(data$age, na.rm = TRUE), ylim = range(data$Bio_Age, na.rm = TRUE))
  
  # Try to add marginal histograms
  p_marginal <- tryCatch({
    ggExtra::ggMarginal(p, type = "histogram", fill = plot_color, alpha = 0.7)
  }, error = function(e) {
    p # Fallback to standard plot if ggMarginal fails
  })
  
  ggsave(file_path, plot = p_marginal, width = 7, height = 7, device = "pdf")
  cat(paste0(sex_label, "-", dataset_type, " Metrics: r=", round(metrics$r, 4), 
             ", R²=", round(metrics$r_squared, 4), ", MAE=", round(metrics$mae, 4), "\n"))
  
  return(metrics)
}

#' Run Complete Evaluation and Plotting Pipeline
#'
#' @param train_data Training dataset.
#' @param test_data Validation dataset.
#' @param disease_test_data Disease test dataset.
#' @param clock_name_prefix Character string indicating the name of the clock.
#' @param output_dir Character string indicating output directory path.
#' @param sex_label Character string for the cohort sex ("Male" or "Female").
#' @param age_mean Mean value used for standardizing age.
#' @param age_sd Standard deviation used for standardizing age.
#' @return A list containing the combined dataset with predictions and the metrics across sets.
run_evaluation_and_plotting <- function(train_data, test_data, disease_test_data, clock_name_prefix, output_dir, sex_label, age_mean = 64.9, age_sd = 7.84) {
  plot_color <- ifelse(sex_label == "Male", "#6B7EB9", "#F28080")
  metrics_list <- list()
  
  metrics_list[["Training"]] <- plot_age_correlation(train_data, "Training", clock_name_prefix, output_dir, sex_label, plot_color)
  metrics_list[["Validation"]] <- plot_age_correlation(test_data, "Validation", clock_name_prefix, output_dir, sex_label, plot_color)
  metrics_list[["Disease.Test"]] <- plot_age_correlation(disease_test_data, "Disease.Test", clock_name_prefix, output_dir, sex_label, plot_color)
  
  combined_data <- rbind(train_data, test_data, disease_test_data)
  names(combined_data)[names(combined_data) == "Bio_Age"] <- clock_name_prefix
  
  return(list(final_data = combined_data, performance_metrics = metrics_list))
}

#' Build Gompertz-Cox Biological Aging Clock
#'
#' @param train_data Training dataset.
#' @param test_data Validation dataset.
#' @param disease_test_data Disease test dataset.
#' @param sex_label Character string for the cohort sex ("Male" or "Female").
#' @param clock_name_prefix Name of the biological age clock.
#' @param output_dir_base Target output directory.
#' @param use_lasso Logical. Whether to use LASSO selection.
#' @param filter_collinearity Logical. Whether to filter out collinear variables.
#' @param age_in_lasso Logical. Whether to include chronological age during LASSO selection.
#' @param age_in_gompertz Logical. Whether to include chronological age in the Gompertz regression model.
#' @param time_variable Character. Name of the survival time column.
#' @param status_variable Character. Name of the survival status column.
#' @return A list containing final datasets and evaluation metrics.
build_aging_clock_main <- function(train_data, test_data, disease_test_data, sex_label, clock_name_prefix, output_dir_base,
                                   use_lasso = TRUE, filter_collinearity = FALSE, age_in_lasso = TRUE, age_in_gompertz = TRUE,
                                   time_variable = "death.ins2py", status_variable = "death") {
  
  # --- 1. Initialization and Denormalization of Age ---
  train_data <- as.data.frame(train_data)
  test_data <- as.data.frame(test_data)
  disease_test_data <- as.data.frame(disease_test_data)
  
  # Scale parameters (based on UKB cohort characteristics)
  age_mean <- 64.9
  age_sd <- 7.84
  
  train_data$age <- train_data$age_at_ins2_std * age_sd + age_mean
  test_data$age <- test_data$age_at_ins2_std * age_sd + age_mean
  disease_test_data$age <- disease_test_data$age_at_ins2_std * age_sd + age_mean
  
  cat(paste("\n\n====== Start Building Clock:", clock_name_prefix, "for", sex_label, " ======\n"))
  if (!dir.exists(output_dir_base)) dir.create(output_dir_base, recursive = TRUE)
  
  # --- 2. Feature Identification ---
  pattern <- switch(clock_name_prefix,
                    "Total.Fat.Age" = "^(Central\\.Fat|Peripheral\\.Fat|Total\\.Fat)\\.Age.*_std$",
                    "Total.Muscle.Age" = "^(Axial\\.Muscle|Peripheral\\.Muscle|Total\\.Muscle)\\.Age.*_std$",
                    "Total.Skeletal.Age" = "^(Axial\\.Skeletal|Peripheral\\.Skeletal|Total\\.Skeletal)\\.Age.*_std$",
                    "Body.Composition.Age" = "_std$",
                    paste0("^", clock_name_prefix, ".*_std$"))
  
  all_std_vars <- grep(pattern, names(train_data), value = TRUE)
  features_base <- setdiff(all_std_vars, "age_at_ins2_std")
  
  # --- 3. Multicollinearity Filtering ---
  model_vars_base <- features_base 
  
  if (filter_collinearity) {
    cat("--- Step: Filtering out highly collinear variables (Threshold = 0.9) ---\n")
    if (length(features_base) > 1) {
      cor_matrix <- cor(train_data[, features_base], use = "pairwise.complete.obs")
      highly_correlated <- caret::findCorrelation(cor_matrix, cutoff = 0.9, verbose = FALSE)
      if (length(highly_correlated) > 0) {
        vars_to_remove <- features_base[highly_correlated]
        cat("Removed due to collinearity:", paste(vars_to_remove, collapse = ", "), "\n")
        model_vars_base <- features_base[-highly_correlated]
      } else {
        cat("No variables removed due to collinearity.\n")
      }
    }
  }
  
  # --- 4. LASSO Feature Selection ---
  if (use_lasso) {
    cat("--- Step: Feature selection using LASSO-Cox ---\n")
    lasso_features <- if (age_in_lasso) c(model_vars_base, "age_at_ins2_std") else model_vars_base
    
    if (length(lasso_features) > 0) {
      # Remove any rows with missing survival outcomes
      valid_train_idx <- !is.na(train_data[[time_variable]]) & !is.na(train_data[[status_variable]])
      X_train <- as.matrix(train_data[valid_train_idx, lasso_features])
      y_train <- Surv(train_data[[time_variable]][valid_train_idx], train_data[[status_variable]][valid_train_idx])
      
      set.seed(123)
      cv.lasso <- cv.glmnet(X_train, y_train, family = "cox", alpha = 1, nfolds = 5)
      final_model <- glmnet(X_train, y_train, family = "cox", alpha = 1, lambda = cv.lasso$lambda.min)
      coef_lasso <- coef(final_model)
      selected_coefs <- data.frame(Feature = rownames(coef_lasso), Coefficient = as.numeric(coef_lasso[, 1]))
      selected_features <- selected_coefs[selected_coefs$Coefficient != 0, ]
      
      cat(sex_label, "LASSO Selected Features:\n"); print(selected_features)
      write.table(selected_features, file = file.path(output_dir_base, paste0(sex_label, "_LASSO_Features.txt")), 
                  sep = "\t", quote = FALSE, row.names = FALSE)
      
      model_vars_base <- selected_features$Feature[selected_features$Feature != "age_at_ins2_std"]
    } else {
      model_vars_base <- character(0) 
    }
  }
  
  # --- 5. Gompertz Model Fitting ---
  cat("--- Step: Fitting Gompertz survival models ---\n")
  formula_crude <- as.formula(paste("Surv(", time_variable, ",", status_variable, ") ~ age_at_ins2_std"))
  fit_gompertz_crude <- flexsurv::flexsurvreg(formula_crude, data = train_data, dist = "gompertz")
  
  beta0 <- log(fit_gompertz_crude$res["rate", "est"])
  beta1 <- fit_gompertz_crude$res["age_at_ins2_std", "est"]
  gamma0 <- fit_gompertz_crude$res["shape", "est"]
  
  all_model_vars <- if (age_in_gompertz) c("age_at_ins2_std", model_vars_base) else model_vars_base
  
  if (length(all_model_vars) == 0) {
    formula_full_str <- paste("Surv(", time_variable, ",", status_variable, ") ~ 1")
  } else {
    formula_full_str <- paste("Surv(", time_variable, ",", status_variable, ") ~", paste(all_model_vars, collapse = " + "))
  }
  
  fit_gompertz_full <- flexsurv::flexsurvreg(as.formula(formula_full_str), data = train_data, dist = "gompertz")
  gamma1 <- fit_gompertz_full$res["shape", "est"]
  log_rate_intercept <- if (length(all_model_vars) > 0) log(fit_gompertz_full$res["rate", "est"]) else fit_gompertz_full$res["rate", "est"]
  beta_coeffs <- if (length(all_model_vars) > 0) fit_gompertz_full$res[all_model_vars, "est", drop = FALSE] else NULL
  
  # --- 6. Biological Age Scoring ---
  # Cumulative hazard M(t) formula for Gompertz: M(t) = 1 - S(t) = 1 - exp(-exp(xb) * (exp(gamma1 * t) - 1) / gamma1)
  calc_survival_M <- function(t, xb, gamma1) { 
    1 - exp(-exp(xb) * (exp(gamma1 * t) - 1) / gamma1) 
  }
  
  # Inverse function to find standard biological age
  calc_bio_age <- function(M, t, beta0, beta1, gamma0) {
    M[M >= 1] <- 0.999999
    M[M <= 0] <- 0.000001
    numerator <- log(-log(1 - M)) - log((exp(gamma0 * t) - 1) / gamma0)
    (numerator - beta0) / beta1
  }
  
  datasets <- list(train = train_data, test = test_data, disease = disease_test_data)
  for (d_name in names(datasets)) {
    d <- datasets[[d_name]]
    if (!is.null(beta_coeffs)) {
      X_matrix <- as.matrix(d[, rownames(beta_coeffs), drop = FALSE])
      d$xb <- X_matrix %*% beta_coeffs[, "est"] + log_rate_intercept
    } else {
      d$xb <- rep(log_rate_intercept, nrow(d))
    }
    
    # Calculate cumulative risk M at t = 10 years
    d$M_value <- calc_survival_M(t = 10, xb = d$xb, gamma1 = gamma1)
    
    # Project back to standardized biological age
    d$Bio_Age_std <- calc_bio_age(M = d$M_value, t = 10, beta0 = beta0, beta1 = beta1, gamma0 = gamma0)
    
    # Denormalize to get Biological Age in years
    d$Bio_Age <- d$Bio_Age_std * age_sd + age_mean 
    datasets[[d_name]] <- d
  }
  
  # --- 7. Evaluation and Visualizations ---
  evaluation_results <- run_evaluation_and_plotting(
    train_data = datasets$train, 
    test_data = datasets$test, 
    disease_test_data = datasets$disease, 
    clock_name_prefix = clock_name_prefix, 
    output_dir = output_dir_base, 
    sex_label = sex_label
  )
  
  return(evaluation_results)
}
