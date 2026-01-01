# Load necessary libraries
library(glmnet)
library(ggplot2)
library(dplyr)

# Read configuration data
config <- readRDS("transfer_learning_brain_data_fixed.rds")[["A.C.cortex"]]

# ============================================
# 1. Helper Functions
# ============================================

# Compute sparsity index (based on the method in the paper)
compute_sparsity_index <- function(target, source_list, t_star = NULL) {
  # target: target domain data, containing X and y
  # source_list: list of source domain data
  # t_star: threshold for SURE screening, default is min(n_k)^0.75
  
  p <- ncol(target$X)
  K <- length(source_list)
  
  # Calculate minimum sample size
  n_k <- sapply(source_list, function(s) nrow(s$X))
  n_0 <- nrow(target$X)
  n_star <- min(c(n_0, n_k))
  
  # Set t_star (according to the paper, use n_star^alpha, alpha=0.75)
  if(is.null(t_star)) {
    t_star <- ceiling(n_star^0.75)
  }
  
  # Calculate marginal statistics
  marginal_stats <- matrix(0, nrow = K, ncol = p)
  
  # Marginal statistics of target domain
  target_marginal <- t(target$X) %*% target$y / n_0
  
  for(k in 1:K) {
    n_k_sample <- nrow(source_list[[k]]$X)
    # Marginal statistics of source domain
    source_marginal <- t(source_list[[k]]$X) %*% source_list[[k]]$y / n_k_sample
    # Difference estimate
    delta_hat <- source_marginal - target_marginal
    
    # SURE screening: select the top t_star largest absolute values
    abs_delta <- abs(delta_hat)
    threshold <- sort(abs_delta, decreasing = TRUE)[min(t_star, p)]
    selected_idx <- which(abs_delta >= threshold)
    
    # Calculate sparsity index (L2 norm squared after screening)
    if(length(selected_idx) > 0) {
      marginal_stats[k, selected_idx] <- delta_hat[selected_idx]
    }
  }
  
  # Calculate sparsity index for each source domain
  sparsity_indices <- apply(marginal_stats, 1, function(x) sum(x^2))
  
  return(list(
    indices = sparsity_indices,
    t_star = t_star,
    n_star = n_star
  ))
}

# Construct candidate sets (based on the method in the paper)
construct_candidate_sets <- function(sparsity_indices, K) {
  # Sort by sparsity index from smallest to largest
  sorted_indices <- order(sparsity_indices)
  
  # Construct candidate sets: start from empty set, gradually add source domains
  candidate_sets <- list()
  candidate_sets[[1]] <- integer(0)  # Empty set
  
  for(l in 1:K) {
    candidate_sets[[l+1]] <- sorted_indices[1:l]
  }
  
  return(candidate_sets)
}

# Oracle Trans-Lasso algorithm (known informative set)
oracle_trans_lasso <- function(X_target, y_target, X_source_list, y_source_list, 
                               informative_set) {
  # Parameters:
  # X_target, y_target: target domain data
  # X_source_list, y_source_list: list of source domain data
  # informative_set: indices of known informative source domains
  
  n0 <- nrow(X_target)
  p <- ncol(X_target)
  
  # Step 1: Combine target domain and informative source domain data
  if(length(informative_set) > 0) {
    # Combine all informative source domain data
    X_informative <- do.call(rbind, X_source_list[informative_set])
    y_informative <- do.call(c, y_source_list[informative_set])
    
    # Calculate weighted sample size
    n_informative <- sum(sapply(X_source_list[informative_set], nrow))
    
    # Combine data
    X_combined <- rbind(X_target, X_informative)
    y_combined <- c(y_target, y_informative)
    
    # Use Lasso to estimate w_A
    fit_w <- cv.glmnet(X_combined, y_combined, alpha = 1, standardize = TRUE)
    w_hat <- as.numeric(coef(fit_w, s = "lambda.min"))[-1]  # Remove intercept
  } else {
    # If no informative source domains, use target domain data only
    w_hat <- rep(0, p)
    n_informative <- 0
  }
  
  # Step 2: Estimate delta
  # Calculate residuals
  r <- y_target - X_target %*% w_hat
  
  # Use Lasso to estimate delta
  fit_delta <- cv.glmnet(X_target, r, alpha = 1, standardize = TRUE)
  delta_hat <- as.numeric(coef(fit_delta, s = "lambda.min"))[-1]  # Remove intercept
  
  # Final estimate
  beta_hat <- w_hat + delta_hat
  
  return(list(
    beta = beta_hat,
    w_hat = w_hat,
    delta_hat = delta_hat,
    n_informative = n_informative
  ))
}

# Q-aggregation function (simplified version, select best candidate based on validation set)
q_aggregation <- function(candidate_models, X_val, y_val) {
  # Calculate MSE for each candidate model on validation set
  mse_values <- sapply(candidate_models, function(model) {
    y_pred <- X_val %*% model$beta
    mean((y_val - y_pred)^2)
  })
  
  # Select model with minimum MSE
  best_idx <- which.min(mse_values)
  
  return(list(
    best_model = candidate_models[[best_idx]],
    best_idx = best_idx,
    mse_values = mse_values
  ))
}

# Trans-Lasso algorithm (adaptive informative set)
trans_lasso <- function(X_target, y_target, X_source_list, y_source_list, 
                        val_ratio = 0.5, t_star = NULL) {
  n0 <- nrow(X_target)
  p <- ncol(X_target)
  K <- length(X_source_list)
  
  # Step 1: Sample splitting (for candidate model construction and aggregation)
  set.seed(123)
  val_idx <- sample(1:n0, size = round(n0 * val_ratio))
  train_idx <- setdiff(1:n0, val_idx)
  
  X_train <- X_target[train_idx, ]
  y_train <- y_target[train_idx]
  X_val <- X_target[val_idx, ]
  y_val <- y_target[val_idx]
  
  # Step 2: Calculate sparsity index and construct candidate sets
  target_data_train <- list(X = X_train, y = y_train)
  source_data_train <- list()
  for(k in 1:K) {
    source_data_train[[k]] <- list(X = X_source_list[[k]], y = y_source_list[[k]])
  }
  
  sparsity_result <- compute_sparsity_index(target_data_train, source_data_train, t_star)
  candidate_sets <- construct_candidate_sets(sparsity_result$indices, K)
  
  # Step 3: Run Oracle Trans-Lasso for each candidate set
  candidate_models <- list()
  
  for(l in 1:length(candidate_sets)) {
    informative_set <- candidate_sets[[l]]
    
    # Estimate model using training data
    model <- oracle_trans_lasso(X_train, y_train, 
                                X_source_list, y_source_list, 
                                informative_set)
    
    candidate_models[[l]] <- model
  }
  
  # Step 4: Q-aggregation (select best model using validation set)
  aggregation_result <- q_aggregation(candidate_models, X_val, y_val)
  
  # Step 5: Return the best model
  best_model <- aggregation_result$best_model
  
  return(list(
    best_model = best_model,
    candidate_models = candidate_models,
    candidate_sets = candidate_sets,
    sparsity_indices = sparsity_result$indices,
    aggregation_result = aggregation_result,
    train_idx = train_idx,
    val_idx = val_idx
  ))
}

# ============================================
# 2. 5-fold Cross-Validation Evaluation
# ============================================

evaluate_trans_lasso_cv <- function(config, n_folds = 5) {
  # Prepare data
  X_target <- config$target$X
  y_target <- config$target$y
  covariates <- config$covariates
  
  # Prepare source domain data
  X_source_list <- list()
  y_source_list <- list()
  
  for(i in 1:length(config$sources)) {
    X_source_list[[i]] <- config$sources[[i]]$X
    y_source_list[[i]] <- config$sources[[i]]$y
  }
  
  n <- nrow(X_target)
  p <- ncol(X_target)
  
  # Create fold indices
  set.seed(123)
  folds <- sample(rep(1:n_folds, length.out = n))
  
  # Store results
  results <- data.frame(
    fold = 1:n_folds,
    mse_lasso = numeric(n_folds),
    mse_trans_lasso = numeric(n_folds),
    mse_naive_trans_lasso = numeric(n_folds),
    n_informative = numeric(n_folds)
  )
  
  cat("Starting 5-fold cross-validation evaluation...\n")
  cat("Target domain: Anterior cingulate cortex (BA24)\n")
  cat(sprintf("Sample size: n=%d, Number of features: p=%d\n", n, p))
  cat(sprintf("Number of source domains: %d\n", length(X_source_list)))
  cat(strrep("=", 50), "\n")
  
  for(fold in 1:n_folds) {
    cat(sprintf("Processing fold %d/%d...\n", fold, n_folds))
    
    # Split data into training and test sets
    test_idx <- which(folds == fold)
    train_idx <- which(folds != fold)
    
    X_train <- X_target[train_idx, ]
    y_train <- y_target[train_idx]
    X_test <- X_target[test_idx, ]
    y_test <- y_target[test_idx]
    
    # Method 1: Lasso using only target domain (baseline)
    cat("  Training Lasso (baseline)...\n")
    lasso_fit <- cv.glmnet(X_train, y_train, alpha = 1, standardize = TRUE)
    lasso_pred <- predict(lasso_fit, X_test, s = "lambda.min")
    mse_lasso <- mean((y_test - lasso_pred)^2)
    
    # Method 2: Trans-Lasso (adaptive informative set)
    cat("  Training Trans-Lasso...\n")
    trans_lasso_result <- trans_lasso(X_train, y_train, X_source_list, y_source_list)
    trans_lasso_pred <- X_test %*% trans_lasso_result$best_model$beta
    mse_trans_lasso <- mean((y_test - trans_lasso_pred)^2)
    
    # Method 3: Naive Trans-Lasso (assume all source domains are informative)
    cat("  Training Naive Trans-Lasso...\n")
    naive_trans_lasso <- oracle_trans_lasso(X_train, y_train, 
                                            X_source_list, y_source_list,
                                            informative_set = 1:length(X_source_list))
    naive_pred <- X_test %*% naive_trans_lasso$beta
    mse_naive_trans_lasso <- mean((y_test - naive_pred)^2)
    
    # Record results
    results[fold, "mse_lasso"] <- mse_lasso
    results[fold, "mse_trans_lasso"] <- mse_trans_lasso
    results[fold, "mse_naive_trans_lasso"] <- mse_naive_trans_lasso
    results[fold, "n_informative"] <- length(trans_lasso_result$candidate_sets[[trans_lasso_result$aggregation_result$best_idx]])
    
    cat(sprintf("  Done. MSE: Lasso=%.4f, Trans-Lasso=%.4f, Naive=%.4f\n", 
                mse_lasso, mse_trans_lasso, mse_naive_trans_lasso))
    cat(strrep("  -", 25), "\n")
  }
  
  # Calculate average results
  avg_results <- data.frame(
    Method = c("Lasso", "Naive Trans-Lasso", "Trans-Lasso"),
    Avg_MSE = c(mean(results$mse_lasso), 
                mean(results$mse_naive_trans_lasso), 
                mean(results$mse_trans_lasso)),
    SE_MSE = c(sd(results$mse_lasso)/sqrt(n_folds),
               sd(results$mse_naive_trans_lasso)/sqrt(n_folds),
               sd(results$mse_trans_lasso)/sqrt(n_folds))
  )
  
  cat("Cross-validation completed!\n")
  print(avg_results)
  
  return(list(
    results = results,
    avg_results = avg_results,
    folds = folds
  ))
}

# ============================================
# 3. Extract Final Model and Important Genes
# ============================================

extract_final_model_and_genes <- function(config, top_n = 15) {
  # Prepare data
  X_target <- config$target$X
  y_target <- config$target$y
  covariates <- config$covariates
  
  # Prepare source domain data
  X_source_list <- list()
  y_source_list <- list()
  
  for(i in 1:length(config$sources)) {
    X_source_list[[i]] <- config$sources[[i]]$X
    y_source_list[[i]] <- config$sources[[i]]$y
  }
  
  cat("Training final Trans-Lasso model...\n")
  
  # Train Trans-Lasso using all data
  final_result <- trans_lasso(X_target, y_target, X_source_list, y_source_list)
  
  # Extract coefficients of the best model
  final_beta <- final_result$best_model$beta
  names(final_beta) <- covariates
  
  # Select most important genes (sorted by absolute coefficient value)
  beta_abs <- abs(final_beta)
  sorted_idx <- order(beta_abs, decreasing = TRUE)
  
  top_genes <- data.frame(
    Rank = 1:min(top_n, length(final_beta)),
    Gene = covariates[sorted_idx[1:min(top_n, length(final_beta))]],
    Coefficient = final_beta[sorted_idx[1:min(top_n, length(final_beta))]],
    Absolute_Coefficient = beta_abs[sorted_idx[1:min(top_n, length(final_beta))]]
  )
  
  # Output distribution of sparsity indices
  sparsity_df <- data.frame(
    Source = 1:length(final_result$sparsity_indices),
    Sparsity_Index = final_result$sparsity_indices,
    Selected = FALSE
  )
  
  # Mark selected source domains
  selected_set <- final_result$candidate_sets[[final_result$aggregation_result$best_idx]]
  sparsity_df$Selected[selected_set] <- TRUE
  
  cat("\nFinal model information:\n")
  cat(sprintf("Number of informative source domains selected: %d\n", length(selected_set)))
  cat(sprintf("Total number of source domains: %d\n", length(X_source_list)))
  cat(sprintf("Top %d most important genes:\n", top_n))
  print(top_genes)
  
  return(list(
    final_model = final_result$best_model,
    top_genes = top_genes,
    sparsity_indices = sparsity_df,
    candidate_sets = final_result$candidate_sets,
    best_set = selected_set
  ))
}

# ============================================
# 4. Visualization Functions
# ============================================

plot_cv_results <- function(cv_results) {
  # Prepare data
  plot_data <- cv_results$avg_results
  plot_data$Method <- factor(plot_data$Method, 
                             levels = c("Lasso", "Naive Trans-Lasso", "Trans-Lasso"))
  
  # Create bar plot
  p <- ggplot(plot_data, aes(x = Method, y = Avg_MSE, fill = Method)) +
    geom_bar(stat = "identity", position = position_dodge(), width = 0.7) +
    geom_errorbar(aes(ymin = Avg_MSE - SE_MSE, ymax = Avg_MSE + SE_MSE), 
                  width = 0.2, position = position_dodge(0.7)) +
    labs(title = "5-Fold Cross-Validation Results: JAM2 Gene Expression Prediction",
         subtitle = "Target Domain: Anterior Cingulate Cortex (BA24)",
         x = "Method", y = "Mean Squared Error (MSE) ± Standard Error") +
    theme_minimal(base_size = 14) +
    theme(legend.position = "none",
          axis.text.x = element_text(angle = 15, hjust = 1),
          plot.title = element_text(hjust = 0.5, face = "bold"),
          plot.subtitle = element_text(hjust = 0.5)) +
    scale_fill_manual(values = c("#1f77b4", "#ff7f0e", "#2ca02c"))
  
  return(p)
}

plot_top_genes <- function(top_genes_result, top_n = 15) {
  # Prepare data
  top_genes <- top_genes_result$top_genes[1:min(top_n, nrow(top_genes_result$top_genes)), ]
  top_genes$Gene <- factor(top_genes$Gene, levels = rev(top_genes$Gene))
  
  # Create bar plot
  p <- ggplot(top_genes, aes(x = Gene, y = Absolute_Coefficient)) +
    geom_bar(stat = "identity", fill = "steelblue", width = 0.7) +
    coord_flip() +
    labs(title = "Most Important Genes Selected by Trans-Lasso",
         subtitle = "Sorted by Absolute Coefficient Value",
         x = "Gene", y = "Absolute Coefficient Value") +
    theme_minimal(base_size = 14) +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"),
          plot.subtitle = element_text(hjust = 0.5),
          axis.text.y = element_text(size = 10))
  
  return(p)
}

plot_sparsity_indices <- function(top_genes_result) {
  # Prepare data
  sparsity_df <- top_genes_result$sparsity_indices
  sparsity_df$Source <- factor(sparsity_df$Source)
  sparsity_df$Selected <- factor(sparsity_df$Selected, levels = c(TRUE, FALSE))
  
  # Create scatter plot
  p <- ggplot(sparsity_df, aes(x = Source, y = Sparsity_Index, color = Selected, shape = Selected)) +
    geom_point(size = 3) +
    labs(title = "Distribution of Sparsity Indices Across Source Domains",
         subtitle = "Red indicates selected as informative source domains",
         x = "Source Domain Index", y = "Sparsity Index") +
    theme_minimal(base_size = 14) +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"),
          plot.subtitle = element_text(hjust = 0.5),
          axis.text.x = element_text(angle = 45, hjust = 1)) +
    scale_color_manual(values = c("TRUE" = "red", "FALSE" = "blue"),
                       labels = c("TRUE" = "Selected", "FALSE" = "Not Selected")) +
    scale_shape_manual(values = c("TRUE" = 17, "FALSE" = 16),
                       labels = c("TRUE" = "Selected", "FALSE" = "Not Selected"))
  
  return(p)
}

# ============================================
# 5. Main Program
# ============================================

main <- function() {
  cat(strrep("=", 60), "\n", sep = "")
  cat("Transfer Learning Experiment: Gene Expression Prediction using Trans-Lasso\n")
  cat("Target Domain: Anterior Cingulate Cortex (BA24)\n")
  cat("Current working directory:", getwd(), "\n")
  cat(strrep("=", 60), "\n\n")
  
  # Check data
  cat("Data Overview:\n")
  cat(sprintf("Target domain sample size: %d\n", nrow(config$target$X)))
  cat(sprintf("Number of features: %d\n", ncol(config$target$X)))
  cat(sprintf("Number of source domains: %d\n", length(config$sources)))
  
  # Step 1: 5-fold cross-validation evaluation
  cat("\nStep 1: Performing 5-fold cross-validation evaluation...\n")
  cv_results <- evaluate_trans_lasso_cv(config, n_folds = 5)
  
  # Save cross-validation results
  write.csv(cv_results$results, "cv_results.csv", row.names = FALSE)
  write.csv(cv_results$avg_results, "cv_avg_results.csv", row.names = FALSE)
  
  # Step 2: Extract final model and important genes
  cat("\n\nStep 2: Extracting final model and important genes...\n")
  final_result <- extract_final_model_and_genes(config, top_n = 15)
  
  # Save important gene list
  write.csv(final_result$top_genes, "top_15_genes.csv", row.names = FALSE)
  write.csv(final_result$sparsity_indices, "sparsity_indices.csv", row.names = FALSE)
  
  # Step 3: Visualization
  cat("\n\nStep 3: Generating visualization results...\n")
  
  # Create visualizations
  p1 <- plot_cv_results(cv_results)
  p2 <- plot_top_genes(final_result)
  p3 <- plot_sparsity_indices(final_result)
  
  # Save plots
  ggsave("cv_results_plot.pdf", p1, width = 10, height = 6, dpi = 300)
  ggsave("top_genes_plot.pdf", p2, width = 10, height = 8, dpi = 300)
  ggsave("sparsity_indices_plot.pdf", p3, width = 12, height = 6, dpi = 300)
  
  # Display results on screen
  print(p1)
  print(p2)
  print(p3)
  
  # Step 4: Generate comprehensive report
  cat("\n\n", strrep("=", 60), "\n", sep = "")
  cat("Experiment Summary Report\n")
  cat(strrep("=", 60), "\n")
  
  cat("\n1. Prediction Performance Comparison (5-fold CV):\n")
  print(cv_results$avg_results)
  
  cat(sprintf("\n2. Relative Improvement (compared to Lasso):\n"))
  lasso_mse <- cv_results$avg_results$Avg_MSE[1]
  trans_lasso_mse <- cv_results$avg_results$Avg_MSE[3]
  improvement <- 100 * (lasso_mse - trans_lasso_mse) / lasso_mse
  cat(sprintf("   MSE reduction of Trans-Lasso relative to Lasso: %.2f%%\n", improvement))
  
  cat("\n3. Important Gene Identification:\n")
  cat(sprintf("   Number of informative source domains selected by Trans-Lasso: %d\n", length(final_result$best_set)))
  cat("   Top 5 most important genes:\n")
  print(final_result$top_genes[1:5, ])
  
  cat("\n4. Output Files:\n")
  cat("   - cv_results.csv: Detailed 5-fold cross-validation results\n")
  cat("   - cv_avg_results.csv: Average cross-validation results\n")
  cat("   - top_15_genes.csv: Top 15 most important genes\n")
  cat("   - sparsity_indices.csv: Sparsity indices of source domains\n")
  cat("   - cv_results_plot.pdf: Cross-validation results plot\n")
  cat("   - top_genes_plot.pdf: Important genes bar plot\n")
  cat("   - sparsity_indices_plot.pdf: Sparsity indices distribution plot\n")
  
  cat("\nExperiment completed!\n")
  
  return(list(
    cv_results = cv_results,
    final_result = final_result
  ))
}

# Run main program
if (interactive()) {
  results <- main()
}