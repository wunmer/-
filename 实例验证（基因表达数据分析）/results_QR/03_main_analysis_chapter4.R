library(hqreg)
library(ggplot2)

config <- readRDS("transfer_learning_brain_data_fixed.rds")[["A.C.cortex"]]

construct_candidates <- function(target_train, source_list) {
  cor_target <- apply(target_train$X, 2, function(x) cor(x, target_train$y))
  diffs <- sapply(source_list, function(src) {
    cor_src <- apply(src$X, 2, function(x) cor(x, src$y))
    mean(abs(cor_target - cor_src), na.rm = TRUE)
  })
  ordered <- order(diffs)
  candidates <- list()
  candidates[[1]] <- integer(0)
  for (k in 1:min(6, length(source_list))) candidates[[k+1]] <- ordered[1:k]
  return(candidates)
}

q_aggregate <- function(beta_list, X_test, y_test) {
  mse <- sapply(beta_list, function(b) mean((y_test - cbind(1, X_test) %*% b)^2))
  beta_list[[which.min(mse)]]
}

# 所有coef()都是长度p+1（直接整体用）
qr_trans_adlasso <- function(target, sources, tau = 0.5) {
  n <- nrow(target$X); set.seed(123)
  train_idx <- sample(n, floor(0.8*n))
  X_tr <- target$X[train_idx, ]; y_tr <- target$y[train_idx]
  X_te <- target$X[-train_idx, ]; y_te <- target$y[-train_idx]
  
  cands <- construct_candidates(list(X=X_tr, y=y_tr), sources)
  beta_list <- list()
  
  for (idx in cands) {
    if (length(idx) == 0) {
      omega <- rep(0, ncol(target$X) + 1)
    } else {
      src_X <- do.call(rbind, lapply(sources[idx], function(s) s$X))
      src_y <- unlist(lapply(sources[idx], function(s) s$y))
      
      fit0 <- cv.hqreg(src_X, src_y, method = "quantile", tau = tau, alpha = 1)
      w <- 1 / (abs(coef(fit0)[-1]) + 1e-8)               # 长度 p
      
      fit_omega <- cv.hqreg(src_X, src_y, method = "quantile", tau = tau, 
                            alpha = 1, penalty.factor = w)
      omega <- coef(fit_omega)                           # 长度 p+1
    }
    
    resid <- y_tr - cbind(1, X_tr) %*% omega
    fit_delta <- cv.hqreg(X_tr, resid, method = "quantile", tau = tau, alpha = 1)
    delta <- coef(fit_delta)                            # 长度 p+1
    
    beta_list[[length(beta_list)+1]] <- omega + delta
  }
  
  final_beta <- q_aggregate(beta_list, X_te, y_te)
  pred <- cbind(1, X_te) %*% final_beta                 # 整体乘
  return(mean((y_te - pred)^2))
}

qr_trans_adenet <- function(target, sources, tau = 0.5, alpha = 0.5) {
  n <- nrow(target$X); set.seed(123)
  train_idx <- sample(n, floor(0.8*n))
  X_tr <- target$X[train_idx, ]; y_tr <- target$y[train_idx]
  X_te <- target$X[-train_idx, ]; y_te <- target$y[-train_idx]
  
  cands <- construct_candidates(list(X=X_tr, y=y_tr), sources)
  beta_list <- list()
  
  for (idx in cands) {
    if (length(idx) == 0) {
      omega <- rep(0, ncol(target$X) + 1)
    } else {
      src_X <- do.call(rbind, lapply(sources[idx], function(s) s$X))
      src_y <- unlist(lapply(sources[idx], function(s) s$y))
      
      fit0 <- cv.hqreg(src_X, src_y, method = "quantile", tau = tau, alpha = alpha)
      w <- 1 / (abs(coef(fit0)[-1]) + 1e-8)
      
      fit_omega <- cv.hqreg(src_X, src_y, method = "quantile", tau = tau, 
                            alpha = alpha, penalty.factor = w)
      omega <- coef(fit_omega)
    }
    
    resid <- y_tr - cbind(1, X_tr) %*% omega
    fit_delta <- cv.hqreg(X_tr, resid, method = "quantile", tau = tau, alpha = alpha)
    delta <- coef(fit_delta)
    
    beta_list[[length(beta_list)+1]] <- omega + delta
  }
  
  final_beta <- q_aggregate(beta_list, X_te, y_te)
  pred <- cbind(1, X_te) %*% final_beta
  return(mean((y_te - pred)^2))
}

qr_naive <- function(target, tau = 0.5) {
  n <- nrow(target$X); set.seed(123)
  train_idx <- sample(n, floor(0.8*n))
  X_tr <- target$X[train_idx, ]; y_tr <- target$y[train_idx]
  X_te <- target$X[-train_idx, ]; y_te <- target$y[-train_idx]
  
  fit <- cv.hqreg(X_tr, y_tr, method = "quantile", tau = tau, alpha = 0.5)
  beta <- coef(fit)
  pred <- cbind(1, X_te) %*% beta
  return(mean((y_te - pred)^2))
}

# 运行
tau_values <- c(0.25, 0.5, 0.75)
results <- data.frame()

for (tau in tau_values) {
  cat("\n正在计算 tau =", tau, "...\n")
  m1 <- qr_naive(config$target, tau = tau)
  m2 <- qr_trans_adlasso(config$target, config$sources, tau = tau)
  m3 <- qr_trans_adenet(config$target, config$sources, tau = tau)
  
  results <- rbind(results, data.frame(
    Tau = tau,
    Method = c("Naïve QR-Enet", "QR Trans-AdLASSO", "QR Trans AdE-net"),
    MSE = c(m1, m2, m3)
  ))
}

print(results)

ggplot(results, aes(x = Method, y = MSE, fill = Method)) +
  geom_col(width = 0.7, color = "black") +
  facet_wrap(~ Tau, scales = "free_y") +
  labs(title = "实证结果：JAM2基因表达量分位数回归预测",
       subtitle = "目标域：Anterior cingulate cortex (BA24), n=147, p=407",
       y = "预测均方误差(MSE)") +
  theme_minimal(base_size = 13) +
  theme(axis.text.x = element_text(angle = 15, hjust = 1),
        legend.position = "none")
ggsave("JAM2基因表达量分位数回归预测图.pdf", width = 11, height = 6, dpi = 300)

cat("\n分位数回归预测完成\n")


#第二步 输出最佳模型的 β + 选出 Top15 基因
library(dplyr)

# 重新拟合一次最终模型（用全部训练数据），得到最终的 β（τ=0.5 时表现最好）
final_fit <- qr_trans_adlasso(config$target, config$sources, tau = 0.5)

# final_fit 现在是一个数值（MSE），需要重新跑一次但返回系数 beta 的版本（只跑一次全数据）

get_final_beta <- function(target, sources, tau = 0.5) {
  n <- nrow(target$X); set.seed(123)
  train_idx <- sample(n, floor(0.8*n))
  X_tr <- target$X[train_idx, ]; y_tr <- target$y[train_idx]
  X_te <- target$X[-train_idx, ]; y_te <- target$y[-train_idx]
  
  cands <- construct_candidates(list(X=X_tr, y=y_tr), sources)
  beta_list <- list()
  
  for (idx in cands) {
    if (length(idx)==0) {
      omega <- rep(0, ncol(target$X)+1)
    } else {
      src_X <- do.call(rbind, lapply(sources[idx], function(s) s$X))
      src_y <- unlist(lapply(sources[idx], function(s) s$y))
      fit0 <- cv.hqreg(src_X, src_y, method="quantile", tau=tau, alpha=1)
      w <- 1 / (abs(coef(fit0)[-1]) + 1e-8)
      fit_omega <- cv.hqreg(src_X, src_y, method="quantile", tau=tau, alpha=1, penalty.factor=w)
      omega <- coef(fit_omega)
    }
    
    resid <- y_tr - cbind(1, X_tr) %*% omega
    fit_delta <- cv.hqreg(X_tr, resid, method="quantile", tau=tau, alpha=1)
    delta <- coef(fit_delta)
    
    beta_list[[length(beta_list)+1]] <- omega + delta
  }
  
  # Q聚合选最好的
  mse_vec <- sapply(beta_list, function(b) mean((y_te - cbind(1, X_te) %*% b)^2))
  best_beta <- beta_list[[which.min(mse_vec)]]
  names(best_beta) <- c("Intercept", config$covariates)   
  return(best_beta)
}

# 得到最终系数向量（含截距）
beta_final <- get_final_beta(config$target, config$sources, tau = 0.5)

# 输出 Top 15 最重要的基因（绝对值最大的）
top15 <- beta_final[-1] %>% abs() %>% sort(decreasing = TRUE) %>% head(15)
top15_df <- data.frame(
  Rank = 1:15,
  Gene = names(top15),
  Coefficient = beta_final[-1][names(top15)],
  Abs_Coefficient = top15
)

print(top15_df)

# 保存表格
write.csv(top15_df, "Top15_重要基因.csv", row.names = FALSE)

# 绘制条形图
library(ggplot2)
ggplot(top15_df, aes(x = reorder(Gene, Abs_Coefficient), y = Coefficient)) +
  geom_col(fill = "steelblue", width = 0.7) +
  coord_flip() +
  labs(title = "QR Trans-AdLASSO 选出的前15个最重要的调控JAM2的基因 (τ=0.5)",
       x = "基因", y = "分位数回归系数") +
  theme_minimal()
ggsave("Top15基因条形图.pdf", width = 10, height = 6, dpi = 300)