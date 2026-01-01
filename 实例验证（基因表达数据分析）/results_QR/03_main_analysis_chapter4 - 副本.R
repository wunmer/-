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

# 创建带CV次数统计的cv.hqreg包装函数
cv_hqreg_with_counter <- function(X, y, counter_type = "unknown", counter_num = 0, 
                                  method = "quantile", tau = 0.5, alpha = 1, 
                                  penalty.factor = rep(1, ncol(X)), ...) {
  cat(sprintf("【CV循环统计】开始第%d次10折CV (类型: %s)\n", counter_num, counter_type))
  
  # 记录开始时间
  start_time <- Sys.time()
  
  # 运行交叉验证
  fit <- cv.hqreg(X, y, method = method, tau = tau, alpha = alpha, 
                  penalty.factor = penalty.factor, ...)
  
  # 计算耗时
  time_elapsed <- difftime(Sys.time(), start_time, units = "secs")
  
  cat(sprintf("【CV循环统计】第%d次10折CV完成 (类型: %s), 耗时: %.2f秒\n\n", 
              counter_num, counter_type, time_elapsed))
  
  return(fit)
}

# 修改qr_trans_adlasso函数，添加CV计数器
qr_trans_adlasso <- function(target, sources, tau = 0.5) {
  n <- nrow(target$X); set.seed(123)
  train_idx <- sample(n, floor(0.8*n))
  X_tr <- target$X[train_idx, ]; y_tr <- target$y[train_idx]
  X_te <- target$X[-train_idx, ]; y_te <- target$y[-train_idx]
  
  cands <- construct_candidates(list(X=X_tr, y=y_tr), sources)
  beta_list <- list()
  
  cat(sprintf("=== 开始QR Trans-AdLASSO (tau=%.2f), 共有%d个候选集 ===\n", tau, length(cands)))
  
  # 初始化CV计数器
  cv_counter <- 0
  
  for (cand_idx in 1:length(cands)) {
    idx <- cands[[cand_idx]]
    cat(sprintf("\n--- 处理候选集 %d/%d (使用源域数量: %d) ---\n", 
                cand_idx, length(cands), length(idx)))
    
    if (length(idx) == 0) {
      omega <- rep(0, ncol(target$X) + 1)
      cat("  候选集为空，omega设为0向量\n")
    } else {
      src_X <- do.call(rbind, lapply(sources[idx], function(s) s$X))
      src_y <- unlist(lapply(sources[idx], function(s) s$y))
      
      # 第一次CV: 初始拟合
      cv_counter <- cv_counter + 1
      cat(sprintf("  步骤1/3: 初始源域拟合 (CV %d)\n", cv_counter))
      fit0 <- cv_hqreg_with_counter(
        src_X, src_y, 
        counter_type = "初始源域拟合", 
        counter_num = cv_counter,
        method = "quantile", 
        tau = tau, 
        alpha = 1
      )
      w <- 1 / (abs(coef(fit0)[-1]) + 1e-8)
      
      # 第二次CV: 带权重的自适应LASSO
      cv_counter <- cv_counter + 1
      cat(sprintf("  步骤2/3: 自适应LASSO拟合 (CV %d)\n", cv_counter))
      fit_omega <- cv_hqreg_with_counter(
        src_X, src_y,
        counter_type = "自适应LASSO",
        counter_num = cv_counter,
        method = "quantile",
        tau = tau,
        alpha = 1,
        penalty.factor = w
      )
      omega <- coef(fit_omega)
    }
    
    # 第三次CV: 目标域残差拟合
    cv_counter <- cv_counter + 1
    cat(sprintf("  步骤3/3: 目标域残差拟合 (CV %d)\n", cv_counter))
    fit_delta <- cv_hqreg_with_counter(
      X_tr, y_tr - cbind(1, X_tr) %*% omega,
      counter_type = "目标域残差拟合",
      counter_num = cv_counter,
      method = "quantile",
      tau = tau,
      alpha = 1
    )
    delta <- coef(fit_delta)
    
    beta_list[[length(beta_list)+1]] <- omega + delta
  }
  
  cat(sprintf("\n=== QR Trans-AdLASSO完成，总计%d次10折CV ===\n\n", cv_counter))
  
  final_beta <- q_aggregate(beta_list, X_te, y_te)
  pred <- cbind(1, X_te) %*% final_beta
  return(mean((y_te - pred)^2))
}

# 修改qr_trans_adenet函数
qr_trans_adenet <- function(target, sources, tau = 0.5, alpha = 0.5) {
  n <- nrow(target$X); set.seed(123)
  train_idx <- sample(n, floor(0.8*n))
  X_tr <- target$X[train_idx, ]; y_tr <- target$y[train_idx]
  X_te <- target$X[-train_idx, ]; y_te <- target$y[-train_idx]
  
  cands <- construct_candidates(list(X=X_tr, y=y_tr), sources)
  beta_list <- list()
  
  cat(sprintf("=== 开始QR Trans-AdE-net (tau=%.2f, alpha=%.2f), 共有%d个候选集 ===\n", 
              tau, alpha, length(cands)))
  
  # 初始化CV计数器
  cv_counter <- 0
  
  for (cand_idx in 1:length(cands)) {
    idx <- cands[[cand_idx]]
    cat(sprintf("\n--- 处理候选集 %d/%d (使用源域数量: %d) ---\n", 
                cand_idx, length(cands), length(idx)))
    
    if (length(idx) == 0) {
      omega <- rep(0, ncol(target$X) + 1)
      cat("  候选集为空，omega设为0向量\n")
    } else {
      src_X <- do.call(rbind, lapply(sources[idx], function(s) s$X))
      src_y <- unlist(lapply(sources[idx], function(s) s$y))
      
      # 第一次CV
      cv_counter <- cv_counter + 1
      cat(sprintf("  步骤1/3: 初始源域拟合 (CV %d)\n", cv_counter))
      fit0 <- cv_hqreg_with_counter(
        src_X, src_y, 
        counter_type = "初始源域拟合(AdE-net)", 
        counter_num = cv_counter,
        method = "quantile", 
        tau = tau, 
        alpha = alpha
      )
      w <- 1 / (abs(coef(fit0)[-1]) + 1e-8)
      
      # 第二次CV
      cv_counter <- cv_counter + 1
      cat(sprintf("  步骤2/3: 自适应弹性网拟合 (CV %d)\n", cv_counter))
      fit_omega <- cv_hqreg_with_counter(
        src_X, src_y,
        counter_type = "自适应弹性网",
        counter_num = cv_counter,
        method = "quantile",
        tau = tau,
        alpha = alpha,
        penalty.factor = w
      )
      omega <- coef(fit_omega)
    }
    
    # 第三次CV
    cv_counter <- cv_counter + 1
    cat(sprintf("  步骤3/3: 目标域残差拟合 (CV %d)\n", cv_counter))
    fit_delta <- cv_hqreg_with_counter(
      X_tr, y_tr - cbind(1, X_tr) %*% omega,
      counter_type = "目标域残差拟合(AdE-net)",
      counter_num = cv_counter,
      method = "quantile",
      tau = tau,
      alpha = alpha
    )
    delta <- coef(fit_delta)
    
    beta_list[[length(beta_list)+1]] <- omega + delta
  }
  
  cat(sprintf("\n=== QR Trans-AdE-net完成，总计%d次10折CV ===\n\n", cv_counter))
  
  final_beta <- q_aggregate(beta_list, X_te, y_te)
  pred <- cbind(1, X_te) %*% final_beta
  return(mean((y_te - pred)^2))
}

# 修改qr_naive函数
qr_naive <- function(target, tau = 0.5) {
  n <- nrow(target$X); set.seed(123)
  train_idx <- sample(n, floor(0.8*n))
  X_tr <- target$X[train_idx, ]; y_tr <- target$y[train_idx]
  X_te <- target$X[-train_idx, ]; y_te <- target$y[-train_idx]
  
  cat(sprintf("=== 开始朴素QR-Enet (tau=%.2f) ===\n", tau))
  
  # 只有一次CV
  cv_counter <- 1
  cat(sprintf("  步骤: 直接拟合目标域 (CV %d)\n", cv_counter))
  
  fit <- cv_hqreg_with_counter(
    X_tr, y_tr,
    counter_type = "朴素QR-Enet",
    counter_num = cv_counter,
    method = "quantile",
    tau = tau,
    alpha = 0.5
  )
  
  cat(sprintf("\n=== 朴素QR-Enet完成，总计%d次10折CV ===\n\n", cv_counter))
  
  beta <- coef(fit)
  pred <- cbind(1, X_te) %*% beta
  return(mean((y_te - pred)^2))
}

# 运行实验
cat("===========================================\n")
cat("开始分位数回归迁移学习实验\n")
cat("目标域: Anterior cingulate cortex (BA24)\n")
cat("样本量: n=147, 特征数: p=407\n")
cat("===========================================\n\n")

tau_values <- c(0.25, 0.5, 0.75)
results <- data.frame()

# 初始化总CV计数器
total_cv_count <- 0
method_cv_counts <- list()

for (tau in tau_values) {
  cat("\n", strrep("=", 50), "\n", sep = "")
  cat(sprintf("正在计算 tau = %.2f ...\n", tau))
  cat(strrep("=", 50), "\n\n", sep = "")
  
  # 记录每个方法开始的CV计数
  cat("1. 运行朴素QR-Enet...\n")
  m1 <- qr_naive(config$target, tau = tau)
  
  cat("\n2. 运行QR Trans-AdLASSO...\n")
  m2 <- qr_trans_adlasso(config$target, config$sources, tau = tau)
  
  cat("\n3. 运行QR Trans AdE-net...\n")
  m3 <- qr_trans_adenet(config$target, config$sources, tau = tau)
  
  results <- rbind(results, data.frame(
    Tau = tau,
    Method = c("Naïve QR-Enet", "QR Trans-AdLASSO", "QR Trans AdE-net"),
    MSE = c(m1, m2, m3)
  ))
  
  cat(sprintf("\nτ=%.2f 完成! MSE结果: Naïve=%.4f, AdLASSO=%.4f, AdE-net=%.4f\n\n",
              tau, m1, m2, m3))
}

cat("===========================================\n")
cat("所有分位数水平计算完成！\n")
print(results)
cat("===========================================\n\n")

# 可视化
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

cat("分位数回归预测完成，结果已保存\n")

# 第二步 输出最佳模型的 β + 选出 Top15 基因
library(dplyr)

cat("\n", strrep("=", 60), "\n", sep = "")
cat("开始提取最佳模型系数和Top15重要基因\n")
cat(strrep("=", 60), "\n\n", sep = "")

# 重新定义一个带CV计数器的最终系数提取函数
get_final_beta <- function(target, sources, tau = 0.5) {
  n <- nrow(target$X); set.seed(123)
  train_idx <- sample(n, floor(0.8*n))
  X_tr <- target$X[train_idx, ]; y_tr <- target$y[train_idx]
  X_te <- target$X[-train_idx, ]; y_te <- target$y[-train_idx]
  
  cands <- construct_candidates(list(X=X_tr, y=y_tr), sources)
  beta_list <- list()
  
  cat(sprintf("=== 提取最终系数 (tau=%.2f), 共有%d个候选集 ===\n", tau, length(cands)))
  
  cv_counter <- 0
  
  for (cand_idx in 1:length(cands)) {
    idx <- cands[[cand_idx]]
    cat(sprintf("\n--- 处理候选集 %d/%d (使用源域数量: %d) ---\n", 
                cand_idx, length(cands), length(idx)))
    
    if (length(idx)==0) {
      omega <- rep(0, ncol(target$X)+1)
      cat("  候选集为空，omega设为0向量\n")
    } else {
      src_X <- do.call(rbind, lapply(sources[idx], function(s) s$X))
      src_y <- unlist(lapply(sources[idx], function(s) s$y))
      
      cv_counter <- cv_counter + 1
      cat(sprintf("  步骤1/3: 初始源域拟合 (CV %d)\n", cv_counter))
      fit0 <- cv_hqreg_with_counter(
        src_X, src_y, 
        counter_type = "最终系数-初始拟合", 
        counter_num = cv_counter,
        method = "quantile", 
        tau = tau, 
        alpha = 1
      )
      w <- 1 / (abs(coef(fit0)[-1]) + 1e-8)
      
      cv_counter <- cv_counter + 1
      cat(sprintf("  步骤2/3: 自适应LASSO拟合 (CV %d)\n", cv_counter))
      fit_omega <- cv_hqreg_with_counter(
        src_X, src_y, 
        counter_type = "最终系数-自适应LASSO", 
        counter_num = cv_counter,
        method = "quantile", 
        tau = tau, 
        alpha = 1, 
        penalty.factor = w
      )
      omega <- coef(fit_omega)
    }
    
    cv_counter <- cv_counter + 1
    cat(sprintf("  步骤3/3: 目标域残差拟合 (CV %d)\n", cv_counter))
    fit_delta <- cv_hqreg_with_counter(
      X_tr, y_tr - cbind(1, X_tr) %*% omega,
      counter_type = "最终系数-目标域拟合",
      counter_num = cv_counter,
      method = "quantile",
      tau = tau,
      alpha = 1
    )
    delta <- coef(fit_delta)
    
    beta_list[[length(beta_list)+1]] <- omega + delta
  }
  
  cat(sprintf("\n=== 最终系数提取完成，总计%d次10折CV ===\n\n", cv_counter))
  
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

cat("\n", strrep("=", 60), "\n", sep = "")
cat("Top15 重要基因结果:\n")
cat(strrep("=", 60), "\n", sep = "")
print(top15_df)

# 保存表格
write.csv(top15_df, "Top15_重要基因.csv", row.names = FALSE)

# 绘制条形图
ggplot(top15_df, aes(x = reorder(Gene, Abs_Coefficient), y = Coefficient)) +
  geom_col(fill = "steelblue", width = 0.7) +
  coord_flip() +
  labs(title = "QR Trans-AdLASSO 选出的前15个最重要的调控JAM2的基因 (τ=0.5)",
       x = "基因", y = "分位数回归系数") +
  theme_minimal()

ggsave("Top15基因条形图.pdf", width = 10, height = 6, dpi = 300)

cat("\n所有分析完成！结果已保存到文件。\n")