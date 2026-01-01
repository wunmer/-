# 加载必要的包
library(data.table)
library(dplyr)

# ==================== 1. 读取基因ID列表并处理版本号 ====================

# 读取之前提取的CNS相关基因的ENSEMBL ID（带版本号）
ensembl_ids <- readLines("cns_ensembl_ids_final.txt")
cat("从文件中读取的基因数量:", length(ensembl_ids), "\n")

# 去掉版本号，只保留主ID
remove_version <- function(gene_ids) {
  sapply(gene_ids, function(x) {
    parts <- strsplit(x, "\\.")[[1]]
    if (length(parts) > 1) {
      return(parts[1])  # 返回主ID
    } else {
      return(x)  # 如果没有版本号，直接返回
    }
  }, USE.NAMES = FALSE)
}

# 去掉版本号
ensembl_ids_no_version <- remove_version(ensembl_ids)
ensembl_ids_no_version <- unique(ensembl_ids_no_version)

cat("去掉版本号后的基因数量（去重）:", length(ensembl_ids_no_version), "\n")

# JAM2基因
jam2_id_with_version <- "ENSG00000154721.14"
jam2_id_no_version <- remove_version(jam2_id_with_version)

cat("JAM2基因ID:\n")
cat("  带版本号:", jam2_id_with_version, "\n")
cat("  不带版本号:", jam2_id_no_version, "\n")

# ==================== 2. 定义13个脑组织文件 ====================

brain_tissues <- c(
  "Brain_Amygdala.v8.normalized_expression.bed",
  "Brain_Anterior_cingulate_cortex_BA24.v8.normalized_expression.bed",
  "Brain_Caudate_basal_ganglia.v8.normalized_expression.bed",
  "Brain_Cerebellar_Hemisphere.v8.normalized_expression.bed",
  "Brain_Cerebellum.v8.normalized_expression.bed",
  "Brain_Cortex.v8.normalized_expression.bed",
  "Brain_Frontal_Cortex_BA9.v8.normalized_expression.bed",
  "Brain_Hippocampus.v8.normalized_expression.bed",
  "Brain_Hypothalamus.v8.normalized_expression.bed",
  "Brain_Nucleus_accumbens_basal_ganglia.v8.normalized_expression.bed",
  "Brain_Putamen_basal_ganglia.v8.normalized_expression.bed",
  "Brain_Spinal_cord_cervical_c-1.v8.normalized_expression.bed",
  "Brain_Substantia_nigra.v8.normalized_expression.bed"
)

tissue_names <- c(
  "Amygdala",
  "A.C.cortex",
  "C.B.ganglia",
  "C.Hemisphere",
  "Cerebellum",
  "Cortex",
  "F.cortex",
  "Hippocampus",
  "Hypothalamus",
  "N.A.B.ganglia",
  "P.B.ganglia",
  "S.C.cervical",
  "S.nigra"
)

# ==================== 3. 读取并处理每个组织的基因表达数据 ====================

cat("\n读取13个脑组织数据...\n")

tissue_matrices <- list()
tissue_gene_ids_with_version <- list()
tissue_gene_ids_no_version <- list()
tissue_sample_counts <- numeric(length(brain_tissues))

for (i in seq_along(brain_tissues)) {
  file <- brain_tissues[i]
  tissue_name <- tissue_names[i]
  
  cat(sprintf("读取: %s...\n", tissue_name))
  
  # 读取.bed文件
  bed_data <- fread(file, header = TRUE)
  
  # 提取基因ID（第4列，带版本号）
  gene_ids_with_version <- bed_data[[4]]
  
  # 去掉版本号
  gene_ids_no_version <- remove_version(gene_ids_with_version)
  
  # 提取表达值矩阵
  expr_matrix <- as.matrix(bed_data[, 5:ncol(bed_data)])
  rownames(expr_matrix) <- gene_ids_no_version  # 使用不带版本号的ID作为行名
  
  # 存储数据
  tissue_matrices[[tissue_name]] <- expr_matrix
  tissue_gene_ids_with_version[[tissue_name]] <- gene_ids_with_version
  tissue_gene_ids_no_version[[tissue_name]] <- gene_ids_no_version
  tissue_sample_counts[i] <- ncol(expr_matrix)
  
  cat(sprintf("  样本数: %d, 基因数: %d\n", ncol(expr_matrix), nrow(expr_matrix)))
  
  # 显示前几个基因ID示例
  cat(sprintf("  前3个基因ID示例（带版本号）: %s\n", 
              paste(head(gene_ids_with_version, 3), collapse = ", ")))
  cat(sprintf("  前3个基因ID示例（不带版本号）: %s\n", 
              paste(head(gene_ids_no_version, 3), collapse = ", ")))
}

# ==================== 4. 诊断基因匹配问题 ====================

cat("\n========== 诊断基因匹配问题 ==========\n")

# 检查每个组织的基因ID数量
for (tissue in tissue_names) {
  cat(sprintf("%s: %d个基因\n", tissue, length(tissue_gene_ids_no_version[[tissue]])))
}

# 查找所有组织的基因交集（不带版本号）
all_genes_list <- lapply(tissue_names, function(tissue) {
  tissue_gene_ids_no_version[[tissue]]
})

common_genes_all_tissues <- Reduce(intersect, all_genes_list)
cat("\n在所有13个组织中均存在的基因数量（不带版本号）:", length(common_genes_all_tissues), "\n")

# 与我们感兴趣的基因列表取交集
common_genes <- intersect(common_genes_all_tissues, ensembl_ids_no_version)
cat("与我们感兴趣基因列表的交集数量:", length(common_genes), "\n")

# 如果没有交集，尝试其他方法
if (length(common_genes) == 0) {
  cat("\n警告：没有找到公共基因！尝试其他方法...\n")
  
  # 方法1：检查是否有部分匹配
  cat("检查部分匹配情况...\n")
  matching_counts <- sapply(ensembl_ids_no_version, function(gene) {
    sum(sapply(tissue_names, function(tissue) {
      gene %in% tissue_gene_ids_no_version[[tissue]]
    }))
  })
  
  cat("在至少一个组织中出现的基因数量:", sum(matching_counts > 0), "\n")
  cat("在至少10个组织中出现的基因数量:", sum(matching_counts >= 10), "\n")
  cat("在至少12个组织中出现的基因数量:", sum(matching_counts >= 12), "\n")
  
  # 选择在大多数组织中出现的基因
  candidate_genes <- names(matching_counts)[matching_counts >= 12]
  cat("在至少12个组织中出现的候选基因数量:", length(candidate_genes), "\n")
  
  # 确保JAM2在候选基因中
  if (!jam2_id_no_version %in% candidate_genes) {
    cat("添加JAM2到候选基因...\n")
    candidate_genes <- c(candidate_genes, jam2_id_no_version)
  }
  
  # 如果候选基因超过308个，选择前308个
  if (length(candidate_genes) > 308) {
    set.seed(123)
    common_genes <- sample(candidate_genes, 308)
  } else {
    common_genes <- candidate_genes
  }
}

cat("\n最终选择的协变量基因数量:", length(common_genes), "\n")

# 确保JAM2不在协变量中（它是响应变量）
common_genes <- setdiff(common_genes, jam2_id_no_version)
cat("排除JAM2后的协变量数量:", length(common_genes), "\n")

# 显示前10个协变量
cat("前10个协变量基因ID:\n")
print(head(common_genes, 10))

# ==================== 5. 提取每个组织的表达子矩阵 ====================

cat("\n提取每个组织的表达子矩阵...\n")

tissue_submatrices <- list()
missing_counts <- list()

for (tissue in tissue_names) {
  cat("处理组织:", tissue, "...\n")
  
  full_matrix <- tissue_matrices[[tissue]]
  
  # 找到协变量基因的行索引
  row_indices <- which(rownames(full_matrix) %in% common_genes)
  
  missing <- setdiff(common_genes, rownames(full_matrix))
  missing_counts[[tissue]] <- length(missing)
  
  cat(sprintf("  找到 %d/%d 个协变量基因，缺失 %d 个\n", 
              length(row_indices), length(common_genes), length(missing)))
  
  if (length(missing) > 0 && length(missing) <= 5) {
    cat("  缺失的基因:", paste(missing, collapse = ", "), "\n")
  }
  
  # 提取子矩阵
  if (length(row_indices) > 0) {
    sub_matrix <- full_matrix[row_indices, ]
    
    # 确保行顺序与common_genes一致
    gene_order <- match(common_genes, rownames(sub_matrix))
    gene_order <- gene_order[!is.na(gene_order)]  # 移除NA
    sub_matrix <- sub_matrix[gene_order, ]
    
    # 转置：行=样本，列=基因
    sub_matrix_t <- t(sub_matrix)
    tissue_submatrices[[tissue]] <- sub_matrix_t
    
    cat(sprintf("  完成: %d 样本 × %d 协变量\n", 
                nrow(sub_matrix_t), ncol(sub_matrix_t)))
  } else {
    cat("  警告: 没有找到任何协变量基因！\n")
    tissue_submatrices[[tissue]] <- NULL
  }
}

# ==================== 6. 提取JAM2表达值 ====================

cat("\n提取JAM2表达值...\n")

jam2_expressions <- list()

for (tissue in tissue_names) {
  full_matrix <- tissue_matrices[[tissue]]
  
  # 查找JAM2（使用不带版本号的ID）
  jam2_row <- which(rownames(full_matrix) == jam2_id_no_version)
  
  if (length(jam2_row) == 0) {
    # 如果没找到，尝试带版本号的查找
    cat(sprintf("  %s: 未找到JAM2 (ID: %s)，尝试其他方法...\n", 
                tissue, jam2_id_no_version))
    
    # 查找以JAM2主ID开头的基因
    jam2_candidates <- which(startsWith(rownames(full_matrix), jam2_id_no_version))
    
    if (length(jam2_candidates) > 0) {
      cat(sprintf("    找到可能的JAM2基因: %s\n", 
                  paste(rownames(full_matrix)[jam2_candidates], collapse = ", ")))
      jam2_expressions[[tissue]] <- as.numeric(full_matrix[jam2_candidates[1], ])
    } else {
      cat(sprintf("    ⚠ 仍未找到JAM2，使用NA填充\n"))
      jam2_expressions[[tissue]] <- rep(NA, tissue_sample_counts[which(tissue_names == tissue)])
    }
  } else {
    jam2_expressions[[tissue]] <- as.numeric(full_matrix[jam2_row, ])
    cat(sprintf("  %s: 成功提取JAM2表达值\n", tissue))
  }
}

# ==================== 7. 创建迁移学习数据结构 ====================

cat("\n创建迁移学习数据结构...\n")

# 首先检查哪些组织有有效数据
valid_tissues <- tissue_names[sapply(tissue_submatrices, function(x) !is.null(x) && ncol(x) > 0)]
cat("有效组织数量:", length(valid_tissues), "\n")
print(valid_tissues)

transfer_data <- list()

for (target_tissue in valid_tissues) {
  cat(sprintf("设置目标域: %s\n", target_tissue))
  
  # 目标域数据
  target_X <- tissue_submatrices[[target_tissue]]
  target_y <- jam2_expressions[[target_tissue]]
  
  # 确保长度一致
  if (nrow(target_X) != length(target_y)) {
    cat("警告: 样本数不一致，进行调整...\n")
    min_len <- min(nrow(target_X), length(target_y))
    target_X <- target_X[1:min_len, ]
    target_y <- target_y[1:min_len]
  }
  
  # 源域数据（其他有效组织）
  source_tissues <- setdiff(valid_tissues, target_tissue)
  source_data <- list()
  
  for (source_tissue in source_tissues) {
    source_X <- tissue_submatrices[[source_tissue]]
    source_y <- jam2_expressions[[source_tissue]]
    
    source_data[[source_tissue]] <- list(
      X = source_X,
      y = source_y,
      n_samples = nrow(source_X),
      n_covariates = ncol(source_X)
    )
  }
  
  # 存储配置
  transfer_data[[target_tissue]] <- list(
    target = list(
      X = target_X,
      y = target_y,
      tissue = target_tissue,
      n_samples = nrow(target_X),
      n_covariates = ncol(target_X)
    ),
    sources = source_data,
    covariates = colnames(target_X),  # 协变量名称
    jam2_id = jam2_id_no_version
  )
  
  cat(sprintf("  配置完成: %d样本 × %d协变量，%d个源域\n", 
              nrow(target_X), ncol(target_X), length(source_tissues)))
}

# ==================== 8. 保存结果 ====================

cat("\n保存结果...\n")

# 保存完整数据
saveRDS(transfer_data, file = "transfer_learning_brain_data_fixed.rds")

# 保存协变量列表
if (length(transfer_data) > 0) {
  first_config <- transfer_data[[1]]
  if (!is.null(first_config$covariates) && length(first_config$covariates) > 0) {
    writeLines(first_config$covariates, "covariates_final.txt")
    cat("协变量数量:", length(first_config$covariates), "\n")
  }
}

# 保存元数据
metadata <- data.frame(
  Tissue = tissue_names,
  Samples = tissue_sample_counts,
  Covariates = sapply(tissue_names, function(t) {
    if (t %in% names(tissue_submatrices) && !is.null(tissue_submatrices[[t]])) {
      ncol(tissue_submatrices[[t]])
    } else {
      0
    }
  }),
  Missing_Genes = sapply(tissue_names, function(t) {
    if (t %in% names(missing_counts)) missing_counts[[t]] else NA
  })
)
write.csv(metadata, "brain_tissues_metadata_fixed.csv", row.names = FALSE)

# ==================== 9. 输出汇总信息 ====================

cat("\n========== 提取完成！汇总信息 ==========\n")

if (length(transfer_data) > 0) {
  first_tissue <- names(transfer_data)[1]
  first_config <- transfer_data[[first_tissue]]
  
  cat(sprintf("有效配置数量: %d\n", length(transfer_data)))
  cat(sprintf("示例配置 (%s):\n", first_tissue))
  cat(sprintf("  目标域样本数: %d\n", first_config$target$n_samples))
  cat(sprintf("  协变量数量: %d\n", first_config$target$n_covariates))
  cat(sprintf("  源域数量: %d\n", length(first_config$sources)))
  
  cat("前5个协变量:\n")
  if (length(first_config$covariates) >= 5) {
    print(head(first_config$covariates, 5))
  } else {
    print(first_config$covariates)
  }
  
  cat("\n数据结构预览:\n")
  str(first_config$target, max.level = 1)
} else {
  cat("警告: 没有创建任何有效配置！\n")
}

cat("\n文件已保存:\n")
cat("  1. transfer_learning_brain_data_fixed.rds\n")
cat("  2. covariates_final.txt\n")
cat("  3. brain_tissues_metadata_fixed.csv\n")

cat("\n下一步: 使用 transfer_learning_brain_data_fixed.rds 进行迁移学习建模\n")
