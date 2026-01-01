# 设置工作目录
setwd("C:\\Users\\qing\\Desktop\\线性模型作业\\脑部基因实例分析")

# 读取
text_lines <- readLines("137.txt")

# 提取以数字开头（代表LocusLink）的行，并分割
gene_data_lines <- grep("^\\d", text_lines, value = TRUE)

# 分割每行，提取第二列（Name）
gene_names <- sapply(strsplit(gene_data_lines, "\\s+"), function(x) x[2])

# 清理并去除可能的空值
gene_names <- gene_names[!is.na(gene_names) & gene_names != ""]

cat("成功提取基因数量:", length(gene_names), "\n")
cat("前10个基因:", head(gene_names, 10), "\n")

# 保存为文件
writeLines(gene_names, "cns_gene_symbols_546.txt")


# convert_ids.R
# 将基因符号转换为ENSEMBL ID

# 解决biomaRt查询错误的问题
library(biomaRt)

# 加载基因符号
gene_symbols <- readLines("cns_gene_symbols_546.txt")
cat("总共基因数:", length(gene_symbols), "\n")

# 连接到Ensembl（使用不同的镜像，如果第一个失败）
cat("连接到Ensembl数据库...\n")

# 尝试不同的方式连接
try_connect <- function() {
  tryCatch({
    # 方法1：直接连接
    ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
    cat("✓ 方法1成功\n")
    return(ensembl)
  }, error = function(e) {
    cat("方法1失败:", e$message, "\n")
    
    # 方法2：使用useEnsembl函数
    tryCatch({
      ensembl <- useEnsembl(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")
      cat("✓ 方法2成功\n")
      return(ensembl)
    }, error = function(e2) {
      cat("方法2失败:", e2$message, "\n")
      
      # 方法3：使用ensembl镜像
      tryCatch({
        ensembl <- useMart(biomart = "ENSEMBL_MART_ENSEMBL", 
                           dataset = "hsapiens_gene_ensembl",
                           host = "https://asia.ensembl.org")  # 亚洲镜像
        cat("✓ 方法3成功（使用亚洲镜像）\n")
        return(ensembl)
      }, error = function(e3) {
        cat("所有连接方法都失败\n")
        return(NULL)
      })
    })
  })
}

ensembl <- try_connect()
if (is.null(ensembl)) stop("无法连接到Ensembl数据库")

# 清理基因符号（移除可能的问题字符）
clean_genes <- function(genes) {
  # 移除空值
  genes <- genes[!is.na(genes) & genes != ""]
  # 移除非字母数字开头的基因
  genes <- genes[grepl("^[A-Za-z0-9]", genes)]
  # 去重
  unique(genes)
}

gene_symbols_clean <- clean_genes(gene_symbols)
cat("清理后基因数:", length(gene_symbols_clean), "\n")

# 显示前20个基因
cat("前20个基因:", paste(head(gene_symbols_clean, 20), collapse = ", "), "\n")

# 分批查询函数（解决"invalid token"错误）
batch_query_genes <- function(gene_list, mart, batch_size = 50) {
  total <- length(gene_list)
  batches <- split(gene_list, ceiling(seq_along(gene_list)/batch_size))
  
  all_results <- data.frame()
  
  for (i in seq_along(batches)) {
    batch <- batches[[i]]
    cat(sprintf("查询批次 %d/%d (%d个基因)... ", i, length(batches), length(batch)))
    
    tryCatch({
      result <- getBM(
        attributes = c("ensembl_gene_id", "external_gene_name", "ensembl_gene_id_version"),
        filters = "external_gene_name", 
        values = batch,
        mart = mart
      )
      
      if (nrow(result) > 0) {
        all_results <- rbind(all_results, result)
        cat(sprintf("成功获取 %d 个映射\n", nrow(result)))
      } else {
        cat("未找到映射\n")
      }
      
      # 小延迟避免请求过快
      if (i %% 5 == 0) Sys.sleep(1)
      
    }, error = function(e) {
      cat(sprintf("错误: %s\n", e$message))
      
      # 如果批量失败，尝试单个查询
      if (length(batch) > 1) {
        cat("尝试逐个查询...\n")
        for (gene in batch) {
          tryCatch({
            single_result <- getBM(
              attributes = c("ensembl_gene_id", "external_gene_name", "ensembl_gene_id_version"),
              filters = "external_gene_name", 
              values = gene,
              mart = mart
            )
            if (nrow(single_result) > 0) {
              all_results <- rbind(all_results, single_result)
            }
          }, error = function(e2) {
            # 单个基因也失败，跳过
          })
        }
      }
    })
  }
  
  return(all_results)
}

# 执行分批查询
cat("\n开始分批查询基因映射...\n")
gene_mapping <- batch_query_genes(gene_symbols_clean, ensembl, batch_size = 30)

# 统计结果
cat("\n=== 查询结果统计 ===\n")
cat("输入基因数:", length(gene_symbols_clean), "\n")
cat("成功映射数:", nrow(gene_mapping), "\n")
cat("映射成功率:", round(nrow(gene_mapping)/length(gene_symbols_clean)*100, 1), "%\n")

# 检查JAM2
cat("\n=== 检查JAM2基因 ===\n")
jam2_result <- gene_mapping[gene_mapping$external_gene_name == "JAM2", ]
if (nrow(jam2_result) > 0) {
  cat("✓ 找到JAM2映射:\n")
  print(jam2_result)
  jam2_id <- jam2_result$ensembl_gene_id_version[1]
} else {
  cat("✗ 未找到JAM2，尝试直接查询...\n")
  
  # 直接查询JAM2
  jam2_direct <- tryCatch({
    getBM(
      attributes = c("ensembl_gene_id", "external_gene_name", "ensembl_gene_id_version"),
      filters = "external_gene_name", 
      values = "JAM2",
      mart = ensembl
    )
  }, error = function(e) {
    cat("直接查询失败:", e$message, "\n")
    NULL
  })
  
  if (!is.null(jam2_direct) && nrow(jam2_direct) > 0) {
    cat("✓ 直接查询找到JAM2:\n")
    print(jam2_direct)
    gene_mapping <- rbind(gene_mapping, jam2_direct)
    jam2_id <- jam2_direct$ensembl_gene_id_version[1]
  } else {
    cat("⚠ 使用已知的JAM2 ID\n")
    jam2_id <- "ENSG00000154721.14"  # 论文中的ID
    # 添加到映射表
    gene_mapping <- rbind(gene_mapping, 
                          data.frame(
                            ensembl_gene_id = "ENSG00000154721",
                            external_gene_name = "JAM2",
                            ensembl_gene_id_version = jam2_id
                          ))
  }
}

cat("\nJAM2基因ID:", jam2_id, "\n")

# 保存结果
write.csv(gene_mapping, "cns_gene_ensembl_mapping_final.csv", row.names = FALSE)

# 提取所有ENSEMBL ID
ensembl_ids <- unique(gene_mapping$ensembl_gene_id_version)
writeLines(ensembl_ids, "cns_ensembl_ids_final.txt")

cat("\n=== 保存完成 ===\n")
cat("1. 完整映射表: cns_gene_ensembl_mapping_final.csv\n")
cat("2. ENSEMBL ID列表: cns_ensembl_ids_final.txt\n")
cat("3. JAM2基因ID:", jam2_id, "\n")
cat("\n基因数统计:\n")
cat("  输入: ", length(gene_symbols), "\n")
cat("  清理后: ", length(gene_symbols_clean), "\n")
cat("  成功映射: ", length(ensembl_ids), "\n")

# 查看未映射的基因（如果有的话）
missing_genes <- setdiff(gene_symbols_clean, gene_mapping$external_gene_name)
if (length(missing_genes) > 0) {
  cat("\n未映射的基因 (", length(missing_genes), "个):\n")
  if (length(missing_genes) <= 20) {
    print(missing_genes)
  } else {
    cat("前20个:", paste(head(missing_genes, 20), collapse = ", "), "\n")
  }
}
