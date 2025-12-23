library(rms)
library(dplyr)
library(purrr)

detect_contaminant_relationships <- function(rawdata, contamination_level, nk = 5) {
  
  # 1. 数据对齐与预处理
  # 确保样本（rawdata的列）和污染指标（contamination_level的行）匹配
  common_samples <- intersect(colnames(rawdata), rownames(contamination_level))
  if (length(common_samples) == 0) stop("Error: 样本名称不匹配，请检查 rawdata 的列名和 contamination_level 的行名。")
  
  # 转置蛋白数据
  dt_prot <- as.data.frame(t(rawdata[, common_samples, drop = FALSE]))
  dt_contam <- as.data.frame(contamination_level[common_samples, , drop = FALSE])
  
  # 合并数据
  df_combined <- cbind(dt_prot, dt_contam)
  
  # 获取蛋白名和污染指标名
  protein_names <- rownames(rawdata)
  exposure_types <- colnames(dt_contam)
  
  # 2. 设置 rms 环境 (核心：必须在函数内定义并赋值)
  # 使用 options(datadist = "dd") 确保 ols 能找到分布信息
  dd <<- datadist(df_combined)
  options(datadist = "dd")
  
  results_list <- list()
  
  # 3. 循环计算
  for (exp_type in exposure_types) {
    message(paste("Checking relationship for:", exp_type))
    
    type_results <- map_df(protein_names, function(prot) {
      # 使用反引号包裹变量名以防特殊字符
      # 使用 rcs() 构建限制性立方样条
      formula_str <- paste0("log2(`", prot, "` + 1) ~ rcs(`", exp_type, "`, nk = ", nk, ")")
      
      tryCatch({
        # 拟合模型
        fit <- ols(as.formula(formula_str), data = df_combined)
        an <- anova(fit)
        
        # 提取 P 值
        # rms 的 anova 表第一行通常是变量总效应，含有 'Nonlinear' 的行是分线性效应
        p_overall <- an[exp_type, "P"]
        
        # 查找非线性行（模糊匹配 "Nonlinear"）
        nl_row <- grep("Nonlinear", rownames(an), ignore.case = TRUE)
        p_nonlinear <- if(length(nl_row) > 0) an[nl_row[1], "P"] else NA
        
        # 逻辑判断
        rel <- "None"
        if (!is.na(p_overall) && p_overall < 0.05) {
          if (!is.na(p_nonlinear) && p_nonlinear < 0.05) {
            rel <- "Nonlinear"
          } else {
            rel <- "Linear"
          }
        }
        
        data.frame(
          Protein = prot,
          Contaminant = exp_type,
          P_overall = p_overall,
          P_nonlinear = p_nonlinear,
          Relationship = rel,
          stringsAsFactors = FALSE
        )
      }, error = function(e) {
        # 如果报错，记录失败原因（可选）
        data.frame(
          Protein = prot, Contaminant = exp_type, 
          P_overall = NA, P_nonlinear = NA, Relationship = "Failed"
        )
      })
    })
    results_list[[exp_type]] <- type_results
  }
  
  return(bind_rows(results_list))
}
