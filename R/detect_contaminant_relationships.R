#' Detect Relationship Type Between Proteins and Contaminants
#'
#' Evaluates whether the relationship between protein expression and contamination 
#' levels is Linear, Nonlinear, or non-existent using Restricted Cubic Splines (RCS).
#'
#' @param rawdata A matrix or data frame of protein expression (Proteins as rows, Samples as columns).
#' @param contamination_level A data frame of contamination scores (Samples as rows, Scores as columns).
#' @param nk Number of knots for RCS. For small samples (N < 15), nk=3 is recommended. Default is 5.
#'
#' @return A data frame with columns:
#' \item{Protein}{Name of the protein}
#' \item{Contaminant}{Name of the contaminant (e.g., x_eryth)}
#' \item{P_overall}{P-value for the total effect of the contaminant}
#' \item{P_nonlinear}{P-value specifically for the nonlinear (spline) component}
#' \item{Relationship}{Classification: "Nonlinear", "Linear", or "None"}
#'
#' @details 
#' The function fits an OLS model: \code{log2(Protein + 1) ~ rcs(Contaminant, nk)}.
#' It uses an ANOVA to partition the variance. If the overall P < 0.05, it checks the 
#' nonlinear component. If P_nonlinear < 0.05, the relationship is "Nonlinear", 
#' otherwise "Linear".
#' 
#' @importFrom rms ols rcs datadist anova.rms
#' @importFrom purrr map_df
#' @importFrom dplyr bind_rows
#' @export
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
  df_combined <- cbind(dt_prot, dt_contam)
  
  protein_names <- rownames(rawdata)
  exposure_types <- colnames(dt_contam)
  # 2. 设置 rms 环境 (局部设置以防止污染全局环境)
  # old_opts <- options(datadist = NULL)
  # on.exit(options(old_opts)) # 函数结束时还原设置
  # dd <<- datadist(df_combined)
  # options(datadist = "dd")
  old_dd_opt <- getOption("datadist")
  
  dd_name <- paste0(".dd_tmp_", paste(sample(c(letters, LETTERS, 0:9), 12, TRUE), collapse=""))
  assign(dd_name, rms::datadist(df_combined), envir = .GlobalEnv)
  options(datadist = dd_name)
  
  on.exit({
    options(datadist = old_dd_opt)
    if (exists(dd_name, envir = .GlobalEnv, inherits = FALSE)) {
      rm(list = dd_name, envir = .GlobalEnv)
    }
  }, add = TRUE)
  results_list <- list()
  
  # 3. 循环计算
  for (exp_type in exposure_types) {
    message(paste("Checking relationship for:", exp_type))
    # 针对 N=10 的保护：如果样本量太小，强制下调 nk
    actual_nk <- if(nrow(df_combined) < (nk + 2)) 3 else nk
    type_results <- purrr::map_df(protein_names, function(prot) {
      # 使用反引号包裹变量名以防特殊字符
      # 使用 rcs() 构建限制性立方样条
      formula_str <- paste0("log2(`", prot, "` + 1) ~ rcs(`", exp_type, "`, nk = ", actual_nk, ")")
      
      tryCatch({
        # 拟合模型
        fit <- rms::ols(as.formula(formula_str), data = df_combined)
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
  
  return(dplyr::bind_rows(results_list))
}
