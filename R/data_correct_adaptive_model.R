# =========================================================
# 1) 通用矫正函数：核心逻辑封装
# =========================================================
.apply_contamination_correction <- function(rawdata,
                                            base_data,
                                            contaminant_relationships,
                                            type_vars,
                                            constraint_global = 1,
                                            constraint_map = list(),
                                            nk = 3,
                                            min_non_na = 3) {
  rawdata   <- as.data.frame(rawdata)
  base_data <- as.data.frame(base_data)
  
  # 0. 基础检查
  if (length(type_vars) == 0) {
    return(log2(as.matrix(rawdata) + 1))
  }
  
  missing_cols <- setdiff(type_vars, colnames(base_data))
  if (length(missing_cols) > 0) {
    stop("base_data is missing required columns: ", paste(missing_cols, collapse = ", "))
  }
  
  # 1. 对齐样本顺序
  if (!is.null(colnames(rawdata))) {
    if (is.null(rownames(base_data))) {
      if (nrow(base_data) != ncol(rawdata)) stop("Dimension mismatch between base_data and rawdata")
      rownames(base_data) <- colnames(rawdata)
    } else {
      common_samples <- intersect(colnames(rawdata), rownames(base_data))
      if (length(common_samples) != ncol(rawdata)) stop("Sample names mismatch")
      base_data <- base_data[colnames(rawdata), , drop = FALSE]
    }
  } else {
    if (nrow(base_data) != ncol(rawdata)) stop("Dimension mismatch")
  }
  
  # 2. 确保 rms 包加载 (不仅是安装，需要 Attach 以便 predict.rms 方法生效)
  if (!require("rms", quietly = TRUE)) {
    stop("Package 'rms' is required. Please install and load it.")
  }
  
  # 3. datadist 环境黑魔法 (关键修复)
  # rms 需要 datadist 在全局环境中可见。这里创建一个临时变量并指向它。
  old_dd_opt <- getOption("datadist")
  # 使用随机名防止冲突
  dd_name <- paste0(".dd_tmp_", paste(sample(letters, 10), collapse="")) 
  
  # 将 datadist 对象放入全局环境
  assign(dd_name, rms::datadist(base_data), envir = .GlobalEnv)
  options(datadist = dd_name)
  
  # 函数退出时清理全局环境，不留垃圾
  on.exit({
    options(datadist = old_dd_opt)
    if (exists(dd_name, envir = .GlobalEnv)) {
      rm(list = dd_name, envir = .GlobalEnv)
    }
  }, add = TRUE)
  
  # 4. 初始化输出矩阵
  ndata1 <- matrix(NA_real_, nrow = nrow(rawdata), ncol = ncol(rawdata),
                   dimnames = list(rownames(rawdata), colnames(rawdata)))
  
  # 辅助函数：安全获取约束因子
  get_constraint_val <- function(map, var, prot) {
    if (is.null(map) || is.null(map[[var]])) return(1)
    val_vec <- map[[var]]
    if (!prot %in% names(val_vec)) return(1)
    val <- val_vec[[prot]]
    if (is.na(val) || !is.finite(val)) return(1)
    return(max(0, min(1, val)))
  }
  
  # 5. 逐蛋白迭代矫正
  for (i in seq_len(nrow(rawdata))) {
    y_name <- rownames(rawdata)[i]
    y_vec  <- as.numeric(rawdata[i, ]) # 原始表达量
    
    # 构造当前数据框
    curr_data <- base_data
    curr_data$y_log <- log2(y_vec + 1)
    
    # 获取该蛋白的关系列表
    c_sub <- contaminant_relationships[contaminant_relationships$Protein == y_name, , drop = FALSE]
    
    # 构建模型公式
    formula_terms <- character(0)
    active_vars   <- character(0)
    
    # 存储该蛋白对应的各污染源的具体约束值 (Constraint Factor)
    current_constraints <- setNames(numeric(length(type_vars)), type_vars)
    
    for (v in type_vars) {
      # 检查该污染分在当前是否可用 (非全NA)
      x <- curr_data[[v]]
      x_non_na <- stats::na.omit(x)
      
      # 获取约束值
      current_constraints[[v]] <- get_constraint_val(constraint_map, v, y_name)
      
      # 忽略几乎没有信息的变量
      if (length(x_non_na) < min_non_na || length(unique(x_non_na)) < 2) next
      
      # 获取关系类型
      rel_val <- c_sub[c_sub$Contaminant == v, "Relationship"]
      v_type <- if (length(rel_val) > 0 && !is.na(rel_val)) as.character(rel_val[[1]]) else "None"
      
      if (v_type == "Linear") {
        formula_terms <- c(formula_terms, v)
        active_vars   <- c(active_vars, v)
      } else if (v_type == "Nonlinear") {
        # 如果点太少，rcs 会报错，降级为线性
        if (length(unique(x_non_na)) <= nk) {
          formula_terms <- c(formula_terms, v)
        } else {
          formula_terms <- c(formula_terms, paste0("rms::rcs(", v, ", nk=", nk, ")"))
        }
        active_vars <- c(active_vars, v)
      }
    }
    
    # 如果没有需要矫正的项，保持原样
    if (length(formula_terms) == 0) {
      ndata1[i, ] <- curr_data$y_log
      next
    }
    
    # 拟合模型 (优先尝试带 RCS 的模型)
    f_str <- paste("y_log ~", paste(formula_terms, collapse = " + "))
    # 注意：这里虽然写了 rms::ols，但为了让 formula 中的 rcs 生效，最好前面已 require(rms)
    fit <- try(rms::ols(stats::as.formula(f_str), data = curr_data, x = TRUE, y = TRUE), silent = TRUE)
    
    # 失败回退：尝试纯线性
    if (inherits(fit, "try-error")) {
      active_vars_u <- unique(active_vars)
      if (length(active_vars_u) == 0) {
        ndata1[i, ] <- curr_data$y_log
        next
      }
      f_linear <- paste("y_log ~", paste(active_vars_u, collapse = " + "))
      fit <- try(rms::ols(stats::as.formula(f_linear), data = curr_data, x = TRUE, y = TRUE), silent = TRUE)
    }
    
    if (inherits(fit, "try-error")) {
      ndata1[i, ] <- curr_data$y_log
      next
    }
    
    # 预测各项效应 (Terms)
    # 显式使用 stats::predict，只要 rms 被 attach，S3 会自动分发给 predict.rms
    term_effects <- try(stats::predict(fit, newdata = curr_data, type = "terms"), silent = TRUE)
    
    if (inherits(term_effects, "try-error")) {
      ndata1[i, ] <- curr_data$y_log
      next
    }
    
    term_effects <- as.matrix(term_effects)
    total_effect <- rep(0, nrow(curr_data))
    
    # 计算总矫正量：Sum(Effect * Constraint)
    for (v in unique(active_vars)) {
      # 匹配列名（rcs项会变成 v, v' 等多列，grep 可以匹配到）
      col_idx <- grep(v, colnames(term_effects), fixed = TRUE)
      if (length(col_idx) > 0) {
        raw_eff <- rowSums(term_effects[, col_idx, drop = FALSE])
        total_effect <- total_effect + raw_eff * current_constraints[[v]]
      }
    }
    
    # 应用矫正：原始Log值 - (全局强度 * 计算出的污染效应)
    ndata1[i, ] <- curr_data$y_log - (constraint_global * total_effect)
  }
  
  return(ndata1)
}


# =========================================================
# 2) 主函数：保留原接口，调用通用逻辑
# =========================================================
data_correct <- function(data,
                         type = c("coagulation", "erythrocyte", "platelet"),
                         constraint = 1,
                         erythrocyte_marker,
                         constraint_erythrocyte = 0.95,
                         coagulation_marker,
                         constraint_coagulation = 0.95,
                         platelet_marker,
                         constraint_platelet = 0.95) {
  
  # 1. 参数清洗
  if (missing(erythrocyte_marker) || is.null(erythrocyte_marker)) erythrocyte_marker <- character(0)
  if (missing(platelet_marker)    || is.null(platelet_marker))    platelet_marker    <- character(0)
  if (missing(coagulation_marker) || is.null(coagulation_marker)) coagulation_marker <- character(0)
  
  # data为输入数据，type为矫正类型，默认为c("coagulation","erythrocyte","platelet")，表示矫正所有污染类型
  # type可以为coagulation, erythrocyte, platelet中任意几种
  # constraint_*:参数为计算约束因子的最小分位数，默认为相关系数从小到大排序，计算后95%-100%相关系数的平均数
  # 强制constraint_*参数在(0,1)范围内（包含边界检查）
  constraint_erythrocyte_original <- max(0, min(1, constraint_erythrocyte))
  constraint_platelet_original <- max(0, min(1, constraint_platelet))
  constraint_coagulation_original <- max(0, min(1, constraint_coagulation)) 
  # 当参数被修正时发出警告
  if(constraint_erythrocyte != constraint_erythrocyte_original) 
    warning("constraint_erythrocyte clamped to [0,1]")
  if(constraint_platelet != constraint_platelet_original) 
    warning("constraint_platelet clamped to [0,1]")
  if(constraint_coagulation != constraint_coagulation_original) 
    warning("constraint_coagulation clamped to [0,1]")
  library(MASS)
  # 定义计算均值的函数
  calculation_conta <- function(data) {
    for_mean <- function(data) {
      # 标准化数据（不再需要转置）
      data <- log2(data + 1)
      # 直接使用原始基因 × 样本的结构
      M <- as.data.frame(data)
      # 计算行（基因）均值并中心化
      a <- M - rowMeans(M)
      # 对列（样本）取均值
      mean <- apply(a,2, mean, na.rm = TRUE)
      # if (min(mean) < 0 ){
      #   mean <- mean - min(mean)
      # }
      return(mean)
    }
    list_erythrocyte <-  if (any(rownames(data$rawdata) %in% erythrocyte_marker)) {
      for_mean(data$rawdata[rownames(data$rawdata) %in% erythrocyte_marker,,drop = FALSE])
    } else {
      rep(NA,dim(data$rawdata)[2])
    }
    list_platelet <- if (any(rownames(data$rawdata) %in% platelet_marker)) {
      for_mean(data$rawdata[rownames(data$rawdata) %in% platelet_marker,, drop = FALSE])
    } else {
      rep(NA,dim(data$rawdata)[2])
    }
    list_coagulation <- if (any(rownames(data$rawdata) %in% coagulation_marker)) {
      for_mean(data$rawdata[rownames(data$rawdata) %in% coagulation_marker,,drop = FALSE])
    } else {
      rep(NA,dim(data$rawdata)[2])
    }
    
    components  <- data.frame(
      x_eryth = list_erythrocyte,
      x_plate = list_platelet,
      x_coagu = list_coagulation
    )
    return(components)
  }
  
  # 向量化计算约束因子
  compute_constraint_factors <- function(cor_matrix, marker_set, constraint_value) {
    if (is.null(marker_set) || length(marker_set) == 0) return(NULL)
    
    # 提取当前标记基因的相关系数子集
    sub_cor <- abs(cor_matrix[, marker_set, drop = FALSE])
    
    # 计算每个基因的约束因子
    constraint_factors <- apply(sub_cor, 1, function(x) {
      n <- length(x)
      if (n < 3) return(mean(x, na.rm = TRUE))
      
      # 计算起始位置
      start_idx <- max(1, min(round(n * constraint_value), n - 2))
      # 排序并取平均值
      sorted_vals <- sort(x)
      mean(sorted_vals[start_idx:(n - 1)], na.rm = TRUE)
    })
    
    return(constraint_factors)
  }
  # 计算污染水平
  smpl2 <- calculation_conta(data)
  rawdata <- as.data.frame(data$rawdata)
  ndata1 <- matrix(nrow = nrow(rawdata), ncol = ncol(rawdata))
  rownames(ndata1) <- rownames(rawdata)
  colnames(ndata1) <- colnames(rawdata)
  b <- matrix(nrow = nrow(rawdata), ncol = ncol(smpl2))
  rownames(b) <- rownames(rawdata)
  colnames(b) <- colnames(smpl2)
  
  # 计算约束系数 ----
  result_cor <- cor(as.matrix(t(data$rawdata)), method = "pearson")
  result_cor <- as.data.frame(result_cor)
  # 初始化约束因子为 NULL
  result_cor_erythrocyte <- result_cor_platelet <- result_cor_coagulation <- NULL
  # 计算约束因子
  result_cor_erythrocyte <- if (length(erythrocyte_marker) > 0) {
    compute_constraint_factors(cor_matrix = result_cor, 
                               erythrocyte_marker, 
                               constraint_erythrocyte_original)
  } else NULL
  
  result_cor_platelet <- if (length(platelet_marker) > 0) {
    compute_constraint_factors(cor_matrix = result_cor, 
                               platelet_marker, 
                               constraint_platelet_original)
  } else NULL
  
  result_cor_coagulation <- if (length(coagulation_marker) > 0) {
    compute_constraint_factors(cor_matrix = result_cor, 
                               coagulation_marker,
                               constraint_coagulation_original)
  } else NULL
  constraint_map <- list(
    x_eryth = result_cor_erythrocyte,
    x_plate = result_cor_platelet,
    x_coagu = result_cor_coagulation
  )
  
  # 将预测变量提取出来，改名为代码中统一使用的名字
  base_data <- data.frame(
    x_eryth = smpl2$x_eryth,
    x_plate = smpl2$x_plate,
    x_coagu = smpl2$x_coagu,
    row.names = rownames(smpl2)
  )
  
  
  # 4. 检测关系 (External Function)
  # 假设 detect_contaminant_relationships 存在于环境中
  # ---- 先根据 type 选出本次要用的污染列（非常关键） ----
  type_map <- c(erythrocyte = "x_eryth", platelet = "x_plate", coagulation = "x_coagu")
  target_types <- intersect(unique(type), names(type_map))
  if (length(target_types) == 0) {
    warning("No valid correction types specified. Returning original data.")
    data$correct_data <- data$rawdata
    data$contamination_level <- smpl2
    return(data)
  }
  type_vars <- unname(type_map[target_types])
  
  # ---- 只保留本次要校正的污染列，避免把全 NA 的 x_plate 传进去 ----
  base_model_data <- base_data[, type_vars, drop = FALSE]
  
  # ---- 再过滤掉信息不足的污染列（<2非NA 或 unique<2），避免 datadist 报错 ----
  keep <- vapply(base_model_data, function(x) {
    sum(!is.na(x)) >= 2 && length(unique(stats::na.omit(x))) >= 2
  }, logical(1))
  
  if (any(!keep)) {
    warning("Dropping contaminants due to insufficient data: ",
            paste(names(keep)[!keep], collapse = ", "))
  }
  
  base_model_data <- base_model_data[, keep, drop = FALSE]
  type_vars <- type_vars[keep]
  
  # 如果过滤完没有可用污染列，就不校正，直接返回
  if (length(type_vars) == 0) {
    warning("No usable contaminant scores available. Returning original data.")
    data$correct_data <- data$rawdata
    data$contamination_level <- smpl2
    return(data)
  }
  
  # ---- 只对本次需要的污染列做关系检测 ----
  contaminant_relationships <- detect_contaminant_relationships(
    rawdata = data$rawdata,
    contamination_level = base_model_data,
    nk = 3
  )
  
  # 5. 准备参数并调用通用矫正
  if (length(target_types) == 0) {
    warning("No valid correction types specified. Returning original data.")
    return(data)
  }
  ndata1_log <- .apply_contamination_correction(
    rawdata = data$rawdata,
    base_data = base_model_data,
    contaminant_relationships = contaminant_relationships,
    type_vars = type_vars,
    constraint_global = constraint,
    constraint_map = constraint_map,
    nk = 3
  )
  
  # 6. 还原为线性空间并返回
  data$correct_data <- 2 ^ ndata1_log - 1
  data$contamination_level <- smpl2
  return(data)
}
