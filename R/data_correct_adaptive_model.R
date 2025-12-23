#' Correct Proteomics Data for Contamination Effects
#'
#' Performs contamination correction on proteomics data using robust linear regression and constraint factors.
#' Supports correction for erythrocyte, platelet, and coagulation contamination either individually or in combination.
#'
#' @param data List containing raw expression data (from `data_check()` output) with components:
#'   - `rawdata`: Expression matrix (proteins as rows, samples as columns)
#' @param type Vector of contamination types to correct for. Options: 
#'   "erythrocyte", "platelet", "coagulation". Can specify multiple types.
#' @param constraint Global constraint factor to control correction strength (default: 1)
#' @param erythrocyte_marker Character vector of erythrocyte marker genes
#' @param constraint_erythrocyte Quantile threshold for erythrocyte constraint factor calculation (0-1, default: 0.95)
#' @param coagulation_marker Character vector of coagulation marker genes
#' @param constraint_coagulation Quantile threshold for coagulation constraint factor calculation (0-1, default: 0.95)
#' @param platelet_marker Character vector of platelet marker genes
#' @param constraint_platelet Quantile threshold for platelet constraint factor calculation (0-1, default: 0.95)
#'
#' @return Updated list containing:
#'   - `rawdata`: Original expression matrix
#'   - `correct_data`: Corrected expression matrix
#'   - `contamination_level`: Contamination score matrix (erythrocyte, platelet, coagulation)
#'   - Other components from input `data`
#'
#' @details
#' Correction workflow:
#' 1. Validates constraint parameters (0-1 range)
#' 2. Computes contamination scores for each type:
#'    - Log2 transforms expression
#'    - Centers marker expression by row means
#'    - Averages centered values by column (sample)
#' 3. Computes protein-wise constraint factors:
#'    - Calculates Pearson correlation matrix
#'    - For each protein, takes mean of top (constraint_*) fraction of correlations to markers
#' 4. Performs robust linear regression (MASS::rlm):
#'    - Response: Protein expression (log2)
#'    - Predictors: Contamination scores
#'    - Incorporates constraint factors as weights
#' 5. Calculates residuals and transforms back to linear space
#' 
#' Mathematical representation:
#' \deqn{\text{Corrected} = 2^{(\text{log2}(y + 1) - \beta \cdot \text{Constraint} \cdot X)} - 1}
#' Where:
#'   - y = original expression
#'   - β = regression coefficient
#'   - Constraint = constraint factor vector
#'   - X = contamination scores
#'
#' @note
#' - Requires MASS package
#' - Handles missing markers gracefully (returns NA for missing components)
#' - Uses robust regression to handle outliers
#' - Constraint factors range 0-1 (1 = full correction)
#'
#' @examples
#' \dontrun{
#' # After running data_check()
#' corrected_data <- data_correct(
#'   data = qc_result,
#'   type = c("erythrocyte", "platelet"),
#'   constraint = 0.8,
#'   erythrocyte_marker = c("HBA1", "HBB"),
#'   platelet_marker = c("PF4", "PPBP"),
#'   constraint_erythrocyte = 0.9,
#'   constraint_platelet = 0.85
#' )
#' 
#' # Access corrected data
#' corrected_matrix <- corrected_data$correct_data
#' }
#'
#' @importFrom MASS rlm
#' @importFrom stats cor
#' @export
data_correct <- function(data, 
                         type = c("coagulation","erythrocyte","platelet"),constraint = 1,
                         erythrocyte_marker,constraint_erythrocyte = 0.95,
                         coagulation_marker,constraint_coagulation = 0.95,
                         platelet_marker,constraint_platelet = 0.95) {
  # 函数说明 ----
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
      if (min(mean) < 0 ){
        mean <- mean - min(mean)
      }
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
      erythrocyte = list_erythrocyte,
      platelet = list_platelet,
      coagulation = list_coagulation
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
  # 计算线性关系 ----
  contaminant_relationships <- detect_contaminant_relationships(rawdata = data$rawdata,
                                                                contamination_level = smpl2)
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
  # fit <-ols(log2(y+1)~rcs(x_plate,nk=5)+ x_eryth + x_coagu,data=rawdata, x=TRUE,y=TRUE)
  # anova(fit)
  # 将预测变量提取出来，改名为代码中统一使用的名字
  base_data <- data.frame(
    x_eryth = smpl2$erythrocyte,
    x_plate = smpl2$platelet,
    x_coagu = smpl2$coagulation,
    row.names = rownames(smpl2)
  )
  dd <- datadist(base_data)
  options(datadist = "dd")
  var_names <- c("erythrocyte", "platelet", "coagulation")
  x_vars    <- c("x_eryth", "x_plate", "x_coagu")
  # correct ----
  if (setequal(type, c("coagulation","erythrocyte","platelet"))) {
    ## all ----
    for (i in 1:nrow(ndata1)) {
      y_vec <- as.numeric(rawdata[i, ])
      # 准备当前循环的数据框：包含 Y 和 所有的 X
      curr_data <- base_data
      curr_data$y_log <- log2(y_vec + 1)
      # 准备约束因子
      y_name <- rownames(rawdata)[i] # 获取当前基因名
      # 安全获取约束因子 (处理可能出现的 NA 或 缺失)
      get_constraint <- function(res_cor, g_name) {
        if (!is.null(res_cor) && g_name %in% names(res_cor)) res_cor[g_name] else 1
      }
      # 将约束因子放入列表，方便索引
      constraints_list <- list(
        erythrocyte = get_constraint(result_cor_erythrocyte, y_name),
        platelet    = get_constraint(result_cor_platelet, y_name),
        coagulation = get_constraint(result_cor_coagulation, y_name)
      )
      # 读取线性关系决定方程式
      c_sub <- contaminant_relationships[contaminant_relationships$Protein == y_name, ]
      # 提取关系并命名，例如 c(erythrocyte="Linear", platelet="None", ...)
      rels <- sapply(var_names, function(v) {
        val <- c_sub[c_sub$Contaminant == v, "Relationship"]
        if(length(val) == 0) return("None") else return(as.character(val))
      })
      # --- 2. 动态构建回归公式 ---
      formula_terms <- c()
      active_vars   <- c() # 记录哪些变量参与了模型（非 None）
      
      for (k in seq_along(var_names)) {
        v_type <- rels[var_names[k]]
        v_x    <- x_vars[k]
        
        if (v_type == "Linear") {
          formula_terms <- c(formula_terms, v_x)
          active_vars   <- c(active_vars, var_names[k])
        } else if (v_type == "Nonlinear") {
          formula_terms <- c(formula_terms, paste0("rcs(", v_x, ", nk=5)"))
          active_vars   <- c(active_vars, var_names[k])
        }
      }
      # 如果所有关系都是 None，保留原数据
      if (length(formula_terms) == 0) {
        ndata1[i, ] <- curr_data$y_log
        next 
      }
      # 组合公式：注意左边使用 curr_data 里的 y_log
      f_str <- paste("y_log ~", paste(formula_terms, collapse = " + "))
      
      # --- 4. 模型拟合与校正 ---
      
      # --- 3. 模型拟合 ---
      is_all_linear <- all(rels[active_vars] == "Linear")
      
      intercept_val <- 0
      total_effect  <- rep(0, nrow(curr_data)) # 初始化效应向量
      # 只有当全部是 Linear 且没有 Nonlinear 时才尝试 rlm (如果你想严格复现之前的逻辑)
      # 或者简单判断：只要不含 rcs 就用 rlm？这里为了稳健，若包含 Nonlinear 或 None 导致项减少，建议用 ols
      # 【分支 A】: 纯线性模型 (使用 rlm)
      if (is_all_linear && length(formula_terms) == 3) { # 假设3个都在且都是Linear才用rlm
        
        fit <- try(rlm(as.formula(f_str), data = curr_data, maxit = 30), silent = TRUE)
        
        if (inherits(fit, "try-error")) {
          # 如果 rlm 失败，回退到 ols
          fit <- ols(as.formula(f_str), data = curr_data, x=TRUE, y=TRUE)
          use_rlm <- FALSE
        } else {
          use_rlm <- TRUE
        }
        if (use_rlm) {
          coefs <- coef(fit)
          intercept_val <- coefs["(Intercept)"]
          
          # 记录系数到 b 矩阵 (如果需要)
          if (exists("b")) {
            # 注意：这里需要根据 names(coefs) 来匹配 b 的列，简单起见略过
          }
          
          # 手动计算效应
          for (v_name in active_vars) {
            v_x <- x_vars[which(var_names == v_name)]
            # rlm 的系数名直接就是变量名
            if (v_x %in% names(coefs)) {
              eff <- curr_data[[v_x]] * coefs[v_x] * constraints_list[[v_name]]
              total_effect <- total_effect + eff
            }
          }
        }
      }
      # 【分支 B】: 包含非线性或部分变量 (使用 ols)
      if (!exists("fit") || !is_all_linear || length(formula_terms) < 3) {
        
        # 关键修正：data = curr_data (不要用 rawdata)
        fit <- ols(as.formula(f_str), data = curr_data, x = TRUE, y = TRUE)
        
        intercept_val <- coef(fit)["Intercept"]
        
        # 使用 predict(type="terms") 获取加权前的各项效应矩阵
        # 这会自动处理 rcs() 产生的多列复杂效应
        term_effects <- predict(fit, type = "terms")
        
        # 遍历每个活跃变量，施加约束因子
        for (v_name in active_vars) {
          v_x <- x_vars[which(var_names == v_name)] # e.g. "x_coagu"
          
          # 在 predict 的结果列名中查找该变量
          # 对于 rcs(x_coagu)，列名通常包含 "x_coagu"
          col_idx <- grep(v_x, colnames(term_effects), fixed = TRUE)
          
          if (length(col_idx) > 0) {
            # 如果是 rcs，可能会占用多列(取决于predict实现)，通常 type="terms" 对每个项返回一列
            # 我们取对应那一列的值
            raw_eff <- term_effects[, col_idx]
            weighted_eff <- raw_eff * constraints_list[[v_name]]
            total_effect <- total_effect + weighted_eff
          }
        }
      }
      # --- 5. 计算新的残差/校正后数据 ---
      # 公式逻辑：原始值 - (截距 + constraint * (各部分加权效应)) + 截距
      # 数学上等同于：原始值 - constraint * (各部分加权效应)
      # 但为了保留你原本的代码结构：
      ndata1[i, ] <- curr_data$y_log - (intercept_val + constraint * total_effect) + intercept_val
    }
  } else if (setequal(type, c("coagulation","platelet"))) {
    ## coagulation platelet ----
    for (i in 1:nrow(ndata1)) {
      y <- rawdata[i, ] # expression across samples
      y <- t(y)
      smpl2 <- smpl2[rownames(y),]
      x_plate <- smpl2$platelet
      x_coagu <- smpl2$coagulation
      Constraint_factor_plate <- if (!is.null(result_cor_platelet)) result_cor_platelet[colnames(y)] else 1
      Constraint_factor_coagu <- if (!is.null(result_cor_coagulation)) result_cor_coagulation[colnames(y)] else 1
      a <- summary(rlm(log2(y + 1) ~ x_plate + x_coagu, maxit = 30))
      for (j in 2:nrow(a$coefficients)) {
        b[i, j - 1] = a$coefficients[j, 1]
      }
      # new residuals:
      ndata1[i, ] <- log2(y + 1) - (a[["coefficients"]]["(Intercept)","Value"] + constraint * (
        x_plate * a[["coefficients"]]["x_plate","Value"]* Constraint_factor_plate +
          x_coagu *a[["coefficients"]]["x_coagu","Value"]* Constraint_factor_coagu) ) + 
        a[["coefficients"]]["(Intercept)","Value"]
    }
  } else if (setequal(type, c("coagulation","erythrocyte")))  {
    ## coagulation erythrocyte ----
    for (i in 1:nrow(ndata1)) {
      y <- rawdata[i, ] # expression across samples
      y <- t(y)
      smpl2 <- smpl2[rownames(y),]
      x_eryth <- smpl2$erythrocyte
      x_coagu <- smpl2$coagulation
      Constraint_factor_eryth <- if (!is.null(result_cor_erythrocyte)) result_cor_erythrocyte[colnames(y)] else 1
      Constraint_factor_coagu <- if (!is.null(result_cor_coagulation)) result_cor_coagulation[colnames(y)] else 1
      a <- summary(rlm(log2(y + 1) ~ x_eryth + x_coagu, maxit = 30))
      for (j in 2:nrow(a$coefficients)) {
        b[i, j - 1] = a$coefficients[j, 1]
      }
      # new residuals:
      ndata1[i, ] <- log2(y + 1) - (a[["coefficients"]]["(Intercept)","Value"] + constraint * (
        x_eryth *a[["coefficients"]]["x_eryth","Value"]* Constraint_factor_eryth +
          x_coagu *a[["coefficients"]]["x_coagu","Value"]* Constraint_factor_coagu) ) + 
        a[["coefficients"]]["(Intercept)","Value"]
    }
  } else if (setequal(type, c("erythrocyte","platelet"))){
    ## erythrocyte platelet ----
    for (i in 1:nrow(ndata1)) {
      y <- rawdata[i, ] # expression across samples
      y <- t(y)
      smpl2 <- smpl2[rownames(y),]
      x_eryth <- smpl2$erythrocyte
      x_plate <- smpl2$platelet
      Constraint_factor_eryth <- if (!is.null(result_cor_erythrocyte)) result_cor_erythrocyte[colnames(y)] else 1
      Constraint_factor_plate <- if (!is.null(result_cor_platelet)) result_cor_platelet[colnames(y)] else 1
      a <- summary(rlm(log2(y + 1) ~ x_eryth + x_plate, maxit = 30))
      for (j in 2:nrow(a$coefficients)) {
        b[i, j - 1] = a$coefficients[j, 1]
      }
      # new residuals:
      ndata1[i, ] <- log2(y + 1) - (a[["coefficients"]]["(Intercept)","Value"] + constraint * (
        x_eryth *a[["coefficients"]]["x_eryth","Value"]* Constraint_factor_eryth + 
          x_plate * a[["coefficients"]]["x_plate","Value"]* Constraint_factor_plate) ) + 
        a[["coefficients"]]["(Intercept)","Value"]
    }
  } else if (identical(type, "erythrocyte")) {
    # 红细胞单独校正 ----
    for (i in 1:nrow(ndata1)) {
      y <- rawdata[i, ] # expression across samples
      y <- t(y)
      smpl2 <- smpl2[rownames(y),]
      x_eryth <- smpl2$erythrocyte
      Constraint_factor_eryth <- if (!is.null(result_cor_erythrocyte)) result_cor_erythrocyte[colnames(y)] else 1
      a <- summary(rlm(log2(y + 1) ~ x_eryth, maxit = 30))
      for (j in 2:nrow(a$coefficients)) {
        b[i, j - 1] = a$coefficients[j, 1]
      }
      # new residuals:
      ndata1[i, ] <- log2(y + 1) - (a[["coefficients"]]["(Intercept)","Value"] + 
                                      constraint * x_eryth * a[["coefficients"]]["x_eryth","Value"] * Constraint_factor_eryth) + 
        a[["coefficients"]]["(Intercept)","Value"]
    }
  } else if (identical(type, "platelet")) {
    # 血小板单独校正 ----
    for (i in 1:nrow(ndata1)) {
      y <- rawdata[i, ] # expression across samples
      y <- t(y)
      smpl2 <- smpl2[rownames(y),]
      x_plate <- smpl2$platelet
      Constraint_factor_plate <- if (!is.null(result_cor_platelet)) result_cor_platelet[colnames(y)] else 1
      # x_coagu <- smpl2$coagulation
      a <- summary(rlm(log2(y + 1) ~ x_plate, maxit = 30))
      for (j in 2:nrow(a$coefficients)) {
        b[i, j - 1] = a$coefficients[j, 1]
      }
      # new residuals:
      ndata1[i, ] <- log2(y + 1) - (a[["coefficients"]]["(Intercept)","Value"] + 
                                      constraint * x_plate * a[["coefficients"]]["x_plate","Value"] * Constraint_factor_plate) + 
        a[["coefficients"]]["(Intercept)","Value"]
    }
  } else if (identical(type, "coagulation"))  {
    # 凝血单独校正 ----
    for (i in 1:nrow(ndata1)) {
      y <- rawdata[i, ] # expression across samples
      y <- t(y)
      smpl2 <- smpl2[rownames(y),]
      x_coagu <- smpl2$coagulation
      Constraint_factor_coagu <- if (!is.null(result_cor_coagulation)) result_cor_coagulation[colnames(y)] else 1
      a <- summary(rlm(log2(y + 1) ~ x_coagu, maxit = 30))
      for (j in 2:nrow(a$coefficients)) {
        b[i, j - 1] = a$coefficients[j, 1]
      }
      # new residuals:
      ndata1[i, ] <- log2(y + 1) - (a[["coefficients"]]["(Intercept)","Value"] + 
                                      constraint * x_coagu * a[["coefficients"]]["x_coagu","Value"] * Constraint_factor_coagu) + 
        a[["coefficients"]]["(Intercept)","Value"]
    }
  } else {
    stop("Unsupported correction type combination: ", paste(type, collapse = ", "))
  }
  result <- 2 ^ ndata1 - 1
  data$correct_data <- result
  data$contamination_level <- smpl2
  return(data)
}
