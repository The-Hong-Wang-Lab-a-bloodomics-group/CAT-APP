#' Plot Distribution of Statistical Metrics
#'
#' Generates a histogram for the distribution of correlation coefficients or p-values, 
#' with options to add kernel density curves and significance thresholds.
#'
#' @param data Data matrix or data frame with features as rows and samples as columns
#' @param marker_list Optional list of marker groups for aggregated correlation analysis.
#'    If provided, calculates average expression for each group and correlates with all features.
#'    If NULL, performs pairwise correlation between all features. Default: NULL
#' @param type Correlation type (currently only "pearson" supported; reserved for future extensions)
#' @param statistic Type of statistic to plot:
#'   - "pvalue": P-value distribution histogram (default)
#'   - "correlation": Correlation coefficient distribution histogram
#' @param alpha Significance threshold (only effective when statistic = "pvalue"). Default: 0.05
#' @param title Plot title (auto-generated if NULL)
#' @param xlab X-axis label (auto-generated if NULL)
#' @param ylab Y-axis label. Default: "频数"
#' @param binwidth Histogram bin width (default: 0.05 for p-values, 0.1 for correlations)
#' @param boundary Starting boundary for histogram bins (default: 0 for p-values, -1 for correlations)
#' @param add_fit Whether to add kernel density curve. Default: TRUE
#' @param fit_color Color for density curve. Default: "#2E8B57" (forest green)
#' @param fit_size Line width for density curve. Default: 1.2
#'
#' @return A ggplot2 object
#' 
#' @details
#' When marker_list is provided:
#'   - Computes average expression for each marker group
#'   - Correlates each group average with all features in data
#'   - Visualizes the combined correlation metrics
#' 
#' @importFrom ggplot2 ggplot aes geom_histogram scale_fill_manual geom_vline
#' @importFrom ggplot2 geom_line labs theme_minimal theme element_text
#' @importFrom ggplot2 after_stat
#' @importFrom Rfast cora
#' @importFrom stats pt
#' 
#' @export
plot_stat_distribution <- function(data, marker_list = NULL, type = "pearson", statistic = "pvalue", alpha = 0.05,
                                   title = NULL, xlab = NULL, ylab = "Number",
                                   binwidth = NULL, boundary = NULL,
                                   add_fit = TRUE, fit_color = "#2E8B57", fit_size = 1.2) {
  # 检查必要包安装
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("please install the package 'ggplot2': install.packages('ggplot2')")
  }
  # 检查statistic参数有效性
  if (!statistic %in% c("pvalue", "correlation")) {
    stop("statistic must be 'pvalue' or 'correlation'")
  }
  
  # 设置默认参数
  if (is.null(title)) {
    title <- if (statistic == "pvalue") "P-value distribution histogram" else "correlation coefficient distribution histogram"
  }
  if (is.null(xlab)) {
    xlab <- if (statistic == "pvalue") "P value" else "correlation coefficient"
  }
  if (is.null(binwidth)) {
    binwidth <- if (statistic == "pvalue") 0.05 else 0.1
  }
  if (is.null(boundary)) {
    boundary <- if (statistic == "pvalue") 0 else -1
  }
  
  # 计算相关矩阵
  if (!is.null(marker_list)) {
    # 验证marker_list结构
    if (!is.list(marker_list) || length(marker_list) == 0) {
      stop("marker_list must be a non-empty list that contains the tag group")
    }
    
    # 计算每个组的平均表达 (转置数据使样本为行)
    t_data <- as.data.frame(t(data))
    group_avgs <- sapply(marker_list, function(markers) {
      # 筛选数据中存在的标记
      valid_markers <- intersect(markers, colnames(t_data))
      if (length(valid_markers) == 0) {
        stop("No markers were found in the data")
      }
      rowMeans(t_data[, valid_markers, drop = FALSE], na.rm = TRUE)
    })
    
    # 合并组平均和原始转置数据
    combined_data <- cbind(group_avgs, t_data)
    
    # 计算相关矩阵
    cor_matrix_full <- Rfast::cora(as.matrix(combined_data))
    n_samples <- nrow(combined_data)
    
    # 提取组平均与所有特征的相关性
    n_groups <- ncol(group_avgs)
    total_features <- ncol(combined_data)
    cor_sub <- cor_matrix_full[1:n_groups, (n_groups + 1):total_features, drop = FALSE]
    
    # 计算t统计量和p值
    t_stats <- cor_sub * sqrt((n_samples - 2) / (1 - cor_sub^2))
    p_values_sub <- 2 * pt(abs(t_stats), df = n_samples - 2, lower.tail = FALSE)
    
    # 选择统计量
    if (statistic == "pvalue") {
      values <- as.vector(p_values_sub)
    } else {
      values <- as.vector(cor_sub)
    }
  } else {
    # 原始方法：所有特征间的成对相关 (使用转置数据)
    t_data <- as.data.frame(t(data))
    cor_matrix <- Rfast::cora(as.matrix(t_data))
    n_samples <- nrow(t_data)
    
    # 计算t统计量和p值
    t_stats <- cor_matrix * sqrt((n_samples - 2) / (1 - cor_matrix^2))
    p_values <- 2 * pt(abs(t_stats), df = n_samples - 2, lower.tail = FALSE)
    
    # 移除对角线
    diag(p_values) <- NA
    diag(cor_matrix) <- NA
    
    # 选择统计量
    if (statistic == "pvalue") {
      values <- p_values[upper.tri(p_values)]
    } else {
      values <- cor_matrix[upper.tri(cor_matrix)]
    }
  }
  
  # 创建绘图数据
  df <- data.frame(value = values)
  
  # 创建基础绘图对象
  p <- ggplot2::ggplot(df, ggplot2::aes(x = value)) 
  
  # 根据统计量类型添加图层
  if (statistic == "pvalue") {
    p <- p +
      ggplot2::geom_histogram(
        ggplot2::aes(fill = value < alpha, y = ggplot2::after_stat(count)),
        binwidth = binwidth,
        color = "black",
        boundary = boundary
      ) +
      ggplot2::scale_fill_manual(
        values = c("FALSE" = "#87CEEB", "TRUE" = "#FF6B6B"), 
        guide = "none"
      ) +
      ggplot2::geom_vline(
        xintercept = alpha, 
        color = "#FF5252", 
        linetype = "dashed", 
        linewidth = 0.8
      )
    
    # 添加P值的核密度估计
    if (add_fit) {
      p <- p + 
        ggplot2::geom_line(
          ggplot2::aes(y = ggplot2::after_stat(count * binwidth * density)),
          stat = "density",
          color = fit_color,
          linewidth = fit_size
        )
    }
  } else {
    p <- p +
      ggplot2::geom_histogram(
        ggplot2::aes(y = ggplot2::after_stat(count)),
        fill = "#87CEEB",
        binwidth = binwidth,
        color = "black",
        boundary = boundary
      )
    
    # 添加相关系数的密度曲线
    if (add_fit) {
      p <- p + 
        ggplot2::geom_line(
          ggplot2::aes(y = ggplot2::after_stat(count * binwidth * density)),
          stat = "density",
          color = fit_color,
          linewidth = fit_size
        )
    }
  }
  
  # 添加通用元素
  p <- p +
    ggplot2::labs(
      title = title,
      x = xlab,
      y = ylab
    ) +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, colour = "black", size = 13),
          axis.text.y = element_text(hjust = 1, colour = "black", size = 13),
          axis.title.y = element_text(colour = "black", size = 15),
          axis.title.x = element_text(colour = "black", size = 15),
          plot.title = element_text(hjust = 0.5, face = "bold",size = 15),
          panel.grid.major.x = ggplot2::element_blank())
  return(p)
}
