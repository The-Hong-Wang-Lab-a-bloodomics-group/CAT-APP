#' Plot Distribution of Statistical Metrics
#'
#' Generates a histogram for the distribution of correlation coefficients or p-values, 
#' with options to add kernel density curves and significance thresholds.
#'
#' @param data Data matrix or data frame with samples as rows and features as columns
#' @param type Correlation type (currently only "pearson" supported; reserved for future extensions)
#' @param statistic Type of statistic to plot:
#'   - "pvalue": P-value distribution histogram (default)
#'   - "correlation": Correlation coefficient distribution histogram
#' @param alpha Significance threshold (only effective when statistic = "pvalue"). Default: 0.05
#' @param title Plot title (auto-generated if NULL)
#' @param xlab X-axis label (auto-generated if NULL)
#' @param ylab Y-axis label. Default: "Frequency"
#' @param binwidth Histogram bin width (default: 0.05 for p-values, 0.1 for correlations)
#' @param boundary Starting boundary for histogram bins (default: 0 for p-values, -1 for correlations)
#' @param add_fit Whether to add kernel density curve. Default: TRUE
#' @param fit_color Color for density curve. Default: "#2E8B57" (forest green)
#' @param fit_size Line width for density curve. Default: 1.2
#'
#' @return A ggplot2 object
#' 
#' @details
#' Function workflow:
#' 1. Checks for required package installations (ggplot2)
#' 2. Computes correlation matrix and associated p-values
#' 3. Selects visualization based on statistic parameter
#' 4. Constructs histogram with optional elements (density curve/significance threshold)
#' 
#' When statistic = "pvalue":
#'   - Red dashed line indicates significance threshold (alpha)
#'   - Red bars represent significant p-values (p < alpha)
#'   - Blue bars represent non-significant p-values
#' 
#' When statistic = "correlation":
#'   - All bars shown in blue
#'   - Density curve shows distribution of correlation coefficients
#'
#' @importFrom ggplot2 ggplot aes geom_histogram scale_fill_manual geom_vline
#' @importFrom ggplot2 geom_line labs theme_minimal theme element_text
#' @importFrom ggplot2 after_stat
#' @importFrom Rfast cora
#' @importFrom stats pt
#' 
#' @examples
#' \dontrun{
#' # Example data
#' set.seed(123)
#' data_matrix <- matrix(rnorm(1000), nrow = 100, ncol = 10)
#' 
#' # Plot p-value distribution
#' plot_stat_distribution(data_matrix)
#' 
#' # Plot correlation distribution with custom parameters
#' plot_stat_distribution(
#'   data = data_matrix,
#'   statistic = "correlation",
#'   title = "Correlation Coefficient Distribution",
#'   binwidth = 0.05,
#'   fit_color = "darkblue"
#' )
#' }
#' 
#' @export
plot_stat_distribution <- function(data, type = "pearson", statistic = "pvalue", alpha = 0.05,
                                   title = NULL, xlab = NULL, ylab = "频数",
                                   binwidth = NULL, boundary = NULL,
                                   add_fit = TRUE, fit_color = "#2E8B57", fit_size = 1.2) {
  # 检查必要包安装
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("请先安装ggplot2包: install.packages('ggplot2')")
  }
  # 检查statistic参数有效性
  if (!statistic %in% c("pvalue", "correlation")) {
    stop("参数statistic必须是'pvalue'或'correlation'")
  }
  
  # 设置默认参数
  if (is.null(title)) {
    title <- if (statistic == "pvalue") "P值分布直方图" else "相关系数分布直方图"
  }
  if (is.null(xlab)) {
    xlab <- if (statistic == "pvalue") "P值" else "相关系数"
  }
  if (is.null(binwidth)) {
    binwidth <- if (statistic == "pvalue") 0.05 else 0.1
  }
  if (is.null(boundary)) {
    boundary <- if (statistic == "pvalue") 0 else -1
  }
  
  # 计算相关矩阵和p值矩阵
  cor_matrix  <- Rfast::cora(as.matrix(t(data)))
  n <- nrow(t(data))
  t_stats <- cor_matrix  * sqrt((n - 2) / (1 - cor_matrix ^2))
  p_values <- 2 * pt(abs(t_stats), df = n - 2, lower.tail = FALSE)
  rownames(p_values) <- rownames(cor_matrix)
  colnames(p_values) <- colnames(cor_matrix)
  diag(p_values) <- NA
  test_result <- list(r = cor_matrix,
                      P = p_values)
  if (statistic == "pvalue") {
    # 计算相关矩阵和p值矩阵
    cor_matrix  <- Rfast::cora(as.matrix(t(data)))
    n <- nrow(t(data))
    t_stats <- cor_matrix  * sqrt((n - 2) / (1 - cor_matrix ^2))
    p_values <- 2 * pt(abs(t_stats), df = n - 2, lower.tail = FALSE)
    rownames(p_values) <- rownames(cor_matrix)
    colnames(p_values) <- colnames(cor_matrix)
    diag(p_values) <- NA
    values <- p_values
  } else {
    # 计算相关矩阵和p值矩阵
    cor_matrix  <- Rfast::cora(as.matrix(t(data)))
    values <- cor_matrix
  }
  
  # 创建绘图数据
  df <- data.frame(value = values[upper.tri(values)])
  
  # 创建基础绘图对象
  p <- ggplot2::ggplot(df, ggplot2::aes(x = value))
  
  # 根据统计量类型添加图层
  if (statistic == "pvalue") {
    p <- p +
      ggplot2::geom_histogram(
        ggplot2::aes(fill = value < alpha, y = after_stat(count)),
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
          ggplot2::aes(y = after_stat(count * binwidth * density)),
          stat = "density",
          color = fit_color,
          linewidth = fit_size
        )
    }
  } else {
    p <- p +
      ggplot2::geom_histogram(
        ggplot2::aes(y = after_stat(count)),
        fill = "#87CEEB",
        binwidth = binwidth,
        color = "black",
        boundary = boundary
      )
    
    # 添加相关系数的密度曲线
    if (add_fit) {
      p <- p + 
        ggplot2::geom_line(
          ggplot2::aes(y = after_stat(count * binwidth * density)),
          stat = "density",
          color = fit_color,
          linewidth = fit_size
        )
    }
  }
  
  # 添加通用元素
  # caption_text <- if (statistic == "pvalue") {
  #   sprintf("显著性阈值: α = %.2f (红色虚线)\n绿色曲线：核密度估计", alpha)
  # } else {
  #   "绿色曲线：核密度估计"
  # }
  
  p <- p +
    ggplot2::labs(
      title = title,
      x = xlab,
      y = ylab
      # caption = caption_text
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      plot.title = ggplot2::element_text(hjust = 0.5, face = "bold"),
      panel.grid.major.x = ggplot2::element_blank(),
      plot.caption = ggplot2::element_text(hjust = 0)
    )
  
  return(p)
}
