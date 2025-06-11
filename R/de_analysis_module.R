#' Differential Expression Analysis using Limma for Proteomics Data
#'
#' Performs differential expression analysis on proteomics data using the limma package.
#' Handles both standard proteomics data and count-based proteomics (with voom transformation).
#'
#' @param expr_matrix Expression matrix with proteins as rows and samples as columns
#' @param group_matrix Data frame containing sample grouping information with two required columns:
#'   - `id`: sample IDs matching colnames of `expr_matrix`
#'   - `group`: biological group labels
#' @param compare Character vector specifying the comparison groups in the format `c("numerator", "denominator")`.
#'   Example: `c("Treatment", "Control")` calculates Treatment vs Control. If NULL, uses the first two group levels in reverse order.
#' @param p_cutoff Significance threshold for p-values (default: 0.05)
#' @param logFC_cutoff Minimum absolute log fold change threshold (default: 0)
#' @param adjust_method P-value adjustment method ("BH", "bonferroni", etc., default: "BH")
#' @param apply_voom Logical indicating whether to apply voom transformation (for count-like proteomics data, default: FALSE)
#' @param p_type Type of p-value to use for significance:
#'   - "raw": unadjusted p-value (default)
#'   - "adjusted": adjusted p-value
#'
#' @return A data frame containing:
#'   - Protein identifiers
#'   - logFC: log2 fold change
#'   - P.Value: raw p-value
#'   - adj.P.Val: adjusted p-value
#'   - significant: logical vector indicating significant proteins
#'   - Other limma statistics (t-statistic, B-statistic, etc.)
#'
#' @details
#' Analysis workflow:
#' 1. Input validation: Checks sample-group matching and parameter integrity
#' 2. Group factor creation: Matches sample groups to expression matrix columns
#' 3. Comparison setup: Uses specified comparison or default level ordering
#' 4. Transformation: Applies voom if requested (for count-like data)
#' 5. Model design: Creates design matrix with 0-intercept groups
#' 6. Contrast specification: Builds specified comparison
#' 7. Model fitting: Linear model -> contrasts -> empirical Bayes moderation
#' 8. Result extraction: Top table with all proteins
#' 9. Significance marking: Based on specified p-value type and cutoffs
#'
#' @note
#' - The comparison direction: `compare = c("A", "B")` means A vs B (logFC > 0 indicates upregulation in A)
#' - For proteomics data, `apply_voom = TRUE` is only needed for spectral count-type data
#' - The design matrix uses 0-intercept (~0 + group) for direct group mean parameterization
#'
#' @examples
#' \dontrun{
#' # Example expression matrix (10 proteins x 6 samples)
#' expr_data <- matrix(rnorm(60), nrow = 10, 
#'                    dimnames = list(paste0("Protein", 1:10),
#'                                   paste0("Sample", 1:6)))
#' 
#' # Group information
#' group_info <- data.frame(
#'   id = paste0("Sample", 1:6),
#'   group = rep(c("Control", "Treatment"), each = 3)
#' )
#'
#' # Run analysis
#' deg_results <- limma_proteomics_analysis(
#'   expr_matrix = expr_data,
#'   group_matrix = group_info,
#'   compare = c("Treatment", "Control"),
#'   p_cutoff = 0.05,
#'   logFC_cutoff = 1,
#'   adjust_method = "BH",
#'   apply_voom = FALSE
#' )
#'
#' # View significant results
#' subset(deg_results, significant)
#' }
#'
#' @importFrom limma voom lmFit makeContrasts contrasts.fit eBayes topTable
#' @importFrom stats model.matrix
#' @export
limma_proteomics_analysis <- function(expr_matrix, 
                                      group_matrix, 
                                      compare = NULL, 
                                      p_cutoff = 0.05, 
                                      logFC_cutoff = 0, 
                                      adjust_method = "BH", 
                                      apply_voom = FALSE,
                                      p_type = "raw") {
  # [Function implementation remains unchanged]
}
limma_proteomics_analysis <- function(expr_matrix, 
                                      group_matrix, 
                                      compare = NULL, 
                                      p_cutoff = 0.05, 
                                      logFC_cutoff = 0, 
                                      adjust_method = "BH", 
                                      apply_voom = FALSE,
                                      p_type = "raw") {
  # 检查依赖包
  if (!requireNamespace("limma", quietly = TRUE)) {
    stop("Package 'limma' required but not installed.")
  }
  
  # 验证输入数据结构
  if (!all(c("id", "group") %in% colnames(group_matrix))) {
    stop("分组矩阵必须包含'id'和'group'两列")
  }
  # 校验p_type参数有效性
  if (!p_type %in% c("adjusted", "raw")) {
    stop("p_type参数必须为'adjusted'或'raw'")
  }
  # 确保样本顺序匹配
  sample_ids <- colnames(expr_matrix)
  if (!all(sample_ids %in% group_matrix$id)) {
    missing_samples <- setdiff(sample_ids, group_matrix$id)
    stop("以下样本缺少分组信息: ", paste(missing_samples, collapse = ", "))
  }
  
  # 提取与表达矩阵匹配的分组信息
  matched_groups <- group_matrix$group[match(sample_ids, group_matrix$id)]
  group <- factor(matched_groups)
  
  # 设置默认比较组
  if (is.null(compare)) {
    compare <- levels(group)[2:1] # 默认第二个组 vs 第一个组  
  } else if (!all(compare %in% levels(group))) {
    stop("对比组必须存在于分组信息中")
  }
  
  # 数据转换（根据需求选择）
  if (apply_voom) {
    expr_matrix <- limma::voom(expr_matrix)$E
  }
  
  # 构建线性模型
  design <- model.matrix(~0 + group)
  colnames(design) <- levels(group)
  
  # 创建对比
  contrast_formula <- paste(compare[1], "-", compare[2])
  contrast_matrix <- limma::makeContrasts(
    contrasts = contrast_formula, 
    levels = design
  )
  
  # 模型拟合
  fit <- limma::lmFit(expr_matrix, design)
  fit <- limma::contrasts.fit(fit, contrast_matrix)
  fit <- limma::eBayes(fit)
  
  # 提取结果
  top_table <- limma::topTable(
    fit, 
    number = Inf, 
    adjust.method = adjust_method
  )
  # 动态选择p值列 --------------------------------------------------
  p_column <- ifelse(p_type == "adjusted", "adj.P.Val", "P.Value")
  
  # 添加显著性标记
  top_table$significant <- (abs(top_table$logFC) >= logFC_cutoff) & 
    (top_table[, p_column] <= p_cutoff)
  
  # 添加蛋白质ID列
  if (!"Protein" %in% colnames(top_table)) {
    top_table <- cbind(Protein = rownames(top_table), top_table)
    rownames(top_table) <- NULL
  }
  
  return(top_table)
}
