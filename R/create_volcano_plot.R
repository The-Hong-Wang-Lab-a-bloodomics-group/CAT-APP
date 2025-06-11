#' Create Volcano Plot for Differential Expression Analysis
#'
#' Generates a customizable volcano plot to visualize differentially expressed proteins/genes 
#' with intelligent labeling and automatic font scaling.
#'
#' @param limma_result Data frame containing limma differential expression results
#' @param p_type Type of p-value to display: "adjusted" (adj.P.Val) or "raw" (P.Value)
#' @param p_cutoff Significance cutoff for p-values (default: 0.05)
#' @param logFC_cutoff Minimum absolute log fold change for significance (default: 1)
#' @param gene_col Column name containing gene/protein identifiers (default: "Protein")
#' @param group_names Character vector of comparison group names (c(numerator, denominator))
#' @param colors Color scheme for expression directions:
#'   - Up: upregulated (default: "#E64B35")
#'   - Down: downregulated (default: "#4DBBD5")
#'   - Not: non-significant (default: "grey80")
#' @param base_size Base font size (default: 14)
#' @param label_size Label font size (default: 4.5)
#' @param auto_font Enable automatic font scaling based on data size (default: TRUE)
#'
#' @return A ggplot2 volcano plot object showing:
#'   - X-axis: log2 fold change
#'   - Y-axis: -log10(p-value)
#'   - Points colored by expression direction
#'   - Labeled top differential proteins
#'   - Significance threshold lines
#'
#' @details
#' Key features:
#' 1. Automatically scales fonts based on data size when auto_font=TRUE:
#'    - Base size: 16 (<5k points), 14 (5-10k), 12 (>10k)
#'    - Label size: 4.5 (<5k), 4 (5-10k), 3.5 (>10k)
#' 2. Labels top 10 most significant proteins by absolute logFC
#' 3. Includes dynamic legend showing counts per category
#' 4. Uses ggrepel for non-overlapping labels
#' 5. Shows proper mathematical notation for axes
#'
#' @examples
#' \dontrun{
#' # Using limma results
#' volcano <- create_volcano_plot(
#'   limma_result = deg_results,
#'   group_names = c("Treatment", "Control"),
#'   p_cutoff = 0.01,
#'   logFC_cutoff = 1.5
#' )
#' 
#' # Customize appearance
#' volcano + 
#'   theme(legend.position = "bottom") +
#'   labs(title = "Differential Expression")
#' }
#'
#' @importFrom ggplot2 ggplot aes geom_point theme_classic theme scale_color_manual labs
#' @importFrom ggplot2 geom_hline geom_vline element_text element_blank
#' @importFrom ggrepel geom_text_repel
#' @importFrom dplyr filter arrange desc
#' @export
create_volcano_plot <- function(limma_result, 
                                p_type = "adjusted", 
                                p_cutoff = 0.05, 
                                logFC_cutoff = 1,
                                gene_col = "Protein",
                                group_names = c("Group1", "Group2"),
                                colors = c(Up = "#E64B35", Down = "#4DBBD5", Not = "grey80"),
                                base_size = 14,        # 新增基础字号参数
                                label_size = 4.5,      # 新增标签字号参数
                                auto_font = TRUE) {    # 新增自动调节开关
  # 检查必要包
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' required but not installed.")
  }
  if (!requireNamespace("ggrepel", quietly = TRUE)) {
    stop("Package 'ggrepel' required but not installed.")
  }
  library(ggplot2)
  library(ggrepel)
  library(dplyr)
  # 校验参数
  if (!p_type %in% c("adjusted", "raw")) {
    stop("p_type参数必须为 'adjusted' 或 'raw'")
  }
  if (!all(c("logFC", "P.Value", "adj.P.Val", gene_col) %in% colnames(limma_result))) {
    stop("输入数据缺少必要列")
  }
  
  # 选择使用的p值列
  p_column <- ifelse(p_type == "adjusted", "adj.P.Val", "P.Value")
  
  # 准备绘图数据
  plot_data <- limma_result
  plot_data$gene <- plot_data[[gene_col]]
  
  # 创建显著性分类
  plot_data$sig <- "Not"
  plot_data$sig[plot_data[[p_column]] < p_cutoff & plot_data$logFC >= logFC_cutoff] <- "Up"
  plot_data$sig[plot_data[[p_column]] < p_cutoff & plot_data$logFC <= -logFC_cutoff] <- "Down"
  plot_data$sig <- factor(plot_data$sig, levels = c("Up", "Down", "Not"))
  
  # 智能字号调节逻辑
  if(auto_font) {
    n_points <- nrow(limma_result)
    base_size <- dplyr::case_when(
      n_points > 10000 ~ 12,
      n_points > 5000 ~ 14,
      TRUE ~ 16
    )
    label_size <- dplyr::case_when(
      n_points > 10000 ~ 3.5,
      n_points > 5000 ~ 4,
      TRUE ~ 4.5
    )
  }
  # 创建基础火山图
  volc <- ggplot2::ggplot(data = plot_data, 
                          aes(x = logFC, 
                              y = -log10(.data[[p_column]]), 
                              color = sig)) +
    ggplot2::geom_point(alpha = 0.9) +
    ggplot2::theme_classic(base_size = base_size) +  # 应用基础字号
    ggplot2::theme(
      panel.grid = element_blank(),
      plot.title = element_text(hjust = 0.5, face = "bold", size = base_size * 1.2),
      axis.title = element_text(size = base_size * 1.1),
      axis.text = element_text(size = base_size * 0.9),
      legend.text = element_text(size = base_size * 0.9),
      legend.title = element_text(size = base_size)
    ) +
    ggplot2::scale_color_manual(
      values = colors,
      name = "Expression",
      labels = c(
        Up = paste0("Up (n=", sum(plot_data$sig == "Up"), ")"),
        Down = paste0("Down (n=", sum(plot_data$sig == "Down"), ")"),
        Not = "Not significant"
      )
    ) +
    ggplot2::labs(
      x = expression("log"[2]*" Fold Change"),
      y = ifelse(p_type == "adjusted",
                 expression("-log"[10]*" (Adj.P)"),
                 expression("-log"[10]*" (P)")),
      title = paste0(group_names[1], " vs ", group_names[2])
    )
  
  # 添加显著性阈值线
  volc <- volc +
    ggplot2::geom_hline(
      yintercept = -log10(p_cutoff),
      linetype = "dashed",
      color = "black",
      linewidth = 0.6,
      alpha = 0.8
    ) +
    ggplot2::geom_vline(
      xintercept = c(-logFC_cutoff, logFC_cutoff),
      linetype = "dashed",
      color = "black",
      linewidth = 0.6,
      alpha = 0.8
    )
  
  # 添加标签（显示前10个显著基因）
  top_genes <- plot_data %>%
    dplyr::filter(sig %in% c("Up", "Down")) %>%
    dplyr::arrange(desc(abs(logFC))) %>%
    head(10)
  
  # 添加标签（动态调节标签字号）
  if(nrow(top_genes) > 0) {
    volc <- volc +
      ggrepel::geom_text_repel(
        data = top_genes,
        aes(label = gene),
        size = label_size,          # 应用标签字号
        box.padding = 0.5,
        max.overlaps = Inf,
        segment.color = "grey50"
      )
  }
  return(volc)
}
