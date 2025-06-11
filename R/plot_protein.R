#' Plot Protein Expression Distribution
#'
#' Creates boxplot visualizations of protein expression patterns. Can display either:
#' 1) Distribution of selected proteins across all samples, or
#' 2) Grouped distribution of selected proteins across experimental conditions.
#'
#' @param data Expression matrix with proteins as rows and samples as columns
#' @param protein Character vector of protein names to plot (NULL plots all proteins)
#' @param group Data frame for grouping information with two columns:
#'   - `id`: sample IDs matching column names of `data`
#'   - `group`: experimental group labels
#'
#' @return A ggplot object showing:
#'   - Without grouping: Boxplots of protein expression distributions (x = proteins)
#'   - With grouping: Grouped boxplots of protein expression (x = proteins, fill = groups)
#'   - Y-axis: log2(expression + 1)
#'
#' @details
#' The function performs the following operations:
#' 1. Filters proteins if specified
#' 2. Transposes data to sample-oriented format
#' 3. Merges grouping information if provided
#' 4. Converts to long format for ggplot
#' 5. Applies log2(value + 1) transformation
#' 6. Generates appropriate boxplot visualization:
#'    - Ungrouped: Single box per protein
#'    - Grouped: Multiple boxes per protein (one per group)
#'
#' @note
#' - The +1 transformation prevents -Inf values from log2(0)
#' - Protein names on x-axis are rotated 90° for readability
#' - When grouping is provided, legend shows group-color mapping
#'
#' @examples
#' \dontrun{
#' # Create example data
#' expr_matrix <- matrix(rnorm(600), nrow = 100, 
#'                      dimnames = list(paste0("Protein", 1:100),
#'                                     paste0("Sample", 1:6)))
#' 
#' # Create grouping information
#' group_info <- data.frame(
#'   id = paste0("Sample", 1:6),
#'   group = rep(c("Control", "Treatment"), each = 3)
#' )
#'
#' # Plot specific proteins without grouping
#' p1 <- plot_protein_by_sample(expr_matrix, 
#'                             protein = c("Protein1", "Protein2", "Protein3"))
#' 
#' # Plot specific proteins with grouping
#' p2 <- plot_protein_by_sample(expr_matrix, 
#'                             protein = c("Protein1", "Protein2"),
#'                             group = group_info)
#' 
#' # Plot all proteins with grouping
#' p3 <- plot_protein_by_sample(expr_matrix, group = group_info)
#' }
#'
#' @importFrom ggplot2 ggplot aes geom_boxplot theme_classic theme element_text
#' @importFrom tidyr gather
#' @importFrom dplyr %>%
#' @export
library(ggplot2)
library(tidyr)
library(dplyr)

plot_protein_by_sample <- function(data, protein = NULL, group = NULL) {
  # data 为 matrix，行为 protein, 列为 sample id
  # protein 默认为 NULL，即对所有 protein 进行绘制
  # 如果想要绘制特定的 protein，则输入要绘制的 protein 向量
  # 如果要进行分组比较，则输入 group，默认不进行分组比较
  # group 要求为两列，id 为 sample id，group 为分组
  
  # 如果有指定的 protein，过滤数据
  if (!is.null(protein)) {
    data <- data[rownames(data)%in%protein,]
  }
  # 将矩阵转换为数据框
  data <- as.data.frame(t(data))
  data$sample_id <- rownames(data)
  # 如果有分组信息，合并分组信息
  if (!is.null(group)) {
    data <- merge(data, group, by.x = "sample_id", by.y = "id")
  }
  
  # 将数据转换为长格式
  data_long <- tidyr::gather(data, key = "protein", value = "value", -sample_id, -group)
  
  # 绘制箱线图
  if (is.null(group)) {
    p <- ggplot(data_long, aes(x = protein, y = log2(value + 1))) + 
      geom_boxplot() +
      theme_classic() +
      theme(axis.text.x = element_text(angle = 90, hjust = 1, colour = "black", size = 10),
            axis.text.y = element_text(hjust = 1, colour = "black", size = 10),
            plot.title = element_text(hjust = 0.5),
            legend.position = "none")
  } else {
    p <- ggplot(data_long, aes(x = protein, y = log2(value + 1), fill = group)) + 
      geom_boxplot() +
      theme_classic() +
      theme(axis.text.x = element_text(angle = 90, hjust = 1, colour = "black", size = 10),
            axis.text.y = element_text(hjust = 1, colour = "black", size = 10),
            plot.title = element_text(hjust = 0.5),
            legend.position = "right") 
  }
  
  return(p)
}
