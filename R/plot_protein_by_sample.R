#' Plot Protein Expression Distribution by Sample
#'
#' Creates a boxplot visualization of log2-transformed protein expression values for each sample.
#' Useful for quality control to identify sample-level expression distribution patterns.
#'
#' @param data Expression matrix with proteins as rows and samples as columns
#'
#' @return A ggplot object showing:
#'   - X-axis: Sample IDs
#'   - Y-axis: log2(expression + 1)
#'   - Boxes: Distribution of protein expression per sample
#'   - Color-coded by sample ID
#'
#' @details
#' The function performs the following transformations:
#' 1. Transposes input matrix (samples → rows, proteins → columns)
#' 2. Converts to long format suitable for ggplot
#' 3. Applies log2(value + 1) transformation to handle zeros and normalize variance
#' 4. Generates boxplots with sample-colored fills
#'
#' @note
#' - The +1 transformation prevents -Inf values from log2(0)
#' - Boxplots show median, quartiles, and outliers
#' - Sample IDs are rotated 90° for readability
#'
#' @examples
#' \dontrun{
#' # Create example data (100 proteins x 6 samples)
#' expr_matrix <- matrix(rnorm(600), nrow = 100, 
#'                      dimnames = list(paste0("Protein", 1:100),
#'                                     paste0("Sample", 1:6)))
#' 
#' # Generate plot
#' p <- plot_protein_by_sample(expr_matrix)
#' print(p)
#' 
#' # Customize plot
#' p + labs(title = "Protein Expression Distribution") + 
#'     theme(legend.position = "none")
#' }
#'
#' @importFrom ggplot2 ggplot aes geom_boxplot theme_classic theme element_text labs
#' @importFrom tidyr gather
#' @importFrom dplyr %>%
#' @export
library(ggplot2)
library(tidyr)
library(dplyr)

plot_protein_by_sample <- function(data) {
  # data 为 matrix，行为 protein, 列为 sample id
  
  # 将矩阵转换为数据框
  data <- as.data.frame(t(data))
  data$sample_id <- rownames(data)
  
  # 将数据转换为长格式
  data_long <- tidyr::gather(data, key = "protein", value = "value", -sample_id)
  # 绘制箱线图
  p <- ggplot(data_long, aes(x = sample_id, y = log2(value + 1), fill = sample_id)) + 
    geom_boxplot() +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, colour = "black", size = 10),
          axis.text.y = element_text(hjust = 1, colour = "black", size = 10),
          plot.title = element_text(hjust = 0.5),
          legend.position = "right") 
  return(p)
}


