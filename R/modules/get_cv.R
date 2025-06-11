#' Calculate Coefficient of Variation (CV)
#'
#' Computes the coefficient of variation (CV) for protein expression data, either overall or by sample groups.
#'
#' @param raw_data Raw data matrix with proteins as rows and samples as columns
#' @param protein Vector of protein names to include in calculation (default: all proteins)
#' @param data_group Data frame specifying sample groups with two required columns:
#'   - id: sample IDs matching column names in raw_data
#'   - group: group labels for samples
#' @param by_group Logical indicating whether to calculate CV by group:
#'   - TRUE: calculates separate CVs for each group
#'   - FALSE: calculates overall CV (default)
#'
#' @return Data frame containing CV results:
#'   - When by_group = FALSE: contains three columns (mean, sd, CV) with proteins as rows
#'   - When by_group = TRUE: contains one column per group (named "groupname_CV") with proteins as rows
#'
#' @details
#' The coefficient of variation is calculated as: \deqn{CV = \frac{\sigma}{\mu}}
#' where σ is standard deviation and μ is mean.
#' 
#' Processing workflow:
#' 1. Filters proteins based on 'protein' parameter
#' 2. For group-wise calculation (by_group = TRUE):
#'    a. Iterates through each unique group
#'    b. Subsets samples belonging to current group
#'    c. Computes mean and standard deviation for each protein
#'    d. Calculates group-specific CVs
#' 3. For overall calculation (by_group = FALSE):
#'    a. Computes mean and standard deviation across all samples
#'    b. Calculates overall CV
#'
#' @note
#' - Grouped calculation requires 'data_group' parameter
#' - Column names for grouped results follow "groupname_CV" format
#' - Returns NA for proteins with zero mean (division by zero)
#'
#' @examples
#' \dontrun{
#' # Sample data matrix
#' data_matrix <- matrix(rnorm(100), nrow = 10, 
#'                      dimnames = list(paste0("Protein", 1:10), 
#'                                     paste0("Sample", 1:10)))
#' 
#' # Group information
#' group_info <- data.frame(
#'   id = paste0("Sample", 1:10),
#'   group = rep(c("Control", "Treatment"), each = 5)
#' )
#'
#' # Overall CV calculation
#' cv_all <- get_cv(data_matrix)
#' 
#' # Grouped CV calculation
#' cv_by_group <- get_cv(data_matrix, data_group = group_info, by_group = TRUE)
#' 
#' # CV calculation for specific proteins
#' cv_specific <- get_cv(data_matrix, protein = c("Protein1", "Protein2"))
#' }
#'
#' @export
# get_Cv ----
# raw_data要求行为蛋白，列为样本
get_cv <- function(raw_data,protein = rownames(raw_data),
                   data_group = NULL,
                   by_group = F){
  raw_data <- raw_data[rownames(raw_data)%in%protein,,drop = F]
  result_cv <- matrix(nrow = dim(raw_data)[1])
  result_cv <- as.data.frame(result_cv)
  if(by_group){
    for (group in unique(data_group$group)){
      data_group_one <- data_group[data_group$group%in% group,]
      data <- raw_data[,colnames(raw_data)%in%data_group_one$id]
      mean_sd <- apply(t(data), 2, function(x) c(mean = mean(x), sd = sd(x)))
      mean_sd <- as.data.frame(mean_sd)
      mean_sd <- as.data.frame(t(mean_sd))
      # 计算每个样本的变异系数（CV）
      result_cv[,paste0(group,"_CV")] <- mean_sd$sd / mean_sd$mean
    }
    result_cv <- result_cv[,-1]
  }else if(by_group == F){
    mean_sd <- apply(t(raw_data), 2, function(x) c(mean = mean(x), sd = sd(x)))
    mean_sd <- as.data.frame(mean_sd)
    mean_sd <- as.data.frame(t(mean_sd))
    # 计算每个样本的变异系数（CV）
    result_cv <- mean_sd
    result_cv[,"CV"] <- result_cv$sd / result_cv$mean
  }
  return(result_cv)
}
