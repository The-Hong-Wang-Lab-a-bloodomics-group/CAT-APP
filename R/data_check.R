#' Quality Control for Contamination Markers in Proteomics Data
#'
#' Performs quality control analysis for erythrocyte, platelet, and coagulation contamination markers 
#' in proteomics data. Identifies non-differentially expressed markers with high correlation to 
#' detect potential sample contamination issues.
#'
#' @param data Expression data matrix with proteins as rows and samples as columns
#' @param data_group Sample group information data frame (must contain 'id' and 'group' columns)
#' @param cutoff Correlation threshold for identifying highly correlated markers (default: 0.9)
#' @param DE_filter Whether to filter differentially expressed markers between groups (default: TRUE)
#' @param group1 First comparison group for DE analysis (optional)
#' @param group2 Second comparison group for DE analysis (optional)
#' @param custom_erythrocyte Custom list of erythrocyte markers (default: predefined list)
#' @param custom_coagulation Custom list of coagulation markers (default: predefined list)
#' @param custom_platelet Custom list of platelet markers (default: predefined list)
#'
#' @return A list containing:
#' \describe{
#'   \item{rawdata}{Original expression matrix}
#'   \item{group}{Sample group information}
#'   \item{data}{Contamination marker expression matrices (erythrocyte, coagulation, platelet)}
#'   \item{plot_marker}{Level plots for filtered markers}
#'   \item{plot_contamination}{Box plots showing marker expression across groups}
#'   \item{missing_genes}{Markers missing from the dataset}
#'   \item{correlation}{Correlation matrices for each marker type}
#'   \item{marker_list}{Filtered markers passing QC criteria}
#'   \item{marker_stats}{Statistical summary for each marker}
#' }
#'
#' @details
#' The function performs the following quality control steps:
#' 1. Checks input validity and loads required packages
#' 2. Defines contamination marker panels (default or custom)
#' 3. Transforms data to long format and merges with group information
#' 4. For each marker panel (erythrocyte, platelet, coagulation):
#'    a. Generates box plots with statistical tests
#'    b. Identifies non-DE markers (p > 0.05)
#'    c. Computes correlation matrix between markers
#'    d. Filters markers that are non-DE and highly correlated (r > cutoff)
#' 5. Generates level plots for filtered markers
#' 6. Compiles comprehensive statistics for each marker
#'
#' @note
#' - When DE_filter = TRUE, data_group must be provided
#' - Default marker lists include known contamination markers from literature
#' - Returns NA for analyses with insufficient markers
#'
#' @examples
#' \dontrun{
#' # Example usage
#' data_matrix <- matrix(rnorm(1000), nrow = 100, 
#'                      dimnames = list(paste0("Protein", 1:100),
#'                                     paste0("Sample", 1:10)))
#' group_info <- data.frame(id = paste0("Sample", 1:10),
#'                          group = rep(c("A", "B"), each = 5))
#' 
#' result <- data_check(
#'   data = data_matrix,
#'   data_group = group_info,
#'   cutoff = 0.85,
#'   custom_erythrocyte = c("HBA1", "HBB", "CA1")
#' )
#' 
#' # View erythrocyte contamination plot
#' print(result$plot_contamination$erythrocyte)
#' 
#' # View marker statistics
#' head(result$marker_stats$stats_erythrocyte)
#' }
#'
#' @importFrom ggplot2 ggplot aes geom_boxplot theme_classic labs theme element_text
#' @importFrom ggpubr stat_compare_means
#' @importFrom tidyr gather spread
#' @importFrom dplyr filter distinct bind_cols select mutate pull arrange group_by
#' @importFrom stringr str_detect
#' @importFrom stats cor
#' @export
data_check <- function(data, data_group = NULL, cutoff = 0.9, DE_filter = T,
                       group1 = NULL, group2 = NULL,
                       custom_erythrocyte = NULL,
                       custom_coagulation = NULL,
                       custom_platelet = NULL){
  # 函数说明 ----
  # data为输入数据矩阵，要求为样本为列，蛋白为行，矩阵中无非数值数据
  # data_group为样本分组，要求为：id列为样本名，和colnames(data)相同，group为输入样本的生物学分组
  # cutoff为相关性阈值，默认0.9
  
  # 加载必要的库
  library(ggplot2)
  library(ggpubr)
  library(tidyr)
  library(dplyr)
  library(stringr)
  # 检查输入数据的有效性
  if (DE_filter == T && is.null(data_group)) {
    stop("当DE_filter为TRUE时，必须提供data_group参数")
  }
  
  if (!is.null(data_group) && !all(data_group$id %in% colnames(data))) {
    stop("Error: data_group中的id列与data的列名不匹配。")
  }
  # 定义分析标记物的函数
  # Function to filter markers and perform analysis
  # 删除生物学分组中有差异的marker
  analyze_markers <- function(data = data_ggplot, marker_list, title) {
    target_proteins <- marker_list$GN
    pattern <- paste0("\\b(", paste(target_proteins, collapse = "|"), ")\\b")
    data_filtered <- data_ggplot %>% filter(str_detect(key, pattern))
    
    # 根据是否有分组采用不同的绘图方式
    if (!is.null(data_group) && DE_filter == TRUE) {
      p <- ggplot(data_filtered, aes(x = key, y = log2(value + 1), fill = group)) + 
        geom_boxplot() +
        theme_classic() +
        stat_compare_means(aes(group = group), label = "p.signif", method = "wilcox.test", size = 3) + 
        theme(axis.text.x = element_text(angle = 90, hjust = 1, colour = "black", size = 10),
              axis.text.y = element_text(hjust = 1, colour = "black", size = 10),
              plot.title = element_text(hjust = 0.5)) +
        labs(title = title)
      
      # 提取统计结果
      ggplot_build_obj <- ggplot_build(p)
      stat_results <- ggplot_build_obj$data[[2]]
      stat_results_with_keys <- data_filtered %>% distinct(key) %>% bind_cols(stat_results)
      keys_with_high_pvalue <- stat_results_with_keys %>% filter(p > 0.05) %>% pull(key)
    } else if (is.null(data_group) || DE_filter == FALSE){
      p <- ggplot(data_filtered, aes(x = key, y = log2(value + 1))) + 
        geom_boxplot(fill = "lightblue") +
        theme_classic() +
        theme(axis.text.x = element_text(angle = 90, hjust = 1, colour = "black", size = 10),
              axis.text.y = element_text(hjust = 1, colour = "black", size = 10),
              plot.title = element_text(hjust = 0.5)) +
        labs(title = paste(title, "(no group matrix)"))
      
      keys_with_high_pvalue <- unique(data_filtered$key)  # 无分组时保留所有marker
    }
    
    return(list(data = data_filtered, keys_with_high_pvalue = keys_with_high_pvalue, plot = p))
  }
  
  correlate_markers <- function(data_filtered) {
    data_spread <- spread(data_filtered, key, value)
    numeric_data <- data_spread %>% dplyr::select(-id, -group)
    
    # 检查是否有足够多的蛋白进行相关性分析
    if (ncol(numeric_data) < 2) {
      warning(paste("Insufficient proteins for correlation analysis (only", 
                    ncol(numeric_data), "protein available)"))
      return(list(r = NA, P = NA, insufficient_proteins = TRUE))
    }
    corr_results <- cor(as.matrix(numeric_data), method = "pearson")
    return(list(r = corr_results, P = NULL, insufficient_proteins = FALSE))
  }
  filter_markers <- function(keys_with_high_pvalue, corr_matrix_r) {
    # 检查相关性分析是否成功
    if (is.na(corr_matrix_r)[1]) {
      return(list(high_corr_keys = character(0),
                  filtered_keys = character(0),
                  insufficient_proteins = TRUE))
    }
    
    # 提取相关性高的标记物
    high_corr_keys <- colnames(corr_matrix_r)[apply(corr_matrix_r, 2, 
                                                    function(x) any(abs(x) > cutoff & abs(x) < 1))]
    filtered_keys <- intersect(keys_with_high_pvalue, high_corr_keys)
    
    result <- list(high_corr_keys = high_corr_keys,
                   filtered_keys = filtered_keys,
                   insufficient_proteins = FALSE)
    return(result)
  }
  # 定义绘制水平图的函数
  # Plot and save level plots for filtered markers
  plot_level <- function(data_filtered, filtered_keys, title) {
    if (length(filtered_keys) == 0) {
      warning(paste("No markers found for", title))
      return()
    }
    data_filtered <- data_filtered %>% filter(key %in% filtered_keys)
    p <- ggplot(data_filtered, aes(x = id, y = log2(value + 1), fill = id)) + 
      geom_boxplot() +
      theme_classic() +
      stat_compare_means(aes(group = id), label = "p.format", method = "anova", size = 3) + 
      theme(axis.text.x = element_text(angle = 90, hjust = 1, colour = "black", size = 10),
            axis.text.y = element_text(hjust = 1, colour = "black", size = 10),
            plot.title = element_text(hjust = 0.5),
            legend.position="none") +
      labs(title = paste(title, "level"))
    return(p)
  }
  analyze_markers_full <- function(data_ggplot, marker_list, analyze_markers_result, correlation_result) {
    target_proteins <- marker_list$GN
    data_filtered <- data_ggplot %>% 
      filter(key %in% target_proteins)
    
    # 更健壮的方式检查基因是否存在
    exists_status <- ifelse(target_proteins %in% unique(data_ggplot$key), "Pass", "NA")
    
    # 更安全的方式提取p值
    p_values <- tryCatch({
      ggplot_build_obj <- ggplot_build(analyze_markers_result$plot)
      stat_data <- ggplot_build_obj$data[[2]]  # 假设统计结果在第二个图层
      stat_data$p
    }, error = function(e) rep(NA, length(target_proteins)))
    
    # 创建统计结果数据框
    stats <- data.frame(
      key = target_proteins,
      exists = exists_status,
      DE = ifelse(target_proteins %in% analyze_markers_result$keys_with_high_pvalue, 
                  "Pass", "Non-removable inter-sample heterogeneity"),
      correlation = ifelse(correlation_result$insufficient_proteins,
                           "Insufficient proteins for analysis",
                           ifelse(target_proteins %in% correlation_result$high_corr_keys,
                                  "Pass", "Low statistical correlation")),
      stringsAsFactors = FALSE
    )
    
    # 对于不存在的基因，将DE和correlation设为NA
    stats[stats$exists == "NA", c("DE", "correlation")] <- NA
    
    return(stats)
  }
  # Mann panel data input ----
  result <- list()
  result$rawdata <- data
  result$group <- data_group
  list_erythrocyte <- if (!is.null(custom_erythrocyte)) {
    data.frame(GN = custom_erythrocyte, ContaminationType = "erythrocyte")
  } else { 
    data.frame(GN = c("HBA1","HBB","CA1","HBD","PRDX2","CA2",
                      "CAT","SLC4A1","BLVRB","SPTA1","SPTB",
                      "ANK1","GAPDH","SOD1","PRDX6","BPGM",
                      "ACTB","ACTG1","SELENBP1","EPB41",
                      "LDHB","EIF5A","EIF5A2","EIF5AL1",
                      "PNP","ALDOA","TPI1","NIF3L1",
                      "UBC","UBB","RPS27A","UBA52",
                      "UBBP4","HBE1","HSPA8","PGK1","ALAD"),
               ContaminationType = "erythrocyte") }
  list_platelet <- if (!is.null(custom_platelet)) {
    data.frame(GN = custom_platelet, ContaminationType = "platelet")
  } else {
    data.frame(GN = c("FLNA","TLN1","MYH9","ACTB","VCL",
                      "ACTN1","TPM4","THBS1","TUBB1",
                      "YWHAZ","GSN","TUBA1B","ITGA2B",
                      "F13A1","PFN1","ITGB3","TAGLN2",
                      "FERMT3","RAP1B","PLEK","PPBP",
                      "GAPDH","MMRN1","MYL6","CFL1",
                      "PARVB","SDPR","TUBB4B","TMSB4X","PKM"),
               ContaminationType = "platelet")}
  list_coagulation <- if (!is.null(custom_coagulation)) {
    data.frame(GN = custom_coagulation, ContaminationType = "coagulation")
  } else {
    data.frame(GN = c("FGB","FGG","FGA","F13A1",
                      "SERPINC1","PPBP","F2",
                      "GP1BA","ECM1","CLU",
                      "DSP","WDR1","ATRN",
                      "GP5","SERPINA5","C1RL",
                      "MAN1A1","F11","KNG1",
                      "BCHE","GPLD1","THBS1",
                      "MASP1","TENM4","HSPG2",
                      "APOC3","CDH5","LTF","CST3",
                      "PROC","PF4","PF4V1"),
               ContaminationType = "coagulation")}
  # 检查缺失基因的函数
  check_missing_genes <- function(data, gene_list) {
    missing_genes <- gene_list[!gene_list %in% rownames(data)]
    return(missing_genes)
  }
  
  # 检查缺失的红细胞基因
  missing_erythrocyte_genes <- check_missing_genes(data, list_erythrocyte$GN)
  # 检查缺失的血小板基因
  missing_platelet_genes <- check_missing_genes(data, list_platelet$GN)
  # 检查缺失的凝血基因
  missing_coagulation_genes <- check_missing_genes(data, list_coagulation$GN)
  ## determine panel ----
  # 转换数据格式
  data <- as.data.frame(t(data))
  data$id <- rownames(data)
  
  # 处理无分组情况
  if (is.null(data_group)) {
    data_merge <- data
    data_merge$group <- "all_samples"  # 添加虚拟分组
  } else {
    data_merge <- merge(data_group, data, by = "id")
  }
  data <- subset(data,select = -c(id))
  # 转换为长格式
  data_ggplot <- tidyr::gather(data_merge,
                               key = "key", 
                               value = "value", 
                               -c("id", "group"))
  
  ## Perform analysis for erythrocyte, platelet, and coagulation markers ----
  ## 删除红细胞、血小板和凝血中在分组间存在差异的标记物 ----
  
  erythrocyte_results <- analyze_markers(data_ggplot, list_erythrocyte,
                                         "Erythrocyte marker")
  platelet_results <- analyze_markers(data_ggplot, list_platelet, 
                                      "Platelet marker")
  coagulation_results <- analyze_markers(data_ggplot, list_coagulation, 
                                         "Coagulation marker")
  
  ## Perform correlation analysis ----
  # 相关性分析
  corr_matrix_erythrocyte <- correlate_markers(erythrocyte_results$data)
  corr_matrix_platelet <- correlate_markers(platelet_results$data)
  corr_matrix_coagulation <- correlate_markers(coagulation_results$data)
  
  result_correlation <- list(
    erythrocyte = corr_matrix_erythrocyte,
    coagulation = corr_matrix_coagulation,
    platelet = corr_matrix_platelet
  )
  
  ## Filter markers based on p-value and correlation ----
  # 过滤标记物
  
  filtered_keys_erythrocyte <- filter_markers(erythrocyte_results$keys_with_high_pvalue,
                                              corr_matrix_erythrocyte$r)
  filtered_keys_platelet <- filter_markers(platelet_results$keys_with_high_pvalue, 
                                           corr_matrix_platelet$r)
  filtered_keys_coagulation <- filter_markers(coagulation_results$keys_with_high_pvalue,
                                              corr_matrix_coagulation$r)
  
  
  marker_list <- list(
    erythrocyte = filtered_keys_erythrocyte$filtered_keys,
    coagulation = filtered_keys_coagulation$filtered_keys,
    platelet = filtered_keys_platelet$filtered_keys
  )
  
  ## Plot and save level plots and result for filtered markers ----
  # 绘制并保存水平图
  Erythrocyte_pannel <- plot_level(erythrocyte_results$data, 
                                   marker_list$erythrocyte, 
                                   "Erythrocyte pannel")
  Platelet_pannel <-  plot_level(platelet_results$data,
                                 marker_list$platelet,
                                 "Platelet pannel")
  Coagulation_pannel <- plot_level(coagulation_results$data, 
                                   marker_list$coagulation, 
                                   "Coagulation pannel")
  stats_erythrocyte <- analyze_markers_full(data_ggplot = data_ggplot,
                                            marker_list = list_erythrocyte, 
                                            analyze_markers_result = erythrocyte_results,
                                            correlation_result = filtered_keys_erythrocyte)
  stats_platelet <- analyze_markers_full(data_ggplot = data_ggplot,
                                         marker_list = list_platelet, 
                                         analyze_markers_result = platelet_results,
                                         correlation_result = filtered_keys_platelet)
  stats_coagulation <- analyze_markers_full(data_ggplot = data_ggplot,
                                            marker_list = list_coagulation, 
                                            analyze_markers_result = coagulation_results,
                                            correlation_result = filtered_keys_coagulation)
  
  # 保存结果
  result <- list(
    rawdata = t(data),
    group = data_group,
    data = list(
      erythrocyte = spread(erythrocyte_results$data, key, value),
      coagulation = spread(coagulation_results$data, key, value),
      platelet = spread(platelet_results$data, key, value)
    ),
    plot_marker = list(erythrocyte = Erythrocyte_pannel,
                       coagulation = Coagulation_pannel,
                       platelet = Platelet_pannel),
    plot_contamination = list(erythrocyte = erythrocyte_results$plot,
                              coagulation = coagulation_results$plot,
                              platelet = platelet_results$plot),
    missing_genes = list(missing_erythrocyte_genes = missing_erythrocyte_genes,
                         missing_coagulation_genes = missing_coagulation_genes,
                         missing_platelet_genes = missing_platelet_genes),
    correlation = result_correlation,
    marker_list = marker_list,
    marker_stats = list(
      stats_erythrocyte = stats_erythrocyte,
      stats_platelet = stats_platelet,
      stats_coagulation = stats_coagulation)
  )
  return(result)
}
