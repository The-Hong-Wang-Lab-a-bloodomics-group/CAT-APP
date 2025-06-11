#' Create Venn Diagram for Two Sets
#'
#' Generates a proportional Venn diagram comparing two sets with automatic error handling 
#' and log file cleanup.
#'
#' @param set1 First vector of elements OR list containing both sets
#' @param set2 Second vector of elements (ignored if set1 is list)
#' @param categories Character vector of category names (length 2)
#' @param title Diagram title (default: "Venn Diagram")
#' @param colors Fill colors for sets (default: c("#ae6b81", "#6982b9"))
#' @param alpha Transparency level (0-1, default: 0.5)
#' @param print.mode Value display mode: "raw" (counts) or "percent" (default: "raw")
#' @param save.plot Save to file? (default: FALSE)
#' @param filename Output filename when save.plot=TRUE (default: "venn_diagram.png")
#'
#' @return Invisible list containing:
#'   - total_unique: Total unique elements
#'   - overlap_count: Number of overlapping elements
#'   - overlap_elements: Vector of overlapping elements
#'   Also draws Venn diagram to active device
#'
#' @details
#' Special features:
#' 1. Automatic cleanup of VennDiagram log files
#' 2. Handles empty sets gracefully
#' 3. Uses proportional scaling (scaled=TRUE) for accurate area representation
#' 4. Includes error handling for failed diagram generation
#' 5. Provides both visual output and structured data return
#'
#' @note
#' - When input is a list, set1 should contain both vectors (set2 ignored)
#' - Output includes both visual plot and structured overlap data
#' - Automatically removes NA values from input sets
#'
#' @examples
#' \dontrun{
#' # Compare two vectors
#' venn_plot(
#'   set1 = c("A", "B", "C"),
#'   set2 = c("B", "C", "D"),
#'   categories = c("Group1", "Group2")
#' )
#' 
#' # Compare using list input
#' venn_plot(
#'   set1 = list(SetA = 1:100, SetB = 80:150),
#'   title = "Number Overlap"
#' )
#' 
#' # Save to file
#' venn_plot(set1 = genes1, set2 = genes2, 
#'           save.plot = TRUE, filename = "gene_overlap.png")
#' }
#'
#' @importFrom VennDiagram draw.pairwise.venn
#' @importFrom grid grid.newpage grid.draw grid.text
#' @importFrom ggplot2 ggsave
#' @export
venn_plot <- function(set1, set2, 
                      categories = c("Set A", "Set B"),
                      title = "Venn Diagram",
                      colors = c("#ae6b81", "#6982b9"),
                      alpha = 0.5,
                      print.mode = c("raw", "percent"),
                      save.plot = FALSE,
                      filename = "venn_diagram.png") {
  
  # 强制退出时清理日志的保险机制
  on.exit({
    log_files <- list.files(
      pattern = "^VennDiagram\\.[0-9]{4}-[0-9]{2}-[0-9]{2}_[0-9]{2}-[0-9]{2}-[0-9]{2}\\.[0-9]+\\.log$",
      ignore.case = TRUE
    )
    if (length(log_files) > 0) {
      suppressWarnings(file.remove(log_files))
    }
  })
  
  # 检查并加载必要包
  if (!requireNamespace("VennDiagram", quietly = TRUE)) {
    install.packages("VennDiagram")
    library(VennDiagram)
  }
  
  # 处理输入并移除NA
  if (is.vector(set1) && is.vector(set2)) {
    set_list <- list(
      set1 = unique(na.omit(set1)),
      set2 = unique(na.omit(set2))
    )
  } else if (is.list(set1) && length(set1) == 2) {
    set_list <- lapply(set1, function(s) unique(na.omit(s)))
  } else {
    stop("输入应为两个向量或包含两个集合的列表")
  }
  
  # 检查空集合
  if (any(lengths(set_list) == 0)) {
    grid::grid.newpage()
    grid::grid.text("No elements in one or both sets", gp = grid::gpar(fontsize = 16))
    return(invisible(list(
      total_unique = length(unique(unlist(set_list))),
      overlap_count = 0,
      overlap_elements = character(0)
    )))
  }
  
  overlap <- intersect(set_list[[1]], set_list[[2]])
  
  # 统一使用draw.pairwise.venn并添加scaled参数
  venn <- tryCatch({
    VennDiagram::draw.pairwise.venn(
      area1 = length(set_list[[1]]),
      area2 = length(set_list[[2]]),
      cross.area = length(overlap),
      category = categories,
      fill = colors,
      alpha = alpha,
      lty = "blank",
      cex = 1.5,
      cat.cex = 1.3,
      cat.pos = c(-30, 30),
      cat.dist = c(0.05, 0.05),
      margin = 0.05,
      print.mode = print.mode,
      ind = FALSE,
      scaled = TRUE  # 关键参数，允许自动调整比例
    )
  }, error = function(e) {
    message("绘图错误: ", e$message)
    return(NULL)
  })
  
  if (is.null(venn)) {
    grid::grid.newpage()
    grid::grid.text("无法生成韦恩图", gp = grid::gpar(fontsize = 16))
  } else {
    grid::grid.newpage()
    grid::grid.draw(venn)
    grid::grid.text(title, y = 0.95, gp = grid::gpar(fontsize = 16))
  }
  
  if (save.plot) {
    ggplot2::ggsave(filename, plot = venn, 
                    width = 8, height = 6, 
                    dpi = 300)
    message("图形已保存为：", filename)
  }
  
  invisible(list(
    total_unique = length(unique(unlist(set_list))),
    overlap_count = length(overlap),
    overlap_elements = overlap
  ))
}
