#' Expression correlation analysis and visualization
#'
#' @param exprMatrix Input matrix or correlation matrix
#' @param input_type Input type ("matrix"/"correlation", default "matrix")
#' @param corMethod Correlation calculation method (only when input_type="matrix")
#' @param colors Color gradient (default blue-white-red)
#' @param displayNumbers Whether to display correlation values
#' @param title Plot title
#' @param cluster_rows Whether to cluster rows
#' @param cluster_cols Whether to cluster columns
#' @return Returns correlation matrix and heatmap object
#' @export
#'
#' @examples
#' # Example 1: Input data matrix
#' expr <- matrix(rnorm(1000), nrow=100, 
#'               dimnames=list(paste0("Gene",1:100), paste0("Sample",1:10)))
#' plot_expression_correlation(expr, input_type = "matrix")
#'
#' # Example 2: Input pre-computed correlation matrix
#' cor_matrix <- cor(expr)
#' plot_expression_correlation(cor_matrix, 
#'                           input_type = "correlation", 
#'                           title = "Pre-computed Correlation")
#' 
#' # Example 3: Small matrix (auto-disable clustering)
#' small_matrix <- matrix(rnorm(4), nrow=2)
#' plot_expression_correlation(small_matrix)
plot_expression_correlation <- function(exprMatrix,
                                        input_type = "matrix",
                                        corMethod = "pearson",
                                        colors = colorRampPalette(c("blue", "white", "red"))(100),
                                        displayNumbers = FALSE,
                                        title = "Expression Correlation Matrix",
                                        cluster_rows = TRUE,
                                        cluster_cols = TRUE) {
  
  # Check required packages
  if (!requireNamespace("pheatmap", quietly = TRUE)) {
    stop("Package 'pheatmap' is required: install.packages('pheatmap')")
  }
  
  # Process data based on input type
  if (input_type == "matrix") {
    # Calculate correlation matrix
    cor_matrix <- cor(exprMatrix, method = corMethod)
  } else if (input_type == "correlation") {
    # Use input as correlation matrix directly
    cor_matrix <- exprMatrix
  } else {
    stop("Invalid input_type. Please use 'matrix' or 'correlation'.")
  }
  
  # Handle 1x1 matrix case
  if (all(dim(cor_matrix) == c(1, 1))) {
    # Create a simple plot for 1x1 matrix
    plot(1, 1, type = "n", axes = FALSE, xlab = "", ylab = "", main = title)
    text(1, 1, labels = paste("Correlation:","(only:",rownames(exprMatrix),")",round(cor_matrix[1,1], 2)), cex = 2)
    return(list(correlation_matrix = cor_matrix, plot = recordPlot()))
  }
  
  # Auto-handle clustering for small matrices
  if (nrow(cor_matrix) < 3) {
    message("Matrix too small for clustering (n < 3), disabling clustering")
    cluster_rows <- FALSE
    cluster_cols <- FALSE
  }
  
  # Set up number display
  annotation <- if (displayNumbers) round(cor_matrix, 2) else FALSE
  
  # Create heatmap
  p <- pheatmap::pheatmap(
    mat = cor_matrix,
    color = colors,
    display_numbers = annotation,
    number_color = "black",
    main = title,
    cluster_rows = cluster_rows,
    cluster_cols = cluster_cols,
    border_color = NA,
    show_rownames = TRUE,
    show_colnames = TRUE
  )
  
  # Return results
  return(list(correlation_matrix = cor_matrix, plot = p))
}
