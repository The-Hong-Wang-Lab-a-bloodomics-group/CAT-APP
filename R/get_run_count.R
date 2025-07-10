# 定义读取日志函数
get_run_count <- function() {
  log_file <- "run_check.log"
  
  if (!file.exists(log_file)) {
    return(0)
  }
  
  log_data <- tryCatch(
    read.csv(log_file, stringsAsFactors = FALSE),
    error = function(e) data.frame(timestamp = character(), count = integer())
  )
  
  if (nrow(log_data) > 0) {
    return(max(log_data$count, na.rm = TRUE))
  } else {
    return(0)
  }
}
