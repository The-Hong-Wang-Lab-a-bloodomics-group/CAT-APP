log_run_event <- function() {
  log_file <- "run_check.log"
  
  # 如果日志文件不存在，创建并初始化
  if (!file.exists(log_file)) {
    writeLines("timestamp,count", log_file)
  }
  
  # 读取现有日志
  log_data <- tryCatch(
    read.csv(log_file, stringsAsFactors = FALSE),
    error = function(e) data.frame(timestamp = character(), count = integer())
  )
  
  # 计算新计数
  new_count <- if (nrow(log_data) > 0) max(log_data$count, na.rm = TRUE) + 1 else 1
  new_timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  
  # 添加新记录
  new_entry <- data.frame(timestamp = new_timestamp, count = new_count)
  updated_log <- rbind(log_data, new_entry)
  
  # 写回日志文件
  write.csv(updated_log, log_file, row.names = FALSE)
  
  # 返回当前计数
  return(new_count)
}
