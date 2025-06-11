
QC_boxplot <- function(data, data_group = NULL) {

  library(ggplot2)
  library(tidyr)

  # Check if data_group is NULL or an empty data frame
  if (is.null(data_group) || (is.data.frame(data_group) && nrow(data_group) == 0)) {
    data_ggplot <- as.data.frame(t(data))
    data_ggplot$id <- rownames(data_ggplot)
    data_ggplot <- tidyr::gather(data_ggplot, key = "key", value = "value", -c("id"))

    ggplot(data_ggplot, aes(x = id, y = log2(value))) +
      geom_boxplot() +
      theme_classic() +
      theme(axis.text.x = element_text(angle = 90, hjust = 1, colour = "black", size = 10),
            axis.text.y = element_text(hjust = 1, colour = "black", size = 10),
            plot.title = element_text(hjust = 0.5))
  } else {
    data_ggplot <- as.data.frame(t(data))
    data_ggplot$id <- rownames(data_ggplot)
    data_ggplot <- merge(data_group, data_ggplot, by = 'id')
    data_ggplot <- tidyr::gather(data_ggplot, key = "key", value = "value", -c("id", "group"))

    data_ggplot <- data_ggplot[order(data_ggplot$group), ]
    data_ggplot$id <- factor(data_ggplot$id, levels = unique(data_ggplot$id))
    data_ggplot <- data_ggplot[order(data_ggplot$id), ]

    if (!is.factor(data_ggplot$group)) {
      data_ggplot$group <- factor(data_ggplot$group)
    }

    ggplot(data_ggplot, aes(x = id, y = log2(value), fill = group)) +
      geom_boxplot() +
      theme_classic() +
      theme(axis.text.x = element_text(angle = 90, hjust = 1, colour = "black", size = 10),
            axis.text.y = element_text(hjust = 1, colour = "black", size = 10),
            plot.title = element_text(hjust = 0.5)) 
  }
}
