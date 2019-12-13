#' Volcano plot
#'
#' @description Generate volcano plots based for a dataset
#'
#' @param data A dataframe with 3 columns: gene number, log2-fold change, and p_value
#' @param p_val_cutoff A p-value cutoff to visualize the data set
#'
#' @return volcano Volcano plot
#'
#' @export
#'

volcano <- function(data, p_val_cutoff){
  data <- na.omit(data)
  names(data)[2] <- "REL_LOG2_FOLD_CHANGE"
  names(data)[3] <- "REL_P_VALUE"
  data <- mutate(data, NEG_LOG_P = -log10(data$`REL_P_VALUE`))
  data <- mutate(data, Significance = ifelse(`REL_P_VALUE`< p_val_cutoff, paste("p-value < ", p_val_cutoff), "Not Significant"))

  volcano <- ggplot(data, aes(`REL_LOG2_FOLD_CHANGE`, `NEG_LOG_P`)) +
                      geom_point(aes(col=Significance),
                                 alpha = 0.3,
                                 na.rm = TRUE) +
                      scale_color_manual(values = c("grey50", "red")) +
                      geom_text_repel(data = top_n(data, 20, `NEG_LOG_P`),
                                      aes(label=`GN`),
                                      segment.size = 0.2,
                                      segment.color = "grey50",
                                      point.padding = 1) +
                      ylab(~-Log[10]~ "( p-value )") +
                      xlab("Relative" ~Log[2]~ "Fold Change") +
                      theme_light()
  return(volcano)
}
