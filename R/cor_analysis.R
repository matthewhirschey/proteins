#' Correlation analysis
#'
#' @description: Returns simple statistics, along with a correlation matrix and heatmap
#' @param tidyData Normalized protein abundance data to perform the correlation analysis on, with columns as protein-coding genes and rows as cell lines
#' @param gene1 The first gene of interest in statistics calculation
#' @param gene2 The second gene of interest in statistics calculation
#' @param method Indicates correlation coefficient to be used. One of "pearson" "kendall" or "spearman"
#' @return Summary of statistics between 2 genes, correlation matrix, and heatmap
#' @author Harshavardhan Srijay
#' @export


cor_analysis <- function(tidyData, gene1, gene2, method){
  tidyData <- na.omit(tidyData)
  colnames(tidyData)[1:2] <- c("ID", "Group")

  tidyData <- tidyData %>% dplyr::mutate(Group = as.factor(Group))

  df <- as.matrix(sapply(tidyData[3: ncol(tidyData)], as.numeric))
  
  stats <- cor.test(as.numeric(df$gene1), as.numeric(df$gene2), 
                    method = method, use = "complete.obs")
  
  drops <- c("ID", "Group")
  df <- df[,!(names(df) %in% drops)]
  c.data <- as.matrix(round(cor(sapply((df), as.numeric)), 3))
  c <- corrplot(c.data, method = "color")
  
  return(list(Statistics = stats, CorrelationMatrix = c.data, Heatmap = c))
  
}
