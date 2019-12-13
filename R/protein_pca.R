#' Principle component analysis
#'
#' @description: Returns an exploratory PCA analysis, along with a score plot and a scree plot
#' @param tidyData Normalized data to perform the PCA analysis on
#' @param components The number of components
#' @param center Whether the data is centered or not
#' @param scale Whether the data is scaled or not
#' @return A list containing the principal components and the ggplots of the score plot and the scree plot
#' @author Matthew Huang

protein_pca <- function(tidyData, components = 4, center = FALSE, scale = FALSE){
  Group <- PC1 <- PC2 <- Groups <- Principal_Component <- Percent_Variance_Explained <- NULL

  colnames(tidyData)[1:2] <- c("ID", "Group")

  tidyData <- tidyData %>% dplyr::mutate(Group = as.factor(Group))

  df <- as.matrix(sapply(tidyData[3: ncol(tidyData)], as.numeric))

  X <- as.matrix(df)
  Y <- tidyData$Group

  pca_res <- mixOmics::pca(X, ncomp = components, center = center, scale = scale)
  PCi <- data.frame(pca_res$x, Groups = Y)

  scoresplot <- ggplot2::ggplot(PCi, aes(x = PC1, y = PC2, col = Groups))+
    ggplot2::geom_point(size = 3, alpha = 0.5) +
    ggplot2::xlab(paste0("PC1 (", round(100*(pca_res$explained_variance)[1], 2), "%)")) +
    ggplot2::ylab(paste0("PC2 (", round(100*(pca_res$explained_variance)[2], 2), "%)")) +
    ggplot2::theme_minimal()

  eigenvalues <- data.frame(Percent_Variance_Explained = round(pca_res$explained_variance*100, 4))
  eigenvalues <- eigenvalues %>% tibble::rownames_to_column("Principal_Component")

  screeplot <- ggplot2::ggplot(eigenvalues, aes(x = Principal_Component, y = Percent_Variance_Explained, fill = NULL)) +
    ggplot2::geom_bar(stat = "identity", fill = rep(c("lightblue"), nrow(eigenvalues))) +
    ggplot2::xlab("Principal Component") +
    ggplot2::ylab("% Variance Explained") +
    ggplot2::geom_label(data = eigenvalues, aes(label =  paste0(round(Percent_Variance_Explained, 3), "%"))) +
    ggplot2::theme_minimal()

  score_data <- PCi %>% dplyr::select(-Groups) %>% round(4)

  return(list(screeplot = screeplot, scoresplot = scoresplot,
              score_data = score_data, eigenvalues = eigenvalues))
}
