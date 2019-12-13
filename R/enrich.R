#' Create a plot of the enriched gene sets
#'
#' @description Generates a plot of the top n significantly enriched gene sets based on your omic data
#'
#' @param data A dataframe with 3 columns: gene number, log2-fold change, and p_value
#' @param topn Number of gene sets that you want to visualize
#' @param p_cutoff Maximum p-value to subset dataframe
#' @param enrich_p_cutoff Maximum p-value of enriched gene sets to be invluded in the plot
#'
#' @return enrich_plot A plot of enriched gene sets
#'
#' @export

enrich <- function(data, topn, p_cutoff = 1, enrich_p_cutoff = 1){
  data <- na.omit(data)
  names(data)[2] <- "LOG2_FOLD_CHANGE"
  names(data)[3] <- "P_VALUE"

  dat <- subset(data, `P_VALUE` < p_cutoff)
  dat <- top_n(dat, topn, `LOG2_FOLD_CHANGE`)
  dat <- subset(dat, `LOG2_FOLD_CHANGE` > 0)

  View(dat)

  if(nrow(dat) > 0){
    dat <- enrichr(genes = dat$GN, databases = "KEGG_2019_Human")$KEGG_2019_Human
    dat <- subset(dat, Adjusted.P.value < enrich_p_cutoff)

    enrich_plot <- ggplot(data = dat, aes(x=reorder(Term, -Adjusted.P.value), y = Adjusted.P.value)) +
            geom_col(aes(fill = Adjusted.P.value)) +
            coord_flip() +
            scale_fill_gradient(low = "green4",
                                high = "grey95") +
            xlab("Pathway") +
            ylab("Adjusted P-Value") +
            theme_light()
  }
}
