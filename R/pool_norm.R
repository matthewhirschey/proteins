#' Pool normalizations
#'
#' @description Calculate the pool normalization for a proteomics data set
#'
#' @param data The dataframe that you want to normalize, with all of the samples you want to normalize in the last columns
#' @param num_sample The number of samples in your dataframe that your are normalizing
#'
#' @return data_pool_normalized
#'
#' @export

pool_normalize <- function(data, num_sample) {
    head <- select(data, (1:(ncol(data) - num_sample)))
    data <- select(data, tail(1:ncol(data), num_sample))
    pool_index <- grep("Group = NA", names(data))
    pool_avg <- rowMeans(data[pool_index])
    data <- data - pool_avg
    names(data) <- paste("Pool Normalized", names(data))
    data_pool_normalized <- bind_cols(head, data)
    return(data_pool_normalized)
}
