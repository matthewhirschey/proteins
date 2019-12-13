#' Univariate statistics
#'
#' @description Calculate the univariate statistics for the data set
#'
#' @param data Dataframe used to calculate univariate statistics
#' @param num_sample Number of samples
#' @param num_exp Number of experimental groups
#'
#' @return data_univariate Dataframe including abundances and univariate statistics
#'
#' @export
#'

univariate <- function(data, num_sample, num_exp, num_pool){
  head <- select(data, (1:(ncol(data) - num_sample)))
  data <- select(data, tail(1:ncol(data), num_sample))

  all_index <- array(1:num_sample)
  group_names <- toupper(letters[1:num_exp])
  group_size <- (num_sample-num_pool)/num_exp

  for(i in 0:(num_exp-1)){
    cur_index <- all_index[(i*group_size+1):((i*group_size)+group_size)]
    avg_name <- paste(group_names[i+1], "AVG")
    sd_name <- paste(group_names[i+1], "STDEV")
    data <- add_column(data,
                     C1 = rowMeans(data[cur_index]),
                     C2 = apply(data[cur_index], 1, sd))
    names(data)[names(data) == "C1"] <- avg_name
    names(data)[names(data) == "C2"] <- sd_name
  }

  for(i in 1:num_exp){
    if(i==num_exp){
      break
    }
    for(j in (i+1):num_exp){
      fc_name <- paste(group_names[j],
                       "/",
                       group_names[i],
                       "Log2 Fold Change")
      data <- add_column(data,
                         C1 = data[grep(paste(group_names[j], "AVG"), names(data))] -
                           data[grep(paste(group_names[i], "AVG"), names(data))])
      data$C1 <- unlist(data$C1)
      names(data)[names(data) == "C1"] <- fc_name

      data <- add_column(data,
                         `P-Value` = NA)

      X_index <- grep(paste("Group =", group_names[i]), names(data))
      Y_index <- grep(paste("Group =", group_names[j]), names(data))

      data$`P-Value` <- unlist(row_t_equalvar(data[X_index],
                                                    data[Y_index])$pvalue)
      data$`Adjusted P-Value` <- p.adjust(data$`P-Value`, "fdr")
      names(data)[names(data) == "P-Value"] <- paste(group_names[j],
                                                 "/",
                                                 group_names[i],
                                                 "P-Value")
      names(data)[names(data) == "Adjusted P-Value"] <- paste(group_names[j],
                                                          "/",
                                                          group_names[i],
                                                          "Adjusted P-Value")

    }
  }

  data_univariate <- bind_cols(head, data)
  return(data_univariate)
}
