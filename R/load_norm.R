#' Load normalize proteomic data
#'
#' @description Load normalizes searched Proteome Discoverer data
#'
#' @param proteins_raw A raw data frame with protein abundances
#' @param peptides_raw A raw data frame with peptide abundances
#' @param num_input Number of input kits in the experiment
#' @param num_group Number of experimental plexes
#' @param num_exp Number of experimental subgroups per plex
#' @param num_pool Number of num_pool groups
#'
#' @return data_load_normalized return the load normalized proteomic data
#'
#' @export

load_normalize <- function(proteins_raw, peptides_raw, num_input, num_group, num_exp, num_pool){
  proteins_raw <- peptides_raw <- num_input <- num_group <- num_exp <- num_pool <- NULL
  # Calculate the normalization coefficients based on input columns of the raw peptide data
  norm_coefs <- vector(mode = "list", length = num_input)
  input_string_ids <- vector(mode = "list", length = num_input)
  for (i in 1:num_input) {
    user_input <-
      svDialogs::dlg_input("Enter a string that is specific to and common among the input column names in the raw PEPTIDE data.")$res
    input_string_ids[[i]] <- user_input
    input_index <- grep(user_input, names(peptides_raw))
    sum <- numeric(length(input_index))
    for (j in 1:length(input_index)) {
      sum[j] <- sum(peptides_raw[input_index[j]], na.rm = TRUE)
    }
    avg <- mean(sum)
    norm_coefs[[i]] <- numeric(length(input_index))
    for (j in 1:length(norm_coefs[[i]])) {
      norm_coefs[[i]][j] <- sum[j] / avg
    }
  # generate report of normalization coefficient
  }

  # Initialize the peptides_load_normalized data structure
  index_list <- vector(mode = "list", length = num_group)
  group_string_ids <- vector(mode = "list", length = num_group)
  for(i in 1:num_group){
    user_input <-
      svDialogs::dlg_input("Enter a string that is specific to and common among all column names for your experimental group PEPTIDE data.")$res
    group_string_ids[[i]] <- user_input
    index_list[[i]] <- grep(user_input, names(peptides_raw))
  }
  peptides_load_normalized <- dplyr::bind_cols(tibble::tibble(unique_id = 1:nrow(peptides_raw)),
                                        dplyr::select(peptides_raw, "Master Protein Accessions"),
                                        dplyr::select(peptides_raw, "Sequence"))
  peptides_load_normalized$Accession <- stringr::str_split_fixed(peptides_load_normalized$`Master Protein Accessions`, ";", 2)[,1]
  peptides_load_normalized <- dplyr::left_join(peptides_load_normalized, dplyr::select(proteins_raw, Accession, Description), by = "Accession")
  peptides_load_normalized$GN <- stringr::str_extract(peptides_load_normalized$Description, "(?<=GN=).*?(?=\\s|$)")
  peptides_load_normalized_list <- vector(mode = "list", length = num_group)

  # Calculate the load normalized peptide abundances and add experimental group identifiers (A, B, C...)
  group_names <- toupper(letters[1:num_exp])
  exp_size <- (length(index_list[[1]]) - num_pool) / num_exp
  cur_exp = 0
  k = 1
  for(i in 1:num_group){
    df <- dplyr::bind_cols(unique_id = 1:nrow(peptides_raw),
                    dplyr::select(peptides_raw, index_list[[i]]))
    for(j in 2:length(df)){
      df[j] <- df[j]/norm_coefs[[i]][j-1]
      names(df)[j] <- paste("Load Normalized ", names(df)[j], "Group =", group_names[k])
      cur_exp <- cur_exp + 1
      if (cur_exp == exp_size){
        cur_exp <- 0
        k <- k + 1
      }
    }
    peptides_load_normalized_list[[i]] <- dplyr::left_join(peptides_load_normalized, df, by = "unique_id")
  }

  # Initialize the proteins_load_normalized data structure
  proteins_load_normalized <- dplyr::bind_cols(tibble::tibble(unique_id = 1:nrow(proteins_raw)),
                                        dplyr::select(proteins_raw, "Accession"),
                                        dplyr::select(proteins_raw, "Sequence"),
                                        dplyr::select(proteins_raw, "Master"))
  proteins_load_normalized_list <- vector(mode = "list", length = num_input)

  # Calculate the load normalized protein abundances and add experimental group identifiers (A, B, C...)
  group_names <- toupper(letters[1:num_exp])
  exp_size <- (length(index_list[[1]]) - num_pool) / num_exp
  cur_exp = 0
  k = 1
  for (i in 1:num_input){
    input_index <- grep(input_string_ids[[i]], names(proteins_raw))
    df <- dplyr::bind_cols(unique_id = 1:nrow(proteins_raw),
                    dplyr::select(proteins_raw, input_index))
    for(j in 2:length(df)){
      df[j] <- df[j]/norm_coefs[[i]][j-1]
      names(df)[j] <- paste("Load Normalized ", names(df)[j], "Group =", group_names[k])
      cur_exp <- cur_exp + 1
      if (cur_exp == exp_size){
        cur_exp <- 0
        k <- k + 1
      }
    }
    proteins_load_normalized_list[[i]] <- dplyr::left_join(proteins_load_normalized, df, by = "unique_id")
  }

  # Combine load normalized protein and peptide abundances and place in a single returnable data
  data_load_normalized <- list("proteins_load_normalized" = proteins_load_normalized_list,
                               "peptides_load_normalized" = peptides_load_normalized_list)
  return(data_load_normalized)
}

