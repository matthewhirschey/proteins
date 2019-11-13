#' Matching protein sequence motifs
#'
#' @description find_motif() offers an intuitive method to look-up a protein motif in proteomic data.
#'
#' @param data A vector containing protein sequence
#' @param string The sequence to match using standard regular expressions
#'
#' @export
#'
#' @return A vector of boolean values showing matched strings
#' @author Matthew Hirschey, \email{matthew@@hirschey.org}

find_motif <- function(data, string) {
  stringr::str_detect(data, as.character(string))
}


