#' Calculates the relative occupancy of peptide data
#'
#' @description Calculates the relative occupancy of peptide data based on the protein_input data
#'
#' @param peptides Peptide abundances that are compared to input protein abundances to generate relative occupancies
#' @param protein_input Input protein abundances used to calculate corresponding peptide occupancies
#'
#' @return peptides_relative_occupancy
#'
#' @export

occupancy <- function(peptides, protein_input){
  header <- select(peptides, grep("Group = ", names(peptides), invert = TRUE))
  peptides <- left_join(peptides,
                        select(protein_input, "Accession", grep("Group = ", names(protein_input))),
                        by = "Accession")

  index <- grep("Input", names(peptides))
  names(peptides)[index] <- paste("Corresponding Protein Abundance:", names(peptides[index]))
  peptides_relative_occupancy <- select(peptides, setdiff(grep("Group = ", names(peptides)),
                                                          grep("Corresponding Protein Abundance:", names(peptides))))
  relative_proteins <- select(peptides, intersect(grep("Group = ", names(peptides)),
                                                  grep("Corresponding Protein Abundance:", names(peptides))))
  peptides_relative_occupancy <- peptides_relative_occupancy - relative_proteins
  names(peptides_relative_occupancy) <- paste("Relative Occupancy", names(peptides_relative_occupancy))
  peptides_relative_occupancy <- bind_cols(header, peptides_relative_occupancy)

  return(peptides_relative_occupancy)
}
