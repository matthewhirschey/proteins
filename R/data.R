#' proteins
#'
#' A dataset containing human proteome information from Uniprot.
#'
#' @usage data(proteins)
#'
#' @docType data
#'
#' @format A data frame with 20,430 rows and 8 variables:
#' \describe{
#'   \item{uniprot_id}{Official uniprot entry id}
#'   \item{gene_name}{Gene name}
#'   \item{gene_name_alt}{Alternative gene names associated with entry}
#'   \item{protein_name}{Protein name}
#'   \item{protein_name_alt}{Alternative protein names associated with entry}
#'   \item{sequence}{Primary amino acid sequencxe}
#'   \item{length}{Number of amino acids}
#'   \item{mass}{Molecular weight of protein}
#'   }
#'
#' @keywords datasets
#'
#' @source \url{https://www.uniprot.org/uniprot/?query=*&fil=organism%3A%22Homo+sapiens+%28Human%29+%5B9606%5D%22+AND+reviewed%3Ayes}
"proteins"
