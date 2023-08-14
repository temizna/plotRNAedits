#' Retrieve base RNA edit types from a VCF object
#' This function is a variation of mut_type from MutationalPatterns
#' A function to extract the base substitutions from a vcf and translate to
#' the 12 base substitution types. This function does not collapse G/A to C/U
#'
#' @param vcf A CollapsedVCF object
#' @return Character vector with base substitution types
#'
#' @examples
#' ## See the 'read_vcfs_as_granges()' example for how we obtained the
#' ## following data:
#' vcfs <- readRDS(system.file("states/read_vcfs_as_granges_output.rds",
#'   package = "MutationalPatterns"
#' ))
#'
#' mut_type(vcfs[[1]])
#' @seealso
#' \code{\link{read_vcfs_as_granges}}
#'
#' @export
edit_type=function (vcf) 
{
  muts <- mutations_from_vcf(vcf)
  types <- unlist(muts)
  types <- gsub("T>A", "U>A", types)
  types <- gsub("T>C", "U>C", types)
  types <- gsub("T>G", "U>G", types)
  types <- gsub("A>T", "A>U", types)
  types <- gsub("C>T", "C>U", types)
  types <- gsub("G>T", "G>U", types)
  
  return(types)
  return(types)
}

