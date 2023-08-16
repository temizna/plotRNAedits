#' Retrieve base edits from vcf
#'This function is a variation of mutations_from_vcf from MutationalPatterns
#' A function to extract base substitutions of each position in vcf
#' @param vcf A CollapsedVCF object
#' @return Character vector with base substitutions
#' @import GenomicRanges
#'
#' @examples
#' ## See the 'read_vcfs_as_granges()' example for how we obtained the
#' ## following data:
#' vcfs <- readRDS(system.file("states/read_vcfs_as_granges_output.rds",
#'   package = "MutationalPatterns"
#' ))
#'
#' muts <- mutations_from_vcf(vcfs[[1]])
#' @seealso
#' \code{\link{read_vcfs_as_granges}}
#'
#' @export
edits_from_vcf<-function (vcf)
{
  #.check_no_indels(vcf)
  ref <- as.character(get_rref(vcf))
  alt <- as.character(unlist(get_ralt(vcf)))
  muts <- paste(ref, alt, sep = ">")
  return(muts)
}

