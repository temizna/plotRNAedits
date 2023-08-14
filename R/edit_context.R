#' Retrieve context of RNA base editions 
#' This function is a variation of mut_context from MutationalPatterns
#' A function to extract the bases 3' upstream and 5' downstream of the base
#' substitutions from the reference genome. Currently only 1 base extention
#' is supported.
#'
#' @param vcf A Granges object
#' @param ref_genome Reference genome
#' @param extension The number of bases, that's extracted upstream and
#' downstream of the base substitutions. (Default: 1).
#' @return Character vector with the context of the base substitutions
#'
#' @examples
#' ## See the 'read_vcfs_as_granges()' example for how we obtained the
#' ## following data:
#' vcfs <- readRDS(system.file("states/read_vcfs_as_granges_output.rds",
#'   package = "MutationalPatterns"
#' ))
#'
#' ## Load the corresponding reference genome.
#' ref_genome <- "BSgenome.Mmusculus.UCSC.mm10"
#' library(ref_genome, character.only = TRUE)
#'
#' ## Get the standard context
#' edit_context <- edit_context(vcfs[[1]], ref_genome)
#'
#' @export
edit_context<-function (vcf, ref_genome, extension = 1) 
{
  if (length(vcf) == 0) {
    warning("Detected empty GRanges object.\n                Returning an empty list for this sample.", 
            call. = FALSE)
    res <- list(types = NULL, context = NULL)
    return(res)
  }
  edit_context.mat <- mut_context(vcf, ref_genome, extension)
  muts <- mutations_from_vcf(vcf)
  types <- edit_type(vcf)
  #x <- which(muts != types)
  y <- edit_context.mat
  y <- chartr("T", "U", y)
  edit_context.mat <- y
  res <- list(types, edit_context.mat)
  names(res) <- c("types", "context")
  return(res)
}
