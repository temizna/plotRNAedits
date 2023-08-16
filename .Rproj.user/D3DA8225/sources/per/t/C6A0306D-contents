#' Make RNA edit count matrix of 192 trinucleotides
#' This function is a variation of mut_matrix from MutationalPatterns
#' @description Make 192 trinucleotide RNA edit count matrix
#' @param vcf_list GRangesList or GRanges object.
#' @param ref_genome BSgenome reference genome object
#' @param extension The number of bases, that's extracted upstream and
#' downstream of the base substitutions. (Default: 1).
#' @return 192 RNA edit count matrix
#'
#' @examples
#' ## See the 'read_vcfs_as_granges()' example for how we obtained the
#' ## following data:
#' grl <- readRDS(system.file("states/read_vcfs_as_granges_output.rds",
#'   package = "MutationalPatterns"
#' ))
#'
#' ## Load the corresponding reference genome.
#' ref_genome <- "BSgenome.Mmusculus.UCSC.mm10"
#' library(ref_genome, character.only = TRUE)
#'
#' ## Construct a mutation matrix from the loaded VCFs in comparison to the
#' ## ref_genome.
#' edit_mat <- edit_matrix(vcf_list = grl, ref_genome = ref_genome)
#'
#'
#' @export
edit_matrix<-function (vcf_list, ref_genome, extension = 1)
{
  if (inherits(vcf_list, "list")) {
    vcf_list <- GenomicRanges::GRangesList(vcf_list)
  }
  if (inherits(vcf_list, "CompressedGRangesList")) {
    gr_sizes <- S4Vectors::elementNROWS(vcf_list)
    gr <- BiocGenerics::unlist(vcf_list)
  }
  else if (inherits(vcf_list, "GRanges")) {
    gr <- vcf_list
    gr_sizes <- length(gr)
    names(gr_sizes) <- "My_sample"
  }
  else {
    .not_gr_or_grl(vcf_list)
  }
  edit_context.mat <- edit_context(gr, ref_genome, extension)
  mut_mat <- edit_192_occurrences(edit_context.mat, gr_sizes)
  return(mut_mat)
}
