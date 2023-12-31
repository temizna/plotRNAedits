% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/edit_matrix.R
\name{edit_matrix}
\alias{edit_matrix}
\title{Make RNA edit count matrix of 192 trinucleotides
This function is a variation of mut_matrix from MutationalPatterns}
\usage{
edit_matrix(vcf_list, ref_genome, extension = 1)
}
\arguments{
\item{vcf_list}{GRangesList or GRanges object.}

\item{ref_genome}{BSgenome reference genome object}

\item{extension}{The number of bases, that's extracted upstream and
downstream of the base substitutions. (Default: 1).}
}
\value{
192 RNA edit count matrix
}
\description{
Make 192 trinucleotide RNA edit count matrix
}
\examples{
## See the 'read_vcfs_as_granges()' example for how we obtained the
## following data:
grl <- readRDS(system.file("states/read_vcfs_as_granges_output.rds",
  package = "MutationalPatterns"
))

## Load the corresponding reference genome.
ref_genome <- "BSgenome.Mmusculus.UCSC.mm10"
library(ref_genome, character.only = TRUE)

## Construct a mutation matrix from the loaded VCFs in comparison to the
## ref_genome.
edit_mat <- edit_matrix(vcf_list = grl, ref_genome = ref_genome)


}
