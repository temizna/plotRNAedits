% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/edit_context.R
\name{edit_context}
\alias{edit_context}
\title{Retrieve context of RNA base editions 
This function is a variation of mut_context from MutationalPatterns
A function to extract the bases 3' upstream and 5' downstream of the base
substitutions from the reference genome. Currently only 1 base extention
is supported.}
\usage{
edit_context(vcf, ref_genome, extension = 1)
}
\arguments{
\item{vcf}{A Granges object}

\item{ref_genome}{Reference genome}

\item{extension}{The number of bases, that's extracted upstream and
downstream of the base substitutions. (Default: 1).}
}
\value{
Character vector with the context of the base substitutions
}
\description{
Retrieve context of RNA base editions 
This function is a variation of mut_context from MutationalPatterns
A function to extract the bases 3' upstream and 5' downstream of the base
substitutions from the reference genome. Currently only 1 base extention
is supported.
}
\examples{
## See the 'read_vcfs_as_granges()' example for how we obtained the
## following data:
vcfs <- readRDS(system.file("states/read_vcfs_as_granges_output.rds",
  package = "MutationalPatterns"
))

## Load the corresponding reference genome.
ref_genome <- "BSgenome.Mmusculus.UCSC.mm10"
library(ref_genome, character.only = TRUE)

## Get the standard context
edit_context <- edit_context(vcfs[[1]], ref_genome)

}
