% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get.LACs.R
\name{get.LACs}
\alias{get.LACs}
\title{Generate Local Ancestry Coordinates (LACs)}
\usage{
get.LACs(
  input.bp,
  psi,
  output.la,
  id.type = c("character", "integer", "numeric"),
  select.PC = NULL,
  select.pos = NULL,
  select.id = NULL,
  if.haplo = F
)
}
\arguments{
\item{input.bp}{filename of extracted local ancestry inference at breakpoints.}

\item{psi}{the designated population anchor matrix.}

\item{output.la}{filename of the output LACs files.}

\item{id.type}{the class of sample IDs. Possible values are "character", "integer", and "numeric".}

\item{select.PC}{a vector indicating PCs included in further analysis, e.g. c(1,2,3). Default is equivalent to the number of columns in the psi matrix.}

\item{select.pos}{a vector indicating genetic positions included in further analysis. Ascending Order. Default is set to include all genetic positions in the input.bp.}

\item{select.id}{a vector indicating sample IDs included in further analysis. Default is set to include all samples in the input.bp.}

\item{if.haplo}{a logical switch indicating whether to generate haplotype-based LACs (default = FALSE).}
}
\value{
a compressed binary file that documents sample information, genetic positions at breakpoints and refined LACs given corresponding local ancestry inferences.
}
\description{
This function uses the extracted local ancestry breakpoints and corresponding local ancestry inferences and converts them into binary LACs files.
}
\concept{Generating and retrieving LACs}
