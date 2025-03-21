% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get.breakpoints.R
\name{get.breakpoints}
\alias{get.breakpoints}
\title{Identify local ancestry breakpoints and generate local-ancestry-based global ancestry proportions if applicable.}
\usage{
get.breakpoints(
  input.la,
  input.pos,
  input.id,
  input.anc,
  input.sep,
  output.bp,
  start.pos = 0,
  if.yield.ga = F,
  output.ga = NULL,
  if.haplo.ga = F
)
}
\arguments{
\item{input.la}{filename for pure local ancestry records. Numeric. Do not include any character description.}

\item{input.pos}{a numeric vector of SNP positions in base pair. Ascending Order.}

\item{input.id}{a vector of sample IDs.}

\item{input.anc}{a character vector of reference ancestry labels.}

\item{input.sep}{delimiters in the input.la.}

\item{output.bp}{filename of identified local ancestry breakpoints.}

\item{start.pos}{An integer indicating the starting position for counting on local ancestry (default = 0).}

\item{if.yield.ga}{a logical switch indicating whether to generate global ancestry proportions or not (default = TRUE).}

\item{output.ga}{filename of output global ancestry proportions. Only effective when if.yield.ga = TRUE (default = NULL).}

\item{if.haplo.ga}{a logical switch indicating whether to generate haplotype-based global ancestry proportions (default = FALSE).}
}
\value{
a compressed binary file that documents sample information, genetic positions at breakpoints and the corresponding local ancestry references, with optional global ancestry proportions.
}
\description{
This function reads in local ancestry files in chunks, identifies local ancestry tracts, and yields local-ancestry-based global ancestry proportions if applicable.
}
\concept{Generating and retrieving LACs}
