% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get.la.aware.score.R
\name{get.la.aware.score}
\alias{get.la.aware.score}
\title{Generate local-ancestry-aware GWAS results}
\usage{
get.la.aware.score(
  null.obj,
  geno.file,
  outfile,
  LAC.file,
  n_PC_used = NULL,
  id.type = c("character", "integer", "numeric"),
  MAF.range = c(1e-07, 0.5),
  miss.cutoff = 1,
  missing.method = "impute2mean",
  is.dosage = F,
  n.batch = 100
)
}
\arguments{
\item{null.obj}{a class glmmkin or class glmmkin.multi object, returned by fitting the null GLMMusing glmmkin.}

\item{geno.file}{the full name of a GDS file (including the suffix .gds), or an object of class SeqVarGDSClass.}

\item{outfile}{the output file name.}

\item{LAC.file}{the full name of a binary LAC file.}

\item{n_PC_used}{the dimension of PCs that will be included in analysis. Default is set to the number of dimension of PCs in the LAC.file.}

\item{id.type}{the class of sample IDs. Possible values are "character", "integer", and "numeric".}

\item{MAF.range}{a numeric vector of length 2 defining the minimum and maximum minor allele frequencies of variants that should be included in the analysis (default = c(1e-7, 0.5)).}

\item{miss.cutoff}{the maximum missing rate allowed for a variant to be included (default = 1, including all variants).}

\item{missing.method}{method of handling missing genotypes. Either "impute2mean" or "omit" (default = "impute2mean").}

\item{is.dosage}{a logical switch for whether imputed dosage should be used from a GDS infile (default = FALSE).}

\item{n.batch}{an integer for how many SNPs should be tested in a batch (default = 100). The computational time can increase dramatically if this value is either small or large. The optimal value for best performance depends on the user’s system.}
}
\value{
a dataframe with following components:
\tabular{ll}{
\strong{Name} \tab \strong{Description} \cr
SNP \tab SNP name, as supplied in snps. \cr
CHR \tab Chromosome, copied from .gds file. \cr
POS \tab physical position in base pair, copied from .gds file. \cr
REF \tab reference allele, copied from .gds file. \cr
ALT \tab alternate allele, copied from .gds file. \cr
MISSRATE \tab number of individuals with non-missing genotypes for each SNP. \cr
AF \tab ALT allele frequency for each SNP. \cr
N \tab total sample size. \cr
SCORE \tab the summary score of the effect allele. \cr
VAR \tab the variance of the summary score. \cr
PVAL \tab local-ancestry-aware GWAS p-values. \cr
}
}
\description{
This function uses a glmmkin class object from the null GLMM to perform score tests for local-ancestry-aware association with genotypes in a GDS file .gds file.
}
\concept{Calculating local-ancestry-aware GWAS}
