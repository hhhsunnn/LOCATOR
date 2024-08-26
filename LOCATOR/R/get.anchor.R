#' Get the population anchor matrix
#'
#' This function generates the population anchor matrix that illustrates the ancestry distribution among designated PCs.
#'
#' @param input.ga A numeric matrix of global ancestry. Please do not include any sample IDs.
#' @param input.gc A numeric matrix of global PCs. Please do not include any sample IDs and ensure it corresponds to input.ga.
#' @return a square matrix of population anchor matrix sharing the same number of columns as input.gc.
#' @export
#' @concept Generating and retrieving LACs

get.anchor<-function(input.ga,input.gc){
  ga_file<-if(!is.matrix(input.ga)) as.matrix(input.ga) else input.ga
  gc_file<-if(!is.matrix(input.gc)) as.matrix(input.gc) else input.gc
  psi<-chol2inv(chol(crossprod(ga_file)))
  psi<-tcrossprod(psi,ga_file) %*% gc_file
  return(psi)
}
