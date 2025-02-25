#' Generate Orthogonalized LACs
#'
#' This function converts original LACs into orthogonalized LACs that are uncorrelated to global PCs.
#' @param input.bp filename of original LACs files.
#' @param input.gc dataframe of global PCs of test samples. Please add sample IDs as rownames.
#' @param output.la filename of the output orthogonalized LACs files.
#' @param id.type the class of sample IDs. Possible values are "character", "integer", and "numeric".
#' @param select.PC a vector indicating PCs included in further analysis, e.g. c(1,2,3). Default is equivalent to the number of columns in the psi matrix.
#' @param select.pos a vector indicating genetic positions included in further analysis. Ascending Order. Default is set to include all genetic positions in the input.bp.
#' @param select.id a vector indicating sample IDs included in further analysis. Default is set to include all samples in the input.bp.
#' @param if.haplo a logical switch indicating whether to generate haplotype-based LACs (default = FALSE).
#' @return a compressed binary file that documents sample information, genetic positions at breakpoints and orthogonalized LACs based on original LACs.
#' @export
#' @concept Generating and retrieving LACs

get.orthogonal.LACs<-function(input.la,input.gc,output.la,id.type=c('character','integer','numeric'),select.PC=NULL,select.pos=NULL,select.id=NULL,if.haplo=F){
  con1<-file(input.la,'rb')
  n_PC<-readBin(con1,what=integer(),n=1,size=4)
  n_indi<-readBin(con1,what=integer(),n=1,size=4)
  if(id.type=='character')
    indi_id<-readBin(con=con1,what=character(),n=n_indi,size=4)
  else if(id.type=='integer')
    indi_id<-readBin(con=con1,what=integer(),n=n_indi,size=4)
  else if(id.type=='numeric')
    indi_id<-readBin(con=con1,what=numeric(),n=n_indi,size=4)
  else{
    close(con1)
    stop('Error: Sample IDs are supposed to be one of the type of characters, integers or numerics. Please double check the sample IDs\' type.')
  }
  n_pos<-readBin(con=con1,what=integer(),n=1,size=4)
  pos<-readBin(con=con1,what=integer(),n=n_pos,size=4)
  if(is.null(select.pos))
    pos_index<-seq_along(pos)
  else{
    pos_index<-match(select.pos,pos)
    if(all(is.na(pos_index))){
      close(con1)
      stop('Error: None positions were included in the relevant local ancestry documents. Please double check the input positions.')
    }else if(any(is.na(pos_index))){
      warning(paste('Positions',paste(select.pos[is.na(pos_index)],collapse=', '),'didn\'t exist. Please double check the input positions.'))
      pos_index<-pos_index[!is.na(pos_index)]
    }
  }
  if(is.null(select.id))
    indi_common<-intersect(indi_id,rownames(input.gc))
  else
    indi_common<-Reduce(intersect,list(indi_id,rownames(input.gc),select.id))
  if(is.null(indi_common)){
    close(con1)
    stop('Error: There is no common individual between selected samples, global PCs, and LAC files. Please double check your input data.')
  }else if(!is.null(select.id) && !setequal(indi_common,select.id))
    warning(paste(paste('Test samples',setdiff(select.id,indi_common),collapse=', '),'didn\'t exist. Please double check the input sample IDs.'))
  indi_gc_index<-match(indi_common,rownames(input.gc))
  indi_la_index<-match(indi_common,indi_id)
  n_unit_read_la<-n_PC*n_indi
  n_unit_skip_la<-ifelse(if.haplo,n_unit_read_la*2,n_unit_read_la)
  if(is.null(select.PC)) select.PC<-seq_len(n_PC)
  input.gc<-as.matrix.data.frame(input.gc[indi_gc_index,select.PC])
  A<-cbind(1,input.gc)
  ATA<-crossprod(A,A)
  ATA_i<-chol2inv(chol(ATA))
  I<-diag(nrow(input.gc))
  I<-Matrix(I,sparse=T)
  P<-I - tcrossprod(A %*% ATA_i,A)
  rm(A,ATA,ATA_i,I)
  con2<-file(output.la,'wb')
  writeBin(length(select.PC),con=con2,size=4)
  writeBin(length(indi_common),con=con2,size=4)
  writeBin(indi_common,con=con2,size=4)
  writeBin(length(pos_index),con=con2,size=4)
  writeBin(pos[pos_index],con=con2,size=4)
  prev_j<-0
  for(j in pos_index){
    j_extra<-j-prev_j-1
    if(j_extra>0)
      seek(con=con1,where=j_extra*n_unit_skip_la*4,origin='current')
    la_coord<-readBin(con1,what=numeric(),n=n_unit_read_la,size=4)
    la_coord<-matrix(la_coord,ncol=n_PC,byrow=T)[indi_la_index,select.PC]
    if(all(is.na(la_coord))){
      close(con1)
      close(con2)
      stop(paste('None local ancestry calls at position',pos[pos_index][j],'are mathemathically available. LACs generation was stopped.'))
    }
    la_proj<-as.vector(P %*% la_coord)
    writeBin(la_proj,con=con2,size=4)
    prev_j<-j
  }
  close(con1)
  close(con2)
  return(invisible(NULL))
}
