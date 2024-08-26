#' Generate Local Ancestry Coordinates (LACs)
#'
#' This function uses the extracted local ancestry breakpoints and corresponding local ancestry inferences and converts them into binary LACs files.
#' @param input.bp filename of extracted local ancestry inference at breakpoints.
#' @param psi the designated population anchor matrix.
#' @param output.la filename of the output LACs files.
#' @param id.type the class of sample IDs. Possible values are "character", "integer", and "numeric".
#' @param select.PC a vector indicating PCs included in further analysis, e.g. c(1,2,3). Default is equivalent to the number of columns in the psi matrix.
#' @param select.pos a vector indicating genetic positions included in further analysis. Ascending Order. Default is set to include all genetic positions in the input.bp.
#' @param select.id a vector indicating sample IDs included in further analysis. Default is set to include all samples in the input.bp.
#' @param if.haplo a logical switch indicating whether to generate haplotype-based LACs (default = FALSE).
#' @return a compressed binary file that documents sample information, genetic positions at breakpoints and refined LACs given corresponding local ancestry inferences.
#' @export
#' @concept Generating and retrieving LACs

get.LACs<-function(input.bp,psi,output.la,id.type=c('character','integer','numeric'),select.PC=NULL,select.pos=NULL,select.id=NULL,if.haplo=F){
  con1<-file(input.bp,'rb')
  n_ancestry<-readBin(con=con1,what=integer(),n=1,size=4)
  ancestry_list<-readBin(con=con1,what=character(),n=n_ancestry,size=4)
  n_indi<-readBin(con=con1,what=integer(),n=1,size=4)
  if(id.type=='character')
    indi_id<-readBin(con=con1,what=character(),n=n_indi,size=4)
  else if(id.type=='integer')
    indi_id<-readBin(con=con1,what=integer(),n=n_indi,size=4)
  else if(id.type=='numeric')
    indi_id<-readBin(con=con1,what=numeric(),n=n_indi,size=4)
  else{
    close(con1)
    stop('Sample IDs are supposed to be one of the type of characters, integers or numerics. Please double check the sample IDs\' type.')
  }
  n_pos<-readBin(con=con1,what=integer(),n=1,size=4)
  pos<-readBin(con=con1,what=integer(),n=n_pos,size=4)
  n_unit_read_la<-n_ancestry*n_indi*2
  if(if.haplo)
    haplo_index<-rep(c(F,T),n_indi)
  if(is.null(select.PC)) select.PC<-seq_len(ncol(psi))
  n_PC<-length(select.PC)
  psi<-psi[,select.PC]
  if(is.null(select.id))
    indi_index<-seq_along(indi_id)
  else{
    indi_index<-match(select.id,indi_id)
    if(all(is.na(indi_index))){
      close(con1)
      stop('None test samples were included in the relevant local ancestry documents. Please double check the input sample IDs.')
    }else if(any(is.na(indi_index))){
      warning(paste('Test samples',paste(select.id[is.na(indi_index)],collapse=', '),'didn\'t exist. Please double check the input sample IDs.'))
      indi_index<-indi_index[!is.na(indi_index)]
    }
  }
  if(is.null(select.pos))
    pos_index<-seq_along(pos)
  else{
    pos_index<-match(select.pos,pos)
    if(all(is.na(pos_index))){
      close(con1)
      stop('None positions were included in the relevant local ancestry documents. Please double check the input positions.')
    }else if(any(is.na(pos_index))){
      warning(paste('Positions',paste(select.pos[is.na(pos_index)],collapse=', '),'didn\'t exist. Please double check the input positions.'))
      pos_index<-pos_index[!is.na(pos_index)]
    }
  }
  con2<-file(output.la,'wb')
  writeBin(n_PC,con=con2,size=4)
  writeBin(length(indi_index),con=con2,size=4)
  writeBin(indi_id[indi_index],con=con2,size=4)
  writeBin(length(pos_index),con=con2,size=4)
  writeBin(pos[pos_index],con=con2,size=4)
  prev_j<-0
  for(j in pos_index){
    j_extra<-j-prev_j-1
    if(j_extra>0)
      seek(con=con1,where=j_extra*n_unit_read_la*4,origin='current')
    la_coord<-readBin(con1,what=numeric(),n=n_unit_read_la,size=4)
    la_coord<-matrix(la_coord,ncol=n_PC,byrow=T)
    if(all(is.na(la_coord))){
      close(con1)
      close(con2)
      stop(paste('None local ancestry calls at position',pos[pos_index][j],'are mathemathically available. LACs generation was stopped.'))
    }else if(any(is.na(la_coord))){
      missing_id<-which(is.na(la_coord),arr.ind=T)[,1]
      missing_hap<-ifelse(missing_id %% 2,'hap1','hap2')
      missing_id<-indi_id[indi_index][missing_id %% 2 + 1]
      warning(paste(paste(missing_id,missing_hap,collapse=' and '),'at position',pos[pos_index][j],'are mathematically available. NAs are yielded in the output. Please double check the input local ancestry calls.'))
    }
    if(if.haplo){
      la_proj<-as.vector(la_coord[!haplo_index,][indi_index,] %*% psi)
      writeBin(la_proj,con=con2,size=4)
      la_proj<-as.vector(la_coord[haplo_index,][indi_index,] %*% psi)
      writeBin(la_proj,con=con2,size=4)
    }else{
      la_coord<-as.data.frame.matrix(la_coord) %>% group_by(indi=gl(n()/2,2)) %>% summarise(across(everything(),\(x) mean(x,na.rm=T))) %>% select(-indi)
      la_proj<-as.vector(as.matrix.data.frame(la_coord[indi_index,]) %*% psi)
      writeBin(la_proj,con=con2,size=4)
    }
    prev_j<-j
  }
  close(con1)
  close(con2)
  return(invisible(NULL))
}
