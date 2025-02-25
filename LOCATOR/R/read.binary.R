#' Retrieve relevant local ancestry information at breakpoints
#'
#' This function reads and extracts compressed local ancestry information, either in local ancestry inferences or in LACs, given designated genetic positions.
#' @param input.file filename of relevant compressed binary local ancestry files.
#' @param type file type of input.file. Either be "breakpoints" or "LACs".
#' @param id.type the class of sample IDs. Possible values are "character", "integer", and "numeric".
#' @param select.pos a vector indicating genetic positions included in further analysis. Ascending Order. Default is set to include all genetic positions in the input.file.
#' @param select.id a vector indicating sample IDs included in further analysis. Default is set to include all samples in the input.file.
#' @param if.haplo a logical switch indicating whether it is haplotype-based local ancestry information (default = FALSE).
#' @param col.names a vector indicating column names for local ancestry in local ancestry inference and PCs in LACs. For LACs files, default is set from PC1 to PCn for all columns in the input.file.
#' @return a list of dataframes including genetic-position-specfic local ancestry information in either local ancestry references or in LACs.
#' @export
#' @concept Generating and retrieving LACs

read.binary<-function(input.file,type=c('breakpoints','LACs'),id.type=c('character','integer','numeric'),select.pos,select.id=NULL,if.haplo=F,col.names=NULL){
  con<-file(input.file,'rb')
  if(type=='breakpoints'){
    n_ancestry<-readBin(con=con,what=integer(),n=1,size=4)
    col.names<-readBin(con=con,what=character(),n=n_ancestry,size=4)
  }else if(type=='LACs'){
    n_PC<-readBin(con=con,what=integer(),n=1,size=4)
    if(is.null(col.names))
      col.names<-paste0('PC',seq_len(n_PC))
  }else{
    close(con)
    stop('Input files are supposed to be either compressed breakpoints or LACs files. Please double check the input file link.')
  }
  n_indi<-readBin(con=con,what=integer(),n=1,size=4)
  if(id.type=='character')
    indi_id<-readBin(con=con,what=character(),n=n_indi,size=4)
  else if(id.type=='integer')
    indi_id<-readBin(con=con,what=integer(),n=n_indi,size=4)
  else if(id.type=='numeric')
    indi_id<-readBin(con=con,what=numeric(),n=n_indi,size=4)
  else{
    close(con)
    stop('Sample IDs are supposed to be one of the type of characters, integers or numerics. Please double check the sample IDs\' type.')
  }
  n_pos<-readBin(con=con,what=integer(),n=1,size=4)
  pos<-readBin(con=con,what=integer(),n=n_pos,size=4)
  pos_index<-match(select.pos,pos)
  if(all(is.na(pos_index))){
    close(con)
    stop('None positions were included in the relevant local ancestry documents. Please double check the input positions.')
  }else if(any(is.na(pos_index))){
    warning(paste('Positions',paste(select.pos[is.na(pos_index)],collapse=', '),'didn\'t exist. Please double check the input positions.'))
    pos_index<-pos_index[!is.na(pos_index)]
  }
  if(is.null(select.id))
    indi_index<-seq_along(indi_id)
  else{
    indi_index<-match(select.id,indi_id)
    if(all(is.na(indi_index))){
      close(con)
      stop('None test samples were included in the relevant local ancestry documents. Please double check the input sample IDs.')
    }else if(any(is.na(indi_index))){
      warning(paste('Test samples',paste(select.id[is.na(indi_index)],collapse=', '),'didn\'t exist. Please double check the input sample IDs.'))
      indi_index<-indi_index[!is.na(indi_index)]
    }
  }
  la_coord_list<-list()
  prev_j<-0
  if(type=='breakpoints'){
    n_unit_read_la<-n_ancestry*n_indi*2
    indi_index_2<-c(sapply(indi_index,function(x) return(c(2*x-1,2*x))))
    for(j in pos_index){
      j_extra<-j-prev_j-1
      if(j_extra>0)
        seek(con=con,where=j_extra*n_unit_read_la*4,origin='current')
      la_coord<-readBin(con,what=numeric(),n=n_unit_read_la,size=4)
      la_coord<-as.data.frame.matrix(matrix(la_coord,ncol=n_ancestry,byrow=T)[indi_index_2,,drop=F])
      colnames(la_coord)<-col.names
      if(if.haplo){
        rownames(la_coord)<-c(t(outer(indi_id[indi_index],c('hap1','hap2'),paste,sep=':::')))
        la_coord<-tibble::rownames_to_column(la_coord,'haplotypes')
      }else{
        la_coord<-la_coord %>% group_by(indi=gl(n()/2,2)) %>% summarise(across(everything(),\(x) mean(x,na.rm=T))) %>% select(-indi)
        rownames(la_coord)<-indi_id[indi_index]
        la_coord<-tibble::rownames_to_column(la_coord,'sampleIDs')
      }
      la_coord_list<-append(la_coord_list,list(la_coord))
      prev_j<-j
    }
  }else{
    n_unit_read_la<-n_PC*n_indi
    for(j in pos_index){
      j_extra<-j-prev_j-1
      if(j_extra>0)
        seek(con=con,where=j_extra*n_unit_read_la*4,origin='current')
      la_coord<-readBin(con,what=numeric(),n=n_unit_read_la,size=4)
      la_coord<-as.data.frame.matrix(matrix(la_coord,ncol=n_PC)[indi_index,,drop=F])
      colnames(la_coord)<-col.names
      rownames(la_coord)<-indi_id[indi_index]
      la_coord<-tibble::rownames_to_column(la_coord,'sampleIDs')
      la_coord_list<-append(la_coord_list,list(la_coord))
      prev_j<-j
    }
  }
  close(con)
  return(la_coord_list)
}
