#' Identify local ancestry breakpoints and generate local-ancestry-based global ancestry proportions if applicable.
#'
#' This function reads in local ancestry files in chunks, identifies local ancestry tracts, and yields local-ancestry-based global ancestry proportions if applicable.
#'
#' @param input.la filename for pure local ancestry records. Numeric. Do not include any character description.
#' @param input.pos a numeric vector of SNP positions in base pair. Ascending Order.
#' @param input.id a vector of sample IDs.
#' @param input.anc a character vector of reference ancestry labels.
#' @param input.sep delimiters in the input.la.
#' @param output.bp filename of identified local ancestry breakpoints.
#' @param start.pos An integer indicating the starting position for counting on local ancestry (default = 0).
#' @param if.yield.ga a logical switch indicating whether to generate global ancestry proportions or not (default = TRUE).
#' @param output.ga filename of output global ancestry proportions. Only effective when if.yield.ga = TRUE (default = NULL).
#' @param if.haplo.ga a logical switch indicating whether to generate haplotype-based global ancestry proportions (default = FALSE).
#'
#' @return a compressed binary file that documents sample information, genetic positions at breakpoints and the corresponding local ancestry references, with optional global ancestry proportions.
#' @export
#' @concept Generating and retrieving LACs

get.breakpoints<-function(input.la,input.pos,input.id,input.anc,input.sep,output.bp,start.pos=0,if.yield.ga=F,output.ga=NULL,if.haplo.ga=F){
  n_ancestry<-length(input.anc)
  file1<-paste(output.bp,'info',sep='.')
  con1<-file(file1,'wb')
  writeBin(n_ancestry,con=con1,size=4)
  writeBin(input.anc,con=con1,size=4)
  writeBin(length(input.id),con=con1,size=4)
  writeBin(input.id,con=con1,size=4)
  # read in large local ancestry calls in chunks
  ### need to specify the delimiter in local ancestry calls
  read_file<-function(cmd,input.sep){
    tryCatch({df<-data.table::fread(cmd=cmd,header=F,sep=input.sep);return(list(df=df,index=T))},error=function(e) return(list(index=F)),warning=function(w) return(list(index=F)))
  }
  total_pos<-max(input.pos)-start.pos
  prev_i<-i<-0
  bp_spot<-c()
  prev_bp<-start.pos
  prev_line<-NULL
  la_sum<-0
  file2<-paste(output.bp,'la',sep='.')
  con2<-file(file2,'wb')
  while(T){
    cmd<-sprintf('sed -n "%i,%ip" %s',prev_i+1,prev_i+1000,input.la)
    la_file<-read_file(cmd,input.sep)
    if(la_file$index==FALSE)
      break
    la_file<-la_file$df
    i<-i+nrow(la_file)
    for(j in seq_len(nrow(la_file))){
      la_slot<-la_file[j,]
      if(!identical(prev_line,la_slot,ignore.environment=T)){
        prev_line<-la_slot
        bp_slot<-j+prev_i
        bp_spot<-c(bp_spot,bp_slot)
        la_sum<-la_sum+(input.pos[bp_slot]-prev_bp)/total_pos*la_slot
        prev_bp<-input.pos[bp_slot]
        writeBin(unlist(la_slot),con=con2,size=4)
      }
    }
    prev_i<-i
  }
  close(con2)
  writeBin(length(bp_spot),con=con1,size=4)
  writeBin(as.integer(input.pos[bp_spot]),con=con1,size=4)
  close(con1)
  system(sprintf('cat %s %s > %s',file1,file2,output.bp))
  system(sprintf('rm %s %s',file1,file2))
  if(if.yield.ga){
    la_sum<-unlist(la_sum)
    ga_file<-matrix(la_sum,ncol=n_ancestry,byrow=T)
    ga_file<-as.data.frame(prop.table(ga_file,margin=1))
    colnames(ga_file)<-input.anc
    if(if.haplo.ga){
      rownames(ga_file)<-c(t(outer(input.id,c('hap1','hap2'),paste,sep=':::')))
      ga_file<-tibble::rownames_to_column(ga_file,'haplotypes')
    }else{
      ga_file<-ga_file %>% group_by(indi=gl(n()/2,2)) %>% summarise(across(everything(),\(x) mean(x,na.rm=T))) %>% select(-indi)
      rownames(ga_file)<-input.id
      ga_file<-tibble::rownames_to_column(ga_file,'sampleIDs')
    }
    data.table::fwrite(ga_file,output.ga,sep='\t',row.names=F,col.names=T,quote=F)
  }
  return(invisible(NULL))
}
