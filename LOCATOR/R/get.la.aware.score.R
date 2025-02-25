#' Generate local-ancestry-aware GWAS results
#'
#' This function uses a glmmkin class object from the null GLMM to perform score tests for local-ancestry-aware association with genotypes in a GDS file .gds file.
#' @param null.obj a class glmmkin or class glmmkin.multi object, returned by fitting the null GLMMusing glmmkin.
#' @param geno.file the full name of a GDS file (including the suffix .gds), or an object of class SeqVarGDSClass.
#' @param outfile the output file name.
#' @param LAC.file the full name of a binary LAC file.
#' @param n_PC_used the dimension of PCs that will be included in analysis. Default is set to the number of dimension of PCs in the LAC.file.
#' @param id.type the class of sample IDs. Possible values are "character", "integer", and "numeric".
#' @param MAF.range a numeric vector of length 2 defining the minimum and maximum minor allele frequencies of variants that should be included in the analysis (default = c(1e-7, 0.5)).
#' @param miss.cutoff the maximum missing rate allowed for a variant to be included (default = 1, including all variants).
#' @param missing.method method of handling missing genotypes. Either "impute2mean" or "omit" (default = "impute2mean").
#' @param is.dosage a logical switch for whether imputed dosage should be used from a GDS infile (default = FALSE).
#' @param n.batch an integer for how many SNPs should be tested in a batch (default = 100). The computational time can increase dramatically if this value is either small or large. The optimal value for best performance depends on the user’s system.
#' @return a dataframe with following components:
#' \tabular{ll}{
#' \strong{Name} \tab \strong{Description} \cr
#' SNP \tab SNP name, as supplied in snps. \cr
#' CHR \tab Chromosome, copied from .gds file. \cr
#' POS \tab physical position in base pair, copied from .gds file. \cr
#' REF \tab reference allele, copied from .gds file. \cr
#' ALT \tab alternate allele, copied from .gds file. \cr
#' MISSRATE \tab number of individuals with non-missing genotypes for each SNP. \cr
#' AF \tab ALT allele frequency for each SNP. \cr
#' N \tab total sample size. \cr
#' SCORE \tab the summary score of the effect allele. \cr
#' VAR \tab the variance of the summary score. \cr
#' PVAL \tab local-ancestry-aware GWAS p-values. \cr
#' }
#' @export
#' @concept Calculating local-ancestry-aware GWAS

# arguments
### null.obj: the null glmmkin model
### geno.file: the full name of GDS file with its suffix '.gds'
### outfile: output filename
### LAC.file: LACs binary file
get.la.aware.score<-function(null.obj,geno.file,outfile,LAC.file,n_PC_used=NULL,id.type=c('character','integer','numeric'),MAF.range=c(1e-7,0.5),miss.cutoff=1,missing.method='impute2mean',is.dosage=F,n.batch=100){
	# Import ancilary functions into the working space
	### Define the function for generalized Schur complement
	cal_schur<-function(mat) tryCatch({rje::schur(mat,x=1,y=1)},error = function(e) {if(grepl('singular',e$message,ignore.case=T)) return(mat[1,1] - mat[1,-1,drop=F] %*% MASS::ginv(mat[-1,-1]) %*% mat[-1,1,drop=F])})
	### Define the function to calculate p-values (dense matrix ver.)
	get_lm_p_value_dense<-function(j_g,Gresid2=Gresid,GPL2=GPL,Lresid2=Lresid,LPL2=LPL,GPG2=GPG){
		cal_p<-function(i){
			M_S_beta<-rbind(c(Gresid2[i],GPL2[i,]),cbind(Lresid2,LPL2))
			S_beta<-cal_schur(M_S_beta)
			M_var_beta<-rbind(cbind(GPG2[i],GPL2[i,,drop=F]),cbind(GPL2[i,],LPL2))
			var_beta<-cal_schur(M_var_beta)
			T_beta<-crossprod(S_beta)/var_beta
			p_val<-pchisq(T_beta,df=1,lower.tail=F)
			return(c(S_beta,var_beta,p_val))
		}
		p_values<-t(sapply(j_g,cal_p))
		return(p_values)
	}
	### Define the function to generate M_var_beta matrix
	create_M_var<-function(Xsigma_iL,Lsigma_iL,n_X,n_PC_used){
		M_var_beta<-matrix(0,nrow=1+n_X+n_PC_used,ncol=1+n_X+n_PC_used)
		XL_block<-rbind(cbind(Xsigma_iX,Xsigma_iL),cbind(t(Xsigma_iL),Lsigma_iL))
		M_var_beta[-1,-1][upper.tri(M_var_beta[-1,-1],diag=T)]<-XL_block[upper.tri(XL_block,diag=T)]
		return(M_var_beta)
	}
	### Define the function to calculate p-values (sparse matrix ver.)
	get_lm_p_value_sparse<-function(j_g,Gresid2=Gresid,GPL2=GPL,Lresid2=Lresid,LPL2=LPL,Gsigma_iG2=Gsigma_iG,Gsigma_iX2=Gsigma_iX,Gsigma_iL2=Gsigma_iL,M_var_beta2=M_var_beta){
		cal_p<-function(i){
			M_S_beta<-rbind(c(Gresid2[i],GPL2[i,]),cbind(Lresid2,LPL2))
			S_beta<-cal_schur(M_S_beta)
			M_var_beta2[1,]<-c(Gsigma_iG2[i],Gsigma_iX2[i,],Gsigma_iL2[i,])
			M_var_beta2<-M_var_beta2+t(M_var_beta2)-diag(diag(M_var_beta2))
			var_beta<-cal_schur(M_var_beta2)
			T_beta<-crossprod(S_beta)/var_beta
			p_val<-pchisq(T_beta,df=1,lower.tail=F)
			return(c(S_beta,var_beta,p_val))
		}
		p_values<-t(sapply(j_g,cal_p))
		return(p_values)
	}
#	if(!require('Matrix'))
#		install.packages('Matrix')
#	if(!require('rje'))
#		install.packages('rje')
	missing.method<-try(match.arg(missing.method,c('impute2mean','impute2zero')))
	if(inherits(missing.method,'try-error')) stop('Error: \"missing.method\" must be \"impute2mean\" or \"impute2zero\".')
	# read in pop.ids from the GDS file
	if(grepl('\\.gds$', geno.file[1])){
		if(!inherits(geno.file,'SeqVarGDSClass'))
			gds<-SeqArray::seqOpen(geno.file)
	}else
		gds<-geno.file
	sample.id<-SeqArray::seqGetData(gds,'sample.id')
	# read in pop.ids from the LAC file
	con_la<-file(LAC.file,'rb')
	n_PC<-readBin(con_la,what=integer(),n=1,size=4)
	if(is.null(n_PC_used))
	  n_PC_used<-n_PC
	if(n_PC<n_PC_used)
		stop('Error: The dimension of LACs is less than requested.')
	n_indi_la<-readBin(con_la,what=integer(),n=1,size=4)
	if(id.type=='character')
	  indi_la<-readBin(con_la,what=character(),n_indi_la,size=4)
	else if(id.type=='integer')
	  indi_la<-readBin(con_la,what=integer(),n=n_indi_la,size=4)
	else if(id.type=='numeric')
	  indi_la<-readBin(con_la,what=numeric(),n=n_indi_la,size=4)
	else{
	  close(con_la)
	  stop('Error: Sample IDs in LACs are supposed to be one of the type of characters, integers or numerics. Please double check the sample IDs\' type.')
	}
	n_pos_la<-readBin(con_la,what=integer(),n=1,size=4)
	pos_la<-readBin(con_la,what=integer(),n=n_pos_la,size=4)
	n_unit_read_la<-n_PC*n_indi_la
	indi_common<-Reduce(intersect,list(indi_la,sample.id,null.obj$id_include))
#	indi_common<-sort(indi_common)  # Need to keep the common sample ids in the ascending order like LACs do, as the sample ids for genotypes are not subject to the change in orders.
	if(length(indi_common)==0)
		stop('Error: There is no common individual between the null model and the LAC file.')
	indi_la_index<-match(indi_common,indi_la)
	indi_null_index<-match(indi_common,null.obj$id_include)
	indi_gds_index<-match(indi_common,sample.id)
	sample.id<-sample.id[indi_gds_index]
	resid<-null.obj$scaled.residuals[indi_null_index]
	X<-null.obj$X[indi_null_index,]
	n_X<-ncol(X)
#	phi<-null.obj$theta['dispersion']
	if(is.null(null.obj$Sigma_i)){
		is.sparse<-F
		P<-null.obj$P[indi_null_index,indi_null_index]
	}else{
		is.sparse<-T
		sigma_i<-null.obj$Sigma_i[indi_null_index,indi_null_index]
		sigma_iX<-as.matrix(null.obj$Sigma_iX[indi_null_index,])
		Xsigma_iX<-crossprod(X,sigma_iX)
		cov<-null.obj$cov
	}
	# calculation
	variant.id<-SeqArray::seqGetData(gds,'variant.id')
	snp_read_slot<-split(variant.id,ceiling(seq_along(variant.id)/1000))
	L_list<-list()
	prev_j<-0
	prev_l_seq_j_g<-0
	for(k in seq_along(snp_read_slot)){
		SeqArray::seqSetFilter(gds,variant.id=snp_read_slot[[k]],sample.id=indi_common,verbose=F)
		miss.rate<-if(is.dosage) SeqArray::seqApply(gds,'annotation/format/DS',function(x) mean(is.na(x)),margin='by.variant',as.is='double') else SeqVarTools::missingGenotypeRate(gds,margin='by.variant')
		af<-if(is.dosage) SeqArray::seqApply(gds,'annotation/format/DS',mean,margin='by.variant',as.is='double',na.rm=T)/2 else 1-SeqVarTools::alleleFrequency(gds)
		include<-(miss.rate<=miss.cutoff & ((af>=MAF.range[1] & af<=MAF.range[2]) | (af>=1-MAF.range[2] & af<=1-MAF.range[1])))
		if(sum(include)==0){
			SeqArray::seqResetFilter(gds,verbose=F)
			next
		}
		SeqArray::seqSetFilter(gds,variant.id=snp_read_slot[[k]][include],sample.id=indi_common,verbose=F)
		pos<-SeqArray::seqGetData(gds,'position')
		j_la<-findInterval(pos,pos_la)
		j_la[j_la==0]<-1  # Need to extend the first LAC to the beginning of the chromosome if necessary, findInterval() automatically extend the last LAC to the end
#		if(k==1 && min(j_la)!=1)
#			seek(con_la,where=(min(j_la)-1)*n_unit_read_la*4,origin='current')  # will jump to the designated position by prev_j and j_extra
		SNP<-SeqArray::seqGetData(gds,'annotation/id')
		SNP[SNP=='']<-NA
		out<-data.frame(SNP=SNP,CHR=SeqArray::seqGetData(gds,'chromosome'),POS=pos)
		alleles.list<-strsplit(SeqArray::seqGetData(gds,'allele'),',')
		out$REF<-unlist(lapply(alleles.list,function(x) x[1]))
		out$ALT<-unlist(lapply(alleles.list, function(x) paste(x[-1],collapse=',')))
		out$MISSRATE<-miss.rate[include]
		out$AF<-af[include]
		geno_file<-if(is.dosage) SeqVarTools::imputedDosage(gds,use.names=FALSE) else SeqVarTools::altDosage(gds,use.names=FALSE)
		out$N<-nrow(geno_file)-colSums(is.na(geno_file))
		miss.idx<-which(is.na(geno_file))
		if(length(miss.idx)>0) {
		  geno_file[miss.idx]<-if(missing.method=='impute2mean') colMeans(geno_file,na.rm=TRUE)[ceiling(miss.idx/nrow(geno_file))] else 0
		}
		n_geno<-ncol(geno_file)
		geno_read_slot<-split(seq_len(n_geno),ceiling(seq_len(n_geno)/n.batch))
		pval_df<-data.frame()
		for(j in seq_along(geno_read_slot)){
			G<-geno_file[,geno_read_slot[[j]]]
			Gresid<-crossprod(G,resid)
			j_g<-j_la[geno_read_slot[[j]]]
			uniq_j_g<-unique(j_g)
			seq_j_g<-seq_along(uniq_j_g)
			l_seq_j_g<-length(seq_j_g)
			for(l in seq_j_g){
				if(prev_j==uniq_j_g[l]){
					L_list[[l]]<-L_list[[length(L_list)]]
				}else{
					j_extra<-uniq_j_g[l]-prev_j-1
					if(j_extra>0)
						seek(con_la,where=j_extra*n_unit_read_la*4,origin='current')
					la_coord<-readBin(con_la,what=numeric(),n=n_unit_read_la,size=4)
					la_coord<-matrix(la_coord,ncol=n_PC)[indi_la_index,]
					colnames(la_coord)<-paste0('PC',seq_len(n_PC),'_la')
					if(anyNA(la_coord) || all(la_coord==0)){
					  if(length(L_list)==0)
					    stop(paste0('Error: The initial LACs tested (',uniq_j_g[l],'th LACs) was supposed to be complete. Please check it out and fill in the missing values.'))
					  is_na_indices<-is.na(la_coord)
					  if(sum(is_na_indices)/length(is_na_indices)>0.1)
					    stop(paste0('Error: The ',uniq_j_g[l],'th LACs have more than 10% missing values. Please check it out and fill it with non-missing values.'))
					  la_coord[is_na_indices]<-L_list[[length(L_list)]][is_na_indices]
					}
					L_list[[l]]<-la_coord[,1:n_PC_used]
				}
				prev_j<-uniq_j_g[l]
			}
			L_list[seq_along(L_list)>l]<-NULL
			L<-Reduce(cbind,L_list)
			Lresid<-crossprod(L,resid)
			j_g<-split(seq_along(j_g),j_g)
			if(prev_l_seq_j_g!=l_seq_j_g){
				j_L<-rep(seq_j_g,each=n_PC_used)
				j_L<-split(seq_along(j_L),j_L)
				prev_l_seq_j_g<-l_seq_j_g
			}
			if(is.sparse){
				Gsigma_iG<-diag(crossprod(G,sigma_i %*% G))
				Gsigma_iX<-crossprod(G,sigma_iX)
				sigma_iL<-as.matrix(sigma_i %*% L)
				Xsigma_iL<-crossprod(X,sigma_iL)
				Lsigma_iL<-crossprod(L,sigma_iL)
				PL<-sigma_iL - sigma_iX %*% cov %*% Xsigma_iL
				LPL<-crossprod(L,PL)
				Gsigma_iL<-crossprod(G,sigma_iL)
				GPL<-crossprod(G,PL)
				M_var_beta<-lapply(j_L,function(x) create_M_var(Xsigma_iL=Xsigma_iL[,x],Lsigma_iL=Lsigma_iL[x,x],n_X=n_X,n_PC_used=n_PC_used))
				pval_slot<-lapply(seq_j_g,function(x) get_lm_p_value_sparse(j_g[[x]],Gresid,GPL[,j_L[[x]],drop=F],Lresid[j_L[[x]]],LPL[j_L[[x]],j_L[[x]]],Gsigma_iG,Gsigma_iX,Gsigma_iL[,j_L[[x]],drop=F],M_var_beta[[x]]))
			}else{
				GPG<-diag(crossprod(G,P %*% G))
				PL<-P %*% L
				GPL<-crossprod(G,PL)
				LPL<-crossprod(L,PL)
				pval_slot<-lapply(seq_j_g,function(x) get_lm_p_value_dense(j_g[[x]],Gresid,GPL[,j_L[[x]],drop=F],Lresid[j_L[[x]]],LPL[j_L[[x]],j_L[[x]]],GPG))
			}
			pval_slot<-Reduce(rbind,pval_slot)
			pval_df<-rbind(pval_df,pval_slot)
		}
		colnames(pval_df)<-c('SCORE','VAR','PVAL')
		out<-cbind(out,pval_df)
		write.table(out,outfile,sep='\t',row.names=F,col.names=(k==1),quote=F,append=(k>1))
		SeqArray::seqResetFilter(gds,verbose=F)
	}
	close(con_la)
	if(!inherits(geno.file,'SeqVarGDSClass'))
		SeqArray::seqClose(gds)
	else
		SeqArray::seqClose(geno.file)
	return(invisible(NULL))
}
