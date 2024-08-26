#' generate local-ancestry-aware phenotype distribution plots
#'
#' This function generates utilizes ggplot2 to generate phenotype or covariate distribution in the context of local ancestry at a specfic genetic position and global PCs as reference. The proportion of phenotype subtypes will be distributed in pie chart ratios, and the size of pie chart indicates the number of samples sharing that specific local ancestry combination. Only works for categorical variables.
#'
#' @param input.file a dataframe that includes information of global PCs, LACs at genetic positions and phenotype or covariate of interest.
#' @param PC1,PC2 a pair of strings indicating the columns contributed to the x-, and y-coordinate in the background, which are usually supposed to be from global PCs, e.g. PC1 and PC2. Please do not include any single quote nor double quote surrounding the string.
#' @param anc a vector indicating reference ancestry populations.
#' @param PC1_la,PC2_la a pair of string indicating the columns corresponding to x-, and y-coordinate in LACs, e.g. PC1_la and PC2_la, on which the combination of local ancestry and phenotype subtypes will be counted. Please do not include any single quote nor double quote surrounding the string.
#' @param pheno a string indicating of the column of phenotype or covariate of interest.
#' @param pheno.name a string as the label of the phenotype or covariate.
#' @param phenotype.name a string vector as labels of phenotype subtypes.
#' @param y.scale an integer or a positive value, usually extreme large, adjusting the radius of pie charts to match the coordinates of PCs (default = 1e6).
#' @param if.jitter a logical switch indicating whether to slightly jitter the y-coordinates of pie charts to reduce overlapping between local ancestry patterns and pie charts (default = TRUE).
#' @param y.jitter a positive number to adjust the y-coordinates of pie charts (default = 0.01). Only effective when if.jitter = TRUE.
#' @param if.round a logical switch indicating whether to round LACs to reduce the impact of deviation in the LACs (default = TRUE).
#' @param round.digit an integer indicating the number decimal places to be used (default = 2). Only effective when if.round = TRUE.
#' @param if.filter a logical switch indicating whether to filter the number of local ancestry patterns to reduce the impact of deviation in the LACs (default = TRUE).
#' @param n.filter an integer indicating the number of threshold for local ancestry patterns preserving for further analysis (default = 5). Only effective when if.filter = TRUE.
#' @param r.x.axis,r.y.axis a pair of numerics indicating the coordinates of legend of pie chart sizes (default = c(-0.025,0.025)).
#' @param r.n.level an integer indicating the layers of circles in the legend of pie chart sizes (default = 4).
#' @return a ggplot2-based plots that illustrates the relationship between global PCs, LACs and phenotype distribution per local ancestry pattern accordingly.
#' @export
#' @concept plotting relevant features

plot.features<-function(input.file,PC1,PC2,anc,PC1_la,PC2_la,pheno,pheno.name,phenotype.name,y.scale=1e6,if.jitter=T,y.jitter=0.01,if.round=T,round.digit=2,if.filter=F,n.filter=5,r.x.axis=-0.025,r.y.axis=0.025,r.n.level=4){
  if(if.round)
    input.file<-input.file %>% mutate({{PC1_la}}:=round({{PC1_la}},round.digit),{{PC2_la}}:=round({{PC2_la}},round.digit))
  pie_file<-input.file %>% filter(!is.na({{pheno}})) %>% group_by({{pheno}},{{PC1_la}},{{PC2_la}}) %>% summarise(n=n()) %>% ungroup() %>% tidyr::complete(!!rlang::ensym(pheno),tidyr::nesting(!!rlang::ensym(PC1_la),!!rlang::ensym(PC2_la)),fill=list(n=0)) %>% tidyr::pivot_wider(names_from=!!rlang::ensym(pheno),names_prefix='n',values_from=n)
  pie_file<-pie_file %>% mutate(n=rowSums(select(.,starts_with('n'))))
  if(if.filter)
    pie_file<-pie_file %>% filter(n>=n.filter)
  if(if.jitter)
    pie_file<-pie_file %>% mutate({{PC2_la}}:={{PC2_la}}+y.jitter)
  colnames(pie_file)[-c(1,2,ncol(pie_file))]<-phenotype.name
  p<-ggplot(input.file,aes(x={{PC1}},y={{PC2}},color={{anc}},shape={{anc}})) + geom_point() + theme_bw()
  p + scatterpie::geom_scatterpie(aes(x={{PC1_la}},y={{PC2_la}},r=n/y.scale),data=pie_file,cols=phenotype.name,legend_name=pheno.name) + coord_equal() + scatterpie::geom_scatterpie_legend(pie_file$n/y.scale,x=r.x.axis,y=r.y.axis,n=r.n.level,labeller=function(x) x*y.scale)
}
