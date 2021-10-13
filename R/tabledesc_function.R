#' Descriptive: Table descritive analysis
#' @description Function for generating a data.frame with averages or other descriptive measures grouped by a categorical variable
#' @author Gabriel Danilo Shimizu, \email{shimizu@uel.br}
#' @author Leandro Simoes Azeredo Goncalves
#' @author Rodrigo Yudi Palhaci Marubayashi
#' @param data data.frame containing the first column with the categorical variable and the remaining response columns
#' @param fun Function of descriptive statistics (default is mean)
#' @return Returns a data.frame with a measure of dispersion or position from a dataset and separated by a factor
#' @keywords descriptive
#' @export
#' @examples
#' data(pomegranate)
#' tabledesc(pomegranate)


tabledesc=function(data,
                   fun=mean){
  dados=data[,-1]
  trat=as.vector(unlist(data[,1]))
  n=nlevels(as.factor(trat))
  nr=ncol(dados)-1
  medias=data.frame(matrix(1:(n*nr),ncol = nr))
  for(i in 1:ncol(dados)){
    medias[,i]=tapply(as.vector(unlist(dados[,i])), trat, fun, na.rm=TRUE)[unique(trat)]}
  colnames(medias)=colnames(dados)
  rownames(medias)=unique(trat)
  print(medias)}
