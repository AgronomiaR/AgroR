#' Descriptive: Boxplot with standardized data
#'
#' @description It makes a graph with the variables and/or treatments with the standardized data.
#' @author Gabriel Danilo Shimizu, \email{shimizu@uel.br}
#' @author Leandro Simoes Azeredo Goncalves
#' @author Rodrigo Yudi Palhaci Marubayashi
#' @param data data.frame containing the response of the experiment.
#' @param trat Numerical or complex vector with treatments
#' @param theme ggplot2 theme (\emph{default} is theme_bw())
#' @param ylab Variable response name (Accepts the \emph{expression}() function)
#' @param xlab Treatments name (Accepts the \emph{expression}() function)
#' @param fill Defines chart color
#' @param textsize Font size
#' @param family Font family
#' @keywords Descriptive
#' @keywords Experimental
#' @export
#' @return Returns a chart of boxes with standardized data
#' @examples
#' library(AgroR)
#' data("pomegranate")
#' dispvar(pomegranate[,-1])
#' trat=pomegranate$trat
#' dispvar(pomegranate[,-1], trat)


dispvar=function(data,
                 trat=NULL,
                 theme=theme_bw(),
                 ylab="Standard mean",
                 xlab="Variable",
                 family="serif",
                 textsize=12,
                 fill="lightblue"){
  requireNamespace("ggplot2")

if(is.null(trat)==TRUE){datap=scale(data)
  datap=data.frame(datap)
  resp=unlist(c(datap))
  trat1=rep(colnames(datap),e=length(datap[,1]))
  trat1=factor(trat1,levels = unique(trat1))
  dados=data.frame(trat1,resp)
  grafico=ggplot(dados,aes(x=trat1,y=resp))
  if(fill=="trat"){grafico=grafico+geom_boxplot(aes(fill="trat"))}
  else{grafico=grafico+geom_boxplot(fill=fill)+
    geom_jitter(fill=fill, width=0.1,alpha=0.2)+
    stat_summary(fill=fill,fun="mean",color="red",geom="point",size=2,shape=8)}
  grafico=grafico+theme+
    ylab(ylab)+
    xlab(xlab)+
    theme(text = element_text(size=textsize, family = family, colour = "black"),
          axis.title = element_text(size=textsize, family = family, colour = "black"),
          axis.text = element_text(size=textsize, family = family, colour = "black"))
  print(grafico)}
if(is.null(trat[1])==FALSE){
  datap=scale(data)
  datap=data.frame(datap)
  resp=unlist(c(datap))
  trat1=as.factor(rep(colnames(datap),e=length(datap[,1])))
  trat=as.factor(rep(trat, length(colnames(datap))))
  trat1=factor(trat1, levels=unique(trat1))
  trat=factor(trat, levels=unique(trat))
  dados=data.frame(trat1,trat,resp)
  grafico=ggplot(dados,aes(x=trat,y=resp))
  grafico=grafico+geom_boxplot(aes(fill=trat))+
    geom_jitter(aes(fill=trat), width=0.1,alpha=0.2)+
    stat_summary(aes(fill=trat),fun="mean",color="red",geom="point",size=2,shape=8)
  grafico=grafico+theme+
    ylab(ylab)+
    xlab(xlab)+facet_wrap(facets = trat1)+
    theme(text = element_text(size=textsize, family = family, colour = "black"),
          axis.title = element_text(size=textsize, family = family, colour = "black"),
          axis.text = element_text(size=textsize, family = family, colour = "black"))
  print(grafico)
}
  grafico=as.list(grafico)
}
