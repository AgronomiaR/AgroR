#' Graph: Plot correlation
#' @description Correlation analysis function (Pearson or Spearman)
#' @param x Numeric vector with independent variable
#' @param y Numeric vector with dependent variable
#' @param method Method correlation (\emph{default} is Pearson)
#' @param ylab Variable response name (Accepts the \emph{expression}() function)
#' @param xlab Treatments name (Accepts the \emph{expression}() function)
#' @param theme ggplot2 theme (\emph{default} is theme_classic())
#' @param pointsize Point size
#' @param shape shape format
#' @param fill Fill point
#' @param color Color point
#' @param axis.size Axis text size
#' @param ic add interval of confidence
#' @param title title
#' @param family Font family
#' @author Gabriel Danilo Shimizu, \email{shimizu@uel.br}
#' @author Leandro Simoes Azeredo Goncalves
#' @author Rodrigo Yudi Palhaci Marubayashi
#' @return The function returns a graph for correlation
#' @export
#' @examples
#' data("pomegranate")
#' with(pomegranate, plot_cor(WL, SS, xlab="WL", ylab="SS"))

plot_cor=function(x,y,
                  method="pearson",
                  ylab="Dependent",
                  xlab="Independent",
                  theme=theme_classic(),
                  pointsize=5,
                  shape=21,
                  fill="gray",
                  color="black",
                  axis.size=12,
                  ic=TRUE,
                  title=NA,
                  family="sans"){
  if(is.na(title)==TRUE){
    if(method=="pearson"){title="Pearson correlation"}
    if(method=="spearman"){title="Spearman correlation"}
  }
  if(method=="pearson"){corre=cor(y,x)
  pvalor=Hmisc::rcorr(y,x)$P[1,2]}
  if(method=="spearman"){corre=cor(y,x,method = "spearman")
  pvalor=Hmisc::rcorr(y,x,type = "spearman")$P[1,2]}
  requireNamespace("ggplot2")
  data=data.frame(y,x)
  ggplot(data,aes(y=y,x=x))+
    geom_point(shape=shape,fill=fill,color=color,size=pointsize)+
    geom_smooth(color=color,method = "lm")+
    theme+labs(title=title,x=xlab,y=ylab)+
    theme(axis.text = element_text(size=axis.size,
                                   color="black",family = family),
          axis.title = element_text(family = family))+
    annotate(geom = "text",x = -Inf,y=Inf,
             label=paste("R = ", round(corre,2), ", p-value =",
                         format(round(pvalor,3),scientific = TRUE),sep = ""),
             hjust = -0.1, vjust = 1.1)
  }
