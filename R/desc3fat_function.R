#' Descriptive: Descriptive analysis (Three factors)
#'
#' @description Performs the descriptive graphical analysis of an experiment with three factors of interest.
#' @author Gabriel Danilo Shimizu, \email{shimizu@uel.br}
#' @author Leandro Simoes Azeredo Goncalves
#' @author Rodrigo Yudi Palhaci Marubayashi
#' @param f1 Numeric or complex vector with factor 1 levels
#' @param f2 Numeric or complex vector with factor 2 levels
#' @param f3 Numeric or complex vector with factor 3 levels
#' @param response Numerical vector containing the response of the experiment.
#' @param legend.title Legend title
#' @param ylab Variable response name (Accepts the \emph{expression}() function)
#' @param xlab x name (Accepts the \emph{expression}() function)
#' @param theme ggplot theme
#' @param plot "interaction" or "box"
#' @keywords Descriptive
#' @keywords Experimental
#' @return The function returns a triple interaction graph.
#' @export
#' @examples
#' library(AgroR)
#' data(enxofre)
#' with(enxofre, desc3fat(f1, f2, f3, resp))

######################################################################################
## Analise descritiva
######################################################################################

desc3fat=function(f1,
                  f2,
                  f3,
                  response,
                  legend.title="Legend",
                  xlab="xlab",
                  ylab="ylab",
                  theme=theme_classic(),
                  plot="interaction"){
  f1=as.factor(f1)
  f2=as.factor(f2)
  f3=as.factor(f3)
  requireNamespace("ggplot2")

  #===========================
  # Fator 1 x Fator 2
  #===========================
  dados=data.frame(f1,f2,f3, response)

  if(plot=="box"){
  grafico=ggplot(dados,aes(x=f1,y=response, fill=f2))+
    stat_boxplot(geom='errorbar', linetype=1,
                 position = position_dodge(width = 0.75),width=0.5)+
    geom_boxplot()+xlab(xlab)+labs(fill=legend.title)+
    ylab(ylab)+theme+facet_wrap(~f3)
  grafico=grafico+
    theme(text = element_text(size=12,color="black"),
          axis.title = element_text(size=12,color="black"),
          axis.text = element_text(size=12,color="black"),
          strip.text = element_text(size=13))
  }
  #===========================
  # Interacao
  #===========================
  if(plot=="interaction"){
  grafico=ggplot(dados,aes(x=f1,y=response, color=f2))+
    stat_summary(fun.data = "mean_se")+
    stat_summary(aes(color=f2, group=f2),
                 geom="line", fun.data = "mean_se")+
    ylab(ylab)+xlab(xlab)+labs(fill=legend.title)+theme+facet_wrap(~f3)
  grafico=grafico+
    theme(text = element_text(size=12,color="black"),
          axis.title = element_text(size=12,color="black"),
          axis.text = element_text(size=12,color="black"),
          strip.text = element_text(size=13))}
  print(grafico)
  grafico=as.list(grafico)
}

