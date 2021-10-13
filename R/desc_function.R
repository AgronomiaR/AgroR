#' Descriptive: Descriptive analysis
#' @description Performs the descriptive analysis of an experiment with a factor of interest.
#' @author Gabriel Danilo Shimizu, \email{shimizu@uel.br}
#' @author Leandro Simoes Azeredo Goncalves
#' @author Rodrigo Yudi Palhaci Marubayashi
#' @param response Numerical vector containing the response of the experiment.
#' @param trat Numerical or complex vector with treatments
#' @param ylab Variable response name (Accepts the \emph{expression}() function)
#' @param xlab x name (Accepts the \emph{expression}() function)
#' @param ylim y-axis scale
#' @keywords Descriptive
#' @keywords Experimental
#' @seealso \link{desc2fat}, \link{tabledesc},\link{dispvar}
#' @return The function returns exploratory measures of position and dispersion, such as mean, median, maximum, minimum, coefficient of variation, etc ...
#' @export
#' @examples
#' library(AgroR)
#' data("pomegranate")
#' with(pomegranate, desc(trat,WL))

######################################################################################
## Analise descritiva
######################################################################################

desc=function(trat,
              response,
              ylab="Response",
              xlab="Treatment",
              ylim=NA){
  requireNamespace("crayon")
  requireNamespace("ggplot2")
  trat=as.factor(trat)
  Media=mean(response)
  Mediana=median(response)
  Minimo=min(response)
  Maximo=max(response)
  Variancia=var(response)
  Desvio=sd(response)
  CV=Desvio/Media*100
  juntos=cbind(Media,Mediana,Minimo,Maximo,Variancia,Desvio,CV)
  rownames(juntos)="General"
  colnames(juntos)=c("Mean","Median","Min","Max","Variance","SD","CV(%)")

  Media=tapply(response, trat, mean, na.rm=TRUE)
  Mediana=tapply(response, trat, median, na.rm=TRUE)
  Minimo=tapply(response, trat, min, na.rm=TRUE)
  Maximo=tapply(response, trat, max, na.rm=TRUE)
  Variancia=tapply(response, trat, var, na.rm=TRUE)
  Desvio=tapply(response, trat, sd, na.rm=TRUE)
  CV=Desvio/Media*100
  juntos1=cbind(Media,Mediana,Minimo,Maximo,Variancia,Desvio,CV)
  colnames(juntos1)=c("Mean","Median","Min","Max","Variance","SD","CV(%)")
  dados=data.frame(trat,response)
  grafico=ggplot(dados,aes(x=trat,y=response))+
    geom_boxplot(aes(fill=trat, group=trat),show.legend = FALSE)+
    geom_jitter(aes(group=trat),show.legend = F, width=0.1,alpha=0.2)+
    ylab(ylab)+xlab(xlab)+theme_classic()
  if(is.na(ylim)==TRUE){grafico=grafico}else{grafico=grafico+ylim(ylim)}
  grafico=grafico+
    theme(text = element_text(size=14,color="black"),
          axis.text = element_text(size=12,color="black"),
          axis.title = element_text(size=14,color="black"))+
    geom_text(aes(label=rownames(dados)),size=4, nudge_x = 0.1)
  print(grafico)
  cat(green(bold("\n-----------------------------------------------------------------\n")))
  green(italic(cat("General description")))
  cat(green(bold("\n-----------------------------------------------------------------\n")))
  print(juntos)
  cat(green(bold("\n-----------------------------------------------------------------\n")))
  green(italic(cat("Treatment")))
  cat(green(bold("\n-----------------------------------------------------------------\n")))
  print(juntos1)
  grafico=as.list(grafico)
  }

