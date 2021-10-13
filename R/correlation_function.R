#' Graph: Correlogram
#' @description Correlation analysis function (Pearson or Spearman)
#' @param data data.frame with responses
#' @param axissize Axes font size (\emph{default} is 12)
#' @param legendsize Legend font size (\emph{default} is 12)
#' @param legendposition Legend position (\emph{default} is c(0.9,0.2))
#' @param legendtitle Legend title (\emph{default} is "Correlation")
#' @param method Method correlation (\emph{default} is Pearson)
#' @author Gabriel Danilo Shimizu, \email{shimizu@uel.br}
#' @author Leandro Simoes Azeredo Goncalves
#' @author Rodrigo Yudi Palhaci Marubayashi
#' @return The function returns a correlation matrix
#' @export
#' @examples
#' data("pomegranate")
#' corgraph(pomegranate[,-1])

corgraph=function(data,
                     axissize=12,
                     legendsize=12,
                     legendposition=c(0.9,0.2),
                     legendtitle="Correlation",
                     method="pearson"){
  dm=data
  requireNamespace("ggplot2")
  requireNamespace("Hmisc")
  requireNamespace("reshape2")
  cr <- cor(dm,method = method)
  cr[upper.tri(cr, diag=TRUE)] <- NA
  dados=melt(cr, na.rm=TRUE, value.name="cor")
  pvalor=rcorr(as.matrix(dm),type = method)
  pvalor=pvalor$P
  pvalor[upper.tri(pvalor, diag=TRUE)] <- NA
  pvalor=melt(pvalor, na.rm=TRUE, value.name="p")
  p=ifelse(unlist(pvalor$p)<0.01,"**", ifelse(unlist(pvalor$p)<0.05,"*"," "))
  dados$p=p
  dados1=data.frame(dados[,1:3],p=pvalor$p)
  print(dados1)
  Var2=dados$Var2
  Var1=dados$Var1
  cor=dados$cor
  p=dados$p
  grafico=ggplot(dados,aes(x=Var2,
                           y=Var1,
                           fill=cor))+
    geom_tile(color="gray50",size=1)+
    scale_x_discrete(position = "top")+
    scale_fill_distiller(palette = "RdBu",direction = 1,limits=c(-1,1))+
    geom_label(aes(label=paste(format(cor,digits=2),p)),
               fill="lightyellow",label.size = 1)+
    ylab("")+xlab("")+
    labs(fill=legendtitle)+
    theme(axis.text = element_text(size=axissize,color="black"),
          legend.text = element_text(size=legendsize),
          legend.position = legendposition,
          axis.ticks = element_blank(),
          panel.background = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank())+
    labs(caption = "*p<0.05; **p<0.01")
  print(grafico)
}
