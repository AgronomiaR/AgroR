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
  pearson=function(data,method="pearson"){
    corre=cor(data,method=method)
    col_combinations = expand.grid(names(data), names(data))
    cor_test_wrapper = function(col_name1, col_name2, data_frame) {
      cor.test(data_frame[[col_name1]], data_frame[[col_name2]],method=method,exact=FALSE)$p.value}
    p_vals = mapply(cor_test_wrapper,
                    col_name1 = col_combinations[[1]],
                    col_name2 = col_combinations[[2]],
                    MoreArgs = list(data_frame = data))
    pvalue=matrix(p_vals, ncol(data), ncol(data), dimnames = list(names(data), names(data)))
    list(r=corre,P=pvalue)}

  cr <- cor(dm,method = method)
  cr[upper.tri(cr, diag=TRUE)] <- NA
  dnovo=expand.grid(rownames(cr),colnames(cr))
  dnovo$cor=c(cr)
  dados=dnovo[!is.na(dnovo$cor),]
  pvalor=pearson(dm,method=method)
  pvalor=pvalor$P
  pvalor[upper.tri(pvalor, diag=TRUE)] <- NA
  dnovo1=expand.grid(rownames(pvalor),colnames(pvalor))
  dnovo1$p=c(pvalor)
  pvalor=dnovo1[!is.na(dnovo1$p),]
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
