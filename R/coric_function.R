#' Graph: Plot Pearson correlation with interval of confidence
#'
#' @description Plot Pearson correlation with interval of confidence
#' @param data data.frame with responses
#' @param background background fill (\emph{default} is TRUE)
#' @param axis.size Axes font size (\emph{default} is 12)
#' @param ylab Variable response name (Accepts the \emph{expression}() function)
#' @param xlab Treatments name (Accepts the \emph{expression}() function)
#' @param theme ggplot theme (\emph{default} is theme_classic())
#' @author Gabriel Danilo Shimizu, \email{shimizu@uel.br}
#' @author Leandro Simoes Azeredo Goncalves
#' @author Rodrigo Yudi Palhaci Marubayashi
#' @return The function returns a new graphical approach to correlation.
#' @importFrom RColorBrewer brewer.pal
#' @importFrom grDevices colorRampPalette
#' @importFrom utils combn
#' @importFrom grDevices blues9
#' @export
#' @examples
#' data("pomegranate")
#' cor_ic(pomegranate[,-1])

cor_ic=function(data,
                background=TRUE,
                axis.size=12,
                ylab="",
                xlab="Correlation (r)",
                theme=theme_classic()){
  method="pearson"
  requireNamespace("RColorBrewer")
  requireNamespace("ggplot2")
  make_gradient <- function(deg = 45, n = 100, cols = blues9) {
    cols <- colorRampPalette(cols)(n + 1)
    rad <- deg / (180 / pi)

    mat <- matrix(
      data = rep(seq(0, 1, length.out = n) * sin(rad), n),
      byrow = FALSE,
      ncol = n
    )+matrix(
      data = rep(seq(0, 1, length.out = n) * cos(rad), n),
      byrow = TRUE,
      ncol = n
    )
    mat <- mat - min(mat)
    mat <- mat / max(mat)
    mat <- 1 + mat * n
    mat <- matrix(data = cols[round(mat)], ncol = n)
    grid::rasterGrob(
      image = mat,
      width = unit(1, "npc"),
      height = unit(1, "npc"),
      interpolate = TRUE
    )
  }
  g <- make_gradient(
    deg = 180, n = 500, cols = brewer.pal(9, "RdBu")[9:1])
  df_list <- lapply(1:(ncol(combn(1:ncol(data), m = 2))),
                    function(y) data[, combn(1:ncol(data), m = 2)[,y]])
  # combs=factorial(length(colnames(data))-1)
  combs=length(df_list)
  combin=1:combs
  combin1=1:combs
  combin2=1:combs
  vari=1:combs
  pvalor=1:combs
  for(i in 1:combs){
    vari[i]=paste(colnames(df_list[[i]])[1],"x",
                  colnames(df_list[[i]])[2])
    combin[i]=cor.test(unlist(df_list[[i]][,1]),
                       unlist(df_list[[i]][,2]),method = method)$estimate
    combin1[i]=cor.test(unlist(df_list[[i]][,1]),
                        unlist(df_list[[i]][,2]),method = method)$conf.int[1]
    combin2[i]=cor.test(unlist(df_list[[i]][,1]),
                        unlist(df_list[[i]][,2]),method = method)$conf.int[2]
    pvalor[i]=cor.test(unlist(df_list[[i]][,1]),
                       unlist(df_list[[i]][,2]),method = method)$p.value
  }
  pvalue=ifelse(pvalor<0.01,"**",ifelse(pvalor<0.05,"*",""))
  data=data.frame(combin,combin1,combin2,vari)
  graph=ggplot(data,aes(x=combin,y=vari))

  if(background==TRUE){graph=graph+
    annotation_custom(
      grob = g, xmin = -1, xmax = 1, ymin = -Inf, ymax = Inf)}
  graph=graph+geom_vline(xintercept = c(-1,0,1),
                         lty=c(2,2,2),color=c("red","black","blue"),size=1)+
    geom_errorbar(aes(xmin=combin2,xmax=combin1),size=1,width=0.1)+
    geom_point(size=5,shape=21,color="black",fill="gray")+theme+
    geom_label(aes(label=paste(round(combin,2),pvalue,sep = "")),
               vjust=-0.5)+
    theme(axis.text = element_text(size=axis.size))+
    labs(y=ylab,
         x=xlab)
  print(graph)
}


