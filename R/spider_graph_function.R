#' Graph: Spider graph for sensorial analysis
#'
#' @description Spider chart or radar chart. Usually used for graphical representation of acceptability in sensory tests
#' @author Gabriel Danilo Shimizu, \email{shimizu@uel.br}
#' @author Leandro Simoes Azeredo Goncalves
#' @author Rodrigo Yudi Palhaci Marubayashi
#' @param resp Vector containing notes
#' @param vari Vector containing the variables
#' @param blend Vector containing treatments
#' @param legend.title Caption title
#' @param xlab x axis title
#' @param ylab y axis title
#' @param ymin Minimum value of y
#' @seealso \link{radargraph}, \link{sk_graph}, \link{plot_TH}, \link{corgraph}, \link{barplot_positive}, \link{line_plot}
#' @return Returns a spider or radar chart. This graph is commonly used in studies of sensory analysis.
#' @export
#' @examples
#' library(AgroR)
#' data(sensorial)
#' with(sensorial, spider_graph(resp, variable, Blend))

spider_graph=function(resp,
                      vari,
                      blend,
                      legend.title="",
                      xlab="",
                      ylab="",
                      ymin=0){
  requireNamespace("ggplot2")
  dados=data.frame(blend,vari,resp)
  dados$vari=factor(dados$vari,levels = unique(dados$vari))
  dados$blend=factor(dados$blend,levels = unique(dados$blend))
  grafico=ggplot(dados,aes(x=vari,y=resp))+
    theme_bw()+
    geom_point(aes(color=blend),size=4)+
    scale_color_manual(values=1:6)+ylim(ymin,ceiling(max(resp)))
  for(i in 1:nlevels(as.factor(blend))){
    grafico=grafico+
      geom_point(data=dados[dados$blend==levels(as.factor(dados$blend))[i],], aes(x=vari,y=resp),color=i, size=4)+
      geom_polygon(data=dados[dados$blend==levels(as.factor(dados$blend))[i],], aes(x=vari,y=resp,group=1), color=i, fill=i,alpha=0.05,size=1)}
  theta = "x"; start = 0; direction = 1
  r <- if (theta == "x") "y" else "x"
  grafico=grafico+
    coord_polar()+
    ggproto("CordRadar", CoordPolar, theta = theta, r = r,
            start = start, direction = sign(direction),
            is_linear = function(coord) TRUE)+
    theme(axis.text = element_text(size=12,color="black"),
          panel.border= element_rect(color = "white"),
          panel.grid = element_line(color="gray70"))+
    labs(y=ylab,color=legend.title,x=xlab)
  graficos=as.list(grafico)
  print(grafico)}
