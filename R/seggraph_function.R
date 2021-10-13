#' Graph: Segment graph for one factor
#'
#' @description This is a function of the bar graph for one factor
#' @author Gabriel Danilo Shimizu, \email{shimizu@uel.br}
#' @author Leandro Simoes Azeredo Goncalves
#' @author Rodrigo Yudi Palhaci Marubayashi
#' @param model DIC, DBC or DQL object
#' @param fill fill bars
#' @param horiz Horizontal Column (\emph{default} is TRUE)
#' @param pointsize Point size
#' @export
#' @return Returns a segment chart for one factor
#' @seealso \link{radargraph}, \link{barplot_positive}, \link{plot_TH}, \link{corgraph}, \link{spider_graph}, \link{line_plot}
#' @examples
#' data("laranja")
#'a=with(laranja, DBC(trat, bloco, resp,
#'      mcomp = "sk",angle=45,
#'      ylab = "Number of fruits/plants"))
#'seg_graph(a,horiz = FALSE)


seg_graph=function(model,
                   fill="lightblue",
                   horiz=TRUE,
                   pointsize=4.5){
  requireNamespace("ggplot2")
  data=model[[1]]$data
  media=data$media
  desvio=data$desvio
  trats=data$trats
  limite=data$limite
  letra=data$letra
  groups=data$groups
  if(horiz==TRUE){
    graph=ggplot(data,aes(y=trats,
                          x=media))+
      model[[1]]$theme+
      geom_errorbar(aes(xmin=media-desvio,
                        xmax=media+desvio),width=0.2,size=0.8)+
      geom_point(size=pointsize,shape=21, fill=fill, color="black")+
      geom_text(aes(x=media+desvio+1/15*as.vector(media),
                    y=trats,
                    label = letra),hjust=0)+
      labs(y=model[[1]]$labels$x,
           x=model[[1]]$labels$y)+
      theme(axis.text = element_text(size=12,color="black"),
            strip.text = element_text(size=12),
            legend.position = "none")+
      scale_y_discrete(limits=trats)+
      xlim(layer_scales(model[[1]])$y$range$range)}
  if(horiz==FALSE){
    graph=ggplot(data,aes(x=trats,
                          y=media))+
      model[[1]]$theme+
      geom_errorbar(aes(ymin=media-desvio,
                        ymax=media+desvio),width=0.2,size=0.8)+
      geom_point(fill=fill,size=pointsize,shape=21,color="black")+
      geom_text(aes(y=media+desvio+1/15*media,
                    x=trats,
                    label = letra),vjust=0)+
      labs(x=model[[1]]$labels$x,
           y=model[[1]]$labels$y)+
      theme(axis.text = element_text(size=12,color="black"),
            strip.text = element_text(size=12),
            legend.position = "none")+
      scale_x_discrete(limits=trats)+
      ylim(layer_scales(model[[1]])$y$range$range)}
  graph
}
