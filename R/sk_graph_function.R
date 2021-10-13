#' Graph: Scott-Knott graphics
#'
#' @description This is a function of the bar graph for the Scott-Knott test
#' @author Gabriel Danilo Shimizu, \email{shimizu@uel.br}
#' @author Leandro Simoes Azeredo Goncalves
#' @author Rodrigo Yudi Palhaci Marubayashi
#' @param model DIC, DBC or DQL object
#' @param horiz Horizontal Column (\emph{default} is TRUE)
#' @export
#' @return Returns a bar chart with columns separated by color according to the Scott-Knott test
#' @seealso \link{radargraph}, \link{barplot_positive}, \link{plot_TH}, \link{corgraph}, \link{spider_graph}, \link{line_plot}
#' @examples
#' data("laranja")
#'a=with(laranja, DBC(trat, bloco, resp,
#'      mcomp = "sk",angle=45,
#'      ylab = "Number of fruits/plants"))
#'sk_graph(a,horiz = FALSE)


sk_graph=function(model,
                  horiz=TRUE){
  requireNamespace("ggplot2")
  data=model[[1]]$data
  family=model[[1]]$plot$family
  media=data$media
  desvio=data$desvio
  trats=data$trats
  limite=data$limite
  letra=data$letra
  groups=data$groups
  sup=model[[1]]$plot$sup
  ploterror=model[[1]]$plot$errorbar
  # if(transf==FALSE){data=data[,c(5,1,2)]}
  # if(transf==TRUE){data=data[,c(6,3,2)]}
  if(horiz==TRUE){
    graph=ggplot(data,aes(y=as.vector(trats),
                        x=as.vector(media)))+
    model[[1]]$theme+
    geom_col(aes(fill=as.vector(groups)),
             size=0.3,
             color="black")
    if(ploterror==TRUE){graph=graph+geom_errorbar(aes(xmax=media+desvio,xmin=media-desvio),
                  width=model[[1]]$layers[3][[1]]$geom_params$width)+
    geom_label(aes(x=as.vector(media)+sup+desvio,
                   y=as.vector(trats),
                   label = letra),family=family,
               fill="lightyellow",hjust=0)}
    if(ploterror==FALSE){graph=graph+geom_label(aes(x=as.vector(media)+sup,
                     y=as.vector(trats),
                     label = letra),
                 fill="lightyellow",hjust=0)}
    graph=graph+
    labs(y=model[[1]]$labels$x,
         x=model[[1]]$labels$y)+
    theme(axis.text = element_text(size=12,color="black"),
          strip.text = element_text(size=12),
          legend.position = "none")+
    scale_y_discrete(limits=data$trats)+
    xlim(layer_scales(model[[1]])$y$range$range)}
  if(horiz==FALSE){
    graph=ggplot(data,aes(x=as.vector(trats),
                          y=as.vector(media)))+
      model[[1]]$theme+
      geom_col(aes(fill=as.vector(groups)),size=0.3,
               color="black")
    if(ploterror==TRUE){graph=graph+geom_errorbar(aes(ymax=media+desvio,ymin=media-desvio),
                                                  width=model[[1]]$layers[3][[1]]$geom_params$width)+
      geom_label(aes(y=as.vector(media)+sup+desvio,
                     x=as.vector(trats),
                     label = letra),family=family,
                 fill="lightyellow",hjust=0)}
    if(ploterror==FALSE){graph=graph+geom_label(aes(y=as.vector(media)+sup,
                                                    x=as.vector(trats),
                                                    label = letra),
                                                fill="lightyellow",hjust=0)}
    graph=graph+
      labs(x=model[[1]]$labels$x,
           y=model[[1]]$labels$y)+
      theme(axis.text = element_text(size=12,color="black"),
            strip.text = element_text(size=12),
            legend.position = "none")+
      scale_x_discrete(limits=data$trats)+
      ylim(layer_scales(model[[1]])$y$range$range)}
  graph
}
