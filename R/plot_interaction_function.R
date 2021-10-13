#' Graph: Interaction plot
#' @author Gabriel Danilo Shimizu, \email{shimizu@uel.br}
#' @author Leandro Simoes Azeredo Goncalves
#' @author Rodrigo Yudi Palhaci Marubayashi
#' @description Performs an interaction graph from an output of the FAT2DIC, FAT2DBC, PSUBDIC or PSUBDBC commands.
#' @param a FAT2DIC, FAT2DBC, PSUBDIC or PSUBDBC object
#' @param repel a boolean, whether to use ggrepel to avoid overplotting text labels or not.
#' @param pointsize Point size
#' @param linesize Line size (Trendline and Error Bar)
#' @param box_label Add box in label
#' @param width.bar width of the error bars.
#' @param add.errorbar Add error bars.
#' @export
#' @import ggplot2
#' @importFrom crayon green
#' @importFrom crayon bold
#' @importFrom crayon italic
#' @importFrom crayon red
#' @importFrom crayon blue
#' @import stats
#' @return Returns an interaction graph with averages and letters from the multiple comparison test
#' @examples
#' data(cloro)
#' a=with(cloro, FAT2DIC(f1, f2, resp))
#' plot_interaction(a)

plot_interaction=function(a,
                          box_label=TRUE,
                          repel=FALSE,
                          pointsize=3,
                          linesize=0.8,
                          width.bar=0.05,
                          add.errorbar=TRUE){
  data=a[[2]]
  requireNamespace("ggplot2")
  graph=ggplot(data$data,
               aes(x=data$data[,1],
                   y=data$data$media,
                   color=data$data[,2],
                   group=data$data[,2]))+
    geom_point(show.legend = TRUE, size=pointsize)

  if(add.errorbar==TRUE){
  graph=graph+
    geom_errorbar(aes(ymax=data$data$media+data$data$desvio,
                      ymin=data$data$media-data$data$desvio),
                  width=width.bar,size=linesize)}
  graph=graph+
    geom_line(size=linesize)

  if(isTRUE(repel)==FALSE & box_label==TRUE){
    graph=graph+
    geom_label(aes(label=data$data$numero),show.legend = FALSE)}
  if(isTRUE(repel)==TRUE & box_label==TRUE){
    requireNamespace("ggrepel")
    graph=graph+
      geom_label_repel(aes(label=data$data$numero),show.legend = FALSE)}
  if(isTRUE(repel)==FALSE & box_label==FALSE){
    graph=graph+
      geom_text(aes(label=data$data$numero),show.legend = FALSE)}
  if(isTRUE(repel)==TRUE & box_label==FALSE){
    requireNamespace("ggrepel")
    graph=graph+
      geom_text_repel(aes(label=data$data$numero),show.legend = FALSE)}

  graph=graph+data$theme+
    labs(caption=data$labels$caption,
         color=data$labels$fill,
         x=data$labels$x,
         y=data$labels$y)
  print(graph)
  graphs=list(graph)
}
