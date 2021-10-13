#' Graph: Circular column chart
#' @author Gabriel Danilo Shimizu, \email{shimizu@uel.br}
#' @author Leandro Simoes Azeredo Goncalves
#' @author Rodrigo Yudi Palhaci Marubayashi
#' @description Circular column chart of an experiment with a factor of interest or isolated effect of a factor
#' @param model DIC, DBC or DQL object
#' @param ylim y-axis limit
#' @param transf If the data has been transformed (\emph{default} is FALSE)
#' @param labelsize Font size of the labels
#' @seealso \link{barplot_positive}, \link{sk_graph}, \link{plot_TH}, \link{corgraph}, \link{spider_graph}, \link{line_plot}
#' @return Returns pie chart with averages and letters from the Scott-Knott cluster test
#' @export
#' @examples
#' data("laranja")
#' a=with(laranja, DBC(trat,bloco,resp, mcomp = "sk"))
#' radargraph(a)

radargraph=function(model,
                    ylim=NA,
                    labelsize=4,
                    transf=FALSE){
  a=model
  requireNamespace("ggplot2")
  data=a[[1]]$data
  size=a[[1]]$theme$axis.text$size
  ylab=a[[1]]$labels$y
  xlab=a[[1]]$labels$x
  trats=data$trats
  media=data$media
  # if(transf==FALSE){
  #   data$resp=data$resp
  #   resp=data$resp}
  # if(transf==TRUE){
  #   data$resp=data$respo
  #   resp=data$respo}
  letra=data$letra
  groups=data$groups
  data$id=id=c(1:length(data$media))
  label_data = data
  number_of_bar = nrow(data)
  angle =  90 - 360 * (label_data$id-0.5)/number_of_bar
  label_data$hjust=hjust=ifelse( angle < -90, 1, 0)
  label_data$angle=ifelse(angle < -90, angle+180, angle)
  limite=label_data$limite
  graph=ggplot(data, aes(x=trats, y=media))+
    geom_bar(aes(fill=groups),
             stat="identity",
             color="black",show.legend = FALSE)
  if(is.na(ylim)==TRUE){graph=graph+ylim(-min(media),1.2*max(media))}
  if(is.na(ylim)==FALSE){graph=graph+ylim(ylim)}
  graph=graph+theme_minimal() +
    theme(axis.text = element_blank(),
          axis.title = element_blank(),
          panel.grid = element_blank(),
          plot.margin = unit(rep(-1,4), "cm")) +
    coord_polar(start = 0) +
    geom_text(data=label_data,
              aes(x=id,
                  y=limite,
                  label=paste(data$trats,"\n",letra),
                  hjust=hjust),
              color="black",
              fontface="bold",
              alpha=0.6, size=4,
              angle= label_data$angle, inherit.aes = FALSE )
  print(graph)
  grafico=list(graph)[[1]]
  }
