#' Graph: Climate chart of temperature and humidity (Model 2)
#'
#' @description The plot_TH1 function allows the user to build a column/line graph with climatic parameters of temperature (maximum, minimum and average) and relative humidity (UR) or precipitation. This chart is widely used in scientific work in agrarian science
#' @author Gabriel Danilo Shimizu, \email{shimizu@uel.br}
#' @author Leandro Simoes Azeredo Goncalves
#' @author Rodrigo Yudi Palhaci Marubayashi
#' @param tempo Vector with times
#' @param Tmed Vector with mean temperature
#' @param Tmax Vector with maximum temperature
#' @param Tmin Vector with minimum temperature
#' @param UR Vector with relative humidity or precipitation
#' @param xlab x axis name
#' @param yname1 y axis name
#' @param yname2 Secondary y-axis name
#' @param legend.T faceted title legend 1
#' @param legend.H faceted title legend 2
#' @param legend.tmed Legend mean temperature
#' @param legend.tmin Legend minimum temperature
#' @param legend.tmax Legend maximum temperature
#' @param facet.fill faceted title fill color (\emph{default} is #FF9933)
#' @param panel.grid remove grid line (\emph{default} is FALSE)
#' @param x x scale type (days or data, default is "days")
#' @param breaks Range for x scale when x = "date" (default is 1 months)
#' @param textsize Axis text size
#' @param legendsize Legend text size
#' @param titlesize Axis title size
#' @param linesize Line size
#' @param date_format Date format for x="data"
#' @param legend.position Legend position
#' @param colormax Maximum line color (\emph{default} is "red")
#' @param colormin Minimum line color (\emph{default} is "blue")
#' @param colormean Midline color (\emph{default} is "darkgreen")
#' @param fillarea area fill color (\emph{default} is "darkblue")
#' @param angle x-axis scale text rotation
#' @return Returns row and column graphs for graphical representation of air temperature and relative humidity. Graph normally used in scientific articles
#' @seealso \link{radargraph}, \link{sk_graph}, \link{barplot_positive}, \link{corgraph}, \link{spider_graph}, \link{line_plot}
#' @export
#' @examples
#' library(AgroR)
#' data(weather)
#' with(weather, plot_TH1(tempo, Tmed, Tmax, Tmin, UR))

plot_TH1=function(tempo,
                 Tmed,
                 Tmax,
                 Tmin,
                 UR,
                 xlab="Time",
                 yname1=expression("Humidity (%)"),
                 yname2=expression("Temperature ("^o*"C)"),
                 legend.T="Temperature",
                 legend.H="Humidity",
                 legend.tmed="Tmed",
                 legend.tmin="Tmin",
                 legend.tmax="Tmax",
                 colormax="red",
                 colormin="blue",
                 colormean="darkgreen",
                 fillarea="darkblue",
                 facet.fill="#FF9933",
                 panel.grid=FALSE,
                 x="days",
                 breaks="1 months",
                 textsize=12,
                 legendsize=12,
                 titlesize=12,
                 linesize=1,
                 date_format="%m-%Y",
                 angle=0,
                 legend.position=c(0.1,0.3)){
  data=data.frame(tempo,Tmed,Tmax,Tmin,UR)
  requireNamespace("ggplot2")
  if(x=="days"){
      a=ggplot(data, aes(x = tempo)) +
        geom_line(aes(y = Tmed,color=legend.tmed),size=linesize)+
        geom_line(aes(y = Tmin,color=legend.tmin),size=linesize)+
        geom_line(aes(y = Tmax,color=legend.tmax),size=linesize)+
        scale_colour_manual(values = c(colormax,colormean,colormin))+
        theme_bw()+
        labs(x=xlab,y=yname2,color="")+
        facet_wrap(~as.character(legend.T))+
        theme(axis.text.y = element_text(size=textsize,
                                       color="black"),
              axis.text.x = element_blank(),
              strip.text = element_text(size=13),
              strip.background = element_rect(fill=facet.fill),
              axis.title.x = element_blank(),
              legend.position = legend.position,
              legend.text = element_text(size=legendsize,hjust = 0))
      b=ggplot(data,aes(x=tempo))+
        geom_area(aes(y=UR),fill=fillarea)+
        theme_bw()+
        labs(x=xlab,y=yname1,color="")+
        facet_wrap(~as.character(legend.H))+
        theme(axis.text = element_text(size=textsize,
                                       color="black"),
              axis.text.x = element_text(angle=angle),
              strip.text = element_text(size=13),
              strip.background = element_rect(fill=facet.fill),
              axis.title.x = element_text(size=titlesize),
              legend.position = legend.position,
              legend.text = element_text(size=legendsize))}
  if(x=="data"){
      a=ggplot(data, aes(x = tempo)) +
        scale_x_datetime(breaks=breaks,labels = date_format(date_format)) +
        geom_line(aes(y = Tmed,color=legend.tmed),size=linesize)+
        geom_line(aes(y = Tmin,color=legend.tmin),size=linesize)+
        geom_line(aes(y = Tmax,color=legend.tmax),size=linesize)+
        scale_colour_manual(values = c(colormax,colormean,colormin))+
        theme_bw()+
        labs(y=yname2,y=ylab,color="")+
        facet_wrap(~as.character(legend.T))+
        theme(axis.text.y = element_text(size=textsize,
                                       color="black"),
              axis.text.x = element_blank(),
              axis.title.x = element_blank(),
              strip.text = element_text(size=13),
              legend.position = legend.position,
              legend.text = element_text(size=legendsize,hjust = 0))
      b=ggplot(data,aes(x=tempo))+
        geom_area(aes(y=UR),fill=fillarea)+
        theme_bw()+
        labs(x=xlab,y=yname1)+
        facet_wrap(~as.character(legend.H))+
        theme(axis.text = element_text(size=textsize,
                                       color="black"),
              strip.text = element_text(size=13),
              axis.text.x = element_text(angle=angle),
              axis.title = element_text(size=titlesize),
              legend.position = legend.position,
              legend.text = element_text(size=legendsize))}
  if(panel.grid==FALSE){
    a=a+theme(panel.grid = element_blank())
    b=b+theme(panel.grid = element_blank())}
  cowplot::plot_grid(a,b,ncol=1,align = "v",rel_heights = c(3/5,2/5))
 # print(grafico)
}
