#' Graph: Climate chart of temperature and humidity
#'
#' @description The plot_TH function allows the user to build a column/line graph with climatic parameters of temperature (maximum, minimum and average) and relative humidity (UR) or precipitation. This chart is widely used in scientific work in agrarian science
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
#' @param legend.H Legend column
#' @param legend.tmed Legend mean temperature
#' @param legend.tmin Legend minimum temperature
#' @param legend.tmax Legend maximum temperature
#' @param x x scale type (days or data, default is "days")
#' @param breaks Range for x scale when x = "date" (default is 1 months)
#' @param textsize Axis text size
#' @param legendsize Legend text size
#' @param titlesize Axis title size
#' @param linesize Line size
#' @param date_format Date format for x="data"
#' @param sc Scale for secondary y-axis in relation to primary y-axis (declare the number of times that y2 is less than or greater than y1, the default being 2.5)
#' @param legend.position Legend position
#' @param theme ggplot2 theme
#' @param colormax Maximum line color (\emph{default} is "red")
#' @param colormin Minimum line color (\emph{default} is "blue")
#' @param colormean Midline color (\emph{default} is "darkgreen")
#' @param fillbar Column fill color (\emph{default} is "gray80")
#' @param limitsy1 Primary y-axis scale (\emph{default} is c(0,100))
#' @param angle x-axis scale text rotation
#' @return Returns row and column graphs for graphical representation of air temperature and relative humidity. Graph normally used in scientific articles
#' @seealso \link{radargraph}, \link{sk_graph}, \link{barplot_positive}, \link{corgraph}, \link{plot_TH1}, \link{spider_graph}, \link{line_plot}
#' @export
#' @examples
#' library(AgroR)
#' data(weather)
#' with(weather, plot_TH(tempo, Tmed, Tmax, Tmin, UR))

plot_TH=function(tempo,
                 Tmed,
                 Tmax,
                 Tmin,
                 UR,
                 xlab="Time",
                 yname1=expression("Humidity (%)"),
                 yname2=expression("Temperature ("^o*"C)"),
                 legend.H="Humidity",
                 legend.tmed="Tmed",
                 legend.tmin="Tmin",
                 legend.tmax="Tmax",
                 colormax="red",
                 colormin="blue",
                 colormean="darkgreen",
                 fillbar="gray80",
                 limitsy1=c(0,100),
                 x="days",
                 breaks="1 months",
                 textsize=12,
                 legendsize=12,
                 titlesize=12,
                 linesize=1,
                 date_format="%m-%Y",
                 sc=2.5,
                 angle=0,
                 legend.position="bottom",
                 theme=theme_classic()){
  data=data.frame(tempo,Tmed,Tmax,Tmin,UR)
  requireNamespace("ggplot2")
  #requireNamespace("scales")
  # Dados completos
  if(is.na(Tmed[1])==FALSE && is.na(Tmin[1])==FALSE && is.na(Tmax[1])==FALSE && is.na(UR[1])==FALSE){
    if(x=="days"){
      a=ggplot(data, aes(x = tempo)) +
        geom_col(aes(y = UR, fill=legend.H))+
        scale_x_continuous() +
        scale_y_continuous(sec.axis = sec_axis(~ . *1 ))+
        geom_line(aes(y = Tmed*sc,color=legend.tmed),size=linesize)+
        geom_line(aes(y = Tmin*sc,color=legend.tmin),size=linesize)+
        geom_line(aes(y = Tmax*sc,color=legend.tmax),size=linesize)+
        scale_y_continuous(sec.axis = sec_axis(~ . /sc,name = yname2),
                           limits = limitsy1,name=yname1)+
        scale_colour_manual(values = c(colormax,colormean,colormin))+
        scale_fill_manual(values = fillbar)+labs(fill="",color="")
      grafico=a}
    if(x=="data"){
      b=ggplot(data, aes(x = tempo)) +
        geom_col(aes(y = UR, fill=legend.H))+
        scale_x_datetime(breaks=breaks,labels = date_format(date_format)) +
        scale_y_continuous(sec.axis = sec_axis(~ . *1 ))+
        geom_line(aes(y = Tmed*sc,color=legend.tmed),size=linesize)+
        geom_line(aes(y = Tmin*sc,color=legend.tmin),size=linesize)+
        geom_line(aes(y = Tmax*sc,color=legend.tmax),size=linesize)+
        scale_y_continuous(sec.axis = sec_axis(~ . /sc,name = yname2),
                           limits = limitsy1,name=yname1)+
        scale_colour_manual(values = c(colormax,colormean,colormin))+
        scale_fill_manual(values = fillbar)+labs(fill="",color="")
      grafico=b}}

  if(is.na(Tmed[1])==FALSE && is.na(Tmin[1])==TRUE && is.na(Tmax[1])==TRUE && is.na(UR[1])==TRUE){
    if(x=="days"){
      a=ggplot(data, aes(x = tempo)) +
        scale_x_continuous() +
        geom_line(aes(y = Tmed*sc,color=legend.tmed),size=linesize,show.legend = F)+
        labs(color="",y=yname2)
      grafico=a}
    if(x=="data"){
      b=ggplot(data, aes(x = tempo)) +
        scale_x_datetime(breaks=breaks,
                         labels = date_format(date_format)) +
        geom_line(aes(y = Tmed*sc,color=legend.tmed),show.legend = FALSE,size=linesize)+
        labs(fill="",y=yname2)
      grafico=b}}

  # Sem UR
  if(is.na(Tmed[1])==FALSE && is.na(Tmin[1])==FALSE && is.na(Tmax[1])==FALSE && is.na(UR[1])==TRUE){
    if(x=="days"){
      a=ggplot(data, aes(x = tempo)) +scale_x_continuous() +
        scale_y_continuous(name=yname2)+
        geom_line(aes(y = Tmed,color=legend.tmed),size=linesize)+
        geom_line(aes(y = Tmin,color=legend.tmin),size=linesize)+
        geom_line(aes(y = Tmax,color=legend.tmax),size=linesize)+
        scale_colour_manual(values = c(colormax,colormean,colormin))+
        labs(color="")
      grafico=a}
    if(x=="data"){
      b=ggplot(data, aes(x = tempo)) +
        scale_x_datetime(breaks=breaks,labels = date_format(date_format)) +
        scale_y_continuous(name=yname2)+
        geom_line(aes(y = Tmed,color=legend.tmed),size=linesize)+
        geom_line(aes(y = Tmin,color=legend.tmin),size=linesize)+
        geom_line(aes(y = Tmax,color=legend.tmax),size=linesize)+
        scale_colour_manual(values = c(colormax,colormean,colormin))+
        labs(color="")
      grafico=b}}

  # Sem T medio
  if(is.na(Tmed[1])==TRUE && is.na(Tmin[1])==FALSE && is.na(Tmax[1])==FALSE && is.na(UR[1])==FALSE){
    if(x=="days"){
      a=ggplot(data, aes(x = tempo)) +
        geom_col(aes(y = UR, fill=legend.H))+scale_x_continuous() +
        scale_y_continuous(sec.axis = sec_axis(~ . *1 ))+
        geom_line(aes(y = Tmin*sc,color=legend.tmin),size=linesize)+
        geom_line(aes(y = Tmax*sc,color=legend.tmax),size=linesize)+
        scale_y_continuous(sec.axis = sec_axis(~ . /sc,name = yname2),
                           limits = limitsy1,name=yname1)+
        scale_colour_manual(values = c(colormax,colormin))+
        scale_fill_manual(values = fillbar)+labs(fill="",color="")
      grafico=a}
    if(x=="data"){
      b=ggplot(data, aes(x = tempo)) +
        geom_col(aes(y = UR, fill=legend.H))+
        scale_x_datetime(breaks=breaks, labels = date_format(date_format)) +
        scale_y_continuous(sec.axis = sec_axis(~ . *1 ))+
        geom_line(aes(y = Tmin*sc,color=legend.tmin),size=linesize)+
        geom_line(aes(y = Tmax*sc,color=legend.tmax),size=linesize)+
        scale_y_continuous(sec.axis = sec_axis(~ . /sc,name = yname2),
                           limits = limitsy1,name=yname1)+
        scale_colour_manual(values = c(colormax,colormin))+
        scale_fill_manual(values = fillbar)+labs(fill="",color="")
      grafico=b}}

  # Com medio e sem maximo/minimo
  if(is.na(Tmed[1])==FALSE && is.na(Tmin[1])==TRUE && is.na(Tmax[1])==TRUE && is.na(UR[1])==FALSE){
    if(x=="days"){
      a=ggplot(data, aes(x = tempo)) +
        geom_col(aes(y = UR, fill=legend.H))+scale_x_continuous() +
        scale_y_continuous(sec.axis = sec_axis(~ . *1 ))+
        geom_line(aes(y = Tmed*sc,color=legend.tmed),size=linesize)+
        scale_y_continuous(sec.axis = sec_axis(~ . /sc,name = yname2),
                           limits = limitsy1,name=yname1)+
        scale_colour_manual(values = colormean)+
        scale_fill_manual(values = fillbar)+labs(fill="",color="")
      grafico=a}
    if(x=="data"){
      b=ggplot(data, aes(x = tempo)) +
        geom_col(aes(y = UR, fill=legend.H))+
        scale_x_datetime(breaks=breaks, labels = date_format(date_format)) +
        scale_y_continuous(sec.axis = sec_axis(~ . *1 ))+
        geom_line(aes(y = Tmed*sc,color=legend.tmed),size=linesize)+
        scale_y_continuous(sec.axis = sec_axis(~ . /sc,name = yname2),
                           limits = limitsy1,name=yname1)+
        scale_colour_manual(values = colormean)+
        scale_fill_manual(values = fillbar)+labs(fill="",color="")
      grafico=b}}

  graficos=list(grafico+
                  theme+labs(x=xlab)+
                  theme(axis.text = element_text(size=textsize,
                                                 color="black"),
                        axis.text.x = element_text(angle=angle),
                        axis.title = element_text(size=titlesize),
                        legend.position = legend.position,
                        legend.text = element_text(size=legendsize)))[[1]]
  print(graficos)
}
