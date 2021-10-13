#' Graph: Line chart
#'
#' @description Performs a descriptive line graph with standard deviation bars
#' @author Gabriel Danilo Shimizu, \email{shimizu@uel.br}
#' @author Leandro Simoes Azeredo Goncalves
#' @author Rodrigo Yudi Palhaci Marubayashi
#' @param time Vector containing the x-axis values
#' @param response Vector containing the y-axis values
#' @param factor Vector containing a categorical factor
#' @param errorbar Error bars (sd or se)
#' @param legend.position Legend position
#' @param ylab y axis title
#' @param xlab x axis title
#' @param theme ggplot2 theme (\emph{default} is theme_classic())
#' @seealso \link{radargraph}, \link{sk_graph}, \link{plot_TH}, \link{corgraph}, \link{spider_graph}
#' @return Returns a line chart with error bars
#' @export
#' @examples
#' dose=rep(c(0,2,4,6,8,10),e=3,2)
#' resp=c(seq(1,18,1),seq(2,19,1))
#' fator=rep(c("A","B"),e=18)
#' line_plot(dose,resp,fator)

line_plot=function(time,
                   response,
                   factor=NA,
                   errorbar="sd",
                   ylab="Response",
                   xlab="Time",
                   legend.position="right",
                   theme=theme_classic()){
  requireNamespace("ggplot2")
  time=as.numeric(as.character(time))
  if(is.na(factor[1])==FALSE){
    factor=as.factor(factor)
    media=tapply(response,list(time,factor), mean, na.rm=TRUE)
    if(errorbar=="sd"){erro=tapply(response,list(time,factor), sd, na.rm=TRUE)}
    if(errorbar=="se"){erro=tapply(response,list(time,factor), sd, na.rm=TRUE)/
      sqrt(tapply(response,list(time,factor), length, na.rm=TRUE))}
    dados=data.frame(time=as.numeric(as.character(rep(rownames(erro),length(levels(factor))))),
                     factor=rep(colnames(erro),e=length(unique(time))),
                     response=c(media),
                     erro=c(erro))
    graph=ggplot(dados,aes(x=time,y=response,color=factor,
                           fill=factor,group=factor))+
      geom_errorbar(aes(ymax=response+erro,
                        ymin=response-erro), size=0.8,
                    width=0.3)+
      geom_point(size=4,color="black",shape=21)+
      geom_line(size=1)+
      theme+theme(axis.text = element_text(size=12,color="black"),
                       axis.title = element_text(size=13),
                       legend.text = element_text(size=12),
                       legend.title = element_text(size=13),
                       legend.position = legend.position)+
      xlab(xlab)+ylab(ylab)
    print(graph)}

  if(is.na(factor[1])==TRUE){
    media=tapply(response,time, mean, na.rm=TRUE)
    if(errorbar=="sd"){erro=tapply(response,time, sd, na.rm=TRUE)}
    if(errorbar=="se"){erro=tapply(response,time, sd, na.rm=TRUE)/sqrt(tapply(response,time, length, na.rm=TRUE))}
    dados=data.frame(time=as.numeric(as.character(names(erro))),
                     response=c(media),
                     erro=c(erro))
    graph=ggplot(dados,aes(x=time, y=response))+
      geom_errorbar(aes(ymax=response+erro,ymin=response-erro),
                    width=0.3,size=0.8)+
      geom_line(size=1)+
      geom_point(size=5,fill="gray",color="black",shape=21)+
      theme_bw()+theme(axis.text = element_text(size=12,color="black"),
                       axis.title = element_text(size=13))+
      xlab(xlab)+ylab(ylab)
    print(graph)}
  graphs=as.list(graph)
}
