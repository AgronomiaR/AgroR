#' Graph: Reverse graph of DICT, DBCT and DQL output when geom="bar"
#'
#' @description The function performs the construction of a reverse graph on the output of DICT, DBCT and DQL when geom="bar".
#' @author Gabriel Danilo Shimizu, \email{shimizu@uel.br}
#' @author Leandro Simoes Azeredo Goncalves
#' @author Rodrigo Yudi Palhaci Marubayashi
#' @param plot.t DICT, DBCT or DQLT output when geom="bar"
#' @note All layout and subtitles are imported from DICT, DBCT and DQLT functions
#' @keywords Experimental
#' @seealso \link{DICT}, \link{DBCT}, \link{DQLT}
#' @return Returns a reverse graph of the output of DICT, DBCT or DQLT when geom="bar".
#' @export
#' @examples
#' data(simulate1)
#' a=with(simulate1, DICT(trat, tempo, resp,geom="bar",sup=40))
#' TBARPLOT.reverse(a)

TBARPLOT.reverse=function(plot.t){
  a=plot.t
  colo=a$plot$fill
  sup=a$plot$sup
  labelsize=a$plot$labelsize
  family=a$plot$family
  a=plot.t
  a$data$trat=factor(a$data$trat,unique(a$data$trat))
  trat=a$data$trat
  media=a$data$media
  desvio=a$data$desvio
  time=a$data$time
  letra=a$data$letra
  graph=ggplot(a$data,
         aes(x=trat,y=media,fill=time))+
    geom_col(position = position_dodge(),
             color="black")+
    a$theme+
    ylab(a$plot$ylab)+
    xlab(a$plot$xlab)+
    labs(fill=a$labels$fill)+
    geom_errorbar(aes(ymin=media-desvio,ymax=media+desvio),
                  position = position_dodge(0.90),width=0.5)+
    geom_text(aes(label=letra,y=media+desvio+sup),
              position = position_dodge(0.90),size=labelsize,
              family=family)
  if(colo=="gray"){graph=graph+scale_fill_grey()}
  graph
  }
