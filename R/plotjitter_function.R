#' Graph: Column, box or segment chart with observations
#'
#' @description The function performs the construction of graphs of boxes, columns or segments with all the observations represented in the graph.
#' @author Gabriel Danilo Shimizu, \email{shimizu@uel.br}
#' @author Leandro Simoes Azeredo Goncalves
#' @author Rodrigo Yudi Palhaci Marubayashi
#' @param model DIC, DBC or DQL object
#' @export
#' @return Returns with graph of boxes, columns or segments with all the observations represented in the graph.
#' @examples
#' data("pomegranate")
#' a=with(pomegranate,DIC(trat,WL,geom="point"))
#' plot_jitter(a)

plot_jitter=function(model){
  a=model
  requireNamespace("ggplot2")
  if(a[[1]]$plot$geom=="bar"){
  data=a[[1]]$data
  if(colnames(data)[3]=="respO"){
    data=data[,-1]
    colnames(data)[2]="resp"}
  sup=a[[1]]$plot$sup
  resp=a[[1]]$plot$response
  trat=a[[1]]$plot$trat
  fill=a[[1]]$plot$fill
  theme=a[[1]]$plot$theme
  ylab=a[[1]]$plot$ylab
  xlab=a[[1]]$plot$xlab
  cap=a[[1]]$labels$caption
  textsize=a[[1]]$plot$textsize
  media=data$media
  desvio=data$desvio
  limite=data$limite
  letra=data$limite
  graph=ggplot(data,aes(y=media,x=rownames(data)))+
    geom_col(fill=fill,color="black")+
    theme+
    geom_text(aes(y=limite+sup,label=letra))+
    labs(x=xlab,y=ylab,caption = cap)+
    geom_jitter(data=data.frame(trat,resp),
                aes(y=resp,x=trat),size=2,
                color="gray10", width = 0.1, alpha = 0.2)+
  geom_errorbar(aes(ymax=media+desvio,
                    ymin=media-desvio),width=0.3)+
    theme(axis.text=element_text(size = textsize,
                                 color="black"))
  print(graph)
  }
  if(a[[1]]$plot$geom=="point"){
    data=a[[1]]$data
    if(colnames(data)[3]=="respO"){
      data=data[,-1]
      colnames(data)[2]="resp"}
    sup=a[[1]]$plot$sup
    resp=a[[1]]$plot$response
    trat=a[[1]]$plot$trat
    fill=a[[1]]$plot$fill
    theme=a[[1]]$plot$theme
    ylab=a[[1]]$plot$ylab
    xlab=a[[1]]$plot$xlab
    cap=a[[1]]$labels$caption
    textsize=a[[1]]$plot$textsize
    media=data$media
    desvio=data$desvio
    limite=data$limite
    letra=data$limite
    graph=ggplot(data,aes(y=media,x=rownames(data)))+
      theme+
      geom_text(aes(y=limite+sup,label=letra))+
      labs(x=xlab,y=ylab,caption = cap)+
      geom_jitter(data=data.frame(trat,resp),
                  aes(y=resp,x=trat),size=2,
                  color="gray10", width = 0.1, alpha = 0.2)+
      geom_errorbar(aes(ymax=media+desvio,
                        ymin=media-desvio),width=0.3)+
      geom_point(fill=fill,color="black",shape=21,size=5)+
      theme(axis.text=element_text(size = textsize,
                                   color="black"))
    print(graph)
  }
  if(a[[1]]$plot$geom=="box"){
    data1=a[[1]]$data
    response=data1$response
    trat=data1$trat
    sup=a[[1]]$plot$sup
    data=a[[1]]$plot$dadosm2
    if(colnames(data)[3]=="respO"){
      data=data[,-1]
      colnames(data)[2]="resp"}
    resp=a[[1]]$plot$response
    trat=a[[1]]$plot$trat
    fill=a[[1]]$plot$fill
    theme=a[[1]]$plot$theme
    ylab=a[[1]]$plot$ylab
    xlab=a[[1]]$plot$xlab
    cap=a[[1]]$labels$caption
    textsize=a[[1]]$plot$textsize
    media=data$media
    desvio=data$desvio
    limite=data$limite
    letra=data$limite
    superior=data$superior
    graph=ggplot(data1,aes(y=response,x=trat))+
      geom_boxplot(fill=fill,color="black",outlier.color = NA)+
      theme+
      geom_text(data=data,aes(y=superior,
                              x=rownames(data),
                              label=letra))+
      labs(x=xlab,y=ylab,caption = cap)+
      geom_jitter(aes(y=resp,x=trat),size=2,
                  color="gray10", width = 0.1, alpha = 0.2)+
      theme(axis.text=element_text(size = textsize,
                                   color="black"))
    print(graph)

    }
  }
