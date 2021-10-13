#' Analysis: Principal components analysis
#'
#' @author Gabriel Danilo Shimizu
#' @description This function performs principal component analysis.
#' @return The eigenvalues and eigenvectors, the explanation percentages of each principal component, the correlations between the vectors with the principal components, as well as graphs are returned.
#' @param data Data.frame with data set. Line name must indicate the treatment
#' @param scale Performs data standardization (\emph{default} is TRUE)
#' @param text Add label (\emph{default} is TRUE)
#' @param pointsize Point size (\emph{default} is 5)
#' @param textsize Text size (\emph{default} is 12)
#' @param labelsize Label size (\emph{default} is 4)
#' @param linesize Line size (\emph{default} is 0.8)
#' @param repel Avoid text overlay (\emph{default} is TRUE)
#' @param ylab Names y-axis
#' @param xlab Names x-axis
#' @param groups Define grouping
#' @param sc Secondary axis scale ratio (\emph{default} is 1)
#' @param font.family Font family (\emph{default} is sans)
#' @param theme Theme ggplot2 (\emph{default} is theme_bw())
#' @param label.legend Legend title (when group is not NA)
#' @param type.graph Type of chart (\emph{default} is biplot)
#' @details
#'
#' The type.graph argument defines the graph that will be returned,
#' in the case of "biplot" the biplot graph is returned with the
#' first two main components and with eigenvalues and eigenvectors.
#' In the case of "scores" only the treatment scores are returned,
#' while for "cor" the correlations are returned. For "corPCA" a
#' correlation between the vectors with the components is returned.
#'
#' @export
#' @examples
#' data(pomegranate)
#' medias=tabledesc(pomegranate)
#' PCA_function(medias)


PCA_function=function(data,
                      scale=TRUE,
                      text=TRUE,
                      pointsize=5,
                      textsize=12,
                      labelsize=4,
                      linesize=0.6,
                      repel=TRUE,
                      ylab=NA,
                      xlab=NA,
                      groups=NA,
                      sc=1,
                      font.family="sans",
                      theme=theme_bw(),
                      label.legend="Cluster",
                      type.graph="biplot"){
  nVar = ncol(data)
  if(scale==TRUE){D=scale(data)}else{D=data}
  comp=length(colnames(D))
  #fat=rep(c(-1,1),comp/2)
  Eig = eigen(var(D))
  Avl = Eig$values#*fat
  Avt = Eig$vectors#*fat
  autov=data.frame(rbind(Eigenvalue=Avl,
      Perc=Avl/sum(Avl),
      CumPer=cumsum(Avl/sum(Avl))))
  colnames(autov)=paste("PC",1:length(Avl),sep = "")
  if(is.na(xlab)==TRUE){xlab=paste("PC1 (",round(autov$PC1[2]*100,2),"%)",sep = "")}
  if(is.na(ylab)==TRUE){ylab=paste("PC2 (",round(autov$PC2[2]*100,2),"%)",sep = "")}
  # Avt[,1]=Avt[,1]*-1
  Escores = as.matrix(D) %*% Avt
  Escores2 = apply(Escores, 2, function(x) (x - mean(x))/sd(x))
  Escores2=data.frame(Escores2)
  colnames(Escores2)=paste("PC",1:ncol(Escores2),sep = "")
  rownames(Escores2)=rownames(data)
  Escores2=Escores2*sc
  setas=cor(D,Escores2)
  setas=data.frame(setas)
  requireNamespace("ggplot2")
  if(is.na(groups[1])==TRUE){groups=NA}else{groups=as.factor(groups)}

  if(type.graph=="biplot"){
  if(is.na(groups[1])==TRUE){
  PC1=Escores2$PC1
  PC2=Escores2$PC2
  graph=ggplot(Escores2,aes(y=PC2,x=PC1))+
    geom_point(data=Escores2,aes(y=PC2,x=PC1),
               shape=21,color="black",fill="gray",size=pointsize)}
  if(is.na(groups[1])==FALSE){
    PC1=Escores2$PC1
    PC2=Escores2$PC2
    graph=ggplot(Escores2,aes(y=PC2,x=PC1))+
    geom_point(data=Escores2,aes(y=PC2,x=PC1,groups=groups,fill=groups),
               shape=21,color="black",size=pointsize)+labs(fill=label.legend)}
  graph=graph+theme+
    geom_vline(xintercept = 0,lty=2,size=linesize)+
    geom_hline(yintercept = 0,lty=2,size=linesize)+
    geom_segment(data=setas,aes(x=0,y=0,
                                yend=PC2,
                                xend=PC1),size=linesize,
                 arrow = arrow(length = unit(0.14,"inches")))+
    theme(axis.text = element_text(size=textsize,color="black"))+
    xlab(xlab)+ylab(ylab)
  if(repel==FALSE){
    if(text==TRUE){
    graph=graph+
      geom_text(data=Escores2,aes(label=rownames(Escores2)),color="black",family=font.family,size=labelsize)}
    graph=graph+geom_text(data=setas,aes(label=rownames(setas)),color="black",family=font.family,size=labelsize)}
  else{requireNamespace("ggrepel")
      if(text==TRUE){graph=graph+
          geom_text_repel(data=Escores2,aes(label=rownames(Escores2)),color="black",family=font.family,size=labelsize)}
        graph=graph+geom_text_repel(data=setas,aes(label=rownames(setas)),color="black",family=font.family,size=labelsize)}}
  if(type.graph=="scores"){
    if(is.na(groups[1])==TRUE){
      PC1=Escores2$PC1
      PC2=Escores2$PC2
      graph=ggplot(Escores2,aes(y=PC2,x=PC1))+
        geom_point(data=Escores2,aes(y=PC2,x=PC1),
                   shape=21,color="black",fill="gray",size=pointsize)}
    if(is.na(groups[1])==FALSE){
      graph=ggplot(Escores2,aes(y=PC2,x=PC1))+
        geom_point(data=Escores2,aes(y=PC2,x=PC1,groups=groups,fill=groups),
                   shape=21,color="black",size=pointsize)+labs(fill=label.legend)}

    graph=graph+theme+
      geom_vline(xintercept = 0,lty=2,size=linesize)+
      geom_hline(yintercept = 0,lty=2,size=linesize)+
      theme(axis.text = element_text(size=textsize,color="black"))+
      xlab(xlab)+ylab(ylab)
    if(repel==FALSE & text==TRUE){
      graph=graph+
        geom_text(data=Escores2,aes(label=rownames(Escores2)),color="black",family=font.family,size=labelsize)}
    else{if(text==TRUE){requireNamespace("ggrepel")
      graph=graph+
        geom_text_repel(data=Escores2,aes(label=rownames(Escores2)),color="black",family=font.family,size=labelsize)}}}
  if(type.graph=="cor"){
    circle <- function(center = c(0, 0), npoints = 100) {
      r = 1
      tt = seq(0, 2 * pi, length = npoints)
      xx = center[1] + r * cos(tt)
      yy = center[1] + r * sin(tt)
      x=xx
      y=yy
      return(data.frame(x = xx, y = yy))
    }
    corcir = circle(c(0, 0), npoints = 100)
    x=corcir$x
    y=corcir$y
    graph=ggplot()+theme+
      geom_vline(xintercept = 0,lty=2,size=linesize)+
      geom_hline(yintercept = 0,lty=2,size=linesize)+
      geom_segment(data=setas,aes(x=0,y=0,
                                  yend=PC2,
                                  xend=PC1),size=linesize,
                   arrow = arrow(length = unit(0.14,"inches")))+
      theme(axis.text = element_text(size=textsize,color="black"))+
      xlab(xlab)+ylab(ylab)+
      geom_path(data = corcir, aes(x = x, y = y), colour = "gray65")

    if(repel==FALSE & text==TRUE){
      PC1=setas$PC1
      PC2=setas$PC2
      graph=graph+
        geom_text(data=setas,aes(y=PC2,x=PC1,label=rownames(setas)),
                  color="black",family=font.family,size=labelsize)}
    else{if(text==TRUE){requireNamespace("ggrepel")
      graph=graph+
        geom_text_repel(data=setas,aes(y=PC2,x=PC1,label=rownames(setas)),
                        color="black",family=font.family,size=labelsize)}}}
  autovetor=Eig$vectors
  colnames(autovetor)=colnames(autov)
  rownames(autovetor)=rownames(setas)
  colnames(Escores)=colnames(autov)
  if(type.graph=="screeplot"){
    dados=data.frame(t(autov))
    Perc=dados$Perc
    graph=ggplot(data=dados,aes(x=rownames(dados),y=Perc,group=1))+
      geom_col(fill="blue",color="black")+
      theme(axis.text = element_text(size=textsize,color="black"))+
      geom_point(size=3,color="red")+geom_line(color="red",size=linesize)+theme+
      xlab("Dimensions")+ylab("Percentage of explained variances")
  }
  if(type.graph=="corPCA"){
    cor=as.vector(unlist(setas))
    Var1=rep(rownames(setas),ncol(setas))
    Var2=rep(colnames(setas),
             e=nrow(setas))
    dados=data.frame(Var1,Var2,cor)
    graph=ggplot(dados,aes(x=Var2,
                             y=Var1,
                             fill=cor))+
      geom_tile(color="gray50",size=1)+
      scale_x_discrete(position = "top")+
      scale_fill_distiller(palette = "RdBu",direction = 1,limits=c(-1,1))+
      geom_label(aes(label=format(cor,digits=2)),
                 fill="lightyellow",label.size = 1)+
      ylab("")+xlab("")+
      labs(fill="Correlation")+
      theme(axis.text = element_text(size=12,color="black"),
            legend.text = element_text(size=12),
            legend.position = "right",
            axis.ticks = element_blank(),
            panel.background = element_blank(),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank())
  }
  #print(graph)
  list(Eigenvalue=autov,
     "Eigenvector"=autovetor,
     "Scores PCs"=Escores,
     "Correlation var x PC"=setas,
     graph=graph)}
