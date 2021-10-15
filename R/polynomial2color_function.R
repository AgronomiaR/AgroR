#' Analysis: Linear regression graph in double factorial with color graph
#' @author Gabriel Danilo Shimizu, \email{shimizu@uel.br}
#' @author Leandro Simoes Azeredo Goncalves
#' @author Rodrigo Yudi Palhaci Marubayashi
#' @description Linear regression analysis for significant interaction of an experiment with two factors, one quantitative and one qualitative
#' @param fator1 Numeric or complex vector with factor 1 levels
#' @param resp Numerical vector containing the response of the experiment.
#' @param fator2 Numeric or complex vector with factor 2 levels
#' @param color Graph color (\emph{default} is NA)
#' @param grau Degree of the polynomial (1,2 or 3)
#' @param ylab Dependent variable name (Accepts the \emph{expression}() function)
#' @param xlab Independent variable name (Accepts the \emph{expression}() function)
#' @param theme ggplot2 theme (\emph{default} is theme_classic())
#' @param se Adds confidence interval (\emph{default} is FALSE)
#' @param legend.title Title legend
#' @param textsize Font size (\emph{default} is 12)
#' @param family Font family (\emph{default} is sans)
#' @param point Defines whether to plot all points ("all"), mean ("mean"), mean with standard deviation ("mean_sd") or mean with standard error (\emph{default} - "mean_se").
#' @param ylim y-axis scale
#' @param posi Legend position
#' @param width.bar width of the error bars of a regression graph.
#' @param pointsize Point size (\emph{default} is 4)
#' @param linesize line size (Trendline and Error Bar)
#' @param separate Separation between treatment and equation (\emph{default} is c("(\"","\")"))
#' @param n Number of decimal places for regression equations
#' @param SSq Sum of squares of the residue
#' @param DFres Residue freedom degrees
#' @keywords regression
#' @keywords Experimental
#' @seealso \link{polynomial}, \link{polynomial2}
#' @return Returns two or more linear, quadratic or cubic regression analyzes.
#' @export
#' @examples
#' dose=rep(c(0,0,0,2,2,2,4,4,4,6,6,6),3)
#' resp=c(8,7,5,23,24,25,30,34,36,80,90,80,
#' 12,14,15,23,24,25,50,54,56,80,90,40,
#' 12,14,15,3,4,5,50,54,56,80,90,40)
#' trat=rep(c("A","B","C"),e=12)
#' polynomial2_color(dose, resp, trat, grau=c(1,2,3))

polynomial2_color=function(fator1,
                           resp,
                           fator2,
                           color=NA,
                           grau=NA,
                           ylab="Response",
                           xlab="independent",
                           theme=theme_classic(),
                           se=FALSE,
                           point="mean_se",
                           legend.title="Tratamentos",
                           posi="top",
                           textsize=12,
                           ylim=NA,
                           family="sans",
                           width.bar=NA,
                           pointsize=5,
                           linesize=0.8,
                           separate=c("(\"","\")"),
                           n=NA,
                           DFres=NA,
                           SSq=NA){
  if(is.na(width.bar)==TRUE){width.bar=0.05*mean(fator1)}
    requireNamespace("ggplot2")
    anovaquali=anova(aov(resp~as.factor(fator1)*fator2))

  Fator2=fator2=factor(fator2,unique(fator2))
  if(is.na(color)[1]==TRUE){color=1:length(levels(Fator2))}
  if(is.na(grau)[1]==TRUE){grau=rep(1,length(levels(Fator2)))}
  curvas=c()
  texto=c()
  desvios=c()
  data=data.frame(fator1,fator2,resp)
  grafico=ggplot(data,aes(y=resp,x=fator1))+
    theme+
    ylab(ylab)+
    xlab(xlab)
  if(point=="all"){grafico=grafico+
    geom_point(aes(color=fator2),size=pointsize)}

  if(point=="mean_sd"){grafico=grafico+
    stat_summary(aes(color=fator2),fun = mean,size=linesize,
                 geom = "errorbar",
                 fun.max = function(x) mean(x) + sd(x),
                 fun.min = function(x) mean(x) - sd(x),
                 width=width.bar,
                 na.rm = TRUE)+
    stat_summary(aes(color=fator2),fun = "mean",  geom="point",
                 size=pointsize,na.rm = TRUE)}
  if(point=="mean_se"){grafico=grafico+
    stat_summary(aes(color=fator2),fun.data=mean_se, geom="errorbar",size=linesize,
                 width=width.bar,na.rm = TRUE)+
    stat_summary(aes(color=fator2),fun="mean",  geom="point", #shape=21,fill="gray80",
                 size=pointsize,na.rm = TRUE)}
  if(point=="mean"){grafico=grafico+
    stat_summary(aes(color=fator2),fun="mean", #shape=21,
                 size=pointsize, geom="point",na.rm = TRUE)}

  if(is.na(ylim[1])==TRUE){grafico=grafico}else{grafico=grafico+ylim(ylim)}
  for(i in 1:length(levels(Fator2))){
    y=resp[Fator2==levels(Fator2)[i]]
    x=fator1[Fator2==levels(Fator2)[i]]
    f2=fator2[Fator2==levels(Fator2)[i]]
    d1=data.frame(y,x,f2)
    adj=grau[i]
    numero=order(levels(Fator2))[levels(Fator2)==levels(Fator2)[i]]
    if(adj==0){mod="ns"}
    if(adj==1){mod=lm(y~x)}
    if(adj==2){mod=lm(y~x+I(x^2))}
    if(adj==3){mod=lm(y~x+I(x^2)+I(x^3))}
    if(adj==1 | adj==2 | adj==3){
      modf1=lm(y~x); modf1ql=anova(modf1)
      modf2=lm(y~x+I(x^2)); modf2ql=anova(modf2)
      modf3=lm(y~x+I(x^2)+I(x^3)); modf3ql=anova(modf3)
      modf1q=aov(y~as.factor(x))
      fadj1=anova(modf1,modf1q)[2,c(3,4,5,6)]
      fadj2=anova(modf2,modf1q)[2,c(3,4,5,6)]
      fadj3=anova(modf3,modf1q)[2,c(3,4,5,6)]
      if(is.na(DFres)==TRUE){DFres=anovaquali[4,1]}
      if(is.na(SSq)==TRUE){SSq=anovaquali[4,2]}
      df1=c(modf3ql[1,1],fadj1[1,1],DFres)
      df2=c(modf3ql[1:2,1],fadj2[1,1],DFres)
      df3=c(modf3ql[1:3,1],fadj3[1,1],DFres)
      sq1=c(modf3ql[1,2],fadj1[1,2],SSq)
      sq2=c(modf3ql[1:2,2],fadj2[1,2],SSq)
      sq3=c(modf3ql[1:3,2],fadj3[1,2],SSq)
      qm1=sq1/df1
      qm2=sq2/df2
      qm3=sq3/df3
      if(adj==1){fa1=data.frame(cbind(df1,sq1,qm1))
      fa1$f1=c(fa1$qm1[1:2]/fa1$qm1[3],NA)
      fa1$p=c(pf(fa1$f1[1:2],fa1$df1[1:2],fa1$df1[3],lower.tail = F),NA)
      rownames(fa1)=c("Linear","Deviation","Residual")
      fa=fa1}

      if(adj==2){fa2=data.frame(cbind(df2,sq2,qm2))
      fa2$f2=c(fa2$qm2[1:3]/fa2$qm2[4],NA)
      fa2$p=c(pf(fa2$f2[1:3],fa2$df2[1:3],fa2$df2[4],lower.tail = F),NA)
      rownames(fa2)=c("Linear","Quadratic","Deviation","Residual")
      fa=fa2}

      if(adj==3){
        fa3=data.frame(cbind(df3,sq3,qm3))
        fa3$f3=c(fa3$qm3[1:4]/fa3$qm3[5],NA)
        fa3$p=c(pf(fa3$f3[1:4],fa3$df3[1:4],fa3$df3[5],lower.tail = F),NA)
        rownames(fa3)=c("Linear","Quadratic","Cubic","Deviation","Residual")
        fa=fa3}
      colnames(fa)=c("Df","SSq","MSQ","F","p-value")
      desvios[[i]]=as.matrix(fa)
      curvas[[i]]=summary(mod)$coefficients
      names(desvios)[i]=paste("Anova",levels(Fator2)[i])
      names(curvas)[i]=levels(Fator2)[i]
    }
    fats=as.character(unique(Fator2)[i])
    if(adj==1){grafico=grafico+geom_smooth(data = data[fator2==levels(Fator2)[i],],aes(color=as.character(unique(fator2))),
                                           method="lm", formula = y~x,size=linesize, se=se, fill="gray70")}
    if(adj==2){grafico=grafico+geom_smooth(data = data[fator2==levels(Fator2)[i],],aes(color=unique(fator2)),
                                           method="lm", formula = y~x+I(x^2), size=linesize,se=se, fill="gray70")}
    if(adj==3){grafico=grafico+geom_smooth(data = data[fator2==levels(Fator2)[i],],aes(color=unique(fator2)),
                                           method="lm", formula = y~x+I(x^2)+I(x^3), size=linesize,se=se, fill="gray70")}
    m1=tapply(y,x,mean, na.rm=TRUE); x1=tapply(x,x,mean, na.rm=TRUE)

    if(adj==0){r2=0}
    if(adj==1){mod1=lm(m1~x1)
    r2=round(summary(mod1)$r.squared,2)}
    if(adj==2){mod1=lm(m1~x1+I(x1^2))
    r2=round(summary(mod1)$r.squared,2)}
    if(adj==3){mod1=lm(m1~x1+I(x1^2)+I(x1^3))
    r2=round(summary(mod1)$r.squared,2)}
    if(adj==0){text=sprintf("ns")}
    if(adj==1){
      if(is.na(n)==FALSE){coef1=round(coef(mod)[1],n)}else{coef1=coef(mod)[1]}
      if(is.na(n)==FALSE){coef2=round(coef(mod)[2],n)}else{coef2=coef(mod)[2]}
      text=sprintf("y == %0.3e %s %0.3e*x ~~~~~ italic(R^2) == %0.2f",
                   coef1,
                   ifelse(coef2 >= 0, "+", "-"),
                   abs(coef2),
                   r2)}
    if(adj==2){
      if(is.na(n)==FALSE){coef1=round(coef(mod)[1],n)}else{coef1=coef(mod)[1]}
      if(is.na(n)==FALSE){coef2=round(coef(mod)[2],n)}else{coef2=coef(mod)[2]}
      if(is.na(n)==FALSE){coef3=round(coef(mod)[3],n)}else{coef3=coef(mod)[3]}
      text=sprintf("y == %0.3e %s %0.3e * x %s %0.3e * x^2 ~~~~~ italic(R^2) ==  %0.2f",
                   coef1,
                   ifelse(coef2 >= 0, "+", "-"),
                   abs(coef2),
                   ifelse(coef3 >= 0, "+", "-"),
                   abs(coef3),
                   r2)}
    if(adj==3){
      if(is.na(n)==FALSE){coef1=round(coef(mod)[1],n)}else{coef1=coef(mod)[1]}
      if(is.na(n)==FALSE){coef2=round(coef(mod)[2],n)}else{coef2=coef(mod)[2]}
      if(is.na(n)==FALSE){coef3=round(coef(mod)[3],n)}else{coef3=coef(mod)[3]}
      if(is.na(n)==FALSE){coef4=round(coef(mod)[4],n)}else{coef4=coef(mod)[4]}

      text=sprintf("y == %0.3e %s %0.3e * x %s %0.3e * x^2 %s %0.3e * x^3 ~~~~~~ italic(R^2) == %0.2f",
                   coef1,
                   ifelse(coef2 >= 0, "+", "-"),
                   abs(coef2),
                   ifelse(coef3 >= 0, "+", "-"),
                   abs(coef3),
                   ifelse(coef4 >= 0, "+", "-"),
                   abs(coef4),
                   r2)}
    texto[[i]]=text
    }
  if(point=="all"){grafico=grafico+
    geom_point(aes(color=fator2),size=pointsize)}

  if(point=="mean_sd"){grafico=grafico+
    stat_summary(aes(color=fator2),fun = mean,size=0.7,
                 geom = "errorbar",
                 fun.max = function(x) mean(x) + sd(x),
                 fun.min = function(x) mean(x) - sd(x),
                 width=width.bar,
                 na.rm = TRUE)+
    stat_summary(aes(color=fator2),fun = "mean",  geom="point",
                 # shape=21,fill="gray80",
                 size=pointsize,na.rm = TRUE)}
  if(point=="mean_se"){grafico=grafico+
    stat_summary(aes(color=fator2),fun.data=mean_se, geom="errorbar",size=0.7,
                 width=width.bar,na.rm = TRUE)+
    stat_summary(aes(color=fator2),fun="mean",  geom="point", #shape=21,fill="gray80",
                 size=pointsize,na.rm = TRUE)}
  if(point=="mean"){grafico=grafico+
    stat_summary(aes(color=fator2),fun="mean", #shape=21,
                 size=pointsize, geom="point",na.rm = TRUE)}

  cat("\n----------------------------------------------------\n")
  cat("Regression Models")
  cat("\n----------------------------------------------------\n")
  print(curvas)
  cat("\n----------------------------------------------------\n")
  cat("Anova")
  cat("\n----------------------------------------------------\n")
  print(desvios,na.print=" ")
  grafico=grafico+scale_colour_discrete(label=parse(text=paste(
    separate[1],
    levels(Fator2),
    separate[2],"~",unlist(texto),sep="")))+
    theme(text = element_text(size=textsize,color="black", family = family),
          axis.text = element_text(size=textsize,color="black", family = family),
          axis.title = element_text(size=textsize,color="black", family = family),
          legend.position = posi,
          legend.text=element_text(size=textsize),
          legend.direction = "vertical",
          legend.text.align = 0,
          legend.justification = 0)+
    labs(color=legend.title, shape=legend.title, lty=legend.title)
  print(grafico)

  grafico=as.list(grafico)
}
