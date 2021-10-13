#' Analysis: Linear regression graph
#' @author Gabriel Danilo Shimizu, \email{shimizu@uel.br}
#' @author Leandro Simoes Azeredo Goncalves
#' @author Rodrigo Yudi Palhaci Marubayashi
#' @description Linear regression analysis of an experiment with a quantitative factor or isolated effect of a quantitative factor
#' @param resp Numerical vector containing the response of the experiment.
#' @param trat Numerical vector with treatments (Declare as numeric)
#' @param ylab Dependent variable name (Accepts the \emph{expression}() function)
#' @param xlab Independent variable name (Accepts the \emph{expression}() function)
#' @param xname.poly X name in equation
#' @param yname.poly Y name in equation
#' @param grau Degree of the polynomial (1, 2 or 3)
#' @param theme ggplot2 theme (\emph{default} is theme_classic())
#' @param color Graph color (\emph{default} is gray80)
#' @param posi Legend position
#' @param textsize Font size
#' @param family Font family
#' @param pointsize Point size
#' @param linesize line size (Trendline and Error Bar)
#' @param ylim y-axis scale
#' @param se Adds confidence interval (\emph{default} is FALSE)
#' @param width.bar width of the error bars of a regression graph.
#' @param point Defines whether to plot mean ("mean"), all repetitions ("all"),mean with standard deviation ("mean_sd") or mean with standard error (\emph{default} - "mean_se").
#' @param n Number of decimal places for regression equations
#' @return Returns linear, quadratic or cubic regression analysis.
#' @keywords Regression
#' @keywords Experimental
#' @seealso \link{polynomial2}, \link{polynomial2_color}
#' @export
#' @examples
#' data("phao")
#' with(phao, polynomial(dose,comp, grau = 2))


polynomial=function(trat,
               resp,
               ylab="Response",
               xlab="Independent",
               yname.poly="y",
               xname.poly="x",
               grau=NA,
               theme=theme_classic(),
               point="mean_sd",
               color="gray80",
               posi="top",
               textsize=12,
               se=FALSE,
               ylim=NA,
               family="sans",
               pointsize=4.5,
               linesize=0.8,
               width.bar=NA,
               n=NA)
{requireNamespace("ggplot2")
  if(is.na(width.bar)==TRUE){width.bar=0.1*mean(trat)}
  if(is.na(grau)==TRUE){grau=1}
  # ================================
  # vetores
  # ================================
  dados=data.frame(trat,resp)
  medias=c()
  dose=tapply(trat, trat, mean, na.rm=TRUE)
  mod=c()
  mod1=c()
  mod2=c()
  modm=c()
  mod1m=c()
  mod2m=c()
  text1=c()
  text2=c()
  text3=c()
  mods=c()
  mod1s=c()
  mod2s=c()
  fparcial1=c()
  fparcial2=c()
  fparcial3=c()
  media=tapply(resp, trat, mean, na.rm=TRUE)
  desvio=tapply(resp, trat, sd, na.rm=TRUE)
  erro=tapply(resp, trat, sd, na.rm=TRUE)/sqrt(table(trat))
  dose=tapply(trat, trat, mean, na.rm=TRUE)
  moda=lm(resp~trat)
  mod1a=lm(resp~trat+I(trat^2))
  mod2a=lm(resp~trat+I(trat^2)+I(trat^3))
  mods=summary(moda)$coefficients
  mod1s=summary(mod1a)$coefficients
  mod2s=summary(mod2a)$coefficients
  modm=lm(media~dose)
  mod1m=lm(media~dose+I(dose^2))
  mod2m=lm(media~dose+I(dose^2)+I(dose^3))

  modquali=lm(resp~as.factor(trat))
  fparcial1=anova(modquali,moda)
  fparcial2=anova(modquali,mod1a)
  fparcial3=anova(modquali,mod2a)

  fparcial1=data.frame(abs(fparcial1[2,c(3,4,5,6)]))
  fparcial2=data.frame(abs(fparcial2[2,c(3,4,5,6)]))
  fparcial3=data.frame(abs(fparcial3[2,c(3,4,5,6)]))

  colnames(fparcial1)=c("GL","SQ","F","p-value");rownames(fparcial1)=""
  colnames(fparcial2)=c("GL","SQ","F","p-value");rownames(fparcial2)=""
  colnames(fparcial3)=c("GL","SQ","F","p-value");rownames(fparcial3)=""

  if(grau=="1"){r2=round(summary(modm)$r.squared, 2)}
  if(grau=="2"){r2=round(summary(mod1m)$r.squared, 2)}
  if(grau=="3"){r2=round(summary(mod2m)$r.squared, 2)}
  if(grau=="1"){
    if(is.na(n)==FALSE){coef1=round(coef(moda)[1],n)}else{coef1=coef(moda)[1]}
    if(is.na(n)==FALSE){coef2=round(coef(moda)[2],n)}else{coef2=coef(moda)[2]}
    s1=s <- sprintf("%s == %e %s %e*%s ~~~~~ italic(R^2) == %0.2f",
                    yname.poly,
                    coef1,
                    ifelse(coef2 >= 0, "+", "-"),
                    abs(coef2),
                    xname.poly,
                    r2)}
  if(grau=="2"){
    if(is.na(n)==FALSE){coef1=round(coef(mod1a)[1],n)}else{coef1=coef(mod1a)[1]}
    if(is.na(n)==FALSE){coef2=round(coef(mod1a)[2],n)}else{coef2=coef(mod1a)[2]}
    if(is.na(n)==FALSE){coef3=round(coef(mod1a)[3],n)}else{coef3=coef(mod1a)[3]}
    s2=s <- sprintf("%s == %e %s %e * %s %s %e * %s^2 ~~~~~ italic(R^2) ==  %0.2f",
                    yname.poly,
                    coef1,
                    ifelse(coef2 >= 0, "+", "-"),
                    abs(coef2),
                    xname.poly,
                    ifelse(coef3 >= 0, "+", "-"),
                    abs(coef3),
                    xname.poly,
                    r2)}
  if(grau=="3"){
    if(is.na(n)==FALSE){coef1=round(coef(mod2a)[1],n)}else{coef1=coef(mod2a)[1]}
    if(is.na(n)==FALSE){coef2=round(coef(mod2a)[2],n)}else{coef2=coef(mod2a)[2]}
    if(is.na(n)==FALSE){coef3=round(coef(mod2a)[3],n)}else{coef3=coef(mod2a)[3]}
    if(is.na(n)==FALSE){coef4=round(coef(mod2a)[4],n)}else{coef4=coef(mod2a)[4]}
    s3=s <- sprintf("%s == %e %s %e * %s %s %e * %s^2 %s %0.e * %s^3 ~~~~~ italic(R^2) == %0.2f",
                    yname.poly,
                    coef1,
                    ifelse(coef2 >= 0, "+", "-"),
                    abs(coef2),
                    xname.poly,
                    ifelse(coef3 >= 0, "+", "-"),
                    abs(coef3),
                    xname.poly,
                    ifelse(coef4 >= 0, "+", "-"),
                    abs(coef4),
                    xname.poly,
                    r2)}
  data1=data.frame(trat,resp)
  data1=data.frame(trat=dose,#as.numeric(as.character(names(media))),
                   resp=media,
                   desvio, erro)
  grafico=ggplot(data1,aes(x=trat,y=resp))
  if(point=="all"){grafico=grafico+
    geom_point(data=dados,
               aes(y=resp,x=trat),shape=21,
               fill=color,color="black")}
  if(point=="mean_sd"){grafico=grafico+
    geom_errorbar(aes(ymin=resp-desvio,ymax=resp+desvio),width=width.bar,size=linesize)}
  if(point=="mean_se"){grafico=grafico+
    geom_errorbar(aes(ymin=resp-erro,ymax=resp+erro),width=width.bar,size=linesize)}
  if(point=="mean"){grafico=grafico}
  grafico=grafico+geom_point(aes(fill=as.factor(rep(1,length(resp)))),na.rm=TRUE,
               size=pointsize,shape=21,
               color="black")+
    theme+ylab(ylab)+xlab(xlab)
  if(is.na(ylim[1])==TRUE){grafico=grafico}else{grafico=grafico+ylim(ylim)}

  if(grau=="0"){grafico=grafico+geom_line(y=mean(resp),size=linesize,lty=2)}
  if(grau=="1"){grafico=grafico+geom_smooth(method = "lm",se=se, na.rm=TRUE, formula = y~x,size=linesize,color="black")}
  if(grau=="2"){grafico=grafico+geom_smooth(method = "lm",se=se, na.rm=TRUE, formula = y~x+I(x^2),size=linesize,color="black")}
  if(grau=="3"){grafico=grafico+geom_smooth(method = "lm",se=se, na.rm=TRUE, formula = y~x+I(x^2)+I(x^3),size=linesize,color="black")}
  if(grau=="0"){grafico=grafico+
    scale_fill_manual(values=color,label=paste("y =",round(mean(resp),3)),name="")}
  if(grau=="1"){grafico=grafico+
    scale_fill_manual(values=color,label=c(parse(text=s1)),name="")}
  if(grau=="2"){grafico=grafico+
    scale_fill_manual(values=color,label=c(parse(text=s2)),name="")}
  if(grau=="3"){grafico=grafico+
    scale_fill_manual(values=color,label=c(parse(text=s3)),name="")}

  if(color=="gray"){if(grau=="1"){grafico=grafico+
    scale_fill_manual(values="black",label=c(parse(text=s1)),name="")}
    if(grau=="2"){grafico=grafico+
      scale_fill_manual(values="black",label=c(parse(text=s2)),name="")}
    if(grau=="3"){grafico=grafico+
      scale_fill_manual(values="black",label=c(parse(text=s3)),name="")}
  }

  grafico=grafico+
    theme(text = element_text(size=textsize,color="black",family=family),
          axis.text = element_text(size=textsize,color="black",family=family),
          axis.title = element_text(size=textsize,color="black",family=family),
          legend.position = posi,
          legend.text=element_text(size=textsize),
          legend.direction = "vertical",
          legend.text.align = 0,
          legend.justification = 0)

  print(grafico)
  if(grau==1){
    print(mods)
    cat("\n----------------------------------------------------\n")
    cat("Deviations from regression")
    cat("\n----------------------------------------------------\n")
    print(fparcial1)
  }
  if(grau==2){
    print(mod1s)
    cat("\n----------------------------------------------------\n")
    cat("Deviations from regression")
    cat("\n----------------------------------------------------\n")
    print(fparcial2)
  }
  if(grau==3){
    print(mod2s)
    cat("\n----------------------------------------------------\n")
    cat("Deviations from regression")
    cat("\n----------------------------------------------------\n")
    print(fparcial3)
  }
  graficos=list(grafico)
}
