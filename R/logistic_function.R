#' Analysis: Logistic regression
#'
#' @description Logistic regression is a very popular analysis in agrarian sciences, such as in fruit growth curves, seed germination, etc...The logistic function performs the analysis using 3 or 4 parameters of the logistic model, being imported from the LL function .3 or LL.4 of the drc package (Ritz & Ritz, 2016).
#' @param trat Numerical or complex vector with treatments
#' @param resp Numerical vector containing the response of the experiment.
#' @param npar Number of model parameters
#' @param error Error bar (It can be SE - \emph{default}, SD or FALSE)
#' @param ylab Variable response name (Accepts the \emph{expression}() function)
#' @param xlab Treatments name (Accepts the \emph{expression}() function)
#' @param theme ggplot2 theme (\emph{default} is theme_bw())
#' @param legend.position Legend position (\emph{default} is c(0.3,0.8))
#' @param r2 Coefficient of determination of the mean or all values (\emph{default} is all)
#' @param scale Sets x scale (\emph{default} is none, can be "log")
#' @param width.bar Bar width
#' @param textsize Font size
#' @param font.family Font family (\emph{default} is sans)
#' @return The function allows the automatic graph and equation construction of the logistic model, provides important statistics, such as the Akaike (AIC) and Bayesian (BIC) inference criteria, coefficient of determination (r2), square root of the mean error ( RMSE).
#' @details The three-parameter log-logistic function with lower limit 0 is
#' \deqn{f(x) = 0 + \frac{d}{1+\exp(b(\log(x)-\log(e)))}}
#' The four-parameter log-logistic function is given by the expression
#' \deqn{f(x) = c + \frac{d-c}{1+\exp(b(\log(x)-\log(e)))}}
#' The function is symmetric about the inflection point (e).
#' @author Model imported from the drc package (Ritz et al., 2016)
#' @author Gabriel Danilo Shimizu
#' @author Leandro Simoes Azeredo Goncalves
#' @references Seber, G. A. F. and Wild, C. J (1989) Nonlinear Regression, New York: Wiley \& Sons (p. 330).
#' @references Ritz, C.; Strebig, J.C.; Ritz, M.C. Package ‘drc’. Creative Commons: Mountain View, CA, USA, 2016.
#' @importFrom drc LL.3
#' @importFrom drc LL.4
#' @importFrom drc drm
#' @export
#'
#' @examples
#' data("emerg")
#' with(emerg, logistic(time, resp,xlab="Time (days)",ylab="Emergence (%)"))
#' with(emerg, logistic(time, resp,npar="LL.4",xlab="Time (days)",ylab="Emergence (%)"))
#'
logistic=function(trat,
                  resp,
                  npar="LL.3",
                  error="SE",
                  ylab="Dependent",
                  xlab=expression("Independent"),
                  theme=theme_classic(),
                  legend.position="top",
                  r2="all",
                  width.bar=NA,
                  scale="none",
                  textsize=12,
                  font.family="sans"){
  requireNamespace("drc")
  requireNamespace("crayon")
  requireNamespace("ggplot2")
  ymean=tapply(resp,trat,mean)
  if(is.na(width.bar)==TRUE){width.bar=0.01*mean(trat)}
  if(error=="SE"){ysd=tapply(resp,trat,sd)/sqrt(tapply(resp,trat,length))}
  if(error=="SD"){ysd=tapply(resp,trat,sd)}
  if(error=="FALSE"){ysd=0}
  desvio=ysd
  xmean=tapply(trat,trat,mean)
  if(npar=="LL.3"){mod=drm(resp~trat,fct=LL.3())
  coef=summary(mod)
  b=coef$coefficients[,1][1]
  d=coef$coefficients[,1][2]
  e=coef$coefficients[,1][3]
  if(r2=="all"){r2=cor(resp, fitted(mod))^2}
  if(r2=="mean"){r2=cor(ymean, predict(mod,newdata=data.frame(trat=unique(trat))))^2}
  r2=floor(r2*100)/100
  equation=sprintf("~~~y==frac(%0.3e, 1+e^(%0.3e*(log(x)-log(%0.3e)))) ~~~~~ italic(R^2) == %0.2f",
                   d,b,e,r2)
  xp=seq(min(trat),max(trat),length.out = 1000)
  preditos=data.frame(x=xp,
                      y=predict(mod,newdata = data.frame(trat=xp)))
  }
  if(npar=="LL.4"){mod=drm(resp~trat,fct=LL.4())
  coef=summary(mod)
  b=coef$coefficients[,1][1]
  c=coef$coefficients[,1][2]
  d=coef$coefficients[,1][3]
  e=coef$coefficients[,1][4]
  if(r2=="all"){r2=cor(resp, fitted(mod))^2}
  if(r2=="mean"){r2=cor(ymean, predict(mod,newdata=data.frame(trat=unique(trat))))^2}
  r2=floor(r2*100)/100
  equation=sprintf("~~~y == %0.3e + frac(%0.3e %s %0.3e, 1+e^(%0.3e*(log(x)-log(%0.3e)))) ~~~~~ italic(R^2) == %0.2f",
                   c,
                   d,
                   ifelse(c >= 0, "+", "-"),
                   abs(c),
                   b,
                   e,
                   r2)
  xp=seq(min(trat),max(trat),length.out = 1000)
  preditos=data.frame(x=xp,
                      y=predict(mod,newdata = data.frame(trat=xp)))}
  predesp=predict(mod)
  predobs=resp
  rmse=sqrt(mean((predesp-predobs)^2))
  x=preditos$x
  y=preditos$y
  s=equation
  data=data.frame(xmean,ymean)
  data1=data.frame(trat=xmean,resp=ymean)
  graph=ggplot(data,aes(x=xmean,y=ymean))
  if(error!="FALSE"){graph=graph+geom_errorbar(aes(ymin=ymean-ysd,ymax=ymean+ysd),
                                               width=width.bar,size=0.8)}
  graph=graph+
    geom_point(aes(color="black"),size=4.5,shape=21,fill="gray")+
    theme+
    geom_line(data=preditos,aes(x=x,
                                y=y,color="black"),size=0.8)+
    scale_color_manual(name="",values=1,label=parse(text = equation))+
    theme(axis.text = element_text(size=textsize,color="black",family = font.family),
          axis.title = element_text(family = font.family),
          legend.position = legend.position,
          legend.text = element_text(size=textsize,family = font.family),
          legend.direction = "vertical",
          legend.text.align = 0,
          legend.justification = 0)+
    ylab(ylab)+xlab(xlab)
  if(scale=="log"){graph=graph+scale_x_log10()}
  aic=AIC(mod)
  bic=BIC(mod)
  graphs=data.frame("Parameter"=c("AIC",
                                  "BIC",
                                  "r-squared",
                                  "RMSE"),
                    "values"=c(aic,
                               bic,
                               r2,
                               rmse))
  models=data.frame(coef$coefficients)
  models$Sig=ifelse(models$p.value>0.05,"ns",ifelse(models$p.value<0.01,"**","*"))
  colnames(models)=c("Estimate","Std Error","t value","P-value","")
  graficos=list("Coefficients"=models,
                "values"=graphs,
                graph=graph)
  print(graficos)
}
