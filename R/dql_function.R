#' Analysis: Latin square design
#'
#' @description This is a function of the AgroR package for statistical analysis of experiments conducted in Latin Square and balanced design with a factor considering the fixed model.
#' @author Gabriel Danilo Shimizu, \email{shimizu@uel.br}
#' @author Leandro Simoes Azeredo Goncalves
#' @author Rodrigo Yudi Palhaci Marubayashi
#' @param trat Numerical or complex vector with treatments
#' @param line Numerical or complex vector with lines
#' @param column Numerical or complex vector with columns
#' @param response Numerical vector containing the response of the experiment.
#' @param norm Error normality test (\emph{default} is Shapiro-Wilk)
#' @param homog Homogeneity test of variances (\emph{default} is Bartlett)
#' @param mcomp Multiple comparison test (Tukey (\emph{default}), LSD, Scott-Knott and Duncan)
#' @param quali Defines whether the factor is quantitative or qualitative (\emph{default} is qualitative)
#' @param alpha.f Level of significance of the F test (\emph{default} is 0.05)
#' @param alpha.t Significance level of the multiple comparison test (\emph{default} is 0.05)
#' @param grau Degree of polynomial in case of quantitative factor (\emph{default} is 1)
#' @param transf Applies data transformation (default is 1; for log consider 0)
#' @param geom Graph type (columns, boxes or segments)
#' @param theme ggplot2 theme (\emph{default} is theme_classic())
#' @param sup Number of units above the standard deviation or average bar on the graph
#' @param CV Plotting the coefficient of variation and p-value of Anova (\emph{default} is TRUE)
#' @param ylab Variable response name (Accepts the \emph{expression}() function)
#' @param xlab Treatments name (Accepts the \emph{expression}() function)
#' @param fill Defines chart color (to generate different colors for different treatments, define fill = "trat")
#' @param angle x-axis scale text rotation
#' @param family Font family
#' @param textsize Font size
#' @param dec Number of cells
#' @param addmean Plot the average value on the graph (\emph{default} is TRUE)
#' @param errorbar Plot the standard deviation bar on the graph (In the case of a segment and column graph) - \emph{default} is TRUE
#' @param posi Legend position
#' @param point Defines whether to plot mean ("mean"), mean with standard deviation ("mean_sd" - \emph{default}) or mean with standard error (\emph{default} - "mean_se").
#' @param angle.label label angle
#' @note The ordering of the graph is according to the sequence in which the factor levels are arranged in the data sheet. The bars of the column and segment graphs are standard deviation.
#' @note CV and p-value of the graph indicate coefficient of variation and p-value of the F test of the analysis of variance.
#' @note In the final output when transformation (transf argument) is different from 1, the columns resp and respo in the mean test are returned, indicating transformed and non-transformed mean, respectively.
#' @keywords DQL
#' @keywords Experimental
#' @references
#'
#' Principles and procedures of statistics a biometrical approach Steel, Torry and Dickey. Third Edition 1997
#'
#' Multiple comparisons theory and methods. Departament of statistics the Ohio State University. USA, 1996. Jason C. Hsu. Chapman Hall/CRC.
#'
#' Ramalho M.A.P., Ferreira D.F., Oliveira A.C. 2000. Experimentacao em Genetica e Melhoramento de Plantas. Editora UFLA.
#'
#' Scott R.J., Knott M. 1974. A cluster analysis method for grouping mans in the analysis of variance. Biometrics, 30, 507-512.
#'
#' Mendiburu, F., and de Mendiburu, M. F. (2019). Package ‘agricolae’. R Package, Version, 1-2.
#'
#' @return The table of analysis of variance, the test of normality of errors (Shapiro-Wilk, Lilliefors, Anderson-Darling, Cramer-von Mises, Pearson and Shapiro-Francia), the test of homogeneity of variances (Bartlett or Levene), the test of independence of Durbin-Watson errors, the test of multiple comparisons (Tukey, LSD, Scott-Knott or Duncan) or adjustment of regression models up to grade 3 polynomial, in the case of quantitative treatments. The column, segment or box chart for qualitative treatments is also returned. The function also returns a standardized residual plot.
#' @seealso \link{DIC}, \link{DBC}
#' @export
#' @examples
#' library(AgroR)
#' data(porco)
#' with(porco, DQL(trat, linhas, colunas, resp, ylab="Weigth (kg)"))

######################################################################################
## Analise de variancia para experimentos em DQL
######################################################################################
DQL=function(trat,
             line,
             column,
             response,
             norm="sw",
             homog="bt",
             mcomp="tukey",
             alpha.f=0.05,
             alpha.t=0.05,
             quali=TRUE,
             transf=1,
             grau=1,
             geom="bar",
             theme=theme_classic(),
             sup=NA,
             CV=TRUE,
             ylab="Response",
             xlab="",
             textsize=12,
             fill="lightblue",
             angle=0,
             family="sans",
             dec=3,
             addmean=TRUE,
             errorbar=TRUE,
             posi="top",
             point="mean_sd",
             angle.label=0)
{if(is.na(sup==TRUE)){sup=0.2*mean(response)}
  if(angle.label==0){hjust=0.5}else{hjust=0}
  requireNamespace("ScottKnott")
  requireNamespace("crayon")
  requireNamespace("ggplot2")
  requireNamespace("nortest")
  if(transf==1){resp=response}else{resp=(response^transf-1)/transf}
  if(transf==0){resp=log(response)}
  if(transf==0.5){resp=sqrt(response)}
  if(transf==-0.5){resp=1/sqrt(response)}
  if(transf==-1){resp=1/response}
  trat1=trat
  trat=as.factor(trat)
  line=as.factor(line)
  column=as.factor(column)

  a = anova(aov(resp ~ trat + line + column))
  b = aov(resp ~ trat + line + column)
  media = tapply(response, trat, mean, na.rm=TRUE)
  anava=a
  colnames(anava)=c("GL","SQ","QM","Fcal","p-value")
  respad=b$residuals/sqrt(a$`Mean Sq`[4])
  out=respad[respad>3 | respad<(-3)]
  out=names(out)
  out=if(length(out)==0)("No discrepant point")else{out}

  ## Normalidade dos erros
  if(norm=="sw"){norm1 = shapiro.test(b$res)}
  if(norm=="li"){norm1=lillie.test(b$residuals)}
  if(norm=="ad"){norm1=ad.test(b$residuals)}
  if(norm=="cvm"){norm1=cvm.test(b$residuals)}
  if(norm=="pearson"){norm1=pearson.test(b$residuals)}
  if(norm=="sf"){norm1=sf.test(b$residuals)}
  if(homog=="bt"){
    homog1 = bartlett.test(b$res ~ trat)
    statistic=homog1$statistic
    phomog=homog1$p.value
    method=paste("Bartlett test","(",names(statistic),")",sep="")
  }
  if(homog=="levene"){
    homog1 = levenehomog(b$res~trat)[1,]
    statistic=homog1$`F value`[1]
    phomog=homog1$`Pr(>F)`[1]
    method="Levene's Test (center = median)(F)"
    names(homog1)=c("Df", "F value","p.value")}

  indep = dwtest(b)

  resids=b$residuals/sqrt(a$`Mean Sq`[4])
  Ids=ifelse(resids>3 | resids<(-3), "darkblue","black")
  residplot=ggplot(data=data.frame(resids,Ids),aes(y=resids,x=1:length(resids)))+
    geom_point(shape=21,color="gray",fill="gray",size=3)+
    labs(x="",y="Standardized residuals")+
    geom_text(x=1:length(resids),label=1:length(resids),color=Ids,size=4)+
    scale_x_continuous(breaks=1:length(resids))+
    theme_classic()+theme(axis.text.y = element_text(size=12),
                          axis.text.x = element_blank())+
    geom_hline(yintercept = c(0,-3,3),lty=c(1,2,2),color="red",size=1)
  print(residplot)
  cat(green(bold("\n-----------------------------------------------------------------\n")))
  cat(green(bold("Normality of errors (Shapiro-Wilk")))
  cat(green(bold("\n-----------------------------------------------------------------\n")))
  normal=data.frame(Method=paste(norm1$method,"(",names(norm1$statistic),")",sep=""),
                    Statistic=norm1$statistic,
                    "p-value"=norm1$p.value)
  rownames(normal)=""
  print(normal)
  cat("\n")
  message(if(norm1$p.value>0.05){
    black("As the calculated p-value is greater than the 5% significance level, hypothesis H0 is not rejected. Therefore, errors can be considered normal")}
    else {"As the calculated p-value is less than the 5% significance level, H0 is rejected. Therefore, errors do not follow a normal distribution"})
  cat(green(bold("\n-----------------------------------------------------------------\n")))
  cat(green(bold("Homogeneity of Variances")))
  cat(green(bold("\n-----------------------------------------------------------------\n")))
  homoge=data.frame(Method=method,
                    Statistic=statistic,
                    "p-value"=phomog)
  rownames(homoge)=""
  print(homoge)
  cat("\n")
  message(if(homog1$p.value>0.05){
    black("As the calculated p-value is greater than the 5% significance level, hypothesis H0 is not rejected. Therefore, the variances can be considered homogeneous")}
    else {"As the calculated p-value is less than the 5% significance level, H0 is rejected. Therefore, the variances are not homogeneous"})
  cat(green(bold("\n-----------------------------------------------------------------\n")))
  cat(green(bold("Independence from errors")))
  cat(green(bold("\n-----------------------------------------------------------------\n")))
  indepe=data.frame(Method=paste(indep$method,"(",
                                 names(indep$statistic),")",sep=""),
                    Statistic=indep$statistic,
                    "p-value"=indep$p.value)
  rownames(indepe)=""
  print(indepe)
  cat("\n")
  message(if(indep$p.value>0.05){
    black("As the calculated p-value is greater than the 5% significance level, hypothesis H0 is not rejected. Therefore, errors can be considered independent")}
    else {"As the calculated p-value is less than the 5% significance level, H0 is rejected.\n Therefore, errors are not independent"})
  cat(green(bold("\n-----------------------------------------------------------------\n")))
  cat(green(bold("Additional Information")))
  cat(green(bold("\n-----------------------------------------------------------------\n")))
  cat(paste("\nCV (%) = ",round(sqrt(a$`Mean Sq`[4])/mean(resp,na.rm=TRUE)*100,2)))
  cat(paste("\nR-squared = ",round(a$`Mean Sq`[1]/(a$`Mean Sq`[4]+a$`Mean Sq`[3]+a$`Mean Sq`[2]+a$`Mean Sq`[1]),2)))
  cat(paste("\nMean = ",round(mean(response,na.rm=TRUE),4)))
  cat(paste("\nMedian = ",round(median(response,na.rm=TRUE),4)))
  cat("\nPossible outliers = ", out)
  cat("\n")
  cat(green(bold("\n-----------------------------------------------------------------\n")))
  cat(green(bold("Analysis of Variance")))
  cat(green(bold("\n-----------------------------------------------------------------\n")))
  anava1=as.matrix(data.frame(anava))
  colnames(anava1)=c("Df","Sum Sq","Mean.Sq","F value","Pr(F)" )
  print(anava1,na.print = "")
  cat("\n")
  message(if (a$`Pr(>F)`[1]<alpha.f){
    black("As the calculated p-value, it is less than the 5% significance level. The hypothesis H0 of equality of means is rejected. Therefore, at least two treatments differ")}
      else {"As the calculated p-value is greater than the 5% significance level, H0 is not rejected"})
  cat(green(bold("\n-----------------------------------------------------------------\n")))
  if(quali==TRUE){cat(green(bold("Multiple Comparison Test")))}else{cat(green(bold("Regression")))}
  cat(green(bold("\n-----------------------------------------------------------------\n")))

  # ================================
  # Comparação multipla
  # ================================
  if(quali==TRUE){
  ## Tukey
  if(mcomp=="tukey"){
    letra <- TUKEY(b, "trat", alpha=alpha.t)
    letra1 <- letra$groups; colnames(letra1)=c("resp","groups")}

  ## Scott-Knott
  if(mcomp=="sk"){
    letra=SK(b,"trat",sig.level=alpha.t)
    letra1=data.frame(resp=letra$m.inf[,1],groups=letters[letra$groups])
    letra1$resp=as.numeric(letra1$resp)}

  ## Duncan
    if(mcomp=="duncan"){
      letra <- duncan(b, "trat", alpha=alpha.t)
      letra1 <- letra$groups; colnames(letra1)=c("resp","groups")}

  ## LSD
    if(mcomp=="lsd"){
      letra <- LSD(b, "trat", alpha=alpha.t)
      letra1 <- letra$groups; colnames(letra1)=c("resp","groups")}

  media = tapply(response, trat, mean, na.rm=TRUE)

  if(transf=="1"){letra1}else{letra1$respO=media[rownames(letra1)]}
  print(if(a$`Pr(>F)`[1]<alpha.f){letra1}else{"H0 is not rejected"})

  cat("\n")
  message(if(transf=="1"){}else{blue("resp = transformed means; respO = averages without transforming")})

  if(transf==1 && norm1$p.value<0.05 | transf==1 && indep$p.value<0.05 | transf==1 &&homog1$p.value<0.05){
    message("\n \nWarning!!! Your analysis is not valid, suggests using a try to transform the data")}else{}
  if(point=="mean_sd"){
    dadosm=data.frame(letra1,
                      media=tapply(response, trat, mean, na.rm=TRUE)[rownames(letra1)],
                      desvio=tapply(response, trat, sd, na.rm=TRUE)[rownames(letra1)])}
  if(point=="mean_se"){
    dadosm=data.frame(letra1,
                      media=tapply(response, trat, mean, na.rm=TRUE)[rownames(letra1)],
                      desvio=tapply(response, trat, sd, na.rm=TRUE)/sqrt(tapply(response, trat, length))[rownames(letra1)])}
  dadosm$trats=factor(rownames(dadosm),levels = unique(trat))
  dadosm$limite=dadosm$media+dadosm$desvio
  dadosm=dadosm[unique(as.character(trat)),]
  if(addmean==TRUE){dadosm$letra=paste(format(dadosm$media,digits = dec),dadosm$groups)}
  if(addmean==FALSE){dadosm$letra=dadosm$groups}
  trats=dadosm$trats
  limite=dadosm$limite
  media=dadosm$media
  desvio=dadosm$desvio
  letra=dadosm$letra

  if(geom=="bar"){grafico=ggplot(dadosm,
                                 aes(x=trats,
                                     y=media))
    if(fill=="trat"){grafico=grafico+
      geom_col(aes(fill=trats),color=1)}
  else{grafico=grafico+
    geom_col(aes(fill=trats),fill=fill,color=1)}
  if(errorbar==TRUE){grafico=grafico+
    geom_text(aes(y=media+sup+if(sup<0){-desvio}else{desvio},
                  label=letra),family=family,angle=angle.label, hjust=hjust)}
  if(errorbar==FALSE){grafico=grafico+
    geom_text(aes(y=media+sup,label=letra),family=family,angle=angle.label, hjust=hjust)}
  if(errorbar==TRUE){grafico=grafico+
    geom_errorbar(data=dadosm,
                  aes(ymin=media-desvio,
                      ymax=media+desvio,color=1),
                  color="black", width=0.3)}}
  if(geom=="point"){grafico=ggplot(dadosm,
                                   aes(x=trats,
                                       y=media))
  if(errorbar==TRUE){grafico=grafico+
    geom_text(aes(y=media+sup+if(sup<0){-desvio}else{desvio},
                  label=letra),family=family,angle=angle.label, hjust=hjust)}
  if(errorbar==FALSE){grafico=grafico+
    geom_text(aes(y=media+sup,
                  label=letra),
              family=family,angle=angle.label, hjust=hjust)}
  if(errorbar==TRUE){grafico=grafico+
    geom_errorbar(data=dadosm,
                  aes(ymin=media-desvio,
                      ymax=media+desvio,color=1),
                  color="black", width=0.3)}
    if(fill=="trat"){grafico=grafico+
      geom_point(aes(color=trats),size=5)}
  else{grafico=grafico+
    geom_point(aes(color=trats),
               color=fill,
               size=5)}}

  if(geom=="box"){
  datam1=data.frame(trats=factor(trat,levels = unique(as.character(trat))),response)
  dadosm2=data.frame(letra1,superior=tapply(response, trat, mean, na.rm=TRUE)[rownames(letra1)])
  dadosm2$trats=rownames(dadosm2)
  dadosm2=dadosm2[unique(as.character(trat)),]
  dadosm2$limite=dadosm$media+dadosm$desvio
  dadosm2$letra=paste(format(dadosm$media,digits = dec),dadosm$groups)
  trats=dadosm2$trats
  limite=dadosm2$limite
  superior=dadosm2$superior
  letra=dadosm2$letra
  stat_box=ggplot(datam1,aes(x=trats,y=response))+geom_boxplot()
  superior=ggplot_build(stat_box)$data[[1]]$ymax
  dadosm2$superior=superior+sup

  grafico=ggplot(datam1,
                 aes(x=trats,
                     y=response))
  if(fill=="trat"){grafico=grafico+geom_boxplot(aes(fill=trats))}
  else{grafico=grafico+geom_boxplot(aes(fill=trats),fill=fill)}
  grafico=grafico+
    geom_text(data=dadosm2,
              aes(y=superior,
                  label=letra),
              family=family,angle=angle.label, hjust=hjust)}

  grafico=grafico+theme+
    ylab(ylab)+
    xlab(xlab)+
    theme(text = element_text(size=textsize,color="black",family=family),
          axis.text = element_text(size=textsize,color="black",family=family),
          axis.title = element_text(size=textsize,color="black",family=family),
          legend.position = "none")
  if(angle !=0){grafico=grafico+theme(axis.text.x=element_text(hjust = 1.01,angle = angle))}
  if(CV==TRUE){grafico=grafico+labs(caption=paste("p-value = ", if(a$`Pr(>F)`[1]<0.0001){paste("<", 0.0001)}
                                                  else{paste("=", round(a$`Pr(>F)`[1],4))},"; CV = ",
                                                  round(abs(sqrt(a$`Mean Sq`[4])/mean(resp))*100,2),"%"))}
  }

  if(quali==FALSE){
    trat=trat1
    if(grau==1){graph=polynomial(trat,response, grau = 1,xlab=xlab,ylab=ylab,textsize=textsize, family=family,posi=posi,point=point)}
    if(grau==2){graph=polynomial(trat,response, grau = 2,xlab=xlab,ylab=ylab,textsize=textsize, family=family,posi=posi,point=point)}
    if(grau==3){graph=polynomial(trat,response, grau = 3,xlab=xlab,ylab=ylab,textsize=textsize, family=family,posi=posi,point=point)}
    grafico=graph[[1]]
  }
  if(quali==TRUE){print(grafico)}
  graficos=list(grafico)#[[1]]
  }
