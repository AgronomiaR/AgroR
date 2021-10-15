#' Analysis: Randomized block design
#'
#' @author Gabriel Danilo Shimizu, \email{shimizu@uel.br}
#' @author Leandro Simoes Azeredo Goncalves
#' @author Rodrigo Yudi Palhaci Marubayashi
#' @description This is a function of the AgroR package for statistical analysis of experiments conducted in a randomized block and balanced design with a factor considering the fixed model.
#' @param trat Numerical or complex vector with treatments
#' @param block Numerical or complex vector with blocks
#' @param response Numerical vector containing the response of the experiment.
#' @param norm Error normality test (\emph{default} is Shapiro-Wilk)
#' @param homog Homogeneity test of variances (\emph{default} is Bartlett)
#' @param mcomp Multiple comparison test (Tukey (\emph{default}), LSD, Scott-Knott and Duncan)
#' @param quali Defines whether the factor is quantitative or qualitative (\emph{default} is qualitative)
#' @param alpha.f Level of significance of the F test (\emph{default} is 0.05)
#' @param alpha.t Significance level of the multiple comparison test (\emph{default} is 0.05)
#' @param grau Degree of polynomial in case of quantitative factor (\emph{default} is 1)
#' @param transf Applies data transformation (default is 1; for log consider 0)
#' @param test "parametric" - Parametric test or "noparametric" - non-parametric test
#' @param geom graph type (columns, boxes or segments)
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
#' @note The ordering of the graph is according to the sequence in which the factor levels are arranged in the data sheet. The bars of the column and segment graphs are standard deviation.
#' @note CV and p-value of the graph indicate coefficient of variation and p-value of the F test of the analysis of variance.
#' @note In the final output when transformation (transf argument) is different from 1, the columns resp and respo in the mean test are returned, indicating transformed and non-transformed mean, respectively.
#' @param point Defines whether to plot mean ("mean"), mean with standard deviation ("mean_sd" - \emph{default}) or mean with standard error (\emph{default} - "mean_se").
#' @param angle.label label angle
#' @keywords DBC
#' @keywords Experimental
#' @references
#'
#' Principles and procedures of statistics a biometrical approach Steel, Torry and Dickey. Third Edition 1997
#'
#' Multiple comparisons theory and methods. Departament of statistics the Ohio State University. USA, 1996. Jason C. Hsu. Chapman Hall/CRC.
#'
#' Practical Nonparametrics Statistics. W.J. Conover, 1999
#'
#' Ramalho M.A.P., Ferreira D.F., Oliveira A.C. 2000. Experimentacao em Genetica e Melhoramento de Plantas. Editora UFLA.
#'
#' Scott R.J., Knott M. 1974. A cluster analysis method for grouping mans in the analysis of variance. Biometrics, 30, 507-512.
#'
#' Mendiburu, F., and de Mendiburu, M. F. (2019). Package ‘agricolae’. R Package, Version, 1-2.
#'
#' @seealso \link{DIC}, \link{DQL}
#' @export
#' @return The table of analysis of variance, the test of normality of errors (Shapiro-Wilk, Lilliefors, Anderson-Darling, Cramer-von Mises, Pearson and Shapiro-Francia), the test of homogeneity of variances (Bartlett or Levene), the test of independence of Durbin-Watson errors, the test of multiple comparisons (Tukey, LSD, Scott-Knott or Duncan) or adjustment of regression models up to grade 3 polynomial, in the case of quantitative treatments. Non-parametric analysis can be used by the Friedman test. The column, segment or box chart for qualitative treatments is also returned. The function also returns a standardized residual plot.
#' @examples
#' library(AgroR)
#'
#' #=============================
#' # Example laranja
#' #=============================
#' data(laranja)
#' attach(laranja)
#' DBC(trat, bloco, resp, mcomp = "sk", angle=45, ylab = "Number of fruits/plants")
#'
#' #=============================
#' # Friedman test
#' #=============================
#' DBC(trat, bloco, resp, test="noparametric", ylab = "Number of fruits/plants")
#'
#' #=============================
#' # Example soybean
#' #=============================
#' data(soybean)
#' with(soybean, DBC(cult, bloc, prod,
#'                   ylab=expression("Grain yield"~(kg~ha^-1))))

DBC=function(trat,
             block,
             response,
             norm="sw",
             homog="bt",
             mcomp="tukey",
             quali=TRUE,
             alpha.f=0.05,
             alpha.t=0.05,
             transf=1,
             test="parametric",
             grau=1,
             geom="bar",
             theme=theme_classic(),
             sup=NA,
             CV=TRUE,
             ylab="response",
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
  {if(is.na(sup==TRUE)){sup=0.2*mean(response, na.rm=TRUE)}
  if(angle.label==0){hjust=0.5}else{hjust=0}
  if(test=="parametric"){
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
  bloco=as.factor(block)
  a = anova(aov(resp ~ trat + bloco))
  b = aov(resp ~ trat + bloco)
  anava=a
  colnames(anava)=c("GL","SQ","QM","Fcal","p-value")
  respad=b$residuals/sqrt(a$`Mean Sq`[3])
  out=respad[respad>3 | respad<(-3)]
  out=names(out)
  out=if(length(out)==0)("No discrepant point")else{out}

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
  resids=b$residuals/sqrt(a$`Mean Sq`[3])
  Ids=ifelse(resids>3 | resids<(-3), "darkblue","black")
  residplot=ggplot(data=data.frame(resids,Ids),aes(y=resids,x=1:length(resids)))+
    geom_point(shape=21,color="gray",fill="gray",size=3)+
    labs(x="",y="Standardized residuals")+
    geom_text(x=1:length(resids),label=1:length(resids),color=Ids,size=4)+
    scale_x_continuous(breaks=1:length(resids))+
    theme_classic()+
    theme(axis.text.y = element_text(size=12),
          axis.text.x = element_blank())+
    geom_hline(yintercept = c(0,-3,3),lty=c(1,2,2),color="red",size=1)
  print(residplot)
  cat(green(bold("\n-----------------------------------------------------------------\n")))
  cat(green(bold("Normality of errors")))
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
    else {"As the calculated p-value is less than the 5% significance level, H0 is rejected. Therefore, errors are not independent"})
  cat(green(bold("\n-----------------------------------------------------------------\n")))
  cat(green(bold("Additional Information")))
  cat(green(bold("\n-----------------------------------------------------------------\n")))
  cat(paste("\nCV (%) = ",round(sqrt(a$`Mean Sq`[3])/mean(resp,na.rm=TRUE)*100,2)))
  cat(paste("\nR-squared = ",round(a$`Mean Sq`[1]/(a$`Mean Sq`[3]+a$`Mean Sq`[2]+a$`Mean Sq`[1]),2)))
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
  if(quali==TRUE){
  if(mcomp=="tukey"){
    letra <- TUKEY(b, "trat", alpha=alpha.t)
    letra1 <- letra$groups; colnames(letra1)=c("resp","groups")}
  if(mcomp=="sk"){
    letra=SK(b,"trat",sig.level=alpha.t)
    letra1=data.frame(resp=letra$m.inf[,1],groups=letters[letra$groups])
    letra1$resp=as.numeric(letra1$resp)}
  if(mcomp=="duncan"){
    letra <- duncan(b, "trat", alpha=alpha.t)
    letra1 <- letra$groups; colnames(letra1)=c("resp","groups")}
  if(mcomp=="lsd"){
      letra <- LSD(b, "trat", alpha=alpha.t)
      letra1 <- letra$groups; colnames(letra1)=c("resp","groups")}
    media = tapply(response, trat, mean, na.rm=TRUE)
  if(transf=="1"){letra1}else{letra1$respO=media[rownames(letra1)]}
  print(if(a$`Pr(>F)`[1]<alpha.f){letra1}else{"H0 is not rejected"})
  cat("\n")
  cat(if(transf=="1"){}else{blue("resp = transformed means; respO = averages without transforming")})
  if(transf==1 && norm1$p.value<0.05 | transf==1 && indep$p.value<0.05 | transf==1 &&homog1$p.value<0.05){
    message("\nYour analysis is not valid, suggests using a non-parametric \ntest and try to transform the data")}
  else{}
  if(transf != 1 && norm1$p.value<0.05 | transf!=1 && indep$p.value<0.05 | transf!=1 && homog1$p.value<0.05){cat(red("\n \nWarning!!! Your analysis is not valid, suggests using a non-parametric \ntest"))}else{}
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
  if(geom=="bar"){grafico=ggplot(dadosm,aes(x=trats,
                                            y=media))
    if(fill=="trat"){grafico=grafico+
      geom_col(aes(fill=trats),color=1)}
  else{grafico=grafico+geom_col(aes(fill=trats),
                                fill=fill,color=1)}
  if(errorbar==TRUE){grafico=grafico+
    geom_text(aes(y=media+sup+if(sup<0){-desvio}else{desvio},label=letra),
              family=family,angle=angle.label, hjust=hjust)}
  if(errorbar==FALSE){grafico=grafico+
    geom_text(aes(y=media+sup,label=letra),
              family=family,angle=angle.label, hjust=hjust)}
  if(errorbar==TRUE){grafico=grafico+
    geom_errorbar(data=dadosm,aes(ymin=media-desvio,
                                  ymax=media+desvio,color=1),
                  color="black", width=0.3)}}
  if(geom=="point"){grafico=ggplot(dadosm,aes(x=trats,
                                              y=media))
  if(errorbar==TRUE){grafico=grafico+
    geom_text(aes(y=media+sup+if(sup<0){-desvio}else{desvio},label=letra),
              family=family,angle=angle.label, hjust=hjust)}
  if(errorbar==FALSE){grafico=grafico+
    geom_text(aes(y=media+sup,label=letra),family=family,angle=angle.label, hjust=hjust)}
  if(errorbar==TRUE){grafico=grafico+
    geom_errorbar(data=dadosm,aes(ymin=media-desvio,
                                  ymax=media+desvio,color=1),
                  color="black", width=0.3)}
    if(fill=="trat"){grafico=grafico+
      geom_point(aes(color=trats),size=5)}
  else{grafico=grafico+
    geom_point(aes(color=trats),fill="gray",pch=21,color="black",size=5)}}

  if(geom=="box"){
  datam1=data.frame(trats=factor(trat,levels = unique(as.character(trat))),response)
  dadosm2=data.frame(letra1,
                     superior=tapply(response, trat, mean, na.rm=TRUE)[rownames(letra1)])
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

  grafico=ggplot(datam1,aes(x=trats,
                            y=response))
    if(fill=="trat"){grafico=grafico+
      geom_boxplot(aes(fill=trats))}
  else{grafico=grafico+
    geom_boxplot(aes(fill=trats),fill=fill)}
  grafico=grafico+
    geom_text(data=dadosm2,aes(y=superior,label=letra),
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
                                                  round(abs(sqrt(a$`Mean Sq`[3])/mean(resp))*100,2),"%"))}
  }
  if(quali==FALSE){
    trat=trat1
    if(grau==1){graph=polynomial(trat,response, grau = 1,xlab=xlab,ylab=ylab,textsize=textsize, family=family,posi=posi,point=point,SSq = a$`Sum Sq`[3],DFres = a$Df[3])}
    if(grau==2){graph=polynomial(trat,response, grau = 2,xlab=xlab,ylab=ylab,textsize=textsize, family=family,posi=posi,point=point,SSq = a$`Sum Sq`[3],DFres = a$Df[3])}
    if(grau==3){graph=polynomial(trat,response, grau = 3,xlab=xlab,ylab=ylab,textsize=textsize, family=family,posi=posi,point=point,SSq = a$`Sum Sq`[3],DFres = a$Df[3])}
    grafico=graph[[1]]
  }}
  if(test=="noparametric"){
    friedman=function(judge,trt,evaluation,alpha=0.05,group=TRUE,main=NULL,console=FALSE){
      name.x <- paste(deparse(substitute(judge)))
      name.y <- paste(deparse(substitute(evaluation)))
      name.t <- paste(deparse(substitute(trt)))
      name.j <- paste(deparse(substitute(judge)))
      if(is.null(main))main<-paste(name.y,"~", name.j,"+",name.t)
      datos <- data.frame(judge, trt, evaluation)
      matriz <- by(datos[,3], datos[,1:2], function(x) mean(x,na.rm=TRUE))
      matriz <-as.data.frame(matriz[,])
      name<-as.character(colnames(matriz))
      ntr <-length(name)
      m<-dim(matriz)
      v<-array(0,m)
      for (i in 1:m[1]){
        v[i,]<-rank(matriz[i,])
      }
      vv<-as.numeric(v)
      junto <- data.frame(evaluation, trt)
      medians<-mean.stat(junto[,1],junto[,2],stat="median")
      for(i in c(1,5,2:4)) {
        x <- mean.stat(junto[,1],junto[,2],function(x)quantile(x)[i])
        medians<-cbind(medians,x[,2])}
      medians<-medians[,3:7]
      names(medians)<-c("Min","Max","Q25","Q50","Q75")
      Means <- mean.stat(junto[,1],junto[,2],stat="mean")
      sds <-   mean.stat(junto[,1],junto[,2],stat="sd")
      nn <-   mean.stat(junto[,1],junto[,2],stat="length")
      nr<-unique(nn[,2])
      s<-array(0,m[2])
      for (j in 1:m[2]){
        s[j]<-sum(v[,j])}
      Means<-data.frame(Means,std=sds[,2],r=nn[,2],medians)
      names(Means)[1:2]<-c(name.t,name.y)
      means<-Means[,c(1:2,4)]
      rownames(Means)<-Means[,1]
      Means<-Means[,-1]
      means[,2]<-s
      rs<-array(0,m[2])
      rs<-s-m[1]*(m[2]+1)/2
      T1<-12*t(rs)%*%rs/(m[1]*m[2]*(m[2]+1))
      T2<-(m[1]-1)*T1/(m[1]*(m[2]-1)-T1)
      if(console){
        cat("\nStudy:",main,"\n\n")
        cat(paste(name.t,",",sep="")," Sum of the ranks\n\n")
        print(data.frame(row.names = means[,1], means[,-1]))
        cat("\nFriedman")
        cat("\n===============")}
      A1<-0
      for (i in 1:m[1]) A1 <- A1 + t(v[i,])%*%v[i,]
      DFerror <-(m[1]-1)*(m[2]-1)
      Tprob<-qt(1-alpha/2,DFerror)
      LSD<-as.numeric(Tprob*sqrt(2*(m[1]*A1-t(s)%*%s)/DFerror))
      C1 <-m[1]*m[2]*(m[2]+1)^2/4
      T1.aj <-(m[2]-1)*(t(s)%*%s-m[1]*C1)/(A1-C1)
      T2.aj <-(m[1]-1)*T1.aj/(m[1]*(m[2]-1)-T1.aj)
      p.value<-1-pchisq(T1.aj,m[2]-1)
      p.noadj<-1-pchisq(T1,m[2]-1)
      PF<-1-pf(T2.aj, ntr-1, (ntr-1)*(nr-1) )
      if(console){
        cat("\nAdjusted for ties")
        cat("\nCritical Value:",T1.aj)
        cat("\nP.Value Chisq:",p.value)
        cat("\nF Value:",T2.aj)
        cat("\nP.Value F:",PF,"\n")
        cat("\nPost Hoc Analysis\n")
      }
      statistics<-data.frame(Chisq=T1.aj,Df=ntr-1,p.chisq=p.value,F=T2.aj,DFerror=DFerror,p.F=PF,t.value=Tprob,LSD)
      if ( group & length(nr) == 1 & console){
        cat("\nAlpha:",alpha,"; DF Error:",DFerror)
        cat("\nt-Student:",Tprob)
        cat("\nLSD:", LSD,"\n")
      }
      if ( group & length(nr) != 1 & console) cat("\nGroups according to probability of treatment differences and alpha level(",alpha,")\n")
      if ( length(nr) != 1) statistics<-data.frame(Chisq=T1.aj,Df=ntr-1,p.chisq=p.value,F=T2.aj,DFerror=DFerror,p.F=PF)
      comb <-utils::combn(ntr,2)
      nn<-ncol(comb)
      dif<-rep(0,nn)
      pvalue<-rep(0,nn)
      LCL<-dif
      UCL<-dif
      sig<-NULL
      LSD<-rep(0,nn)
      stat<-rep("ns",nn)
      for (k in 1:nn) {
        i<-comb[1,k]
        j<-comb[2,k]
        dif[k]<-s[comb[1,k]]-s[comb[2,k]]
        sdtdif<- sqrt(2*(m[1]*A1-t(s)%*%s)/DFerror)
        pvalue[k]<- round(2*(1-pt(abs(dif[k])/sdtdif,DFerror)),4)
        LSD[k]<-round(Tprob*sdtdif,2)
        LCL[k] <- dif[k] - LSD[k]
        UCL[k] <- dif[k] + LSD[k]
        sig[k]<-" "
        if (pvalue[k] <= 0.001) sig[k]<-"***"
        else  if (pvalue[k] <= 0.01) sig[k]<-"**"
        else  if (pvalue[k] <= 0.05) sig[k]<-"*"
        else  if (pvalue[k] <= 0.1) sig[k]<-"."
      }
      if(!group){
        tr.i <- means[comb[1, ],1]
        tr.j <- means[comb[2, ],1]
        comparison<-data.frame("difference" = dif, pvalue=pvalue,"signif."=sig,LCL,UCL)
        rownames(comparison)<-paste(tr.i,tr.j,sep=" - ")
        if(console){cat("\nComparison between treatments\nSum of the ranks\n\n")
          print(comparison)}
        groups=NULL
      }
      if (group) {
        Q<-matrix(1,ncol=ntr,nrow=ntr)
        p<-pvalue
        k<-0
        for(i in 1:(ntr-1)){
          for(j in (i+1):ntr){
            k<-k+1
            Q[i,j]<-p[k]
            Q[j,i]<-p[k]
          }
        }
        groups <- ordenacao(means[, 1], means[, 2],alpha, Q,console)
        names(groups)[1]<-"Sum of ranks"
        if(console) {
          cat("\nTreatments with the same letter are not significantly different.\n\n")
          print(groups)
        }
        comparison<-NULL
      }
      parameters<-data.frame(test="Friedman",name.t=name.t,ntr = ntr,alpha=alpha)
      rownames(parameters)<-" "
      rownames(statistics)<-" "
      Means<-data.frame(rankSum=means[,2],Means)
      Means<-Means[,c(2,1,3:9)]
      output<-list(statistics=statistics,parameters=parameters,
                   means=Means,comparison=comparison,groups=groups)
      class(output)<-"group"
      invisible(output)
    }
    fried=friedman(bloco,trat,response,alpha=alpha.t)
    cat(green(bold("\n\n-----------------------------------------------------------------\n")))
    cat(green(italic("Statistics")))
    cat(green(bold("\n-----------------------------------------------------------------\n")))
    print(fried$statistics)
    cat(green(bold("\n\n-----------------------------------------------------------------\n")))
    cat(green(italic("Parameters")))
    cat(green(bold("\n-----------------------------------------------------------------\n")))
    print(fried$parameters)
    cat(green(bold("\n\n-----------------------------------------------------------------\n")))
    cat(green(italic("Multiple Comparison Test")))
    cat(green(bold("\n-----------------------------------------------------------------\n")))
    saida=cbind(fried$means[,c(1,3)],fried$groups[rownames(fried$means),])
    colnames(saida)=c("Mean","SD","Rank","Groups")
    print(saida)
    dadosm=data.frame(fried$means,fried$groups[rownames(fried$means),])
    dadosm$trats=factor(rownames(dadosm),levels = unique(trat))
    dadosm$media=tapply(response,trat,mean, na.rm=TRUE)[rownames(fried$means)]

    if(point=="mean_sd"){dadosm$std=tapply(response,trat,sd, na.rm=TRUE)[rownames(fried$means)]}
    if(point=="mean_se"){dadosm$std=tapply(response,trat,sd, na.rm=TRUE)/
      sqrt(tapply(response,trat,length))[rownames(fried$means)]}

    dadosm$limite=dadosm$response+dadosm$std
    dadosm$letra=paste(format(dadosm$response,digits = dec),dadosm$groups)
    if(addmean==TRUE){dadosm$letra=paste(format(dadosm$response,digits = dec),dadosm$groups)}
    if(addmean==FALSE){dadosm$letra=dadosm$groups}
    dadosm=dadosm[unique(trat),]
    trats=dadosm$trats
    limite=dadosm$limite
    media=dadosm$media
    std=dadosm$std
    letra=dadosm$letra
    if(geom=="bar"){grafico=ggplot(dadosm,aes(x=trats,
                                              y=response))
      if(fill=="trat"){grafico=grafico+
        geom_col(aes(fill=trats),color=1)}
    else{grafico=grafico+
      geom_col(aes(fill=trats),fill=fill,color=1)}
    if(errorbar==TRUE){grafico=grafico+
      geom_text(aes(y=media+sup+if(sup<0){-std}else{std},label=letra),
                family=family,angle=angle.label, hjust=hjust)}
    if(errorbar==FALSE){grafico=grafico+
      geom_text(aes(y=media+sup,label=letra),family=family,angle=angle.label, hjust=hjust)}
    if(errorbar==TRUE){grafico=grafico+
      geom_errorbar(aes(ymin=response-std,
                        ymax=response+std),
                    color="black", width=0.3)}}
    if(geom=="point"){grafico=ggplot(dadosm,
                                     aes(x=trats,
                                         y=response))
    if(errorbar==TRUE){grafico=grafico+
      geom_text(aes(y=media+sup+if(sup<0){-std}else{std},label=letra),
                family=family,angle=angle.label, hjust=hjust)}
    if(errorbar==FALSE){grafico=grafico+
      geom_text(aes(y=media+sup,label=letra),family=family,angle=angle.label, hjust=hjust)}
    if(errorbar==TRUE){grafico=grafico+
      geom_errorbar(aes(ymin=response-std,
                        ymax=response+std),
                    color="black", width=0.3)}
      if(fill=="trat"){grafico=grafico+
        geom_point(aes(color=trats),size=5)}
    else{grafico=grafico+
      geom_point(aes(color=trats),fill="gray",pch=21,color="black",size=5)}}
    if(geom=="box"){
    datam1=data.frame(trats=factor(trat,levels = unique(as.character(trat))),response)
    dadosm2=data.frame(fried$means)
    dadosm2$trats=factor(rownames(dadosm),levels = unique(trat))
    dadosm2$limite=dadosm2$response+dadosm2$std
    dadosm2$letra=paste(format(dadosm$response,digits = dec),
                        dadosm$groups)
    dadosm2=dadosm2[unique(as.character(trat)),]
    trats=dadosm2$trats
    limite=dadosm2$limite
    letra=dadosm2$letra
    stat_box=ggplot(datam1,aes(x=trats,y=response))+geom_boxplot()
    superior=ggplot_build(stat_box)$data[[1]]$ymax
    dadosm2$superior=superior+sup

    grafico=ggplot(datam1,aes(x=trats,
                              y=response))
      if(fill=="trat"){grafico=grafico+
        geom_boxplot(aes(fill=trats))}
    else{grafico=grafico+
      geom_boxplot(aes(fill=trats),fill=fill)}
    grafico=grafico+
      geom_text(data=dadosm2,
                aes(y=superior,
                    label=letra),family=family,angle=angle.label, hjust=hjust)}

    grafico=grafico+theme+
      ylab(ylab)+
      xlab(xlab)+
      theme(text = element_text(size=textsize,color="black",family=family),
            axis.text = element_text(size=textsize,color="black",family=family),
            axis.title = element_text(size=textsize,color="black",family=family),
            legend.position = "none")
    if(angle !=0){grafico=grafico+theme(axis.text.x=element_text(hjust = 1.01,
                                                                 angle = angle))}
  }
  if(quali==TRUE){print(grafico)}
  grafico=list(grafico)
  }
