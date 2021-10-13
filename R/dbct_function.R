#' Analysis: Randomized block design evaluated over time
#'
#' @description Function of the AgroR package for analysis of experiments conducted in a balanced qualitative, single-factorial randomized block design with multiple assessments over time, however without considering time as a factor.
#' @author Gabriel Danilo Shimizu, \email{shimizu@uel.br}
#' @author Leandro Simoes Azeredo Gonçalves
#' @author Rodrigo Yudi Palhaci Marubayashi
#' @param trat Numerical or complex vector with treatments
#' @param block Numerical or complex vector with blocks
#' @param time Numerical or complex vector with times
#' @param response Numerical vector containing the response of the experiment.
#' @param alpha.f Level of significance of the F test (\emph{default} is 0.05)
#' @param alpha.t Significance level of the multiple comparison test (\emph{default} is 0.05)
#' @param mcomp Multiple comparison test (Tukey (\emph{default}), LSD ("lsd"), Scott-Knott ("sk"), Duncan ("duncan") and Friedman ("fd"))
#' @param ylab Variable response name (Accepts the \emph{expression}() function)
#' @param xlab Treatments name (Accepts the \emph{expression}() function)
#' @param fill Defines chart color (to generate different colors for different treatments, define fill = "trat")
#' @param theme ggplot2 theme (\emph{default} is theme_classic())
#' @param error Add error bar
#' @param sup Number of units above the standard deviation or average bar on the graph
#' @param addmean Plot the average value on the graph (\emph{default} is TRUE)
#' @param textsize Font size of the texts and titles of the axes
#' @param labelsize Font size of the labels
#' @param family Font family
#' @param dec Number of cells
#' @param geom Graph type (columns - "bar" or segments "point")
#' @param legend Legend title
#' @param posi Legend position
#' @param ylim y-axis scale
#' @param width.bar width errorbar
#' @param xnumeric Declare x as numeric (\emph{default} is FALSE)
#' @param all.letters Adds all label letters regardless of whether it is significant or not.
#' @note The ordering of the graph is according to the sequence in which the factor levels are arranged in the data sheet. The bars of the column and segment graphs are standard deviation.
#' @keywords dbct
#' @keywords Experimental
#' @seealso \link{DBC}, \link{DICT}, \link{DQLT}
#' @references
#' @references
#'
#' Principles and procedures of statistics a biometrical approach Steel & Torry & Dickey. Third Edition 1997
#'
#' Multiple comparisons theory and methods. Departament of statistics the Ohio State University. USA, 1996. Jason C. Hsu. Chapman Hall/CRC.
#'
#' Practical Nonparametrics Statistics. W.J. Conover, 1999
#'
#' Ramalho M.A.P., Ferreira D.F., Oliveira A.C. 2000. Experimentacao em Genetica e Melhoramento de Plantas. Editora UFLA.
#'
#' Scott R.J., Knott M. 1974. A cluster analysis method for grouping mans in the analysis of variance. Biometrics, 30, 507-512.
#' @details The p-value of the analysis of variance, the normality test for Shapiro-Wilk errors, the Bartlett homogeneity test of variances, the independence of Durbin-Watson errors and the multiple comparison test (Tukey, Scott-Knott, LSD or Duncan).
#' @export
#' @return The function returns the p-value of Anova, the assumptions of normality of errors, homogeneity of variances and independence of errors, multiple comparison test, as well as a line graph
#' @examples
#' rm(list=ls())
#' data(simulate2)
#' attach(simulate2)
#'
#' #===================================
#' # default
#' #===================================
#' DBCT(trat, bloco, tempo, resp)
#'
#' #===================================
#' # segment chart
#' #===================================
#' DBCT(trat, bloco, tempo, resp, geom="point")

DBCT=function(trat,
              block,
              time,
              response,
              alpha.f=0.05,
              alpha.t=0.05,
              theme=theme_classic(),
              geom="bar",
              fill="gray",
              ylab="Response",
              xlab="Independent",
              mcomp="tukey",
              textsize=12,
              labelsize=5,
              error=TRUE,
              family="sans",
              sup=0,
              addmean=FALSE,
              posi=c(0.1,0.8),
              legend="Legend",
              ylim=NA,
              width.bar=0.1,
              dec=3,
              xnumeric=FALSE,
              all.letters=FALSE){
  requireNamespace("ScottKnott")
  requireNamespace("crayon")
  requireNamespace("ggplot2")
  requireNamespace("nortest")
  requireNamespace("ggrepel")
  trat=as.factor(trat)
  resp=response
  block=as.factor(block)
  time=factor(time,unique(time))
  dados=data.frame(resp,trat,block,time)
  if(mcomp=="tukey"){
    tukeyg=c()
    ordem=c()
    normg=c()
    homog=c()
    indepg=c()
    anovag=c()
    cv=c()
    for(i in 1:length(levels(time))){
      mod=aov(resp~trat+block, data=dados[dados$time==levels(dados$time)[i],])
      anovag[[i]]=anova(mod)$`Pr(>F)`[1]
      cv[[i]]=sqrt(anova(mod)$`Mean Sq`[3])/mean(mod$model$resp)*100
      tukey=TUKEY(mod,"trat",alpha = alpha.t)
      tukey$groups=tukey$groups[unique(as.character(trat)),2]
      if(all.letters==FALSE){
      if(anova(mod)$`Pr(>F)`[1]>alpha.f){tukey$groups=c("ns",rep(" ",length(unique(trat))-1))}}
      tukeyg[[i]]=as.character(tukey$groups)
      ordem[[i]]=rownames(tukey$groups)
      norm=shapiro.test(mod$residuals)
      homo=with(dados[dados$time==levels(dados$time)[i],], bartlett.test(mod$residuals~trat))
      indep=dwtest(mod)
      normg[[i]]=norm$p.value
      homog[[i]]=homo$p.value
      indepg[[i]]=indep$p.value
    }
    m=unlist(tukeyg)
    nor=unlist(normg)
    hom=unlist(homog)
    ind=unlist(indepg)
    an=unlist(anovag)
    cv=unlist(cv)
    press=data.frame(an,nor,hom,ind,cv)
    colnames(press)=c("p-value ANOVA","Shapiro-Wilk","Bartlett","Durbin-Watson","CV (%)")}

  if(mcomp=="lsd"){
    lsdg=c()
    ordem=c()
    normg=c()
    homog=c()
    indepg=c()
    anovag=c()
    cv=c()
    for(i in 1:length(levels(time))){
      mod=aov(resp~trat+block, data=dados[dados$time==levels(dados$time)[i],])
      anovag[[i]]=anova(mod)$`Pr(>F)`[1]
      cv[[i]]=sqrt(anova(mod)$`Mean Sq`[3])/mean(mod$model$resp)*100
      lsd=LSD(mod,"trat",alpha = alpha.t)
      lsd$groups=lsd$groups[unique(as.character(trat)),2]
      if(all.letters==FALSE){
      if(anova(mod)$`Pr(>F)`[1]>alpha.f){lsd$groups=c("ns",rep(" ",length(unique(trat))-1))}}
      lsdg[[i]]=as.character(lsd$groups)
      ordem[[i]]=rownames(lsd$groups)
      norm=shapiro.test(mod$residuals)
      homo=with(dados[dados$time==levels(dados$time)[i],], bartlett.test(mod$residuals~trat))
      indep=dwtest(mod)
      normg[[i]]=norm$p.value
      homog[[i]]=homo$p.value
      indepg[[i]]=indep$p.value
    }
    m=unlist(lsdg)
    nor=unlist(normg)
    hom=unlist(homog)
    ind=unlist(indepg)
    an=unlist(anovag)
    cv=unlist(cv)
    press=data.frame(an,nor,hom,ind,cv)
    colnames(press)=c("p-value ANOVA","Shapiro-Wilk","Bartlett","Durbin-Watson","CV (%)")}

  if(mcomp=="duncan"){
    duncang=c()
    ordem=c()
    normg=c()
    homog=c()
    indepg=c()
    anovag=c()
    cv=c()
    for(i in 1:length(levels(time))){
      mod=aov(resp~trat+block, data=dados[dados$time==levels(dados$time)[i],])
      anovag[[i]]=anova(mod)$`Pr(>F)`[1]
      cv[[i]]=sqrt(anova(mod)$`Mean Sq`[3])/mean(mod$model$resp)*100
      duncan=duncan(mod,"trat",alpha = alpha.t)
      duncan$groups=duncan$groups[unique(as.character(trat)),2]
      if(all.letters==FALSE){
      if(anova(mod)$`Pr(>F)`[1]>alpha.f){duncan$groups=c("ns",rep(" ",length(unique(trat))-1))}}
      duncang[[i]]=as.character(duncan$groups)
      ordem[[i]]=rownames(duncan$groups)
      norm=shapiro.test(mod$residuals)
      homo=with(dados[dados$time==levels(dados$time)[i],], bartlett.test(mod$residuals~trat))
      indep=dwtest(mod)
      normg[[i]]=norm$p.value
      homog[[i]]=homo$p.value
      indepg[[i]]=indep$p.value
    }
    m=unlist(duncang)
    nor=unlist(normg)
    hom=unlist(homog)
    ind=unlist(indepg)
    an=unlist(anovag)
    cv=unlist(cv)
    press=data.frame(an,nor,hom,ind,cv)
    colnames(press)=c("p-value ANOVA","Shapiro-Wilk","Bartlett","Durbin-Watson","CV (%)")}


  if(mcomp=="sk"){
    scott=c()
    normg=c()
    homog=c()
    indepg=c()
    anovag=c()
    cv=c()
    for(i in 1:length(levels(time))){
      mod=aov(resp~trat+block, data=dados[dados$time==levels(dados$time)[i],])
      anovag[[i]]=anova(mod)$`Pr(>F)`[1]
      cv[[i]]=sqrt(anova(mod)$`Mean Sq`[3])/mean(mod$model$resp)*100
      letra=SK(mod,"trat",sig.level = alpha.t)
      data=data.frame(sk=letters[letra$groups])
      rownames(data)=rownames(letra$m.inf)
      data=data[unique(as.character(trat)),]
      if(all.letters==FALSE){
      if(anova(mod)$`Pr(>F)`[1]>alpha.f){data=c("ns",rep(" ",length(unique(trat))-1))}}
      scott[[i]]=as.character(data)
      norm=shapiro.test(mod$residuals)
      homo=with(dados[dados$time==levels(dados$time)[i],], bartlett.test(mod$residuals~trat))
      indep=dwtest(mod)
      normg[[i]]=norm$p.value
      homog[[i]]=homo$p.value
      indepg[[i]]=indep$p.value
    }
    m=unlist(scott)
    nor=unlist(normg)
    hom=unlist(homog)
    ind=unlist(indepg)
    an=unlist(anovag)
    cv=unlist(cv)
    press=data.frame(an,nor,hom,ind,cv)
    colnames(press)=c("p-value ANOVA","Shapiro-Wilk","Bartlett","Durbin-Watson","CV (%)")}

  if(mcomp=="fd"){
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
    fdg=c()
    ordem=c()
    anovag=c()
    for(i in 1:length(levels(time))){
      fd=friedman(block,trat,resp,alpha=alpha.t)
      anovag[[i]]=mod$statistics$p.chisq
      fd$groups=fd$groups[unique(as.character(trat)),2]
      if(all.letters==TRUE){
        if(anova(mod)$`Pr(>F)`[1]>alpha.f){fd$groups=c("ns",rep(" ",length(unique(trat))-1))}}
      fdg[[i]]=as.character(duncan$groups)
      ordem[[i]]=rownames(fd$groups)
    }
    m=unlist(fdg)
    an=unlist(anovag)
    press=data.frame(fd);colnames(press)=c("p-value Chisq Friedman")}


  cat(green(bold("\n-----------------------------------------------------------------\n")))
  cat(green(bold("ANOVA and assumptions")))
  cat(green(bold("\n-----------------------------------------------------------------\n")))
  print(press)

  dadosm=data.frame(time=as.character(rep(unique(time),e=length(unique(as.character(trat))))),
                    trat=rep(unique(as.character(trat)),length(unique(time))),
                    media=c(tapply(resp,list(trat,time),mean, na.rm=TRUE)[unique(as.character(trat)),]),
                    desvio=c(tapply(resp,list(trat,time),sd, na.rm=TRUE)[unique(as.character(trat)),]),
                    letra=m)
  if(xnumeric==TRUE){dadosm$time=as.numeric(as.character(dadosm$time))}
  time=dadosm$time
  trat=dadosm$trat
  media=dadosm$media
  desvio=dadosm$desvio
  letra=dadosm$letra

  if(geom=="point"){
    grafico=ggplot(dadosm,aes(y=media,
                              x=time))+
      geom_point(aes(shape=factor(trat, levels=unique(as.character(trat))),
                     group=factor(trat, levels=unique(as.character(trat)))),size=3)+
      geom_line(aes(lty=factor(trat, levels=unique(as.character(trat))),
                    group=factor(trat, levels=unique(as.character(trat)))),size=0.8)+
      ylab(ylab)+
      xlab(xlab)+theme+
      theme(text = element_text(size=textsize,color="black", family = family),
            axis.title = element_text(size=textsize,color="black", family = family),
            axis.text = element_text(size=textsize,color="black", family = family),
            legend.position = posi,
            legend.text = element_text(size = textsize))+labs(shape=legend, lty=legend)
    if(error==TRUE){grafico=grafico+
      geom_errorbar(aes(ymin=media-desvio,
                        ymax=media+desvio), width=width.bar)}
    if(addmean==FALSE && error==FALSE){grafico=grafico+
      geom_text_repel(aes(y=media+sup,label=letra),family=family,size=labelsize)}
    if(addmean==TRUE && error==FALSE){grafico=grafico+
      geom_text_repel(aes(y=media+sup,
                          label=paste(format(media,digits = dec),
                                      letra)),family=family,size=labelsize)}
    if(addmean==FALSE && error==TRUE){grafico=grafico+
      geom_text_repel(aes(y=desvio+media+sup,
                          label=letra),family=family,size=labelsize)}
    if(addmean==TRUE && error==TRUE){grafico=grafico+
      geom_text_repel(aes(y=desvio+media+sup,
                          label=paste(format(media,digits = dec),
                                      letra)),family=family,size=labelsize)}
  }

  if(geom=="bar"){
    if(sup==0){sup=0.1*mean(dadosm$media)}
    grafico=ggplot(dadosm,aes(y=media,
                              x=as.factor(time),
                              fill=factor(trat,levels = unique(trat))))+
      geom_col(position = "dodge",color="black")+
      ylab(ylab)+
      xlab(xlab)+theme+
      theme(text = element_text(size=textsize,color="black", family = family),
            axis.title = element_text(size=textsize,color="black", family = family),
            axis.text = element_text(size=textsize,color="black", family = family),
            legend.position = posi,
            legend.text = element_text(size = textsize))+labs(fill=legend)
    if(error==TRUE){grafico=grafico+
      geom_errorbar(aes(ymin=media-desvio,
                        ymax=media+desvio),
                    width=width.bar, position = position_dodge(width=0.9))}
    if(addmean==FALSE && error==FALSE){grafico=grafico+
      geom_text(aes(y=media+sup,label=letra),
                size=labelsize,family=family,
                position = position_dodge(width=0.9))}
    if(addmean==TRUE && error==FALSE){grafico=grafico+
      geom_text(aes(y=media+sup,
                    label=paste(format(media,digits = dec),letra)),
                size=labelsize,family=family,
                position = position_dodge(width=0.9))}
    if(addmean==FALSE && error==TRUE){grafico=grafico+
      geom_text(aes(y=desvio+media+sup,label=letra),
                size=labelsize,family=family,
                position = position_dodge(width=0.9))}
    if(addmean==TRUE && error==TRUE){grafico=grafico+
      geom_text(aes(y=desvio+media+sup,
                    label=paste(format(media,digits = dec),letra)),
                size=labelsize,family=family,
                position = position_dodge(width=0.9))}
  }
  if(fill=="gray"){grafico=grafico+scale_fill_grey(start = 1, end = 0.1)}
  if(is.na(ylim)==FALSE){grafico=grafico+scale_y_continuous(breaks = ylim)}
  graficos=as.list(grafico)
  print(grafico)
}
