#' Analysis: Completely randomized design evaluated over time
#'
#' @description Function of the AgroR package for the analysis of experiments conducted in a completely randomized, qualitative, uniform qualitative design with multiple assessments over time, however without considering time as a factor.
#' @author Gabriel Danilo Shimizu, \email{shimizu@uel.br}
#' @author Leandro Simoes Azeredo Goncalves
#' @author Rodrigo Yudi Palhaci Marubayashi
#' @param trat Numerical or complex vector with treatments
#' @param time Numerical or complex vector with times
#' @param response Numerical vector containing the response of the experiment.
#' @param alpha.f Level of significance of the F test (\emph{default} is 0.05)
#' @param alpha.t Significance level of the multiple comparison test (\emph{default} is 0.05)
#' @param mcomp Multiple comparison test (Tukey (\emph{default}), LSD ("lsd"), Scott-Knott ("sk"), Duncan ("duncan") and Kruskal-Wallis ("kw"))
#' @param ylab Variable response name (Accepts the \emph{expression}() function)
#' @param xlab treatments name (Accepts the \emph{expression}() function)
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
#' @param p.adj Method for adjusting p values for Kruskal-Wallis ("none","holm","hommel", "hochberg", "bonferroni", "BH", "BY", "fdr")
#' @param all.letters Adds all label letters regardless of whether it is significant or not.
#' @note The ordering of the graph is according to the sequence in which the factor levels are arranged in the data sheet. The bars of the column and segment graphs are standard deviation.
#' @keywords dict
#' @keywords Experimental
#' @seealso \link{DIC}, \link{DBCT}, \link{DQLT}
#' @return The function returns the p-value of Anova, the assumptions of normality of errors, homogeneity of variances and independence of errors, multiple comparison test, as well as a line graph
#' @references
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
#' @export
#' @examples
#' rm(list=ls())
#' data(simulate1)
#' attach(simulate1)
#' with(simulate1, DICT(trat, tempo, resp))
#' with(simulate1, DICT(trat, tempo, resp,geom="bar",sup=40))

DICT=function(trat,
              time,
              response,
              alpha.f=0.05,
              alpha.t=0.05,
              mcomp="tukey",
              theme=theme_classic(),
              geom="bar",
              xlab="Independent",
              ylab="Response",
              p.adj="holm",
              dec=3,
              fill="gray",
              error=TRUE,
              textsize=12,
              labelsize=5,
              family="sans",
              sup=0,
              addmean=FALSE,
              legend="Legend",
              ylim=NA,
              width.bar=0.2,
              posi=c(0.1,0.8),
              xnumeric=FALSE,
              all.letters=FALSE){
  requireNamespace("crayon")
  requireNamespace("ggplot2")
  requireNamespace("nortest")
  requireNamespace("ggrepel")
resp=response
trat=as.factor(trat)
time=factor(time,unique(time))
dados=data.frame(resp,trat,time)

if(mcomp=="tukey"){
  tukeyg=c()
  ordem=c()
  normg=c()
  homog=c()
  indepg=c()
  anovag=c()
  cv=c()
  for(i in 1:length(levels(time))){
  mod=aov(resp~trat, data=dados[dados$time==levels(dados$time)[i],])
  anovag[[i]]=anova(mod)$`Pr(>F)`[1]
  cv[[i]]=sqrt(anova(mod)$`Mean Sq`[2])/mean(mod$model$resp)*100
  norm=shapiro.test(mod$residuals)
  homo=with(dados[dados$time==levels(dados$time)[i],], bartlett.test(mod$residuals~trat))
  indep=dwtest(mod)
  tukey=TUKEY(mod,"trat",alpha = alpha.t)
  tukey$groups=tukey$groups[unique(as.character(trat)),2]

  if(all.letters==FALSE){
    if(anova(mod)$`Pr(>F)`[1]>alpha.f){tukey$groups=c("ns",rep(" ",length(unique(trat))-1))}}
  tukeyg[[i]]=as.character(tukey$groups)
  normg[[i]]=norm$p.value
  homog[[i]]=homo$p.value
  indepg[[i]]=indep$p.value
  ordem[[i]]=rownames(tukey$groups)
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
    mod=aov(resp~trat, data=dados[dados$time==levels(dados$time)[i],])
    anovag[[i]]=anova(mod)$`Pr(>F)`[1]
    cv[[i]]=sqrt(anova(mod)$`Mean Sq`[2])/mean(mod$model$resp)*100
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
    mod=aov(resp~trat, data=dados[dados$time==levels(dados$time)[i],])
    anovag[[i]]=anova(mod)$`Pr(>F)`[1]
    cv[[i]]=sqrt(anova(mod)$`Mean Sq`[2])/mean(mod$model$resp)*100
    duncan=duncan(mod,"trat",alpha = alpha.t)
    duncan$groups=duncan$groups[unique(as.character(trat)),2]
    if(all.letters==FALSE){
      if(anova(mod)$`Pr(>F)`[1]>alpha.f){duncan$groups=c("ns",rep(" ",length(unique(trat))-1))}}
    norm=shapiro.test(mod$residuals)
    homo=with(dados[dados$time==levels(dados$time)[i],], bartlett.test(mod$residuals~trat))
    indep=dwtest(mod)
    normg[[i]]=norm$p.value
    homog[[i]]=homo$p.value
    indepg[[i]]=indep$p.value
    duncang[[i]]=as.character(duncan$groups)
    ordem[[i]]=rownames(duncan$groups)
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
  mod=aov(resp~trat, data=dados[dados$time==levels(dados$time)[i],])
  anovag[[i]]=anova(mod)$`Pr(>F)`[1]
  cv[[i]]=sqrt(anova(mod)$`Mean Sq`[2])/mean(mod$model$resp)*100
  norm=shapiro.test(mod$residuals)
  homo=with(dados[dados$time==levels(dados$time)[i],], bartlett.test(mod$residuals~trat))
  indep=dwtest(mod)
  nrep=with(dados[dados$time==levels(dados$time)[i],], table(trat)[1])
  ao=anova(mod)
  medias=sort(tapply(resp,trat,mean),decreasing = TRUE)
  letra=scottknott(means = medias,
                   df1 = ao$Df[2],
                   nrep = nrep,
                   QME = ao$`Mean Sq`[2],
                   alpha = alpha.t)
  letra1=data.frame(resp=medias,groups=letra)
  letra1=letra1[unique(as.character(trat)),]
  data=letra1$groups
  if(all.letters==FALSE){
    if(anova(mod)$`Pr(>F)`[1]>alpha.f){data=c("ns",rep(" ",length(unique(trat))-1))}}
  data=data
  scott[[i]]=data
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
if(mcomp=="kw"){
  kruskal=function (y,
                    trt,
                    alpha = 0.05,
                    p.adj = c("none", "holm",
                              "hommel", "hochberg", "bonferroni", "BH", "BY", "fdr"),
                    group = TRUE, main = NULL,console=FALSE){
    name.y <- paste(deparse(substitute(y)))
    name.t <- paste(deparse(substitute(trt)))
    if(is.null(main))main<-paste(name.y,"~", name.t)
    p.adj <- match.arg(p.adj)
    junto <- subset(data.frame(y, trt), is.na(y) == FALSE)
    N <- nrow(junto)
    medians<-mean.stat(junto[,1],junto[,2],stat="median")
    for(i in c(1,5,2:4)) {
      x <- mean.stat(junto[,1],junto[,2],function(x)quantile(x)[i])
      medians<-cbind(medians,x[,2])}
    medians<-medians[,3:7]
    names(medians)<-c("Min","Max","Q25","Q50","Q75")
    Means <- mean.stat(junto[,1],junto[,2],stat="mean")
    sds <-   mean.stat(junto[,1],junto[,2], stat="sd")
    nn <-   mean.stat(junto[,1],junto[,2],stat="length")
    Means<-data.frame(Means,std=sds[,2],r=nn[,2],medians)
    rownames(Means)<-Means[,1]
    Means<-Means[,-1]
    names(Means)[1]<-name.y
    junto[, 1] <- rank(junto[, 1])
    means <- mean.stat(junto[, 1], junto[, 2], stat = "sum")
    sds <- mean.stat(junto[, 1], junto[, 2], stat = "sd")
    nn <- mean.stat(junto[, 1], junto[, 2], stat = "length")
    means <- data.frame(means, r = nn[, 2])
    names(means)[1:2] <- c(name.t, name.y)
    ntr <- nrow(means)
    nk <- choose(ntr, 2)
    DFerror <- N - ntr
    rs <- 0
    U <- 0
    for (i in 1:ntr) {
      rs <- rs + means[i, 2]^2/means[i, 3]
      U <- U + 1/means[i, 3]
    }
    S <- (sum(junto[, 1]^2) - (N * (N + 1)^2)/4)/(N - 1)
    H <- (rs - (N * (N + 1)^2)/4)/S
    p.chisq <- 1 - pchisq(H, ntr - 1)
    if(console){
      cat("\nStudy:", main)
      cat("\nKruskal-Wallis test's\nTies or no Ties\n")
      cat("\nCritical Value:", H)
      cat("\nDegrees of freedom:", ntr - 1)
      cat("\nPvalue Chisq  :", p.chisq, "\n\n")}
    DFerror <- N - ntr
    Tprob <- qt(1 - alpha/2, DFerror)
    MSerror <- S * ((N - 1 - H)/(N - ntr))
    means[, 2] <- means[, 2]/means[, 3]
    if(console){cat(paste(name.t, ",", sep = ""), " means of the ranks\n\n")
      print(data.frame(row.names = means[, 1], means[, -1]))
      cat("\nPost Hoc Analysis\n")}
    if (p.adj != "none") {
      if(console)cat("\nP value adjustment method:", p.adj)
      a <- 1e-06
      b <- 1
      for (i in 1:100) {
        x <- (b + a)/2
        xr <- rep(x, nk)
        d <- p.adjust(xr, p.adj)[1] - alpha
        ar <- rep(a, nk)
        fa <- p.adjust(ar, p.adj)[1] - alpha
        if (d * fa < 0)
          b <- x
        if (d * fa > 0)
          a <- x}
      Tprob <- qt(1 - x/2, DFerror)
    }
    nr <- unique(means[, 3])
    if (group & console){
      cat("\nt-Student:", Tprob)
      cat("\nAlpha    :", alpha)}
    if (length(nr) == 1)  LSD <- Tprob * sqrt(2 * MSerror/nr)
    statistics<-data.frame(Chisq=H,Df=ntr-1,p.chisq=p.chisq)
    if ( group & length(nr) == 1 & console) cat("\nMinimum Significant Difference:",LSD,"\n")
    if ( group & length(nr) != 1 & console) cat("\nGroups according to probability of treatment differences and alpha level.\n")
    if ( length(nr) == 1) statistics<-data.frame(statistics,t.value=Tprob,MSD=LSD)
    comb <- utils::combn(ntr, 2)
    nn <- ncol(comb)
    dif <- rep(0, nn)
    LCL <- dif
    UCL <- dif
    pvalue <- dif
    sdtdif <- dif
    for (k in 1:nn) {
      i <- comb[1, k]
      j <- comb[2, k]
      dif[k] <- means[i, 2] - means[j, 2]
      sdtdif[k] <- sqrt(MSerror * (1/means[i,3] + 1/means[j, 3]))
      pvalue[k] <- 2*(1 - pt(abs(dif[k])/sdtdif[k],DFerror))
    }
    if (p.adj != "none") pvalue <- p.adjust(pvalue, p.adj)
    pvalue <- round(pvalue,4)
    sig <- rep(" ", nn)
    for (k in 1:nn) {
      if (pvalue[k] <= 0.001)
        sig[k] <- "***"
      else if (pvalue[k] <= 0.01)
        sig[k] <- "**"
      else if (pvalue[k] <= 0.05)
        sig[k] <- "*"
      else if (pvalue[k] <= 0.1)
        sig[k] <- "."
    }
    tr.i <- means[comb[1, ], 1]
    tr.j <- means[comb[2, ], 1]
    LCL <- dif - Tprob * sdtdif
    UCL <- dif + Tprob * sdtdif
    comparison <- data.frame(Difference = dif, pvalue = pvalue, "Signif."=sig, LCL, UCL)
    if (p.adj !="bonferroni" & p.adj !="none"){
      comparison<-comparison[,1:3]
      statistics<-data.frame(Chisq=H,p.chisq=p.chisq)}
    rownames(comparison) <- paste(tr.i, tr.j, sep = " - ")
    if (!group) {
      groups<-NULL
      if(console){
        cat("\nComparison between treatments mean of the ranks.\n\n")
        print(comparison)
      }
    }
    if (group) {
      comparison=NULL
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
      groups <- ordenacao(means[, 1], means[, 2],alpha, Q, console)
      names(groups)[1]<-name.y
      if(console) {
        cat("\nTreatments with the same letter are not significantly different.\n\n")
        print(groups)
      }
    }
    ranks=means
    Means<-data.frame(rank=ranks[,2],Means)
    Means<-Means[,c(2,1,3:9)]
    parameters<-data.frame(test="Kruskal-Wallis",p.ajusted=p.adj,name.t=name.t,ntr = ntr,alpha=alpha)
    rownames(parameters)<-" "
    rownames(statistics)<-" "
    output<-list(statistics=statistics,parameters=parameters,
                 means=Means,comparison=comparison,groups=groups)
    class(output)<-"group"
    invisible(output)
  }
  kwg=c()
  ordem=c()
  normg=c()
  homog=c()
  indepg=c()
  anovag=c()
  for(i in 1:length(levels(time))){
    data=dados[dados$time==levels(dados$time)[i],]
    mod=with(data,kruskal(resp,trat,p.adj = p.adj,alpha = alpha.t))
    anovag[[i]]=mod$statistics$p.chisq
    norm=""
    homo=""
    indep=""
    kw=mod
    kw$groups=kw$groups[unique(as.character(trat)),2]
    if(all.letters==FALSE){
      if(mod$statistics$p.chisq>alpha.f){kw$groups=c("ns",
                                                   rep(" ",length(unique(trat))-1))}}
    normg[[i]]=norm
    homog[[i]]=homo
    indepg[[i]]=indep
    kwg[[i]]=as.character(kw$groups)
    ordem[[i]]=rownames(kw$groups)
  }
  m=unlist(kwg)
  nor=unlist(normg)
  hom=unlist(homog)
  ind=unlist(indepg)
  an=unlist(anovag)
  press=data.frame(an)
  colnames(press)=c("p-value Kruskal")}


cat(green(bold("\n-----------------------------------------------------------------\n")))
cat(green(bold("ANOVA and assumptions")))
cat(green(bold("\n-----------------------------------------------------------------\n")))
print(press)

dadosm=data.frame(#time=as.numeric(as.character(rep(unique(time),e=length(unique(as.character(trat)))))),
                  time=as.character(rep(unique(time),e=length(unique(as.character(trat))))),
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
  geom_point(aes(shape=factor(trat,
                              levels=unique(as.character(trat))),
                 group=factor(trat,levels=unique(as.character(trat)))),size=3)+
  geom_line(aes(lty=factor(trat,levels=unique(as.character(trat))),
                group=factor(trat,levels=unique(as.character(trat)))),size=0.8)+
  ylab(ylab)+
  xlab(xlab)+theme+
  theme(text = element_text(size=textsize,color="black", family = family),
        axis.title = element_text(size=textsize,color="black", family = family),
        axis.text = element_text(size=textsize,color="black", family = family),
        legend.position = posi,
        legend.text = element_text(size = textsize))+
  labs(shape=legend, lty=legend)
if(error==TRUE){grafico=grafico+
  geom_errorbar(aes(ymin=media-desvio,
                    ymax=media+desvio),
                width=width.bar)}
if(addmean==FALSE && error==FALSE){grafico=grafico+
  geom_text_repel(aes(y=media+sup,label=letra),family=family,size=labelsize)}
if(addmean==TRUE && error==FALSE){grafico=grafico+
  geom_text_repel(aes(y=media+sup,
                      label=paste(format(media,digits = dec),letra)),family=family,size=labelsize)}
if(addmean==FALSE && error==TRUE){grafico=grafico+
  geom_text_repel(aes(y=desvio+media+sup,
                      label=letra),family=family,size=labelsize)}
if(addmean==TRUE && error==TRUE){grafico=grafico+
  geom_text_repel(aes(y=desvio+media+sup,
                      label=paste(format(media,digits = dec),letra)),family=family,size=labelsize)}
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
                width=width.bar,
                position = position_dodge(width=0.9))}
if(addmean==FALSE && error==FALSE){grafico=grafico+
  geom_text(aes(y=media+sup,label=letra),
            position = position_dodge(width=0.9),family=family,size=labelsize)}
if(addmean==TRUE && error==FALSE){grafico=grafico+
  geom_text(aes(y=media+sup,
                label=paste(format(media,digits = dec),letra)),
            position = position_dodge(width=0.9),family=family,size=labelsize)}
if(addmean==FALSE && error==TRUE){grafico=grafico+
  geom_text(aes(y=desvio+media+sup,label=letra),
            position = position_dodge(width=0.9),family=family,size=labelsize)}
if(addmean==TRUE && error==TRUE){grafico=grafico+
  geom_text(aes(y=desvio+media+sup,
                label=paste(format(media,digits = dec),letra)),
            position = position_dodge(width=0.9),family=family,size=labelsize)}
}
if(fill=="gray"){grafico=grafico+scale_fill_grey(start = 1, end = 0.1)}
if(is.na(ylim)==FALSE){grafico=grafico+scale_y_continuous(breaks = ylim)}
graficos=as.list(grafico)
print(grafico)
}

