#' Analysis: Completely randomized design
#'
#' @description Statistical analysis of experiments conducted in a completely randomized and balanced design with a factor considering the fixed model.
#' @author Gabriel Danilo Shimizu, \email{shimizu@uel.br}
#' @author Leandro Simoes Azeredo Goncalves
#' @author Rodrigo Yudi Palhaci Marubayashi
#' @param trat Numerical or complex vector with treatments
#' @param response Numerical vector containing the response of the experiment.
#' @param norm Error normality test (\emph{default} is Shapiro-Wilk)
#' @param homog Homogeneity test of variances (\emph{default} is Bartlett)
#' @param mcomp Multiple comparison test (Tukey (\emph{default}), LSD, Scott-Knott and Duncan)
#' @param quali Defines whether the factor is quantitative or qualitative (\emph{default} is qualitative)
#' @param alpha.f Level of significance of the F test (\emph{default} is 0.05)
#' @param alpha.t Significance level of the multiple comparison test (\emph{default} is 0.05)
#' @param grau Degree of polynomial in case of quantitative factor (\emph{default} is 1)
#' @param transf Applies data transformation (\emph{default} is 1; for log consider 0)
#' @param test "parametric" - Parametric test or "noparametric" - non-parametric test
#' @param p.adj Method for adjusting p values for Kruskal-Wallis ("none","holm","hommel", "hochberg", "bonferroni", "BH", "BY", "fdr")
#' @param geom Graph type (columns, boxes or segments)
#' @param theme ggplot2 theme (\emph{default} is theme_classic())
#' @param sup Number of units above the standard deviation or average bar on the graph
#' @param CV Plotting the coefficient of variation and p-value of Anova (\emph{default} is TRUE)
#' @param ylab Variable response name (Accepts the \emph{expression}() function)
#' @param xlab Treatments name (Accepts the \emph{expression}() function)
#' @param textsize Font size
#' @param fill Defines chart color (to generate different colors for different treatments, define fill = "trat")
#' @param angle x-axis scale text rotation
#' @param family Font family
#' @param dec Number of cells
#' @param addmean Plot the average value on the graph (\emph{default} is TRUE)
#' @param errorbar Plot the standard deviation bar on the graph (In the case of a segment and column graph) - \emph{default} is TRUE
#' @param posi Legend position
#' @param point Defines whether to plot mean ("mean"), mean with standard deviation ("mean_sd" - \emph{default}) or mean with standard error (\emph{default} - "mean_se"). Only for quali=F.
#' @param angle.label label angle
#' @import ggplot2
#' @import stats
#' @import multcompView
#' @import ScottKnott
#' @importFrom crayon green
#' @importFrom crayon bold
#' @importFrom crayon italic
#' @importFrom crayon red
#' @importFrom crayon blue
#' @importFrom crayon black
#' @importFrom nortest lillie.test
#' @importFrom nortest ad.test
#' @importFrom nortest cvm.test
#' @importFrom nortest pearson.test
#' @importFrom nortest sf.test
#' @importFrom utils setTxtProgressBar
#' @importFrom utils txtProgressBar
#' @importFrom graphics abline
#' @importFrom ggrepel geom_label_repel
#' @importFrom ggrepel geom_text_repel
#' @importFrom emmeans emmeans
#' @importFrom multcomp cld
#' @importFrom ARTool artlm
#' @importFrom ARTool art
#' @importFrom lme4 lmer
#' @importFrom graphics par
#' @importFrom utils read.table
#' @importFrom stringr str_trim
#' @importFrom reshape2 melt
#' @importFrom Hmisc rcorr
#' @importFrom cowplot plot_grid
#' @importFrom lmtest dwtest
#' @note The ordering of the graph is according to the sequence in which the factor levels are arranged in the data sheet. The bars of the column and segment graphs are standard deviation.
#' @note Post hoc test in nonparametric is using the criterium Fisher's least significant difference (p-adj="holm").
#' @note CV and p-value of the graph indicate coefficient of variation and p-value of the F test of the analysis of variance.
#' @note In the final output when transformation (transf argument) is different from 1, the columns resp and respo in the mean test are returned, indicating transformed and non-transformed mean, respectively.
#' @references
#'
#' Principles and procedures of statistics a biometrical approach Steel, Torry and Dickey. Third Edition 1997
#'
#' Multiple comparisons theory and methods. Departament of statistics the Ohio State University. USA, 1996. Jason C. Hsu. Chapman Hall/CRC.
#'
#' W.J. Conover, Practical Nonparametrics Statistics. 1999
#'
#' Ramalho M.A.P., Ferreira D.F., Oliveira A.C. 2000. Experimentacao em Genetica e Melhoramento de Plantas. Editora UFLA.
#'
#' Scott R.J., Knott M. 1974. A cluster analysis method for grouping mans in the analysis of variance. Biometrics, 30, 507-512.
#'
#' Mendiburu, F., and de Mendiburu, M. F. (2019). Package ‘agricolae’. R Package, Version, 1-2.
#'
#' HOTHORN, Torsten et al. Package ‘lmtest’. Testing linear regression models. https://cran. r-project. org/web/packages/lmtest/lmtest. pdf. Accessed, v. 6, 2015.
#'
#' @return The table of analysis of variance, the test of normality of errors (Shapiro-Wilk, Lilliefors, Anderson-Darling, Cramer-von Mises, Pearson and Shapiro-Francia), the test of homogeneity of variances (Bartlett or Levene), the test of independence of Durbin-Watson errors, the test of multiple comparisons (Tukey, LSD, Scott-Knott or Duncan) or adjustment of regression models up to grade 3 polynomial, in the case of quantitative treatments. Non-parametric analysis can be used by the Kruskal-Wallis test. The column, segment or box chart for qualitative treatments is also returned. The function also returns a standardized residual plot.
#' @keywords DIC
#' @keywords Experimental
#' @seealso \link{DBC} \link{DQL}
#' @export
#' @examples
#' library(AgroR)
#' data(pomegranate)
#'
#' with(pomegranate, DIC(trat, WL, ylab = "Weight loss (%)")) # tukey
#' with(pomegranate, DIC(trat, WL, mcomp = "sk", ylab = "Weight loss (%)"))
#' with(pomegranate, DIC(trat, WL, mcomp = "duncan", ylab = "Weight loss (%)"))
#'
#' #=============================
#' # Kruskal-Wallis
#' #=============================
#' with(pomegranate, DIC(trat, WL, test = "noparametric", ylab = "Weight loss (%)"))
#'
#'
#' #=============================
#' # chart type
#' #=============================
#' with(pomegranate, DIC(trat, WL, geom="point", ylab = "Weight loss (%)"))
#' with(pomegranate, DIC(trat, WL, ylab = "Weight loss (%)", xlab="Treatments"))
#'
#' #=============================
#' # quantitative factor
#' #=============================
#' data("phao")
#' with(phao, DIC(dose,comp,quali=FALSE,grau=2,
#'                xlab = expression("Dose"~(g~vase^-1)),
#'                ylab="Leaf length (cm)"))
#'
#' #=============================
#' # data transformation
#' #=============================
#' data("pepper")
#' with(pepper, DIC(Acesso, VitC, transf = 0,ylab="Vitamin C"))

DIC <- function(trat,
                response,
                norm="sw",
                homog="bt",
                mcomp="tukey",
                quali=TRUE,
                alpha.f=0.05,
                alpha.t=0.05,
                grau=1,
                transf=1,
                test="parametric",
                p.adj="holm",
                geom="bar",
                theme=theme_classic(),
                sup=NA,
                CV=TRUE,
                ylab="Response",
                xlab="",
                fill="lightblue",
                angle=0,
                family="sans",
                textsize=12,
                dec=3,
                addmean=TRUE,
                errorbar=TRUE,
                posi="top",
                point="mean_sd",
                angle.label=0){
  if(is.na(sup==TRUE)){sup=0.1*mean(response)}
  if(angle.label==0){hjust=0.5}else{hjust=0}
  requireNamespace("nortest")
  requireNamespace("ScottKnott")
  requireNamespace("crayon")
  requireNamespace("ggplot2")

  if(test=="parametric"){
  if(transf==1){resp=response}else{resp=(response^transf-1)/transf}
  if(transf==0){resp=log(response)}
  if(transf==0.5){resp=sqrt(response)}
  if(transf==-0.5){resp=1/sqrt(response)}
  if(transf==-1){resp=1/response}
  trat1=trat
  trat=as.factor(trat)
  a = anova(aov(resp ~ trat))
  aa = summary(aov(resp ~ trat))
  b = aov(resp ~ trat)
  anava=a
  colnames(anava)=c("GL","SQ","QM","Fcal","p-value")
  respad=b$residuals/sqrt(a$`Mean Sq`[2])
  out=respad[respad>3 | respad<(-3)]
  out=names(out)
  out=if(length(out)==0)("No discrepant point")else{out}
  if(norm=="sw"){norm1 = shapiro.test(b$res)}
  if(norm=="li"){norm1=nortest::lillie.test(b$residuals)}
  if(norm=="ad"){norm1=nortest::ad.test(b$residuals)}
  if(norm=="cvm"){norm1=nortest::cvm.test(b$residuals)}
  if(norm=="pearson"){norm1=nortest::pearson.test(b$residuals)}
  if(norm=="sf"){norm1=nortest::sf.test(b$residuals)}
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
  resids=b$residuals/sqrt(a$`Mean Sq`[2])
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
    black("As the calculated p-value is greater than the 5% significance level,hypothesis H0 is not rejected. Therefore, the variances can be considered homogeneous")}
      else {"As the calculated p-value is less than the 5% significance level, H0 is rejected.Therefore, the variances are not homogeneous"})
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
      else {"As the calculated p-value is less than the 5% significance level, H0 is rejected.Therefore, errors are not independent"})
  cat(green(bold("\n-----------------------------------------------------------------\n")))
  cat(green(bold("Additional Information")))
  cat(green(bold("\n-----------------------------------------------------------------\n")))
  cat(paste("\nCV (%) = ",round(sqrt(a$`Mean Sq`[2])/mean(resp,na.rm=TRUE)*100,2)))
  cat(paste("\nR-squared = ",round(a$`Mean Sq`[1]/(a$`Mean Sq`[2]+a$`Mean Sq`[1]),2)))
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
  cat("\n\n")
  message(if (a$`Pr(>F)`[1]<alpha.f){
    black("As the calculated p-value, it is less than the 5% significance level.The hypothesis H0 of equality of means is rejected. Therefore, at least two treatments differ")}
      else {"As the calculated p-value is greater than the 5% significance level, H0 is not rejected"})
  cat(green(bold("\n\n-----------------------------------------------------------------\n")))
  if(quali==TRUE){cat(green(bold("Multiple Comparison Test")))}else{cat(green(bold("Regression")))}
  cat(green(bold("\n-----------------------------------------------------------------\n")))
  if(quali==TRUE){
  if(mcomp=="tukey"){
    letra <- TUKEY(b, "trat", alpha=alpha.t)
    letra1 <- letra$groups; colnames(letra1)=c("resp","groups")}
  if(mcomp=="sk"){
    # letra=SK(b,"trat",sig.level=alpha.t)
    # letra1=data.frame(resp=letra$m.inf[,1],groups=letters[letra$groups])
    nrep=table(trat)[1]
    medias=sort(tapply(resp,trat,mean),decreasing = TRUE)
    letra=scottknott(means = medias,
             df1 = a$Df[2],
             nrep = nrep,
             QME = a$`Mean Sq`[2],
             alpha = alpha.t)
    letra1=data.frame(resp=medias,groups=letra)}
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
  message(if(transf=="1"){}else{blue("\nNOTE: resp = transformed means; respO = averages without transforming\n")})
  if(transf==1 && norm1$p.value<0.05 | transf==1 && indep$p.value<0.05 | transf==1 &&homog1$p.value<0.05){
    message("\nYour analysis is not valid, suggests using a non-parametric test and try to transform the data")
    }
  else{}
  if(transf != 1 && norm1$p.value<0.05 | transf!=1 && indep$p.value<0.05 | transf!=1 && homog1$p.value<0.05){cat(red("\nWarning!!! Your analysis is not valid, suggests using a non-parametric test"))}else{}
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
  if(geom=="bar"){grafico=ggplot(dadosm,aes(x=trats,y=media))
    if(fill=="trat"){grafico=grafico+
      geom_col(aes(fill=trats),color=1)}else{grafico=grafico+
      geom_col(aes(fill=trats),fill=fill,color=1)}
  if(errorbar==TRUE){grafico=grafico+
    geom_text(aes(y=media+sup+if(sup<0){-desvio}else{desvio},
                  label=letra),family=family,angle=angle.label, hjust=hjust)}
  if(errorbar==FALSE){grafico=grafico+
    geom_text(aes(y=media+sup,label=letra),family=family,angle=angle.label, hjust=hjust)}
  if(errorbar==TRUE){grafico=grafico+
    geom_errorbar(data=dadosm,aes(ymin=media-desvio,
                                  ymax=media+desvio,color=1),
                  color="black",width=0.3)}}
  if(geom=="point"){grafico=ggplot(dadosm,aes(x=trats,
                                              y=media))
  if(errorbar==TRUE){grafico=grafico+
    geom_text(aes(y=media+sup+if(sup<0){-desvio}else{desvio},
                  label=letra),family=family,angle=angle.label, hjust=hjust)}
  if(errorbar==FALSE){grafico=grafico+
    geom_text(aes(y=media+sup,
                  label=letra),family=family,angle=angle.label, hjust=hjust)}
  if(errorbar==TRUE){grafico=grafico+
    geom_errorbar(data=dadosm,
                  aes(ymin=media-desvio,
                      ymax=media+desvio,color=1),
                  color="black",width=0.3)}
  if(fill=="trat"){grafico=grafico+
    geom_point(aes(color=trats),size=5)}
  else{grafico=grafico+
    geom_point(aes(color=trats),
               color="black",
               fill=fill,shape=21,size=5)}}
  if(geom=="box"){
  datam1=data.frame(trats=factor(trat,levels = unique(as.character(trat))),
                    response)
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
  grafico=ggplot(datam1,aes(x=trats,y=response))
  if(fill=="trat"){grafico=grafico+geom_boxplot(aes(fill=trats))}
  else{grafico=grafico+
    geom_boxplot(aes(fill=trats),fill=fill)}
  grafico=grafico+
    geom_text(data=dadosm2,
              aes(y=superior,
                  label=letra),
              family = family,angle=angle.label, hjust=hjust)}
  grafico=grafico+
    theme+
    ylab(ylab)+
    xlab(xlab)+
    theme(text = element_text(size=textsize,color="black", family = family),
          axis.text = element_text(size=textsize,color="black", family = family),
          axis.title = element_text(size=textsize,color="black", family = family),
          legend.position = "none")
  if(angle !=0){grafico=grafico+
    theme(axis.text.x=element_text(hjust = 1.01,angle = angle))}
  if(CV==TRUE){grafico=grafico+
    labs(caption=paste("p-value = ", if(a$`Pr(>F)`[1]<0.0001){paste("<", 0.0001)}
                                                  else{paste("=", round(a$`Pr(>F)`[1],4))},"; CV = ",
                                                  round(abs(sqrt(a$`Mean Sq`[2])/mean(resp))*100,2),"%"))}
  grafico=as.list(grafico)
  }
  if(quali==FALSE){
  trat=trat1
  # trat=as.numeric(as.character(trat))
  if(grau==1){graph=polynomial(trat,response, grau = 1,textsize=textsize,xlab=xlab,ylab=ylab, family=family,posi=posi,point=point)}
  if(grau==2){graph=polynomial(trat,response, grau = 2,textsize=textsize,xlab=xlab,ylab=ylab, family=family,posi=posi,point=point)}
  if(grau==3){graph=polynomial(trat,response, grau = 3,textsize=textsize,xlab=xlab,ylab=ylab, family=family,posi=posi,point=point)}
  grafico=graph[[1]]
  }}
  if(test=="noparametric"){
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
    krusk=kruskal(response,trat,p.adj = p.adj,alpha=alpha.t)
    cat(green(bold("\n\n-----------------------------------------------------------------\n")))
    cat(green(italic("Statistics")))
    cat(green(bold("\n-----------------------------------------------------------------\n")))
    print(krusk$statistics)
    cat(green(bold("\n\n-----------------------------------------------------------------\n")))
    cat(green(italic("Parameters")))
    cat(green(bold("\n-----------------------------------------------------------------\n")))
    print(krusk$parameters)
    cat(green(bold("\n\n-----------------------------------------------------------------\n")))
    cat(green(italic("Multiple Comparison Test")))
    cat(green(bold("\n-----------------------------------------------------------------\n")))
    saida=cbind(krusk$means[,c(1,3)],krusk$groups[rownames(krusk$means),])
    colnames(saida)=c("Mean","SD","Rank","Groups")
    print(saida)
    dadosm=data.frame(krusk$means,krusk$groups[rownames(krusk$means),])
    dadosm$trats=factor(rownames(dadosm),levels = unique(trat))
    dadosm$media=tapply(response,trat,mean, na.rm=TRUE)[rownames(krusk$means)]
    if(point=="mean_sd"){dadosm$std=tapply(response,trat,sd, na.rm=TRUE)[rownames(krusk$means)]}
    if(point=="mean_se"){dadosm$std=tapply(response, trat, sd, na.rm=TRUE)/
      sqrt(tapply(response, trat, length))[rownames(krusk$means)]}
    if(addmean==TRUE){dadosm$letra=paste(format(dadosm$response,digits = dec),dadosm$groups)}
    if(addmean==FALSE){dadosm$letra=dadosm$groups}
    trats=dadosm$trats
    limite=dadosm$limite
    media=dadosm$media
    std=dadosm$std
    letra=dadosm$letra
    if(geom=="bar"){grafico=ggplot(dadosm,
                                   aes(x=trats,y=response))
      if(fill=="trat"){grafico=grafico+
        geom_col(aes(fill=trats),color=1)}
      else{grafico=grafico+
        geom_col(aes(fill=trats),fill=fill,color=1)}
    if(errorbar==TRUE){grafico=grafico+
      geom_text(aes(y=media+sup+if(sup<0){-std}else{std},
                    label=letra),family=family,angle=angle.label, hjust=hjust)}
    if(errorbar==FALSE){grafico=grafico+
      geom_text(aes(y=media+sup,label=letra),family=family,angle=angle.label, hjust=hjust)}
    if(errorbar==TRUE){grafico=grafico+
      geom_errorbar(data=dadosm,aes(ymin=response-std,
                                    ymax=response+std,
                                    color=1),
                    color="black",width=0.3)}}
    if(geom=="point"){grafico=ggplot(dadosm,
                                     aes(x=trats,
                                         y=response))
    if(errorbar==TRUE){grafico=grafico+
      geom_text(aes(y=media+sup+if(sup<0){-std}else{std},
                    label=letra),
                family=family,angle=angle.label, hjust=hjust)}
    if(errorbar==FALSE){grafico=grafico+
      geom_text(aes(y=media+sup,
                    label=letra),
                family=family,angle=angle.label, hjust=hjust)}
    if(errorbar==TRUE){grafico=grafico+
      geom_errorbar(data=dadosm,
                    aes(ymin=response-std,
                        ymax=response+std,
                        color=1),
                    color="black",width=0.3)}
    if(fill=="trat"){grafico=grafico+
      geom_point(aes(color=trats),size=5)}
    else{grafico=grafico+
      geom_point(aes(color=trats),
                 color="black",
                 fill=fill,shape=21,size=5)}}
    if(geom=="box"){
    datam1=data.frame(trats=factor(trat,levels = unique(as.character(trat))),response)
    dadosm2=data.frame(krusk$means)
    dadosm2$trats=rownames(dadosm2)
    dadosm2$limite=dadosm2$response+dadosm2$std
    dadosm2$letra=paste(format(dadosm2$response,digits = dec),
                        dadosm$groups)
    dadosm2=dadosm2[unique(as.character(trat)),]
    trats=dadosm2$trats
    limite=dadosm2$limite
    letra=dadosm2$letra
    stat_box=ggplot(datam1,aes(x=trats,y=response))+geom_boxplot()
    superior=ggplot_build(stat_box)$data[[1]]$ymax
    dadosm2$superior=superior+sup

    grafico=ggplot(datam1,
                   aes(x=trats,
                       y=response))
      if(fill=="trat"){grafico=grafico+
        geom_boxplot(aes(fill=1))}
    else{grafico=grafico+
      geom_boxplot(aes(fill=trats),fill=fill)}
    grafico=grafico+
      geom_text(data=dadosm2,
                aes(y=superior,
                    label=letra),
                family = family,angle=angle.label, hjust=hjust)}
    grafico=grafico+theme+
      ylab(ylab)+
      xlab(xlab)+
      theme(text = element_text(size=textsize,color="black", family = family),
            axis.title = element_text(size=textsize,color="black", family = family),
            axis.text = element_text(size=textsize,color="black", family = family),
            legend.position = "none")
    if(angle !=0){grafico=grafico+theme(axis.text.x=element_text(hjust = 1.01,angle = angle))}
}
  if(quali==TRUE){print(grafico)}
  graficos=list(grafico)#[[1]]
}
