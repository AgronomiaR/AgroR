#' Analysis: DBC experiment in double factorial design with an additional treatment
#' @description Analysis of an experiment conducted in a randomized block design in a double factorial scheme using analysis of variance of fixed effects.
#' @author Gabriel Danilo Shimizu, \email{shimizu@uel.br}
#' @author Leandro Simoes Azeredo Goncalves
#' @author Rodrigo Yudi Palhaci Marubayashi
#' @param f1 Numeric or complex vector with factor 1 levels
#' @param f2 Numeric or complex vector with factor 2 levels
#' @param block Numeric or complex vector with repetitions
#' @param response Numerical vector containing the response of the experiment.
#' @param responseAd Numerical vector with additional treatment responses
#' @param norm Error normality test (\emph{default} is Shapiro-Wilk)
#' @param homog Homogeneity test of variances (\emph{default} is Bartlett)
#' @param mcomp Multiple comparison test (Tukey (\emph{default}), LSD and Duncan)
#' @param quali Defines whether the factor is quantitative or qualitative (\emph{qualitative})
#' @param alpha.f Level of significance of the F test (\emph{default} is 0.05)
#' @param alpha.t Significance level of the multiple comparison test (\emph{default} is 0.05)
#' @param grau Degree of polynomial in case of quantitative factor (\emph{default} is 1)
#' @param transf Applies data transformation (default is 1; for log consider 0)
#' @param geom Graph type (columns or segments (For simple effect only))
#' @param theme ggplot2 theme (\emph{default} is theme_classic())
#' @param ylab Variable response name (Accepts the \emph{expression}() function)
#' @param xlab Treatments name (Accepts the \emph{expression}() function)
#' @param legend Legend title name
#' @param ad.label Aditional label
#' @param fill Defines chart color (to generate different colors for different treatments, define fill = "trat")
#' @param angle x-axis scale text rotation
#' @param textsize Font size
#' @param dec Number of cells
#' @param family Font family
#' @param addmean Plot the average value on the graph (\emph{default} is TRUE)
#' @param errorbar Plot the standard deviation bar on the graph (In the case of a segment and column graph) - \emph{default} is TRUE
#' @param CV Plotting the coefficient of variation and p-value of Anova (\emph{default} is TRUE)
#' @param sup Number of units above the standard deviation or average bar on the graph
#' @param color Column chart color (\emph{default} is "rainbow")
#' @param posi legend position
#' @param ylim y-axis scale
#' @param point if quali=F, defines whether to plot all points ("all"), mean ("mean"), standard deviation ("mean_sd") or mean with standard error (\emph{default} - "mean_se").
#' @param angle.label label angle
#' @note The ordering of the graph is according to the sequence in which the factor levels are arranged in the data sheet. The bars of the column and segment graphs are standard deviation.
#' @note The function does not perform multiple regression in the case of two quantitative factors.
#' @note The assumptions of variance analysis disregard additional treatment
#' @note In the final output when transformation (transf argument) is different from 1, the columns resp and respo in the mean test are returned, indicating transformed and non-transformed mean, respectively.
#' @return The table of analysis of variance, the test of normality of errors (Shapiro-Wilk, Lilliefors, Anderson-Darling, Cramer-von Mises, Pearson and Shapiro-Francia), the test of homogeneity of variances (Bartlett or Levene), the test of independence of Durbin-Watson errors, the test of multiple comparisons (Tukey, LSD, Scott-Knott or Duncan) or adjustment of regression models up to grade 3 polynomial, in the case of quantitative treatments. The column chart for qualitative treatments is also returned.
#' @keywords DBC
#' @keywords Factorial
#' @keywords Aditional
#' @seealso \link{FAT2DBC.art}
#' @seealso \link{FAT2DBC}
#' @seealso \link{dunnett}
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
#' @export
#' @examples
#' library(AgroR)
#' data(cloro)
#' respAd=c(268, 322, 275, 350, 320)
#' with(cloro, FAT2DBC.ad(f1, f2, bloco, resp, respAd, ylab="Number of nodules", legend = "Stages"))

FAT2DBC.ad=function(f1,
                    f2,
                    block,
                    response,
                    responseAd,
                    norm="sw",
                    homog="bt",
                    mcomp="tukey",
                    alpha.f=0.05,
                    alpha.t=0.05,
                    quali=c(TRUE,TRUE),
                    grau=NA,
                    transf=1,
                    geom="bar",
                    theme=theme_classic(),
                    ylab="Response",
                    xlab="",
                    legend="Legend",
                    ad.label="Additional",
                    color="rainbow",
                    fill="lightblue",
                    textsize=12,
                    addmean=TRUE,
                    errorbar=TRUE,
                    CV=TRUE,
                    dec=3,
                    angle=0,
                    posi="right",
                    family="sans",
                    point="mean_sd",
                    sup=NA,
                    ylim=NA,
                    angle.label=0){
  if(angle.label==0){hjust=0.5}else{hjust=0}
  requireNamespace("ScottKnott")
  requireNamespace("crayon")
  requireNamespace("ggplot2")
  requireNamespace("nortest")
  fator1=f1
  fator2=f2
  fator1a=fator1
  fator2a=fator2
  block=as.factor(block)
  # ================================
  # Transformacao de dados
  # ================================
  if(transf==1){resp=response}else{resp=(response^transf-1)/transf}
  if(transf==0){resp=log(response)}
  if(transf==0.5){resp=sqrt(response)}
  if(transf==-0.5){resp=1/sqrt(response)}
  if(transf==-1){resp=1/response}
  if(transf==1){respAd=responseAd}else{respAd=(responseAd^transf-1)/transf}
  if(transf==0){respAd=log(responseAd)}
  if(transf==0.5){respAd=sqrt(responseAd)}
  if(transf==-0.5){respAd=1/sqrt(responseAd)}
  if(transf==-1){respAd=1/responseAd}

  if(is.na(sup==TRUE)){sup=0.1*mean(response)}
  Fator1=factor(fator1, levels = unique(fator1))
  Fator2=factor(fator2, levels = unique(fator2))
  nv1 <- length(summary(Fator1))
  nv2 <- length(summary(Fator2))
  lf1 <- levels(Fator1)
  lf2 <- levels(Fator2)
  fac.names = c("F1", "F2")
  fatores <- cbind(fator1, fator2)
  J = length(respAd)
  n.trat2 <- nv1 * nv2
  anavaF2 <- summary(aov(resp ~ Fator1 * Fator2 + block))
  anava=anavaF2[[1]][c(1:4),]
  col1 <- numeric(0)
  for (i in 1:n.trat2) {
    col1 <- c(col1, rep(i, J))
  }
  col1 <- c(col1, rep("ad", J))
  col2 <- c(block, rep(1:J))
  col3 <- c(resp, respAd)
  tabF2ad <- data.frame(TRAT2 = col1, REP = col2, RESP2 = col3)
  TRAT2 <- factor(tabF2ad[, 1])
  anavaf1 <- aov(tabF2ad[, 3] ~ TRAT2)
  anavaTr <- summary(anavaf1)[[1]]
  anava1=rbind(anava,anavaTr)
  anava1$Df[5]=1
  anava1$`Sum Sq`[5]=anava1$`Sum Sq`[5]-sum(anava1$`Sum Sq`[c(1:3)])
  anava1$`Mean Sq`[5]=anava1$`Sum Sq`[5]/anava1$Df[5]
  anava1$`F value`[1:5]=anava1$`Mean Sq`[1:5]/anava1$`Mean Sq`[6]
  for(i in 1:nrow(anava1)-1){
    anava1$`Pr(>F)`[i]=1-pf(anava1$`F value`[i],anava1$Df[i],anava1$Df[6])
  }
  rownames(anava1)[5]="Ad x Factorial"
  anava=anava1
  b=aov(resp ~ Fator1 * Fator2+block)
  an=anova(b)
  respad=b$residuals/sqrt(an$`Mean Sq`[5])
  out=respad[respad>3 | respad<(-3)]
  out=names(out)
  out=if(length(out)==0)("No discrepant point")else{out}

  if(norm=="sw"){norm1 = shapiro.test(b$res)}
  if(norm=="li"){norm1=lillie.test(b$residuals)}
  if(norm=="ad"){norm1=ad.test(b$residuals)}
  if(norm=="cvm"){norm1=cvm.test(b$residuals)}
  if(norm=="pearson"){norm1=pearson.test(b$residuals)}
  if(norm=="sf"){norm1=sf.test(b$residuals)}
  trat=as.factor(paste(Fator1,Fator2))
  c=aov(resp~trat)
  if(homog=="bt"){
    homog1 = bartlett.test(b$res ~ trat)
    statistic=homog1$statistic
    phomog=homog1$p.value
    method=paste("Bartlett test","(",names(statistic),")",sep="")
  }
  if(homog=="levene"){
    homog1 = levenehomog(c$res~trat)
    statistic=homog1$`F value`[1]
    phomog=homog1$`Pr(>F)`[1]
    method="Levene's Test (center = median)(F)"
    names(homog1)=c("Df", "F value","p.value")}

  indep = dwtest(b)
  resids=b$residuals/sqrt(an$`Mean Sq`[5])
  Ids=ifelse(resids>3 | resids<(-3), "darkblue","black")
  residplot=ggplot(data=data.frame(resids,Ids),
                   aes(y=resids,x=1:length(resids)))+
    geom_point(shape=21,color="gray",fill="gray",size=3)+
    labs(x="",y="Standardized residuals")+
    geom_text(x=1:length(resids),label=1:length(resids),
              color=Ids,size=4)+
    scale_x_continuous(breaks=1:length(resids))+
    theme_classic()+theme(axis.text.y = element_text(size=12),
                          axis.text.x = element_blank())+
    geom_hline(yintercept = c(0,-3,3),lty=c(1,2,2),color="red",size=1)

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

  message(if(homog1$p.value[1]>0.05){
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
  cat(paste("\nCV (%) = ",round(sqrt(anava$`Mean Sq`[6])/mean(c(resp,respAd),na.rm=TRUE)*100,2)))
  cat(paste("\nMean Factorial = ",round(mean(response,na.rm=TRUE),4)))
  cat(paste("\nMedian Factorial = ",round(median(response,na.rm=TRUE),4)))
  cat(paste("\nMean Aditional = ",round(mean(responseAd,na.rm=TRUE),4)))
  cat(paste("\nMedian Aditional = ",round(median(responseAd,na.rm=TRUE),4)))
  cat("\nPossible outliers = ", out)
  cat("\n")
  cat(green(bold("\n-----------------------------------------------------------------\n")))
  cat(green(bold("Analysis of Variance")))
  cat(green(bold("\n-----------------------------------------------------------------\n")))
  anava1=as.matrix(data.frame(anava))
  colnames(anava1)=c("Df","Sum Sq","Mean.Sq","F value","Pr(F)" )
  print(anava1,na.print = "")
  cat("\n")
  if(anava$`Pr(>F)`[5]<alpha.f){"The additional treatment does differ from the factorial by the F test"}else{"The additional treatment does not differ from the factorial by the F test "}
  if(transf==1 && norm1$p.value<0.05 | transf==1 && indep$p.value<0.05 | transf==1 &&homog1$p.value<0.05){
    message("\nYour analysis is not valid, suggests using a non-parametric test and try to transform the data\n")}else{}
  if(transf != 1 && norm1$p.value<0.05 | transf!=1 && indep$p.value<0.05 | transf!=1 && homog1$p.value<0.05){
    message("\nYour analysis is not valid, suggests using the function FATDIC.art\n")}else{}
  message(if(transf !=1){blue("NOTE: resp = transformed means; respO = averages without transforming\n")})


  #------------------------------------
  # Fatores isolados
  #------------------------------------

  if (anava$`Pr(>F)`[4] > alpha.f)
  { cat(green(bold("-----------------------------------------------------------------\n")))
    cat(green(bold("No significant interaction")))
    cat(green(bold("\n-----------------------------------------------------------------\n")))
    fatores <- data.frame(Fator1 = factor(fator1), Fator2 = factor(fator2))
    fatoresa <- data.frame(Fator1 = fator1a, Fator2 = fator2a)
    graficos=list(1,2,3)

    for (i in 1:2) {if (anava$`Pr(>F)`[i] <= alpha.f)
    {cat(green(bold("\n-----------------------------------------------------------------\n")))
      cat(bold(fac.names[i]))
      cat(green(bold("\n-----------------------------------------------------------------\n")))
      if(quali[i]==TRUE){
        ## Tukey
        if(mcomp=="tukey"){
          letra <- TUKEY(resp, fatores[,i], anava$Df[6],
                            anava$`Mean Sq`[6], alpha=alpha.t)
          letra1 <- letra$groups; colnames(letra1)=c("resp","groups")
          if(transf !=1){letra1$respo=tapply(response,fatores[,i],mean, na.rm=TRUE)[rownames(letra1)]}}
        if(mcomp=="lsd"){
          letra <- LSD(resp, fatores[,i], anava$Df[6],anava$`Mean Sq`[6],alpha=alpha.t)
          letra1 <- letra$groups; colnames(letra1)=c("resp","groups")
          if(transf !=1){letra1$respo=tapply(response,fatores[,i],mean, na.rm=TRUE)[rownames(letra1)]}}
        if (mcomp == "sk"){
          letra=sk(resp, fatores[,i], anava$Df[6],anava$`Sum Sq`[6],alpha.t)
          letra1 <- letra; colnames(letra1)=c("resp","groups")
          if(transf !=1){letra1$respo=tapply(response,fatores[,i],mean, na.rm=TRUE)[rownames(letra1)]}}
        if(mcomp=="duncan"){
          letra <- duncan(resp, fatores[,i],anava$Df[6],anava$`Mean Sq`[6], alpha=alpha.t)
          letra1 <- letra$groups; colnames(letra1)=c("resp","groups")
          if(transf !=1){letra1$respo=tapply(response,fatores[,i],mean, na.rm=TRUE)[rownames(letra1)]}}
        print(letra1)
        ordem=unique(as.vector(unlist(fatores[i])))
        dadosm=data.frame(letra1[ordem,],
                          media=tapply(response, c(fatores[i]), mean, na.rm=TRUE)[ordem],
                          desvio=tapply(response, c(fatores[i]), sd, na.rm=TRUE)[ordem])
        dadosm$trats=factor(rownames(dadosm),levels = ordem)
        dadosm$limite=dadosm$media+dadosm$desvio
        lim.y=dadosm$limite[which.max(abs(dadosm$limite))]
        if(is.na(ylim[1])==TRUE && lim.y<0){ylim=c(1.5*lim.y,0)}
        if(is.na(ylim[1])==TRUE && lim.y>0){ylim=c(0,1.5*lim.y)}
        if(addmean==TRUE){dadosm$letra=paste(format(dadosm$media,digits = dec),dadosm$groups)}
        if(addmean==FALSE){dadosm$letra=dadosm$groups}
        media=dadosm$media
        desvio=dadosm$desvio
        trats=dadosm$trats
        letra=dadosm$letra
        if(geom=="bar"){grafico=ggplot(dadosm,
                                       aes(x=trats,
                                           y=media))
        if(fill=="trat"){grafico=grafico+
          geom_col(aes(fill=trats),color=1)}
        else{grafico=grafico+
          geom_col(aes(fill=trats),
                   fill=fill,color=1)}
        grafico=grafico+theme+ylab(ylab)+xlab(xlab)+ylim(ylim)
        if(errorbar==TRUE){grafico=grafico+
          geom_text(aes(y=media+
                          sup+if(sup<0){-desvio}else{desvio},
                        label=letra),family=family,angle=angle.label, hjust=hjust)}
        if(errorbar==FALSE){grafico=grafico+
          geom_text(aes(y=media+sup,
                        label=letra),family=family,angle=angle.label, hjust=hjust)}
        if(errorbar==TRUE){grafico=grafico+
          geom_errorbar(data=dadosm,
                        aes(ymin=media-desvio,
                            ymax=media+desvio,color=1),
                        color="black",width=0.3)}
        if(angle !=0){grafico=grafico+theme(axis.text.x=element_text(hjust = 1.01,angle = angle))}
        grafico=grafico+
          theme(text = element_text(size=textsize,color="black",family=family),
                axis.text = element_text(size=textsize,color="black",family=family),
                axis.title = element_text(size=textsize,color="black",family=family))}

        # ================================
        # grafico de pontos
        # ================================
        if(geom=="point"){grafico=ggplot(dadosm,
                                         aes(x=trats,
                                             y=media))
        if(fill=="trat"){grafico=grafico+
          geom_point(aes(color=trats),size=5)}
        else{grafico=grafico+
          geom_point(aes(color=trats),fill="gray",pch=21,color="black",size=5)}
        grafico=grafico+theme+ylab(ylab)+xlab(xlab)+ylim(ylim)
        if(errorbar==TRUE){grafico=grafico+
          geom_text(aes(y=media+sup+
                          if(sup<0){-desvio}else{desvio},
                        label=letra),family=family,angle=angle.label, hjust=hjust)}
        if(errorbar==FALSE){grafico=grafico+
          geom_text(aes(y=media+sup,label=letra),family=family,angle=angle.label, hjust=hjust)}
        if(errorbar==TRUE){grafico=grafico+
          geom_errorbar(data=dadosm,
                        aes(ymin=media-desvio,
                            ymax=media+desvio,color=1),
                        color="black", width=0.3)}
        if(angle !=0){grafico=grafico+theme(axis.text.x=element_text(hjust = 1.01,angle = angle))}
        grafico=grafico+
          theme(text = element_text(size=textsize,color="black",family=family),
                axis.text = element_text(size=textsize,color="black",family=family),
                axis.title = element_text(size=textsize,color="black",family=family))}
        grafico=grafico+
          geom_hline(aes(color=ad.label,group=ad.label,yintercept=mean(responseAd,na.rm=T)),lty=2)+
          scale_color_manual(values = "black")+labs(color="")
        if(CV==TRUE){grafico=grafico+labs(caption=paste("p-value = ", if(anava$`Pr(>F)`[i]<0.0001){paste("<", 0.0001)}
                                                        else{paste("=", round(anava$`Pr(>F)`[i],4))},"; CV = ",
                                                        round(abs(sqrt(anava$`Mean Sq`[6])/mean(c(resp,respAd),na.rm=TRUE))*100,2),"%"))}
        if(color=="gray"){grafico=grafico+scale_fill_grey()}
        print(grafico)
        cat("\n\n")
      }

      # Regression
      if(quali[i]==FALSE){
        # dose=as.numeric(as.character(as.vector(unlist(fatores[i]))))
        dose=as.vector(unlist(fatoresa[i]))
        grafico=polynomial(dose,
                           response,
                           grau = grau,
                           ylab=ylab,
                           xlab=xlab,
                           posi=posi,
                           theme=theme,
                           textsize=textsize,
                           se=errorbar,
                           point=point,
                           family=family)
        grafico=grafico[[1]]}

      # Ns
      #if (a$`Pr(>F)`[i] > alpha.f) {cat("\nAccording to the F test, the means do not differ\n")}
      graficos[[i+1]]=grafico}}
    graficos[[1]]=residplot
    if(anava$`Pr(>F)`[1]>=alpha.f && anava$`Pr(>F)`[2] <alpha.f){
      cat(green(bold("\n-----------------------------------------------------------------\n")))
      cat(green("Isolated factors 1 not significant"))
      cat(green(bold("\n-----------------------------------------------------------------\n")))
      d1=data.frame(tapply(response,fator1,mean, na.rm=TRUE))
      colnames(d1)="Mean"
      print(d1)
    }
    if(anava$`Pr(>F)`[1]<alpha.f && anava$`Pr(>F)`[2] >=alpha.f){
      cat(green(bold("\n-----------------------------------------------------------------\n")))
      cat(green("Isolated factors 2 not significant"))
      cat(green(bold("\n-----------------------------------------------------------------\n")))
      d1=data.frame(tapply(response,fator2,mean, na.rm=TRUE))
      colnames(d1)="Mean"
      print(d1)}
    if(anava$`Pr(>F)`[1]>=alpha.f && anava$`Pr(>F)`[2] >=alpha.f){
      cat(green(bold("\n-----------------------------------------------------------------\n")))
      cat(green("Isolated factors not significant"))
      cat(green(bold("\n-----------------------------------------------------------------\n")))
      print(tapply(response,list(fator1,fator2),mean, na.rm=TRUE))}
  }

  #-----------------------------------
  # Interação
  #-----------------------------------

  # Desdobramento de F1 dentro de F2

  if (anava$`Pr(>F)`[4]  <= alpha.f) {
    fatores <- data.frame(Fator1, Fator2)
    cat(green(bold("-----------------------------------------------------------------\n")))
    cat(green(bold("Significant interaction: analyzing the interaction")))
    cat(green(bold("\n-----------------------------------------------------------------\n")))
    # cat("\n-----------------------------------------------------------------\n")
    # cat("Analyzing ", fac.names[1], " inside of each level of ",fac.names[2])
    # cat("\n-----------------------------------------------------------------\n")
    # cat("\n")
    des1<-aov(resp~Fator2/Fator1+block)

    l1<-vector('list',nv2)
    names(l1)<-names(summary(Fator2))
    v<-numeric(0)
    for(j in 1:nv2) {
      for(i in 0:(nv1-2)) v<-cbind(v,i*nv2+j)
      l1[[j]]<-v
      v<-numeric(0)
    }
    des1.tab<-summary(des1,split=list('Fator2:Fator1'=l1))[[1]]
    nlinhas=nrow(des1.tab)
    des1.tab=des1.tab[-c(nlinhas),]
    des1.tab$`F value`=des1.tab$`Mean Sq`/anava$`Mean Sq`[6]
    des1.tab$`Pr(>F)`=1-pf(des1.tab$`F value`,des1.tab$Df,anava$Df[6])
    print(des1.tab)
    desdobramento1=des1.tab
    if(quali[1]==TRUE & quali[2]==TRUE){
      #-------------------------------------
      # Teste de comparação
      #-------------------------------------
      if (mcomp == "tukey"){
        tukeygrafico=c()
        ordem=c()
        for (i in 1:nv2) {
          trati=fatores[, 1][Fator2 == lf2[i]]
          respi=resp[Fator2 == lf2[i]]
          tukey=TUKEY(respi,trati,anava$Df[6],anava$`Mean Sq`[6],alpha.t)
          if(transf !="1"){tukey$groups$respo=tapply(response[Fator2 == lf2[i]],trati,
                                                     mean, na.rm=TRUE)[rownames(tukey$groups)]}
          tukeygrafico[[i]]=tukey$groups[as.character(unique(trati)),2]
          ordem[[i]]=rownames(tukey$groups[as.character(unique(trati)),])
        }
        letra=unlist(tukeygrafico)
        datag=data.frame(letra, ordem=unlist(ordem))
        datag=datag[order(factor(datag$ordem,levels=unique(Fator1))),]
        letra=datag$letra
      }
      if (mcomp == "duncan"){
        duncangrafico=c()
        ordem=c()
        for (i in 1:nv2) {
          trati=fatores[, 1][Fator2 == lf2[i]]
          respi=resp[Fator2 == lf2[i]]
          duncan=duncan(respi,trati,anava$Df[6],anava$`Mean Sq`[6],alpha.t)
          if(transf !="1"){duncan$groups$respo=tapply(response[Fator2 == lf2[i]],
                                                      trati,mean, na.rm=TRUE)[rownames(duncan$groups)]}
          duncangrafico[[i]]=duncan$groups[as.character(unique(trati)),2]
          ordem[[i]]=rownames(duncan$groups[as.character(unique(trati)),])
        }
        letra=unlist(duncangrafico)
        datag=data.frame(letra, ordem=unlist(ordem))
        datag=datag[order(factor(datag$ordem,levels=unique(Fator1))),]
        letra=datag$letra
      }
      if (mcomp == "lsd"){
        duncangrafico=c()
        ordem=c()
        for (i in 1:nv2) {
          trati=fatores[, 1][Fator2 == lf2[i]]
          respi=resp[Fator2 == lf2[i]]
          lsd=LSD(respi,trati,anava$Df[6],anava$`Mean Sq`[6],alpha.t)
          if(transf !="1"){lsd$groups$respo=tapply(response[Fator2 == lf2[i]],trati,
                                                   mean, na.rm=TRUE)[rownames(lsd$groups)]}
          duncangrafico[[i]]=lsd$groups[as.character(unique(trati)),2]
          ordem[[i]]=rownames(lsd$groups[as.character(unique(trati)),])
        }
        letra=unlist(duncangrafico)
        datag=data.frame(letra, ordem=unlist(ordem))
        datag=datag[order(factor(datag$ordem,levels=unique(Fator1))),]
        letra=datag$letra
      }
      if (mcomp == "sk"){
        skgrafico=c()
        ordem=c()
        for (i in 1:nv2) {
          trati=fatores[, 1][Fator2 == lf2[i]]
          trati=factor(trati,levels = unique(trati))
          respi=resp[Fator2 == lf2[i]]
          sk=sk(respi,trati,anava$Df[6],anava$`Sum Sq`[6],alpha.t)
          if(transf !="1"){sk$groups$respo=tapply(respi,trati,mean, na.rm=TRUE)[rownames(sk$groups)]}
          skgrafico[[i]]=sk[levels(trati),2]
          ordem[[i]]=rownames(sk[levels(trati),])
        }
        letra=unlist(skgrafico)
        datag=data.frame(letra,ordem=unlist(ordem))
        datag$ordem=factor(datag$ordem,levels = unique(datag$ordem))
        datag=datag[order(datag$ordem),]
        letra=datag$letra}
    }

    # Desdobramento de F2 dentro de F1

    cat("\n-----------------------------------------------------------------\n")
    cat("Analyzing ", fac.names[2], " inside of the level of ",fac.names[1])
    cat("\n-----------------------------------------------------------------\n")
    cat("\n")
    des1<-aov(resp~Fator1/Fator2+block)

    l1<-vector('list',nv1)
    names(l1)<-names(summary(Fator1))
    v<-numeric(0)
    for(j in 1:nv1) {
      for(i in 0:(nv2-2)) v<-cbind(v,i*nv1+j)
      l1[[j]]<-v
      v<-numeric(0)
    }
    des1.tab<-summary(des1,split=list('Fator1:Fator2'=l1))[[1]]
    nlinhas=nrow(des1.tab)
    des1.tab=des1.tab[-c(nlinhas),]
    des1.tab$`F value`=des1.tab$`Mean Sq`/anava$`Mean Sq`[6]
    des1.tab$`Pr(>F)`=1-pf(des1.tab$`F value`,des1.tab$Df,anava$Df[6])
    print(des1.tab)
    desdobramento2=des1.tab

    #-------------------------------------
    # Teste de comparação
    #-------------------------------------
    if(quali[1]==TRUE & quali[2]==TRUE){
      if (mcomp == "tukey"){
        tukeygrafico1=c()
        for (i in 1:nv1) {
          trati=fatores[, 2][Fator1 == lf1[i]]
          respi=resp[Fator1 == lf1[i]]
          tukey=TUKEY(respi,trati,anava$Df[6],anava$`Mean Sq`[6],alpha.t)
          if(transf !="1"){tukey$groups$respo=tapply(response[Fator1 == lf1[i]],trati,mean, na.rm=TRUE)[rownames(tukey$groups)]}
          tukeygrafico1[[i]]=tukey$groups[as.character(unique(trati)),2]
        }
        letra1=unlist(tukeygrafico1)
        letra1=toupper(letra1)}
      if (mcomp == "duncan"){
        duncangrafico1=c()
        for (i in 1:nv1) {
          trati=fatores[, 2][Fator1 == lf1[i]]
          respi=resp[Fator1 == lf1[i]]
          duncan=duncan(respi,trati,anava$Df[6],anava$`Mean Sq`[6],alpha.t)
          if(transf !="1"){duncan$groups$respo=tapply(response[Fator1 == lf1[i]],trati,mean, na.rm=TRUE)[rownames(duncan$groups)]}
          duncangrafico1[[i]]=duncan$groups[as.character(unique(trati)),2]
        }
        letra1=unlist(duncangrafico1)
        letra1=toupper(letra1)}
      if (mcomp == "lsd"){
        lsdgrafico1=c()
        for (i in 1:nv1) {
          trati=fatores[, 2][Fator1 == lf1[i]]
          respi=resp[Fator1 == lf1[i]]
          lsd=LSD(respi,trati,anava$Df[6],anava$`Mean Sq`[6],alpha.t)
          if(transf !="1"){lsd$groups$respo=tapply(response[Fator1 == lf1[i]],trati,mean, na.rm=TRUE)[rownames(lsd$groups)]}
          lsdgrafico1[[i]]=lsd$groups[as.character(unique(trati)),2]
        }
        letra1=unlist(lsdgrafico1)
        letra1=toupper(letra1)}
      if (mcomp == "sk"){
        skgrafico1=c()
        for (i in 1:nv1) {
          trati=fatores[, 2][Fator1 == lf1[i]]
          trati=factor(trati,levels = unique(trati))
          respi=resp[Fator1 == lf1[i]]
          sk=sk(respi,trati,anava$Df[6],anava$`Sum Sq`[6],alpha.t)
          if(transf !="1"){sk$respo=tapply(respi,trati,mean, na.rm=TRUE)[rownames(sk)]}
          skgrafico1[[i]]=sk[levels(trati),2]
        }
        letra1=unlist(skgrafico1)
        letra1=toupper(letra1)}
    }

    if(quali[1]==FALSE && color=="gray"| quali[2]==FALSE && color=="gray"){
      if(quali[2]==FALSE){
        Fator2=fator2a#as.numeric(as.character(Fator2))
        grafico=polynomial2(Fator2,response,Fator1,
                            grau = grau,
                            ylab=ylab,
                            xlab=xlab,
                            theme=theme,
                            posi=posi,
                            point=point,
                            textsize=textsize,
                            se=errorbar,
                            family=family,
                            ylim=ylim)}
      if(quali[2]==TRUE){
        Fator1=fator1a#as.numeric(as.character(Fator1))
        grafico=polynomial2(Fator1,
                            response,
                            Fator2,
                            grau = grau,
                            ylab=ylab,
                            xlab=xlab,
                            theme=theme,
                            posi=posi,
                            point=point,
                            textsize=textsize,
                            se=errorbar,
                            family=family,
                            ylim=ylim)}
    }
    if(quali[1]==FALSE && color=="rainbow"| quali[2]==FALSE && color=="rainbow"){
      if(quali[2]==FALSE){
        Fator2=fator2a#as.numeric(as.character(Fator2))
        grafico=polynomial2_color(Fator2,
                                  response,
                                  Fator1,
                                  grau = grau,
                                  ylab=ylab,
                                  xlab=xlab,
                                  theme=theme,
                                  posi=posi,
                                  point=point,
                                  textsize=textsize,
                                  se=errorbar,
                                  family=family,
                                  ylim=ylim)}
      if(quali[2]==TRUE){
        Fator1=fator1a#as.numeric(as.character(Fator1))
        grafico=polynomial2_color(Fator1,
                                  response,
                                  Fator2,
                                  grau = grau,
                                  ylab=ylab,
                                  xlab=xlab,
                                  theme=theme,
                                  posi=posi,
                                  point=point,
                                  textsize=textsize,
                                  se=errorbar,
                                  family=family,
                                  ylim=ylim)}
    }
    # -----------------------------
    # Gráfico de colunas
    #------------------------------
    if(quali[1] & quali[2]==TRUE){
      media=tapply(response,list(Fator1,Fator2), mean, na.rm=TRUE)
      desvio=tapply(response,list(Fator1,Fator2), sd, na.rm=TRUE)
      graph=data.frame(f1=rep(rownames(media),length(colnames(media))),
                       f2=rep(colnames(media),e=length(rownames(media))),
                       media=as.vector(media),
                       desvio=as.vector(desvio))
      limite=graph$media+graph$desvio
      lim.y=limite[which.max(abs(limite))]
      if(is.na(ylim[1])==TRUE && lim.y<0){ylim=c(1.5*lim.y,0)}
      if(is.na(ylim[1])==TRUE && lim.y>0){ylim=c(0,1.5*lim.y)}
      rownames(graph)=paste(graph$f1,graph$f2)
      graph=graph[paste(rep(unique(Fator1),
                            e=length(colnames(media))),
                        rep(unique(Fator2),length(rownames(media)))),]
      graph$letra=letra
      graph$letra1=letra1
      graph$f1=factor(graph$f1,levels = unique(Fator1))
      graph$f2=factor(graph$f2,levels = unique(Fator2))
      if(addmean==TRUE){graph$numero=paste(format(graph$media,digits = dec), graph$letra,graph$letra1, sep="")}
      if(addmean==FALSE){graph$numero=paste(graph$letra,graph$letra1, sep="")}
      f1=graph$f1
      f2=graph$f2
      media=graph$media
      desvio=graph$desvio
      numero=graph$numero

      colint=ggplot(graph,
                    aes(x=f1,
                        y=media,
                        fill=f2))+
        geom_col(position = "dodge",color="black")+
        ylab(ylab)+xlab(xlab)+ylim(ylim)+
        theme
      if(errorbar==TRUE){colint=colint+
        geom_errorbar(data=graph,
                      aes(ymin=media-desvio,
                          ymax=media+desvio),
                      width=0.3,color="black",
                      position = position_dodge(width=0.9))}
      if(errorbar==TRUE){colint=colint+
        geom_text(aes(y=media+sup+if(sup<0){-desvio}else{desvio},
                      label=numero),
                  position = position_dodge(width=0.9),
                  family = family,angle=angle.label,hjust=hjust)}
      if(errorbar==FALSE){colint=colint+
        geom_text(aes(y=media+sup,label=numero),
                  position = position_dodge(width=0.9),
                  family = family,angle=angle.label, hjust=hjust)}
      colint=colint+theme(text=element_text(size=textsize,family = family),
                          axis.text = element_text(size=textsize,color="black",family = family),
                          axis.title = element_text(size=textsize,color="black",family = family),
                          legend.text = element_text(family = family,size = textsize),
                          legend.title = element_text(family = family,size = textsize),
                          legend.position = posi)+labs(fill=legend)+
        geom_hline(aes(color=ad.label,yintercept=mean(responseAd,na.rm=T)),lty=2)+
        scale_color_manual(values = "black")+labs(color="")
      if(CV==TRUE){colint=colint+labs(caption=paste("p-value ", if(anava$`Pr(>F)`[4]<0.0001){paste("<", 0.0001)}
                                                    else{paste("=", round(anava$`Pr(>F)`[4],4))},"; CV = ",
                                                    round(abs(sqrt(anava$`Mean Sq`[6])/mean(c(resp,respAd),na.rm=TRUE))*100,2),"%"))}
      if(angle !=0){colint=colint+
        theme(axis.text.x=element_text(hjust = 1.01,angle = angle))}
      if(color=="gray"){colint=colint+scale_fill_grey()}
      print(colint)
      grafico=colint
      letras=paste(graph$letra,
                   graph$letra1,
                   sep="")
      matriz=data.frame(t(matrix(paste(format(graph$media,digits = dec),letras),ncol = length(levels(Fator1)))))
      rownames(matriz)=levels(Fator1)
      colnames(matriz)=levels(Fator2)
      cat(green(bold("\n-----------------------------------------------------------------\n")))
      cat(green(bold("Final table")))
      cat(green(bold("\n-----------------------------------------------------------------\n")))
      print(matriz)
      message(black("\n\nAverages followed by the same lowercase letter in the column and \nuppercase in the row do not differ by the",mcomp,"(p<",alpha.t,")"))
    }
  }
  if(anava$`Pr(>F)`[4]>alpha.f){
    names(graficos)=c("residplot","graph1","graph2")
    graficos}else{colints=list(residplot,grafico)}
}
