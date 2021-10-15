#' Analysis: DBC experiments in split-plot
#' @description Analysis of an experiment conducted in a randomized block design in a split-plot scheme using fixed effects analysis of variance.
#' @author Gabriel Danilo Shimizu
#' @param f1 Numeric or complex vector with plot levels
#' @param f2 Numeric or complex vector with subplot levels
#' @param block Numeric or complex vector with blocks
#' @param response Numeric vector with responses
#' @param transf Applies data transformation (default is 1; for log consider 0)
#' @param norm Error normality test (\emph{default} is Shapiro-Wilk)
#' @param homog Homogeneity test of variances (\emph{default} is Bartlett)
#' @param mcomp Multiple comparison test (Tukey (\emph{default}), LSD, Scott-Knott and Duncan)
#' @param alpha.f Level of significance of the F test (\emph{default} is 0.05)
#' @param alpha.t Significance level of the multiple comparison test (\emph{default} is 0.05)
#' @param quali Defines whether the factor is quantitative or qualitative (\emph{qualitative})
#' @param grau Degree of polynomial in case of quantitative factor (\emph{default} is 1)
#' @param geom Graph type (columns or segments (For simple effect only))
#' @param theme ggplot2 theme (\emph{default} is theme_classic())
#' @param ylab Variable response name (Accepts the \emph{expression}() function)
#' @param xlab Treatments name (Accepts the \emph{expression}() function)
#' @param fill Defines chart color (to generate different colors for different treatments, define fill = "trat")
#' @param angle x-axis scale text rotation
#' @param family Font family (\emph{default} is sans)
#' @param color When the columns are different colors (Set fill-in argument as "trat")
#' @param legend Legend title name
#' @param errorbar Plot the standard deviation bar on the graph (In the case of a segment and column graph) - \emph{default} is TRUE
#' @param addmean Plot the average value on the graph (\emph{default} is TRUE)
#' @param textsize Font size (\emph{default} is 12)
#' @param dec Number of cells (\emph{default} is 3)
#' @param ylim y-axis limit
#' @param posi Legend position
#' @param point Point type for regression ("mean_se","mean_sd","mean" or "all")
#' @param angle.label Label angle
#' @note The ordering of the graph is according to the sequence in which the factor levels are arranged in the data sheet. The bars of the column and segment graphs are standard deviation.
#' @note In the final output when transformation (transf argument) is different from 1, the columns resp and respo in the mean test are returned, indicating transformed and non-transformed mean, respectively.
#' @import ggplot2
#' @importFrom crayon green
#' @importFrom crayon bold
#' @importFrom crayon italic
#' @importFrom crayon red
#' @importFrom crayon blue
#' @import stats
#' @keywords DBC
#' @keywords split-plot
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
#' @return The table of analysis of variance, the test of normality of errors (Shapiro-Wilk, Lilliefors, Anderson-Darling, Cramer-von Mises, Pearson and Shapiro-Francia), the test of homogeneity of variances (Bartlett), the test of independence of Durbin-Watson errors, the test of multiple comparisons (Tukey, LSD, Scott-Knott or Duncan) or adjustment of regression models up to grade 3 polynomial, in the case of quantitative treatments. Non-parametric analysis can be used by the Friedman test. The column chart for qualitative treatments is also returned. The function also returns a standardized residual plot.
#' @examples
#'
#' #==============================
#' # Example tomate
#' #==============================
#' library(AgroR)
#' data(tomate)
#' with(tomate, PSUBDBC(parc, subp, bloco, resp, ylab="Dry mass (g)"))
#'
#' #==============================
#' # Example orchard
#' #==============================
#' library(AgroR)
#' data(orchard)
#' with(orchard, PSUBDBC(A, B, Bloco, Resp, ylab="CBM"))

PSUBDBC=function(f1,
                 f2,
                 block,
                 response,
                 norm="sw",
                 homog="bt",
                 alpha.f=0.05,
                 alpha.t=0.05,
                 mcomp = "tukey",
                 quali=c(TRUE,TRUE),
                 grau=NA,
                 transf=1,
                 geom="bar",
                 theme=theme_classic(),
                 ylab="Response",
                 xlab="",
                 color="rainbow",
                 textsize=12,
                 dec=3,
                 legend="Legend",
                 errorbar=TRUE,
                 addmean=TRUE,
                 ylim=NA,
                 point="mean_se",
                 fill="lightblue",
                 angle=0,
                 family="sans",
                 posi="right",
                 angle.label=0){
  if(angle.label==0){hjust=0.5}else{hjust=0}
  requireNamespace("crayon")
  requireNamespace("ggplot2")
  requireNamespace("ScottKnott")
  requireNamespace("nortest")
  fator1=f1
  fator2=f2
  fator1a=fator1
  fator2a=fator2
  fac = c("F1", "F2")
  cont <- c(1, 4)
  Fator1 <- factor(fator1, levels = unique(fator1))
  Fator2 <- factor(fator2, levels = unique(fator2))
  bloco <- factor(block)
  lf1 <- levels(Fator1)
  lf2 <- levels(Fator2)
  nv1 <- length(summary(Fator1))
  nv2 <- length(summary(Fator2))
  num=function(x){as.numeric(x)}
  sup=0.1*mean(response)

  # ================================
  # Transformacao de dados
  # ================================
  if(transf==1){resp=response}else{resp=(response^transf-1)/transf}
  if(transf==0){resp=log(response)}
  if(transf==0.5){resp=sqrt(response)}
  if(transf==-0.5){resp=1/sqrt(response)}
  if(transf==-1){resp=1/response}
  graph=data.frame(Fator1,Fator2,resp)
  # -----------------------------
  # Analise de variancia
  # -----------------------------

  mod=aov(resp~Fator1*Fator2+Fator1:bloco+bloco)
  anova=summary(mod)[[1]]
  anova=anova[c(1,3,5,2,4,6),]
  anova$`F value`[1]=anova$`Mean Sq`[1]/anova$`Mean Sq`[3]
  anova$`F value`[2]=anova$`Mean Sq`[2]/anova$`Mean Sq`[3]
  anova$`F value`[3]=NA
  anova$`Pr(>F)`[3]=NA
  anova$`Pr(>F)`[1]=1-pf(anova[1,4],anova[1,1],anova[3,1])
  anova$`Pr(>F)`[2]=1-pf(anova[2,4],anova[2,1],anova[3,1])
  anova1=anova
  anova=data.frame(anova)
  colnames(anova)=colnames(anova1)
  rownames(anova)=c("F1","Block","Error A", "F2", "F1 x F2", "Error B")
  tab=anova

  # -----------------------------
  # Pressupostos
  # -----------------------------
  modp=lme4::lmer(resp~Fator1*Fator2+(1|bloco/Fator1)+bloco)
  resids=residuals(modp,scaled=TRUE)
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

  # Normalidade dos erros
  if(norm=="sw"){norm1 = shapiro.test(resid(modp))}
  if(norm=="li"){norm1=lillie.test(resid(modp))}
  if(norm=="ad"){norm1=ad.test(resid(modp))}
  if(norm=="cvm"){norm1=cvm.test(resid(modp))}
  if(norm=="pearson"){norm1=pearson.test(resid(modp))}
  if(norm=="sf"){norm1=sf.test(resid(modp))}

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

  homog1=bartlett.test(resid(modp)~Fator1)
  homog2=bartlett.test(resid(modp)~Fator2)
  homog3=bartlett.test(resid(modp)~paste(Fator1,Fator2))
  cat(green(bold("\n\n-----------------------------------------------------------------\n")))
  cat(green(bold("Homogeneity of Variances")))
  cat(green(bold("\n-----------------------------------------------------------------\n")))
  cat(green(bold("Plot\n")))
  statistic1=homog1$statistic
  phomog1=homog1$p.value
  method1=paste("Bartlett test","(",names(statistic1),")",sep="")
  homoge1=data.frame(Method=method1,
                     Statistic=statistic1,
                     "p-value"=phomog1)
  rownames(homoge1)=""
  print(homoge1)
  cat("\n")
  message(if(homog1$p.value[1]>0.05){
    black("As the calculated p-value is greater than the 5% significance level, hypothesis H0 is not rejected. Therefore, the variances can be considered homogeneous")}
      else {"As the calculated p-value is less than the 5% significance level, H0 is rejected. Therefore, the variances are not homogeneous"})
  cat("\n----------------------------------------------------\n")
  cat(green(bold("Split-plot\n")))
  statistic2=homog2$statistic
  phomog2=homog2$p.value
  method2=paste("Bartlett test","(",names(statistic2),")",sep="")
  homoge2=data.frame(Method=method2,
                     Statistic=statistic2,
                     "p-value"=phomog2)
  rownames(homoge2)=""
  print(homoge2)
  cat("\n")
  message(if(homog2$p.value[1]>0.05){
    black("As the calculated p-value is greater than the 5% significance level, hypothesis H0 is not rejected. Therefore, the variances can be considered homogeneous")}
      else {"As the calculated p-value is less than the 5% significance level, H0 is rejected. Therefore, the variances are not homogeneous"})
  cat("\n----------------------------------------------------\n")
  cat(green(bold("Interaction\n")))
  statistic3=homog3$statistic
  phomog3=homog3$p.value
  method3=paste("Bartlett test","(",names(statistic3),")",sep="")
  homoge3=data.frame(Method=method3,
                     Statistic=statistic3,
                     "p-value"=phomog3)
  rownames(homoge3)=""
  print(homoge3)
  cat("\n")
  message(if(homog3$p.value[1]>0.05){
    black("As the calculated p-value is greater than the 5% significance level, hypothesis H0 is not rejected. Therefore, the variances can be considered homogeneous")}
      else {"As the calculated p-value is less than the 5% significance level, H0 is rejected. Therefore, the variances are not homogeneous"})

  cat(green(bold("\n-----------------------------------------------------------------\n")))
  cat(green(bold("Additional Information")))
  cat(green(bold("\n-----------------------------------------------------------------\n")))
  cat(paste("\nCV1 (%) = ",round(sqrt(tab$`Mean Sq`[3])/mean(resp,na.rm=TRUE)*100,2)))
  cat(paste("\nCV2 (%) = ",round(sqrt(tab$`Mean Sq`[6])/mean(resp,na.rm=TRUE)*100,2)))
  cat(paste("\nMean = ",round(mean(response,na.rm=TRUE),4)))
  cat(paste("\nMedian = ",round(median(response,na.rm=TRUE),4)))
  #cat("\nPossible outliers = ", out)
  cat("\n")

  cat(green(bold("\n-----------------------------------------------------------------\n")))
  cat(green(bold("Analysis of Variance")))
  cat(green(bold("\n-----------------------------------------------------------------\n")))
  anova$`Pr(>F)`=ifelse(anova$`Pr(>F)`>0.001,round(anova$`Pr(>F)`,3),"p<0.001")
  print(as.matrix(anova),na.print="",quote = FALSE)

  if(transf==1 && norm1$p.value<0.05 |  transf==1 &&homog1$p.value<0.05){
    message("\n Your analysis is not valid, suggests using a try to transform the data\n")}else{}
  message(if(transf !=1){blue("\nNOTE: resp = transformed means; respO = averages without transforming\n")})
  fat <- data.frame(`fator 1`=factor(fator1), `fator 2`=factor(fator2))
  fata <- data.frame(`fator 1`=fator1a, `fator 2`=fator2a)
  #------------------------------------
  # Fatores isolados
  #------------------------------------
  if (as.numeric(tab[5, 5]) > alpha.f)
  {cat(green(bold("-----------------------------------------------------------------\n")))
    cat("No significant interaction")
    cat(green(bold("\n-----------------------------------------------------------------\n")))
    graficos=list(1,2,3)

    for (i in 1:2) {if (num(tab[cont[i], 5]) <= alpha.f)
    {cat(green(bold("\n-----------------------------------------------------------------\n")))
      cat(bold(fac[i]))
      cat(green(bold("\n-----------------------------------------------------------------\n")))
      if(quali[i]==TRUE){
        ## Tukey
        if(mcomp=="tukey"){
          letra <- TUKEY(resp, fat[, i],num(tab[3*i,1]),
                            num(tab[3*i,2])/num(tab[3*i,1]), alpha.t)
          letra1 <- letra$groups; colnames(letra1)=c("resp","groups")
          if(transf !=1){letra1$respo=tapply(response,fat[,i],mean, na.rm=TRUE)[rownames(letra1)]}}
        if(mcomp=="duncan"){
          letra <- duncan(resp, fat[, i],num(tab[3*i,1]),
                               num(tab[3*i,2])/num(tab[3*i,1]), alpha.t)
          letra1 <- letra$groups; colnames(letra1)=c("resp","groups")
          if(transf !=1){letra1$respo=tapply(response,fat[,i],mean, na.rm=TRUE)[rownames(letra1)]}}
        if(mcomp=="lsd"){
          letra <- LSD(resp, fat[, i],num(tab[3*i,1]),
                            num(tab[3*i,2])/num(tab[3*i,1]), alpha.t)
          letra1 <- letra$groups; colnames(letra1)=c("resp","groups")
          if(transf !=1){letra1$respo=tapply(response,fat[,i],mean, na.rm=TRUE)[rownames(letra1)]}}

        if(mcomp=="sk"){
          letra1 <- sk(resp, fat[, i],num(tab[3*i,1]),
                              num(tab[3*i,2]), alpha.t)
          colnames(letra1)=c("resp","groups")
          if(transf !=1){letra1$respo=tapply(response,fat[,i],mean, na.rm=TRUE)[rownames(letra1)]}}
        print(letra1)
        ordem=unique(as.vector(unlist(fat[i])))
        dadosm=data.frame(letra1,
                          desvio=tapply(response, c(fat[i]), sd, na.rm=TRUE)[rownames(letra1)])
        dadosm$media=tapply(response, c(fat[i]), mean, na.rm=TRUE)[rownames(letra1)]
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
                   fill=fill,
                   color=1)}
        grafico=grafico+
          theme+ylim(ylim)+
          ylab(ylab)+
          xlab(xlab)
        if(errorbar==TRUE){grafico=grafico+
          geom_text(aes(y=media+sup+
                          if(sup<0){-desvio}else{desvio},
                        label=letra),family=family,angle=angle.label, hjust=hjust)}
        if(errorbar==FALSE){grafico=grafico+
          geom_text(aes(y=media+sup,
                        label=letra),family=family,angle=angle.label, hjust=hjust)}
        if(errorbar==TRUE){grafico=grafico+
          geom_errorbar(data=dadosm,
                        aes(ymin=media-desvio,
                            ymax=media+desvio,
                            color=1),
                        color="black",
                        width=0.3)}
        if(angle !=0){grafico=grafico+
          theme(axis.text.x=element_text(hjust = 1.01,
                                         angle = angle))}
        grafico=grafico+
          theme(text = element_text(size=textsize,color="black",family=family),
                axis.text = element_text(size=textsize,color="black",family=family),
                axis.title = element_text(size=textsize,color="black",family=family),
                legend.position = "none")}

        # ================================
        # grafico de segmentos
        # ================================
        if(geom=="point"){grafico=ggplot(dadosm,
                                         aes(x=trats,
                                             y=media))
        if(fill=="trat"){grafico=grafico+
          geom_point(aes(color=trats),size=4)}
        else{grafico=grafico+
          geom_point(aes(color=trats),color=fill,size=4)}
        grafico=grafico+
          theme+
          ylab(ylab)+
          xlab(xlab)
        if(errorbar==TRUE){grafico=grafico+
          geom_text(aes(y=media+sup+
                          if(sup<0){-desvio}else{desvio},
                        label=letra),family=family,angle=angle.label, hjust=hjust)}
        if(errorbar==FALSE){grafico=grafico+
          geom_text(aes(y=media+sup,label=letra),
                    family=family,angle=angle.label, hjust=hjust)}
        if(errorbar==TRUE){grafico=grafico+
          geom_errorbar(data=dadosm,
                        aes(ymin=media-desvio,
                            ymax=media+desvio,color=1),
                        color="black",width=0.3)}
        if(angle !=0){grafico=grafico+
          theme(axis.text.x=element_text(hjust = 1.01,
                                         angle = angle))}
        grafico=grafico+ylim(ylim)+
          theme(text = element_text(size=textsize,color="black",family=family),
                axis.text = element_text(size=textsize,color="black",family=family),
                axis.title = element_text(size=textsize,color="black",family=family),
                legend.position = "none")}
        grafico=grafico
        if(color=="gray"){grafico=grafico+scale_fill_grey()}
        # print(grafico)
        cat("\n\n")
      }


      # # Regression
      if(quali[i]==FALSE){
        # dose=as.numeric(as.character(as.vector(unlist(fat[i]))))
        dose=as.vector(unlist(fata[i]))
        grafico=polynomial(dose,
                           resp,
                           grau = grau, ylab = ylab,
                           xlab = xlab, posi = posi, point = point,
                           theme = theme, textsize = textsize,
                           family=family,DFres = num(tab[3*i,1]),
                           SSq=num(tab[3*i,2]))
        grafico=grafico[[1]]}
      graficos[[i+1]]=grafico
    }
      graficos[[1]]=residplot
    }

    if(as.numeric(tab[1,5])>=alpha.f && as.numeric(tab[4,5])<alpha.f){
      cat(green(bold("-----------------------------------------------------------------\n")))
      cat("Isolated factors 1 not significant")
      cat(green(bold("\n-----------------------------------------------------------------\n")))
      d1=data.frame(tapply(response,fator1,mean, na.rm=TRUE))
      colnames(d1)="Mean"
      print(d1)}
    if(as.numeric(tab[4,5])>=alpha.f && as.numeric(tab[1,5])<alpha.f){
      cat(green(bold("-----------------------------------------------------------------\n")))
      cat("Isolated factors 2 not significant")
      cat(green(bold("\n-----------------------------------------------------------------\n")))
      d1=data.frame(tapply(response,fator2,mean, na.rm=TRUE))
      colnames(d1)="Mean"
      print(d1)
      }
    if(as.numeric(tab[1,5])>=alpha.f && as.numeric(tab[4,5])>=alpha.f){
      cat(green(bold("-----------------------------------------------------------------\n")))
      cat("Isolated factors not significant")
      cat(green(bold("\n-----------------------------------------------------------------\n")))
      print(tapply(response,list(fator1,fator2),mean, na.rm=TRUE))}
  }

  #-------------------------------------
  # Desdobramento de F1 dentro de F2
  #-------------------------------------
  if (as.numeric(tab[5, 5]) <= alpha.f) {
    cat(green(bold("-----------------------------------------------------------------\n")))
    cat("Significant interaction: analyzing the interaction")
    cat(green(bold("\n-----------------------------------------------------------------\n")))
    cat(green(bold("Analyzing ", fac[1], " inside of each level of ", fac[2])))
    cat(green(bold("\n-----------------------------------------------------------------\n")))
    l2 <- names(summary(Fator2))
    sq <- numeric(0)

    for (k in 1:nv2) {
      soma <- numeric(0)
      for (j in 1:nv1) {
        sub <- resp[Fator1 == levels(Fator1)[j] & Fator2 == levels(Fator2)[k]]
        q.som <- length(sub)
        soma <- c(soma, sum(sub))}
      sq <- c(sq, sum(soma^2)/q.som - sum(soma)^2/(q.som * length(soma)))}
    gl.sattert <- (num(tab[3,3])+(nv2-1)*num(tab[6,3]))^2/((num(tab[3,3])^2/num(tab[3,1]))+(((nv2-1)*num(tab[6,3]))^2/num(tab[6,1])))
    gl.f1f2 <- c(rep(nv1 - 1, nv2), gl.sattert)
    sq <- c(sq, NA)
    qm.f1f2 <- sq[1:nv2]/gl.f1f2[1:nv2]
    qm.ecomb <- (num(tab[3,3])+(nv2-1)*num(tab[6,3]))/nv2
    qm.f1f2 <- c(qm.f1f2, qm.ecomb)
    fc.f1f2 <- c(qm.f1f2[1:nv2]/qm.f1f2[nv2 + 1], NA)
    p.f1f2 <- c(1 - pf(fc.f1f2, gl.f1f2, gl.sattert))
    tab.f1f2 <- data.frame(GL = gl.f1f2, SQ = sq, QM = qm.f1f2, Fc = fc.f1f2, "p-value" = p.f1f2)
    nome.f1f2 <- numeric(0)
    for (j in 1:nv2) {nome.f1f2 <- c(nome.f1f2, paste(fac[1], " : ", fac[2], " ", l2[j], " ", sep = ""))}
    nome.f1f2 <- c(nome.f1f2, "Combined error")
    rownames(tab.f1f2) <- nome.f1f2
    tab.f1f2 <- round(tab.f1f2, 6)
    tab.f1f2[nv2 + 1, 2] <- tab.f1f2[nv2 + 1, 3] * tab.f1f2[nv2 + 1, 1]
    tab.f1f2[nv2 + 1, 5] <- tab.f1f2[nv2 + 1, 4] <- ""
    print(tab.f1f2)
    desdobramento1=tab.f1f2
    #-------------------------------------
    # Teste de Tukey
    #-------------------------------------
    if(quali[1]==TRUE & quali[2]==TRUE){
      if (mcomp == "tukey"){
        tukeygrafico=c()
        ordem=c()
        for (i in 1:nv2) {
          tukey=TUKEY(resp[fat[,2]==l2[i]], fat[,1][fat[,2]==l2[i]], num(tab.f1f2[nv2+1,1]),
                         num(tab.f1f2[nv2+1,2])/num(tab.f1f2[nv2+1,1]), alpha.t)
          colnames(tukey$groups)=c("resp","groups")
          tukeygrafico[[i]]=tukey$groups[as.character(unique(fat[,1][fat[,2]==l2[i]])),2]
          ordem[[i]]=rownames(tukey$groups[as.character(unique(fat[,1][fat[,2]==l2[i]])),])
          }
        letra=unlist(tukeygrafico)
        datag=data.frame(letra,ordem=unlist(ordem))
        datag=datag[order(factor(datag$ordem,levels=unique(Fator1))),]
        letra=datag$letra}


      if (mcomp == "duncan"){
        duncangrafico=c()
        ordem=c()
        for (i in 1:nv2) {
          duncan=duncan(resp[fat[,2]==l2[i]], fat[,1][fat[,2]==l2[i]], num(tab.f1f2[nv2+1,1]),
                             num(tab.f1f2[nv2+1,2])/num(tab.f1f2[nv2+1,1]), alpha.t)
          colnames(duncan$groups)=c("resp","groups")
          duncangrafico[[i]]=duncan$groups[as.character(unique(fat[,1][fat[,2]==l2[i]])),2]
          ordem[[i]]=rownames(duncan$groups[as.character(unique(fat[,1][fat[,2]==l2[i]])),])
          }
        letra=unlist(duncangrafico)
        datag=data.frame(letra,ordem=unlist(ordem))
        datag=datag[order(factor(datag$ordem,levels=unique(Fator1))),]
        letra=datag$letra}

      if (mcomp == "lsd"){
        lsdgrafico=c()
        ordem=c()
        for (i in 1:nv2) {
          lsd=lsd(resp[fat[,2]==l2[i]], fat[,1][fat[,2]==l2[i]], num(tab.f1f2[nv2+1,1]),
                             num(tab.f1f2[nv2+1,2])/num(tab.f1f2[nv2+1,1]), alpha.t)
          colnames(lsd$groups)=c("resp","groups")
          lsdgrafico[[i]]=lsd$groups[as.character(unique(fat[,1][fat[,2]==l2[i]])),2]
          ordem[[i]]=rownames(lsd$groups[as.character(unique(fat[,1][fat[,2]==l2[i]])),])
          }
        letra=unlist(lsdgrafico)
        datag=data.frame(letra,ordem=unlist(ordem))
        datag=datag[order(factor(datag$ordem,levels=unique(Fator1))),]
        letra=datag$letra}

      if (mcomp == "sk"){
        skgrafico=c()
        ordem=c()
        for (i in 1:nv2) {
          respi=resp[fat[,2]==l2[i]]
          trati=fat[,1][fat[,2]==l2[i]]
          # trati=fatores[, 1][Fator2 == lf2[i]]
          trati=factor(trati,levels = unique(trati))
          # respi=resp[Fator2 == lf2[i]]
          sk=sk(respi,trati,
                       num(tab.f1f2[nv2+1,1]),
                       num(tab.f1f2[nv2+1,2]),alpha.t)
          if(transf !="1"){sk$respo=tapply(response[Fator2 == lf2[i]],
                                           trati,mean, na.rm=TRUE)[rownames(sk$groups)]}
          skgrafico[[i]]=sk[levels(trati),2]
          ordem[[i]]=rownames(sk[levels(trati),])
        }
        letra=unlist(skgrafico)
        datag=data.frame(letra,ordem=unlist(ordem))
        datag$ordem=factor(datag$ordem,levels = unique(datag$ordem))
        datag=datag[order(datag$ordem),]
        letra=datag$letra}

      }

    #-------------------------------------
    # Desdobramento de F2 dentro de F1
    #-------------------------------------
    cat(green(bold("\n-----------------------------------------------------------------\n")))
    cat(green(bold("Analyzing ", fac[2], " inside of the level of ",fac[1])))
    cat(green(bold("\n-----------------------------------------------------------------\n")))
    l1 <- names(summary(Fator1))
    sq <- numeric(0)
    for (k in 1:nv1) {
      soma <- numeric(0)
      for (j in 1:nv2) {
        parc <- resp[Fator1 == levels(Fator1)[k] & Fator2 == levels(Fator2)[j]]
        q.som <- length(parc)
        soma <- c(soma, sum(parc))}
      sq <- c(sq, sum(soma^2)/q.som-sum(soma)^2/(q.som*length(soma)))}

    gl.f2f1 <- c(rep(nv2 - 1, nv1), tab[6, 1])
    sq <- c(sq, as.numeric(tab[6, 2]))
    qm.f2f1 <- sq/gl.f2f1
    fc.f2f1 <- c(qm.f2f1[1:nv1]/num(tab[6, 3]), NA)
    p.f2f1 <- c(1 - pf(fc.f2f1, gl.f2f1, num(tab[6,1])))
    tab.f2f1 <- data.frame(GL=gl.f2f1, SQ=sq, QM=qm.f2f1, Fc=fc.f2f1, "p-value"=p.f2f1)
    nome.f2f1 <- numeric(0)
    for (j in 1:nv1) {nome.f2f1 <- c(nome.f2f1, paste(fac[2], " : ",fac[1], " ", l1[j], " ", sep = ""))}
    nome.f2f1 <- c(nome.f2f1, "Error b")
    rownames(tab.f2f1) <- nome.f2f1
    tab.f2f1 <- round(tab.f2f1, 6)
    tab.f2f1[nv1 + 1, 5] <- tab.f2f1[nv1 + 1, 4] <- ""
    print(tab.f2f1)
    desdobramento2=tab.f2f1
    #-------------------------------------
    # Teste de Tukey
    #-------------------------------------
    if(quali[1]==TRUE && quali[2]==TRUE){
      if (mcomp == "tukey"){
        tukeygrafico1=c()
        for (i in 1:nv1) {
          tukey=TUKEY(resp[fat[, 1] == l1[i]],
                                    fat[,2][fat[, 1] == l1[i]],
                                    num(tab.f2f1[nv1 +1, 1]),
                                    num(tab.f2f1[nv1 + 1, 2])/num(tab.f2f1[nv1 +1, 1]),alpha.t)
          colnames(tukey$groups)=c("resp","groups")
          tukeygrafico1[[i]]=tukey$groups[as.character(unique(fat[,2][fat[, 1] == l1[i]])),2]
          if(transf !="1"){tukey$groups$respo=tapply(response[fat[, 1] == l1[i]],
                                                     fat[,2][fat[, 1] == l1[i]],mean, na.rm=TRUE)[rownames(tukey$groups)]}
          }
        letra1=unlist(tukeygrafico1)
        letra1=toupper(letra1)}

      if (mcomp == "duncan"){
        duncangrafico1=c()
        for (i in 1:nv1) {
          duncan=duncan(resp[fat[, 1] == l1[i]], fat[,2][fat[, 1] == l1[i]], num(tab.f2f1[nv1 +1, 1]),
                             num(tab.f2f1[nv1 + 1, 2])/num(tab.f2f1[nv1 +1, 1]),alpha.t)
          colnames(duncan$groups)=c("resp","groups")
          duncangrafico1[[i]]=duncan$groups[levels(fat[,2][fat[, 1] == l1[i]]),2]
          if(transf !="1"){duncan$groups$respo=tapply(response[fat[, 1] == l1[i]],
                                                     fat[,2][fat[, 1] == l1[i]],mean, na.rm=TRUE)[rownames(duncan$groups)]}
          }
        letra1=unlist(duncangrafico1)
        letra1=toupper(letra1)}

      if (mcomp == "lsd"){
        lsdgrafico1=c()
        for (i in 1:nv1) {
          lsd=LSD(resp[fat[, 1] == l1[i]], fat[,2][fat[, 1] == l1[i]], num(tab.f2f1[nv1 +1, 1]),
                             num(tab.f2f1[nv1 + 1, 2])/num(tab.f2f1[nv1 +1, 1]),alpha.t)
          colnames(lsd$groups)=c("resp","groups")
          duncangrafico1[[i]]=lsd$groups[levels(fat[,2][fat[, 1] == l1[i]]),2]
          if(transf !="1"){lsd$groups$respo=tapply(response[fat[, 1] == l1[i]],
                                                      fat[,2][fat[, 1] == l1[i]],mean, na.rm=TRUE)[rownames(lsd$groups)]}
          }
        letra1=unlist(lsdgrafico1)
        letra1=toupper(letra1)}
      if (mcomp == "sk"){
        skgrafico1=c()
        for (i in 1:nv1) {
          respi=resp[fat[, 1] == l1[i]]
          trati=fat[,2][fat[, 1] == l1[i]]
          trati=factor(trati,unique(trati))
          sk=sk(respi,trati,
                       num(tab.f2f1[nv1 +1, 1]),
                       num(tab.f2f1[nv1 + 1, 2]),alpha.t)
          if(transf !=1){sk$respo=tapply(response[Fator1 == lf1[i]],trati,
                                         mean, na.rm=TRUE)[rownames(sk)]}
          skgrafico1[[i]]=sk[levels(trati),2]
        }
        letra1=unlist(skgrafico1)
        letra1=toupper(letra1)}
    }


    # -----------------------------
    # Gráfico de colunas
    #------------------------------
    if(quali[1]==TRUE & quali[2]==TRUE){
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
        graph=graph[paste(rep(unique(Fator1),e=length(colnames(media))),
                          rep(unique(Fator2),length(rownames(media)))),]
        graph$letra=letra
        graph$letra1=letra1
        graph$f1=factor(graph$f1,levels = unique(Fator1))
        graph$f2=factor(graph$f2,levels = unique(Fator2))
        if(addmean==TRUE){graph$numero=paste(format(graph$media,digits = dec), graph$letra,graph$letra1, sep="")}
        if(addmean==FALSE){graph$numero=paste(graph$letra,graph$letra1)}
        f1=graph$f1
        f2=graph$f2
        media=graph$media
        desvio=graph$desvio
        numero=graph$numero

        colint=ggplot(graph, aes(x=f1,
                                 y=media,
                                 fill=f2))+
          geom_col(position = "dodge",color="black")+
          ylab(ylab)+xlab(xlab)+ylim(ylim)+
          theme
        if(errorbar==TRUE){colint=colint+
          geom_errorbar(data=graph,
                        aes(ymin=media-desvio,
                            ymax=media+desvio),
                        width=0.3,
                        position = position_dodge(width=0.9))}
        if(errorbar==TRUE){colint=colint+
          geom_text(aes(y=media+sup+
                          if(sup<0){-desvio}else{desvio},
                        label=numero),
                    position = position_dodge(width=0.9),angle=angle.label, hjust=hjust)}
        if(errorbar==FALSE){colint=colint+
          geom_text(aes(y=media+sup,label=numero),
                    position = position_dodge(width=0.9),angle=angle.label, hjust=hjust)}
        colint=colint+theme(text=element_text(size=12),
                            legend.position = posi,
                            axis.text = element_text(size=12,
                                                     color="black"),
                            axis.title = element_text(size=12,
                                                      color="black"))+
          labs(fill=legend)
        if(angle !=0){colint=colint+theme(axis.text.x=element_text(hjust = 1.01,angle = angle))}
        if(color=="gray"){colint=colint+scale_fill_grey()}
        print(colint)
        grafico=colint
        letras=paste(graph$letra,graph$letra1,sep="")
        matriz=data.frame(t(matrix(paste(format(graph$media,digits = dec),letras),ncol = length(levels(Fator1)))))
        rownames(matriz)=levels(Fator1)
        colnames(matriz)=levels(Fator2)
        cat(green(bold("\n-----------------------------------------------------------------\n")))
        cat(green(bold("Final table")))
        cat(green(bold("\n-----------------------------------------------------------------\n")))
        print(matriz)
        message(black("\n\nAverages followed by the same lowercase letter in the column and uppercase in the row do not differ by the",mcomp,"(p<",alpha.t,")"))

      }


    if(quali[1]==FALSE  && color=="gray"| quali[2]==FALSE && color=="gray"){
      if(quali[2]==FALSE){
        if (mcomp == "tukey"){
          for (i in 1:nv2) {
            tukey=TUKEY(resp[fat[,2]==l2[i]], fat[,1][fat[,2]==l2[i]], num(tab.f1f2[nv2+1,1]),
                           num(tab.f1f2[nv2+1,2])/num(tab.f1f2[nv2+1,1]), alpha.t)
            colnames(tukey$groups)=c("resp","groups")
            cat("\n----------------------\n")
            cat("Multiple comparison of F1 within level",lf2[i],"of F2")
            cat("\n----------------------\n")
            print(tukey$groups)}}
        if (mcomp == "duncan"){
        for (i in 1:nv2) {
            duncan=duncan(resp[fat[,2]==l2[i]], fat[,1][fat[,2]==l2[i]], num(tab.f1f2[nv2+1,1]),
                               num(tab.f1f2[nv2+1,2])/num(tab.f1f2[nv2+1,1]), alpha.t)
            colnames(duncan$groups)=c("resp","groups")
            cat("\n----------------------\n")
            cat("Multiple comparison of F1 within level",lf2[i],"of F2")
            cat("\n----------------------\n")
            print(duncan$groups)

          }}
        if (mcomp == "lsd"){
          for (i in 1:nv2) {
            lsd=lsd(resp[fat[,2]==l2[i]], fat[,1][fat[,2]==l2[i]], num(tab.f1f2[nv2+1,1]),
                    num(tab.f1f2[nv2+1,2])/num(tab.f1f2[nv2+1,1]), alpha.t)
            colnames(lsd$groups)=c("resp","groups")
            cat("\n----------------------\n")
            cat("Multiple comparison of F1 within level",lf2[i],"of F2")
            cat("\n----------------------\n")
            print(lsd$groups)

          }}
        if (mcomp == "sk"){
          for (i in 1:nv2) {
            respi=resp[fat[,2]==l2[i]]
            trati=fat[,1][fat[,2]==l2[i]]
            trati=factor(trati,levels = unique(trati))
            sk=sk(respi,trati,
                         num(tab.f1f2[nv2+1,1]),
                         num(tab.f1f2[nv2+1,2]),alpha.t)
            if(transf !="1"){sk$respo=tapply(response[Fator2 == lf2[i]],
                                             trati,mean, na.rm=TRUE)[rownames(sk$groups)]}
            cat("\n----------------------\n")
            cat("Multiple comparison of F1 within level",lf2[i],"of F2")
            cat("\n----------------------\n")
            print(sk)}}
      }
      if(quali[2]==FALSE){
        fator2a=fator2a#as.numeric(as.character(fator2))
        grafico=polynomial2(fator2a,
                            resp,
                            fator1,
                            grau = grau,
                            ylab=ylab,
                            xlab=xlab,
                            theme=theme,
                            point=point,
                            posi= posi,
                            ylim=ylim,
                            textsize=textsize,
                            family=family,
                            DFres = num(tab.f1f2[nv2+1,1]),
                            SSq = num(tab.f1f2[nv2+1,2]))
        if(quali[1]==FALSE & quali[2]==FALSE){
          graf=list(grafico,NA)}
      }
      if(quali[1]==FALSE){
        if (mcomp == "tukey"){
          for (i in 1:nv1) {
            tukey=TUKEY(resp[fat[, 1] == l1[i]],
                           fat[,2][fat[, 1] == l1[i]],
                           num(tab.f2f1[nv1 +1, 1]),
                           num(tab.f2f1[nv1 + 1, 2])/num(tab.f2f1[nv1 +1, 1]),alpha.t)
            colnames(tukey$groups)=c("resp","groups")
            if(transf !="1"){tukey$groups$respo=tapply(response[fat[, 1] == l1[i]],
                                                       fat[,2][fat[, 1] == l1[i]],mean, na.rm=TRUE)[rownames(tukey$groups)]}
            cat("\n----------------------\n")
            cat("Multiple comparison of F2 within level",lf1[i],"of F1")
            cat("\n----------------------\n")
            print(tukey$groups)
          }}
        if (mcomp == "duncan"){
          for (i in 1:nv1) {
            duncan=duncan(resp[fat[, 1] == l1[i]], fat[,2][fat[, 1] == l1[i]], num(tab.f2f1[nv1 +1, 1]),
                               num(tab.f2f1[nv1 + 1, 2])/num(tab.f2f1[nv1 +1, 1]),alpha.t)
            colnames(duncan$groups)=c("resp","groups")
            if(transf !="1"){duncan$groups$respo=tapply(response[fat[, 1] == l1[i]],
                                                        fat[,2][fat[, 1] == l1[i]],mean, na.rm=TRUE)[rownames(duncan$groups)]}
            cat("\n----------------------\n")
            cat("Multiple comparison of F2 within level",lf1[i],"of F1")
            cat("\n----------------------\n")
            print(duncan$groups)
          }}
        if (mcomp == "lsd"){
          for (i in 1:nv1) {
            lsd=LSD(resp[fat[, 1] == l1[i]], fat[,2][fat[, 1] == l1[i]],
                    num(tab.f2f1[nv1 +1, 1]),
                         num(tab.f2f1[nv1 + 1, 2])/num(tab.f2f1[nv1 +1, 1]),alpha.t)
            colnames(lsd$groups)=c("resp","groups")
            if(transf !="1"){lsd$groups$respo=tapply(response[fat[, 1] == l1[i]],
                                                     fat[,2][fat[, 1] == l1[i]],mean, na.rm=TRUE)[rownames(lsd$groups)]}
            cat("\n----------------------\n")
            cat("Multiple comparison of F2 within level",lf1[i],"of F1")
            cat("\n----------------------\n")
            print(lsd$groups)

            }}
        if (mcomp == "sk"){
          skgrafico1=c()
          for (i in 1:nv1) {
            respi=resp[fat[, 1] == l1[i]]
            trati=fat[,2][fat[, 1] == l1[i]]
            trati=factor(trati,unique(trati))
            sk=sk(respi,trati,
                         num(tab.f2f1[nv1 +1, 1]),
                         num(tab.f2f1[nv1 + 1, 2]),alpha.t)
            if(transf !=1){sk$respo=tapply(response[Fator1 == lf1[i]],trati,
                                           mean, na.rm=TRUE)[rownames(sk)]}
            cat("\n----------------------\n")
            cat("Multiple comparison of F2 within level",lf1[i],"of F1")
            cat("\n----------------------\n")
            print(sk)
          }}
      }
      if(quali[1]==FALSE){
        fator1a=fator1a#as.numeric(as.character(fator1))
        grafico=polynomial2(fator1a,
                            resp,
                            fator2,
                            grau = grau,
                            ylab=ylab,
                            xlab=xlab,
                            theme=theme,
                            point=point,
                            posi = posi,
                            ylim=ylim,
                            textsize=textsize,
                            family=family,
                            DFres = num(tab.f2f1[nv1 +1, 1]),
                            SSq = num(tab.f2f1[nv1 + 1, 2]))
        if(quali[1]==FALSE & quali[2]==FALSE){
          graf[[2]]=grafico
          grafico=graf}
      }
    }
    if(quali[1]==FALSE  && color=="rainbow"| quali[2]==FALSE && color=="rainbow"){
      if(quali[2]==FALSE){
        if (mcomp == "tukey"){
          for (i in 1:nv2) {
            tukey=TUKEY(resp[fat[,2]==l2[i]], fat[,1][fat[,2]==l2[i]], num(tab.f1f2[nv2+1,1]),
                           num(tab.f1f2[nv2+1,2])/num(tab.f1f2[nv2+1,1]), alpha.t)
            colnames(tukey$groups)=c("resp","groups")
            cat("\n----------------------\n")
            cat("Multiple comparison of F1 within level",lf2[i],"of F2")
            cat("\n----------------------\n")
            print(tukey$groups)}}
        if (mcomp == "duncan"){
          for (i in 1:nv2) {
            duncan=duncan(resp[fat[,2]==l2[i]], fat[,1][fat[,2]==l2[i]], num(tab.f1f2[nv2+1,1]),
                               num(tab.f1f2[nv2+1,2])/num(tab.f1f2[nv2+1,1]), alpha.t)
            colnames(duncan$groups)=c("resp","groups")
            cat("\n----------------------\n")
            cat("Multiple comparison of F1 within level",lf2[i],"of F2")
            cat("\n----------------------\n")
            print(duncan$groups)

          }}
        if (mcomp == "lsd"){
          for (i in 1:nv2) {
            lsd=lsd(resp[fat[,2]==l2[i]], fat[,1][fat[,2]==l2[i]], num(tab.f1f2[nv2+1,1]),
                    num(tab.f1f2[nv2+1,2])/num(tab.f1f2[nv2+1,1]), alpha.t)
            colnames(lsd$groups)=c("resp","groups")
            cat("\n----------------------\n")
            cat("Multiple comparison of F1 within level",lf2[i],"of F2")
            cat("\n----------------------\n")
            print(lsd$groups)

          }}
        if (mcomp == "sk"){
          for (i in 1:nv2) {
            respi=resp[fat[,2]==l2[i]]
            trati=fat[,1][fat[,2]==l2[i]]
            # trati=fatores[, 1][Fator2 == lf2[i]]
            trati=factor(trati,levels = unique(trati))
            # respi=resp[Fator2 == lf2[i]]
            sk=sk(respi,trati,
                         num(tab.f1f2[nv2+1,1]),
                         num(tab.f1f2[nv2+1,2]),alpha.t)
            if(transf !="1"){sk$respo=tapply(response[Fator2 == lf2[i]],
                                             trati,mean, na.rm=TRUE)[rownames(sk$groups)]}
            cat("\n----------------------\n")
            cat("Multiple comparison of F1 within level",lf2[i],"of F2")
            cat("\n----------------------\n")
            print(sk)}}
      }
      if(quali[2]==FALSE){
        fator2a=fator2a#as.numeric(as.character(fator2))
        grafico=polynomial2_color(fator2a,
                                  resp,
                                  fator1,
                                  grau = grau,
                                  ylab=ylab,
                                  xlab=xlab,
                                  theme=theme,
                                  point=point,
                                  posi=posi,
                                  ylim=ylim,
                                  textsize=textsize,
                                  family=family,
                                  DFres = num(tab.f1f2[nv2+1,1]),
                                  SSq = num(tab.f1f2[nv2+1,2]))
        if(quali[1]==FALSE & quali[2]==FALSE){
          graf=list(grafico,NA)}
      }
      if(quali[1]==FALSE){
        if (mcomp == "tukey"){
          for (i in 1:nv1) {
            tukey=TUKEY(resp[fat[, 1] == l1[i]],
                           fat[,2][fat[, 1] == l1[i]],
                           num(tab.f2f1[nv1 +1, 1]),
                           num(tab.f2f1[nv1 + 1, 2])/num(tab.f2f1[nv1 +1, 1]),alpha.t)
            colnames(tukey$groups)=c("resp","groups")
            if(transf !="1"){tukey$groups$respo=tapply(response[fat[, 1] == l1[i]],
                                                       fat[,2][fat[, 1] == l1[i]],mean, na.rm=TRUE)[rownames(tukey$groups)]}
            cat("\n----------------------\n")
            cat("Multiple comparison of F2 within level",lf1[i],"of F1")
            cat("\n----------------------\n")
            print(tukey$groups)
          }}
        if (mcomp == "duncan"){
          for (i in 1:nv1) {
            duncan=duncan(resp[fat[, 1] == l1[i]], fat[,2][fat[, 1] == l1[i]], num(tab.f2f1[nv1 +1, 1]),
                               num(tab.f2f1[nv1 + 1, 2])/num(tab.f2f1[nv1 +1, 1]),alpha.t)
            colnames(duncan$groups)=c("resp","groups")
            if(transf !="1"){duncan$groups$respo=tapply(response[fat[, 1] == l1[i]],
                                                        fat[,2][fat[, 1] == l1[i]],mean, na.rm=TRUE)[rownames(duncan$groups)]}
            cat("\n----------------------\n")
            cat("Multiple comparison of F2 within level",lf1[i],"of F1")
            cat("\n----------------------\n")
            print(duncan$groups)
          }}
        if (mcomp == "lsd"){
          for (i in 1:nv1) {
            lsd=LSD(resp[fat[, 1] == l1[i]], fat[,2][fat[, 1] == l1[i]], num(tab.f2f1[nv1 +1, 1]),
                         num(tab.f2f1[nv1 + 1, 2])/num(tab.f2f1[nv1 +1, 1]),alpha.t)
            colnames(lsd$groups)=c("resp","groups")
            if(transf !="1"){lsd$groups$respo=tapply(response[fat[, 1] == l1[i]],
                                                     fat[,2][fat[, 1] == l1[i]],mean, na.rm=TRUE)[rownames(lsd$groups)]}
            cat("\n----------------------\n")
            cat("Multiple comparison of F2 within level",lf1[i],"of F1")
            cat("\n----------------------\n")
            print(lsd$groups)

          }}
        if (mcomp == "sk"){
          skgrafico1=c()
          for (i in 1:nv1) {
            respi=resp[fat[, 1] == l1[i]]
            trati=fat[,2][fat[, 1] == l1[i]]
            trati=factor(trati,unique(trati))
            sk=sk(respi,trati,
                         num(tab.f2f1[nv1 +1, 1]),
                         num(tab.f2f1[nv1 + 1, 2]),alpha.t)
            if(transf !=1){sk$respo=tapply(response[Fator1 == lf1[i]],trati,
                                           mean, na.rm=TRUE)[rownames(sk)]}
            cat("\n----------------------\n")
            cat("Multiple comparison of F2 within level",lf1[i],"of F1")
            cat("\n----------------------\n")
            print(sk)
          }}
      }
      if(quali[1]==FALSE){
        fator1a=fator1a#as.numeric(as.character(fator1))
        grafico=polynomial2_color(fator1a,resp,fator2, grau = grau,
                                  ylab=ylab,xlab=xlab,
                                  theme=theme,point=point,
                                  posi = posi,ylim=ylim,
                                  textsize=textsize,
                                  family=family,
                                  DFres = num(tab.f2f1[nv1 +1, 1]),
                                  SSq = num(tab.f2f1[nv1 + 1, 2]))
        if(quali[1]==FALSE & quali[2]==FALSE){
          graf[[2]]=grafico
          grafico=graf}
      }
    }
  }
  if(as.numeric(tab[5, 5])>alpha.f){
    names(graficos)=c("residplot","graph1","graph2")
    graficos}else{colints=list(residplot,grafico)}
}


