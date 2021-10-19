#' Analysis: Joint analysis of experiments in randomized block design
#'
#' @description Function of the AgroR package for joint analysis of experiments conducted in a randomized qualitative or quantitative single-block design with balanced data.
#' @author Gabriel Danilo Shimizu, \email{shimizu@uel.br}
#' @author Leandro Simoes Azeredo Goncalves
#' @author Rodrigo Yudi Palhaci Marubayashi
#' @param trat Numerical or complex vector with treatments
#' @param block Numerical or complex vector with blocks
#' @param local Numeric or complex vector with locations or times
#' @param response Numerical vector containing the response of the experiment.
#' @param transf Applies data transformation (default is 1; for log consider 0)
#' @param norm Error normality test (\emph{default} is Shapiro-Wilk)
#' @param homog Homogeneity test of variances (\emph{default} is Bartlett)
#' @param mcomp Multiple comparison test (Tukey (\emph{default}), LSD, Scott-Knott and Duncan)
#' @param quali Defines whether the factor is quantitative or qualitative (\emph{default} is qualitative)
#' @param alpha.f Level of significance of the F test (\emph{default} is 0.05)
#' @param alpha.t Significance level of the multiple comparison test (\emph{default} is 0.05)
#' @param grau Degree of polynomial in case of quantitative factor (\emph{default} is 1)
#' @param ylab Variable response name (Accepts the \emph{expression}() function)
#' @param xlab Treatments name (Accepts the \emph{expression}() function)
#' @param title Graph title
#' @param theme ggplot2 theme (\emph{default} is theme_classic())
#' @param dec Number of cells
#' @param fill Defines chart color (to generate different colors for different treatments, define fill = "trat")
#' @param angulo x-axis scale text rotation
#' @param textsize Font size
#' @param family Font family
#' @param errorbar Plot the standard deviation bar on the graph (In the case of a segment and column graph) - \emph{default} is TRUE
#' @note The ordering of the graph is according to the sequence in which the factor levels are arranged in the data sheet. The bars of the column and segment graphs are standard deviation.
#' @note In the final output when transformation (transf argument) is different from 1, the columns resp and respo in the mean test are returned, indicating transformed and non-transformed mean, respectively.
#' @return Returns the assumptions of the analysis of variance, the assumption of the joint analysis by means of a QMres ratio matrix, the analysis of variance, the multiple comparison test or regression.
#' @references
#'
#' Ferreira, P. V. Estatistica experimental aplicada a agronomia. Edufal, 2018.
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
#' @keywords DBC
#' @keywords Joint Analysis
#' @export
#' @examples
#' library(AgroR)
#' data(mirtilo)
#' with(mirtilo, conjdbc(trat, bloco, exp, resp))

conjdbc=function(trat,
                 block,
                 local,
                 response,
                 transf=1,
                 norm="sw",
                 homog="bt",
                 theme=theme_classic(),
                 mcomp="tukey",
                 quali=TRUE,
                 alpha.f=0.05,
                 alpha.t=0.05,
                 grau=NA,
                 ylab="response",
                 title="",
                 xlab="",
                 fill="lightblue",
                 angulo=0,
                 textsize=12,
                 dec=3,
                 family="sans",
                 errorbar=TRUE){
  sup=0.2*mean(response, na.rm=TRUE)
  requireNamespace("crayon")
  requireNamespace("ggplot2")

  if(transf==1){resp=response}else{resp=(response^transf-1)/transf}
  if(transf==0){resp=log(response)}
  if(transf==0.5){resp=sqrt(response)}
  if(transf==-0.5){resp=1/sqrt(response)}
  if(transf==-1){resp=1/response}
  tratnum=trat
  tratamento=factor(trat,levels=unique(trat))
  bloco=as.factor(block)
  local=as.factor(local)

  a = anova(aov(resp ~ local + local:bloco + tratamento + local:tratamento))[c(4:5), ]
  b = summary(aov(resp ~ bloco+local + local:bloco + tratamento +
                    Error(local:(bloco + tratamento))))
  c = aov(resp ~ local + local:bloco + tratamento + local:tratamento)
  dados=data.frame(resp,response,tratamento,local,bloco,tratnum)
  anova=c()
  tukey=c()
  graficos=list()
  qmres=data.frame(QM=NA)

  for(i in 1:length(levels(local))){
    anova[[i]]=anova(aov(resp~tratamento+bloco,
                         data=dados[dados$local==levels(dados$local)[i],]))
    qm=anova[[i]]$`Mean Sq`[3]
    qmres[i,]=c(qm)
    qmres=as.vector(qmres)
    names(anova)[i]=levels(local)[i]
    aov1=aov(resp~tratamento+bloco, data=dados[dados$local==levels(dados$local)[i],])}
  matriza=matrix(rep(qmres[[1]],e=length(qmres[[1]])),
           ncol=length(qmres[[1]]))
  matrizb=matrix(rep(qmres[[1]],length(qmres[[1]])),
           ncol=length(qmres[[1]]))
  ratio=matriza/matrizb
  rownames(ratio)=levels(local)
  colnames(ratio)=levels(local)
  razao=data.frame(resp1=c(ratio),
                   var1=rep(rownames(ratio),e=length(rownames(ratio))),
                   var2=rep(colnames(ratio),length(colnames(ratio))))
  var1=razao$var1
  var2=razao$var2
  resp1=razao$resp1

  ratioplot=ggplot(razao,
                   aes(x=var2,
                       y=var1,
                       fill=resp1))+
    geom_tile(color="gray50",size=1)+
    scale_x_discrete(position = "top")+
    scale_fill_distiller(palette = "RdBu",direction = 1)+
    ylab("Numerator")+
    xlab("Denominator")+
    geom_label(aes(label=format(resp1,digits=2)),fill="white")+
    labs(fill="ratio")+
    theme(axis.text = element_text(size=12,color="black"),
          legend.text = element_text(size=12),
          axis.ticks = element_blank(),
          panel.background = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank())+
    labs(caption = "The ratio must be less than 7 (Ferreira et al., 2018)",
         title="Matrix of average square of the residue")
  print(ratioplot)
  QMRES=as.vector(qmres$QM)
  qmresmedio=max(QMRES)/min(QMRES)
  b1=matrix(unlist(b$`Error: local:tratamento`),
            ncol=5,2)
  b2=matrix(c(unlist(b$`Error: local:bloco`),NA,NA,NA,NA,NA,NA),
            ncol=5,3)[2:3,]
  datas=rbind(b1[1,],b2);colnames(datas)=colnames(a)
  datas=rbind(datas,a[1,])
  nexp=length(unique(local))
  ntrat=length(unique(trat))
  nrep=table(trat)/nexp
  GL=nexp*(ntrat*nrep[1]-(ntrat-1)-(nrep[1]-1))
  resmed=data.frame(rbind(c(GL,NA,mean(QMRES),NA,NA)))
  colnames(resmed)=colnames(datas)
  datas=rbind(datas,resmed)
  rownames(datas)=c("Trat","Exp","Block/Local","Exp:Trat","Average residue")
  datas[2,4]=datas[2,3]/datas[4,3]
  datas[3,4]=datas[3,3]/datas[4,3]
  datas[2,5]=1-pf(datas[2,4],datas[2,1],datas[4,1])
  datas[3,5]=1-pf(datas[3,4],datas[3,1],datas[4,1])
  datas[5,2]=datas[5,3]*datas[5,1]

  d=aov(resp~tratamento*local+bloco+local/bloco)
  if(norm=="sw"){norm1 = shapiro.test(d$res)}
  if(norm=="li"){norm1=nortest::lillie.test(d$residuals)}
  if(norm=="ad"){norm1=nortest::ad.test(d$residuals)}
  if(norm=="cvm"){norm1=nortest::cvm.test(d$residuals)}
  if(norm=="pearson"){norm1=nortest::pearson.test(d$residuals)}
  if(norm=="sf"){norm1=nortest::sf.test(d$residuals)}

  if(homog=="bt"){
    homog1 = bartlett.test(d$res ~ trat)
    statistic=homog1$statistic
    phomog=homog1$p.value
    method=paste("Bartlett test","(",names(statistic),")",sep="")
  }
  if(homog=="levene"){
    homog1 = levenehomog(d$res~trat)
    statistic=homog1$`F value`[1]
    phomog=homog1$`Pr(>F)`[1]
    method="Levene's Test (center = median)(F)"
    names(homog1)=c("Df", "F value","p.value")}

  indep = dwtest(d)

  cat(green(bold("\n-----------------------------------------------------------------\n")))
  cat(green(bold("Normality of errors")))
  cat(green(bold("\n-----------------------------------------------------------------\n")))
  normal=data.frame(Method=paste(norm1$method,"(",names(norm1$statistic),")",sep=""),
                    Statistic=norm1$statistic,
                    "p-value"=norm1$p.value)
  rownames(normal)=""
  print(normal)
  cat("\n")
  cat(if(norm1$p.value>0.05){black("As the calculated p-value is greater than the 5% significance level, hypothesis H0 is not rejected. Therefore, errors can be considered normal")}
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
  cat(if(homog1$p.value>0.05){black("As the calculated p-value is greater than the 5% significance level, hypothesis H0 is not rejected. Therefore, the variances can be considered homogeneous")}
      else {"As the calculated p-value is less than the 5% significance level, H0 is rejected. Therefore, the variances are not homogeneous\n"})
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
  cat(black(if(indep$p.value>0.05){black("As the calculated p-value is greater than the 5% significance level, hypothesis H0 is not rejected. Therefore, errors can be considered independent")}
      else {"As the calculated p-value is less than the 5% significance level, H0 is rejected. Therefore, errors are not independent"}))
  cat(green(bold("\n-----------------------------------------------------------------\n")))
  cat(green(bold("Test Homogeneity of experiments")))
  cat(green(bold("\n-----------------------------------------------------------------\n")))
  print(qmresmedio)
  cat("\nBased on the analysis of variance and homogeneity of experiments, it can be concluded that: ")
  if(qmresmedio<7 && a$`Pr(>F)`[1]>alpha.f){
    message(black("The experiments can be analyzed together"))}else{
      message("Experiments cannot be analyzed together (Separate by experiment)")}
  cat("\n\n")
  modres=anova(d)
  respad=d$res/sqrt(modres$`Mean Sq`[6])
  out=respad[respad>3 | respad<(-3)]
  out=names(out)
  out=if(length(out)==0)("No discrepant point")else{out}
  resids=d$res/sqrt(modres$`Mean Sq`[6])
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
  cat(green(bold("Analysis of variance")))
  cat(green(bold("\n-----------------------------------------------------------------\n")))
  #print(b$`Error: local:tratamento`)
  print(as.matrix(datas),na.print = "")
  cat(green(bold("\n-----------------------------------------------------------------\n")))

  anova=as.list(1:length(levels(local)))
  tukey=as.list(1:length(levels(local)))

  if(a$`Pr(>F)`[1] < alpha.f | qmresmedio > 7){
    if(quali==TRUE){
    for(i in 1:length(levels(local))){
        anova1=anova(aov(resp~tratamento+bloco,
                             data=dados[dados$local==levels(dados$local)[i],]))
        qm=anova1$`Mean Sq`[3]
        qmres[i,]=c(qm)
        qmres=as.vector(qmres)
        #names(anova1)[i]=levels(local)[i]
        aov1=aov(resp~tratamento+bloco, data=dados[dados$local==levels(dados$local)[i],])
        if(quali==TRUE){
          if(mcomp=="tukey"){tukey[[i]]=TUKEY(aov1,"tratamento",
                                                 alpha = alpha.t)$groups[unique(as.character(trat)),]
          comp=TUKEY(aov1,"tratamento")$groups}
          if(mcomp=="duncan"){tukey[[i]]=duncan(aov1,"tratamento",
                                                     alpha = alpha.t)$groups[unique(as.character(trat)),]
          comp=duncan(aov1,"tratamento")$groups}
          if(mcomp=="lsd"){tukey[[i]]=LSD(aov1,"tratamento",
                                               alpha = alpha.t)$groups[unique(as.character(trat)),]
          comp=LSD(aov1,"tratamento")$groups}
          if(mcomp=="sk"){
            anova=anova(aov1)
            data=dados[dados$local==levels(dados$local)[i],]
            nrep=table(data$trat)[1]
            medias=sort(tapply(data$resp,data$trat,mean),decreasing = TRUE)
            letra=scottknott(means = medias,
                             df1 = anova$Df[2],
                             nrep = nrep,
                             QME = anova$`Mean Sq`[2],
                             alpha = alpha.t)
            letra1=data.frame(resp=medias,groups=letra)[unique(as.character(trat)),]
            tukey[[i]]=letra1
            comp=letra1

            # anova=anova(aov1)
            # data=dados[dados$local==levels(dados$local)[i],]
            # tukey[[i]]=sk(data$resp,
            #                      data$tratamento,
            #                      anova$Df[3],
            #                      anova$`Sum Sq`[3],alpha = alpha.t)[unique(as.character(trat)),]
            # comp=sk(data$resp,data$tratamento,anova$Df[3],
            #                anova$`Sum Sq`[3],alpha = alpha.t)
            }

          if(transf=="1"){}else{tukey[[i]]$respo=with(dados[dados$local==levels(dados$local)[i],],
                                                      tapply(response, tratamento, mean, na.rm=TRUE))[rownames(comp)]}
          names(tukey)[i]=levels(local)[i]

          dadosm=data.frame(comp,
                            media=with(dados[dados$local==levels(dados$local)[i],],
                                       tapply(response, tratamento, mean, na.rm=TRUE))[rownames(comp)],
                            desvio=with(dados[dados$local==levels(dados$local)[i],],
                                        tapply(response, tratamento, sd, na.rm=TRUE))[rownames(comp)])
          dadosm$trats=factor(rownames(dadosm),unique(trat))
          dadosm$limite=dadosm$media+dadosm$desvio
          dadosm$letra=paste(format(dadosm$media,digits = dec),
                             dadosm$groups)
          dadosm=dadosm[unique(as.character(trat)),]
          trats=dadosm$trats
          limite=dadosm$limite
          letra=dadosm$letra
          media=dadosm$media
          desvio=dadosm$desvio
          grafico=ggplot(dadosm,
                         aes(x=trats,
                             y=media))+
            geom_col(aes(fill=trats),fill=fill,color=1)+
            theme_classic()+
            ylab(ylab)+
            xlab(xlab)+ylim(0,1.5*max(limite))+
            geom_errorbar(aes(ymin=media-desvio,
                              ymax=media+desvio),color="black", width=0.3)+
            geom_text(aes(y=media+desvio+sup,
                          label=letra))+
            theme(text = element_text(size=textsize,color="black", family = family),
                  axis.title = element_text(size=textsize,color="black", family = family),
                  axis.text = element_text(size=textsize,color="black", family = family),
                  legend.position = "none")
          graficos[[i]]=grafico}
        print(graficos)
    }
      print(tukey)}
      if(quali==FALSE){
        for(i in 1:length(levels(local))){
        data=dados[dados$local==levels(dados$local)[i],]
        dose1=data$tratnum#as.numeric(as.character(data$tratamento))
        resp=data$response
        grafico=polynomial(dose1,
                           resp,grau = grau,
                           textsize=textsize,
                           family=family,
                           ylab=ylab,
                           xlab=xlab,
                           theme=theme,
                           posi="top",
                           se=errorbar)[[1]]
        graficos[[i]]=grafico}
      }
  }
  if(a$`Pr(>F)`[1] > alpha.f && qmresmedio < 7){
    if(quali==TRUE){
      if(mcomp=="tukey"){
        tukeyjuntos=(TUKEY(resp,tratamento,a$Df[1], a$`Mean Sq`[1], alpha = alpha.t))
        if(transf!="1"){tukeyjuntos$groups$respo=tapply(response,
                                                        tratamento,
                                                        mean, na.rm=TRUE)[rownames(tukeyjuntos$groups)]}
        tukeyjuntos=tukeyjuntos$groups}
      if(mcomp=="duncan"){
        tukeyjuntos=duncan(resp,tratamento,a$Df[1], a$`Mean Sq`[1], alpha = alpha.t)
        if(transf!="1"){tukeyjuntos$groups$respo=tapply(response, tratamento,
                                                        mean, na.rm=TRUE)[rownames(tukeyjuntos$groups)]}
        tukeyjuntos=tukeyjuntos$groups}
      if(mcomp=="lsd"){
        tukeyjuntos=LSD(resp,tratamento,a$Df[1], a$`Mean Sq`[1], alpha = alpha.t)
        if(transf!="1"){tukeyjuntos$groups$respo=tapply(response, tratamento,
                                                        mean, na.rm=TRUE)[rownames(tukeyjuntos$groups)]}
        tukeyjuntos=tukeyjuntos$groups}
      if(mcomp=="sk"){
        nrep=table(tratamento)[1]
        medias=sort(tapply(resp,tratamento,mean),decreasing = TRUE)
        letra=scottknott(means = medias,
                         df1 = a$Df[1],
                         nrep = nrep,
                         QME = a$`Mean Sq`[1],
                         alpha = alpha.t)
        tukeyjuntos=data.frame(resp=medias,groups=letra)
        #
        # tukeyjuntos=sk(resp,tratamento,a$Df[1], a$`Sum Sq`[1],
        #                       alpha = alpha.t)
        # colnames(tukeyjuntos)=c("resp","groups")
        if(transf!="1"){tukeyjuntos$respo=tapply(response, tratamento,
                                                 mean, na.rm=TRUE)[rownames(tukeyjuntos)]}}

      dadosm=data.frame(tukeyjuntos,
                        media=tapply(response, tratamento, mean, na.rm=TRUE)[rownames(tukeyjuntos)],
                        desvio=tapply(response, tratamento, sd, na.rm=TRUE)[rownames(tukeyjuntos)])
      dadosm$trats=factor(rownames(dadosm),unique(trat))
      dadosm$limite=dadosm$media+dadosm$desvio
      dadosm$letra=paste(format(dadosm$media,digits = dec),dadosm$groups)
      dadosm=dadosm[unique(as.character(trat)),]
      media=dadosm$media
      desvio=dadosm$desvio
      limite=dadosm$limite
      trats=dadosm$trats
      letra=dadosm$letra
      grafico1=ggplot(dadosm,aes(x=trats,y=media))
      if(fill=="trat"){grafico1=grafico+
        geom_col(aes(fill=trats),color=1)}
      else{grafico1=grafico1+
        geom_col(aes(fill=trats),fill=fill,color=1)}
      if(errorbar==TRUE){grafico1=grafico1+
        geom_text(aes(y=media+sup+if(sup<0){-desvio}else{desvio},
                      label=letra),family=family)}
      if(errorbar==FALSE){grafico1=grafico1+
        geom_text(aes(y=media+sup,label=letra),family=family)}
      if(errorbar==TRUE){grafico1=grafico1+
        geom_errorbar(data=dadosm,aes(ymin=media-desvio,
                                      ymax=media+desvio,color=1),
                      color="black", width=0.3)}
      grafico1=grafico1+theme+
        ylab(ylab)+
        xlab(xlab)+
        theme(text = element_text(size=textsize,color="black",family=family),
              axis.text = element_text(size=textsize,color="black",family=family),
              axis.title = element_text(size=textsize,color="black",family=family),
              legend.position = "none")
      if(angulo !=0){grafico1=grafico1+theme(axis.text.x=element_text(hjust = 1.01,angle = angulo))}
      cat("Multiple comparison test (",mcomp,")")
      cat(green(bold("\n-----------------------------------------------------------------\n")))
      print(tukeyjuntos)
      print(grafico1)
      graficos=list(grafico1)}
    if(quali==FALSE){grafico1=polynomial(tratnum,#as.numeric(as.character(tratamento)),
                                     response,grau = grau,
                                     textsize=textsize,
                                     family=family,
                                     ylab=ylab,
                                     xlab=xlab,
                                     theme=theme,
                                     posi="top",
                                     se=errorbar)
    graficos=list(grafico1)
  }
  }
  cat(if(transf=="1"){}else{blue("\nNOTE: resp = transformed means; respO = averages without transforming\n")})

  if(transf==1 && norm1$p.value<0.05 | transf==1 && indep$p.value<0.05 | transf==1 && homog1$p.value<0.05){cat(red("\n \nWarning!!! Your analysis is not valid, suggests using a non-parametric test and try to transform the data"))}
  if(transf != 1 && norm1$p.value<0.05 | transf!=1 && indep$p.value<0.05 | transf!=1 && homog1$p.value<0.05){cat(red("\n \nWarning!!! Your analysis is not valid, suggests using a non-parametric test"))}
  graph=as.list(graficos)
}
