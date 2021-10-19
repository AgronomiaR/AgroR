#' Analysis: DBC experiments in split-split-plot
#' @description Analysis of an experiment conducted in a randomized block design in a split-split-plot scheme using analysis of variance of fixed effects.
#' @author Gabriel Danilo Shimizu, \email{shimizu@uel.br}
#' @author Leandro Simoes Azeredo Goncalves
#' @author Rodrigo Yudi Palhaci Marubayashi
#' @param f1 Numeric or complex vector with plot levels
#' @param f2 Numeric or complex vector with splitplot levels
#' @param f3 Numeric or complex vector with splitsplitplot levels
#' @param block Numeric or complex vector with blocks
#' @param mcomp Multiple comparison test (Tukey (\emph{default}), LSD and Duncan)
#' @param response Numeric vector with responses
#' @param alpha.f Level of significance of the F test (\emph{default} is 0.05)
#' @param alpha.t Significance level of the multiple comparison test (\emph{default} is 0.05)
#' @param dec Number of cells
#' @note The PSUBSUBDBC function does not present residual analysis, interaction breakdown, graphs and implementations of various multiple comparison or regression tests. The function only returns the analysis of variance and multiple comparison test of Tukey, LSD or Duncan.
#' @return Analysis of variance of fixed effects and multiple comparison test of Tukey, LSD or Duncan.
#' @keywords DBC
#' @export
#' @examples
#' library(AgroR)
#' data(enxofre)
#' with(enxofre, PSUBSUBDBC(f1, f2, f3, bloco, resp))

PSUBSUBDBC=function(f1,
                    f2,
                    f3,
                    block,
                    response,
                    alpha.f=0.05,
                    alpha.t=0.05,
                    dec=3,
                    mcomp="tukey"){
  fac.names=c("F1","F2","F3")
  fator1=as.factor(f1)
  fator2=as.factor(f2)
  fator3=as.factor(f3)
  bloco=as.factor(block)
  resp=response
  fatores<-data.frame(fator1,fator2,fator3)
  Fator1<-factor(fator1,levels=unique(fator1))
  Fator2<-factor(fator2,levels=unique(fator2))
  Fator3<-factor(fator3,levels=unique(fator3))
  nv1<-length(summary(Fator1))
  nv2<-length(summary(Fator2))
  nv3<-length(summary(Fator3))
  nbl<-length(summary(bloco))
  J<-(length(response))/(nv1*nv2*nv3)
  lf1<-levels(Fator1); lf2<-levels(Fator2); lf3<-levels(Fator3)

  # pressupostos
  m1=aov(response~Fator1*Fator2*Fator3+bloco/Fator1+bloco/Fator1/Fator2)
  summary(m1)
  norm=shapiro.test(m1$residuals)
  homog=bartlett.test(m1$residuals~paste(Fator1,Fator2,Fator3))

  cat(green(bold("\n-----------------------------------------------------------------\n")))
  cat(green(bold("Normality of errors")))
  cat(green(bold("\n-----------------------------------------------------------------\n")))
  normal=data.frame(Method=paste(norm$method,"(",names(norm$statistic),")",sep=""),
                    Statistic=norm$statistic,
                    "p-value"=norm$p.value)
  rownames(normal)=""
  print(normal)
  cat("\n")

  message(if(norm$p.value>0.05){
    black("As the calculated p-value is greater than the 5% significance level, hypothesis H0 is not rejected. Therefore, errors can be considered normal")}
    else {"As the calculated p-value is less than the 5% significance level, H0 is rejected. Therefore, errors do not follow a normal distribution"})
  cat(green(bold("\n\n-----------------------------------------------------------------\n")))
  cat(green(bold("Homogeneity of Variances")))
  cat(green(bold("\n-----------------------------------------------------------------\n")))
  statistic1=homog$statistic
  phomog1=homog$p.value
  method1=paste("Bartlett test","(",names(statistic1),")",sep="")
  homoge1=data.frame(Method=method1,
                     Statistic=statistic1,
                     "p-value"=phomog1)
  rownames(homoge1)=""
  print(homoge1)
  cat("\n")
  message(if(homog$p.value[1]>0.05){
    black("As the calculated p-value is greater than the 5% significance level, hypothesis H0 is not rejected. Therefore, the variances can be considered homogeneous")}
    else {"As the calculated p-value is less than the 5% significance level, H0 is rejected. Therefore, the variances are not homogeneous"})

  mod=aov(response~Fator1*Fator2*Fator3+
            Error(bloco/Fator1/paste(Fator1,Fator2)))
  a=summary(mod)
  anava=rbind(data.frame(a$`Error: bloco:Fator1`[[1]]),
              data.frame(a$`Error: bloco:Fator1:paste(Fator1, Fator2)`[[1]]),
              data.frame(a$`Error: Within`[[1]]))
  anavap=anava
  anava$F.value=ifelse(is.na(anava$F.value)==TRUE,"",round(anava$F.value,5))
  anava$Pr..F.=ifelse(is.na(anava$Pr..F.)==TRUE,"",round(anava$Pr..F.,5))
  rownames(anava)=c("F1","Error A","F2","F1 x F2", "Error B", "F3", "F1 x F3", "F2 x F3", "F1 x F2 x F3","Residuals")
  colnames(anava)=c("Df","Sum Sq","Mean Sq","F value","Pr(>F)")

  cat(green(bold("\n-----------------------------------------------------------------\n")))
  cat(green(bold("Additional Information")))
  cat(green(bold("\n-----------------------------------------------------------------\n")))
  cat(paste("\nCV plot (%) = ",round(sqrt(anava$`Mean Sq`[2])/mean(resp,na.rm=TRUE)*100,2)))
  cat(paste("\nCV split plot (%) = ",round(sqrt(anava$`Mean Sq`[5])/mean(resp,na.rm=TRUE)*100,2)))
  cat(paste("\nCV split split plot (%) = ",round(sqrt(anava$`Mean Sq`[10])/mean(resp,na.rm=TRUE)*100,2)))
  cat(paste("\nMean = ",round(mean(response,na.rm=TRUE),4)))
  cat(paste("\nMedian = ",round(median(response,na.rm=TRUE),4)))
  cat("\n")

  cat(green(bold("\n-----------------------------------------------------------------\n")))
  cat(green(bold("Analysis of Variance")))
  cat(green(bold("\n-----------------------------------------------------------------\n")))
  print(anava)

  fatores<-data.frame('fator 1'=fator1,
                      'fator 2' = fator2,
                      'fator 3' = fator3)
  qmres=c(as.numeric(anavap[2,3]),
          as.numeric(anavap[5,3]),
          as.numeric(anavap[10,3]))
  GL=c(as.numeric(anavap[2,1]),
       as.numeric(anavap[5,1]),
       as.numeric(anavap[10,1]))
  pvalor=c(as.numeric(anavap[1,5]),
           as.numeric(anavap[3,5]),
           as.numeric(anavap[6,5]))

  ################################################################################################
  # Efeitos simples
  ################################################################################################
  if(as.numeric(anavap[9,5])>alpha.f &&
     as.numeric(anavap[8,5])>alpha.f &&
     as.numeric(anavap[7,5])>alpha.f &&
     as.numeric(anavap[4,5])>alpha.f) {
    graficos=list(1,2,3)
    cat(green(bold("\n------------------------------------------\n")))
    cat(green(bold('Non-significant interaction: analyzing the simple effects')))
    cat(green(bold("\n------------------------------------------\n")))

    for(i in 1:3){
      if(pvalor[i]<=alpha.f) {
        cat(green(bold("\n------------------------------------------\n")))
        cat(fac.names[i])
        cat(green(bold("\n------------------------------------------\n")))
        if(mcomp=="tukey"){
        letra=TUKEY(response,
                       fatores[,i],
                       GL[i],
                       qmres[i],alpha.t)
        letra1 <- letra$groups; colnames(letra1)=c("resp","groups")
        print(letra1)}
        if(mcomp=="lsd"){
          letra=LSD(response,
                         fatores[,i],
                         GL[i],
                         qmres[i],alpha.t)
          letra1 <- letra$groups; colnames(letra1)=c("resp","groups")
          print(letra1)}
        if(mcomp=="duncan"){
          letra=duncan(response,
                         fatores[,i],
                         GL[i],
                         qmres[i],alpha.t)
          letra1 <- letra$groups; colnames(letra1)=c("resp","groups")
          print(letra1)}
        if(mcomp=="sk"){
          nrep=table(fatores[, i])[1]
          medias=sort(tapply(resp,fatores[, i],mean),decreasing = TRUE)
          letra=scottknott(means = medias,
                           df1 = GL[i],
                           nrep = nrep,
                           QME = qmres[i],
                           alpha = alpha.t)
          letra1=data.frame(resp=medias,groups=letra)
          #
          # letra1=sk(response,
          #              fatores[,i],
          #              GL[i],
          #              qmres[i]/GL[i],alpha.t)
          # colnames(letra1)=c("resp","groups")
          print(letra1)}
        }
      if(pvalor[i]>alpha.f) {
        cat(fac.names[i])
        cat(green(bold("\n------------------------------------------\n")))
        mean.table<-mean.stat(response,fatores[,i],mean)
        colnames(mean.table)<-c('Levels','Mean')
        print(mean.table)
        grafico=NA
        cat(green(bold("\n------------------------------------------")))}
      cat('\n')
    }
  }

  #####################################################################
  #Interacao Fator1*Fator2      +     Fator3
  #####################################################################
  # corrigir para variancia complexa
  qmresf1f2=(qmres[1]+(nv2-1)*qmres[2])/nv2

  # gl de Satterthwaite
  (nf1f2=((qmres[1]+(nv2-1)*qmres[2])^2)/
      ((qmres[1]^2)/GL[1]+(((nv2-1)*qmres[2])^2)/GL[2]))
  nf1f2=round(nf1f2)

  if(as.numeric(anavap[9,5])>alpha.f &&
     as.numeric(anavap[4,5])<=alpha.f){
    cat(green(bold("\n------------------------------------------\n")))
    cat(green(bold("Interaction",paste(fac.names[1],'*',fac.names[2],sep='')," significant: unfolding the interaction")))
    cat(green(bold("\n------------------------------------------\n")))

    cat(green(bold("\n------------------------------------------\n")))
    cat("Analyzing ", fac.names[2], " inside of the level of ",fac.names[1])
    cat(green(bold("\n------------------------------------------\n")))

    #################################################################
    #### desdobramento f2 dentro de f1
    #################################################################
    mod=aov(response~Fator1/Fator2+Fator2:Fator3+Fator1:Fator3+Fator1:Fator2:Fator3+
              Error(bloco/Fator1/Fator2))
    summary(mod)
    l2<-vector('list',nv1)
    names(l2)<-names(summary(Fator1))
    v<-numeric(0)
    for(j in 1:nv1) {
      for(i in 0:(nv2-2)) v<-cbind(v,i*nv1+j)
      l2[[j]]<-v
      v<-numeric(0)
    }
    des1.tab<-summary(mod,split=list('Fator1:Fator2'=l2))
    desdf2f1=data.frame(des1.tab$`Error: bloco:Fator1:Fator2`[[1]])
    colnames(desdf2f1)=c("Df","Sum sq","Mean Sq", "F value","Pr(>F)")
    print(as.matrix(desdf2f1),na.print="")

    ###############################################
    #### desdobramento f1 dentro de f2
    ###############################################

    cat(green(bold("\n------------------------------------------\n")))
    cat("Analyzing ", fac.names[1], " inside of the level of ",fac.names[2])
    cat(green(bold("\n------------------------------------------\n")))

    mod=aov(response~Fator2/Fator1+Fator2:Fator3+
              Fator1:Fator3+Fator1:Fator2:Fator3+
              Error(bloco/Fator2))
    summary(mod)
    l1<-vector('list',nv2)
    names(l1)<-names(summary(Fator2))
    v<-numeric(0)
    for(j in 1:nv2) {
      for(i in 0:(nv1-2)) v<-cbind(v,i*nv2+j)
      l1[[j]]<-v
      v<-numeric(0)
    }
    desd1.tab<-summary(mod,split=list('Fator2:Fator1'=l1))
    desd1.tab
    desd=data.frame(desd1.tab$`Error: Within`[[1]])
    desd
    nlinhas=nrow(desd)
    desd=desd[-c(1,nlinhas-3,nlinhas-2,nlinhas-1,nlinhas),]
    qmresf1f2=(qmres[1]+(nv2-1)*qmres[2])/nv2
    nf1f2=((qmres[1]+(nv2-1)*qmres[2])^2)/
      ((qmres[1]^2)/GL[1]+(((nv2-1)*qmres[2])^2)/GL[2])
    nf1f2=round(nf1f2)
    desd$F.value=desd$Mean.Sq/qmresf1f2
    nline=nrow(desd)
    for(i in 1:nline){
      desd$Pr..F.[i]=1-pf(desd$F.value[i],desd$Df[i],nf1f2)}
    f1f2=data.frame(desd1.tab$`Error: Within`[[1]])[1,]
    desdf1f2=rbind(f1f2,desd,c(nf1f2,qmresf1f2/nf1f2,qmresf1f2,NA,NA))
    nline1=nrow(desdf1f2)
    rownames(desdf1f2)[nline1]="Residuals combined"
    colnames(desdf1f2)=c("Df","Sum sq","Mean Sq", "F value","Pr(>F)")
    print(as.matrix(desdf1f2),na.print = "")

    cat(green(bold("\n------------------------------------------\n")))
    cat(green(bold("Final table")))
    cat(green(bold("\n------------------------------------------\n")))

    if(mcomp=="tukey"){
    tukeygrafico=c()
    ordem=c()
    for (i in 1:nv2) {
      trati=fatores[, 1][Fator2 == lf2[i]]
      trati=factor(trati,levels = unique(trati))
      respi=response[Fator2 == lf2[i]]
      tukey=TUKEY(respi,trati,nf1f2,qmresf1f2,alpha.t)
      tukeygrafico[[i]]=tukey$groups[levels(trati),2]
      ordem[[i]]=rownames(tukey$groups[levels(trati),])
    }
    letra=unlist(tukeygrafico)
    datag=data.frame(letra,ordem=unlist(ordem))
    datag$ordem=factor(datag$ordem,levels = unique(datag$ordem))
    datag=datag[order(datag$ordem),]
    letra=datag$letra

    tukeygrafico1=c()
    for (i in 1:nv1) {
      trati=fatores[, 2][Fator1 == lf1[i]]
      trati=factor(trati,levels = unique(trati))
      respi=response[Fator1 == lf1[i]]
      tukey=TUKEY(respi,trati,GL[2],qmres[2],alpha.t)
      tukeygrafico1[[i]]=tukey$groups[levels(trati),2]
    }
    letra1=unlist(tukeygrafico1)
    letra1=toupper(letra1)}
    if(mcomp=="lsd"){
      lsdgrafico=c()
      ordem=c()
      for (i in 1:nv2) {
        trati=fatores[, 1][Fator2 == lf2[i]]
        trati=factor(trati,levels = unique(trati))
        respi=response[Fator2 == lf2[i]]
        lsd=LSD(respi,trati,nf1f2,qmresf1f2,alpha.t)
        lsdgrafico[[i]]=lsd$groups[levels(trati),2]
        ordem[[i]]=rownames(lsd$groups[levels(trati),])
      }
      letra=unlist(lsdgrafico)
      datag=data.frame(letra,ordem=unlist(ordem))
      datag$ordem=factor(datag$ordem,levels = unique(datag$ordem))
      datag=datag[order(datag$ordem),]
      letra=datag$letra

      lsdgrafico1=c()
      for (i in 1:nv1) {
        trati=fatores[, 2][Fator1 == lf1[i]]
        trati=factor(trati,levels = unique(trati))
        respi=response[Fator1 == lf1[i]]
        lsd=LSD(respi,trati,GL[2],qmres[2],alpha.t)
        lsdgrafico1[[i]]=lsd$groups[levels(trati),2]
      }
      letra1=unlist(lsdgrafico1)
      letra1=toupper(letra1)}
    if(mcomp=="duncan"){
      duncangrafico=c()
      ordem=c()
      for (i in 1:nv2) {
        trati=fatores[, 1][Fator2 == lf2[i]]
        trati=factor(trati,levels = unique(trati))
        respi=response[Fator2 == lf2[i]]
        duncan=duncan(respi,trati,nf1f2,qmresf1f2,alpha.t)
        duncangrafico[[i]]=duncan$groups[levels(trati),2]
        ordem[[i]]=rownames(duncan$groups[levels(trati),])
      }
      letra=unlist(duncangrafico)
      datag=data.frame(letra,ordem=unlist(ordem))
      datag$ordem=factor(datag$ordem,levels = unique(datag$ordem))
      datag=datag[order(datag$ordem),]
      letra=datag$letra

      duncangrafico1=c()
      for (i in 1:nv1) {
        trati=fatores[, 2][Fator1 == lf1[i]]
        trati=factor(trati,levels = unique(trati))
        respi=response[Fator1 == lf1[i]]
        duncan=duncan(respi,trati,GL[2],qmres[2],alpha.t)
        duncangrafico1[[i]]=duncan$groups[levels(trati),2]
      }
      letra1=unlist(duncangrafico1)
      letra1=toupper(letra1)}

    if(mcomp=="sk"){
      skgrafico=c()
      ordem=c()
      for (i in 1:nv2) {
        trati=fatores[, 1][Fator2 == lf2[i]]
        trati=factor(trati,levels = unique(trati))
        respi=response[Fator2 == lf2[i]]
        nrep=table(trati)[1]
        medias=sort(tapply(respi,trati,mean),decreasing = TRUE)
        sk=scottknott(means = medias,
                      df1 = nf1f2,
                      nrep = nrep,
                      QME = qmresf1f2,
                      alpha = alpha.t)
        sk=data.frame(respi=medias,groups=sk)
        # sk=sk(respi,trati,nf1f2,qmresf1f2/nf1f2,alpha.t)
        skgrafico[[i]]=sk[levels(trati),2]
        ordem[[i]]=rownames(sk[levels(trati),])
      }
      letra=unlist(skgrafico)
      datag=data.frame(letra,ordem=unlist(ordem))
      datag$ordem=factor(datag$ordem,levels = unique(datag$ordem))
      datag=datag[order(datag$ordem),]
      letra=datag$letra

      skgrafico1=c()
      for (i in 1:nv1) {
        trati=fatores[, 2][Fator1 == lf1[i]]
        trati=factor(trati,levels = unique(trati))
        respi=response[Fator1 == lf1[i]]
        nrep=table(trati)[1]
        medias=sort(tapply(respi,trati,mean),decreasing = TRUE)
        sk=scottknott(means = medias,
                      df1 = GL[2],
                      nrep = nrep,
                      QME = qmres[2],
                      alpha = alpha.t)
        sk=data.frame(respi=medias,groups=sk)
        # sk=sk(respi,trati,GL[2],qmres[2]/GL[2],alpha.t)
        skgrafico1[[i]]=sk[levels(trati),2]
      }
      letra1=unlist(skgrafico1)
      letra1=toupper(letra1)}
    f1=rep(levels(Fator1),e=length(levels(Fator2)))
    f2=rep(unique(as.character(Fator2)),length(levels(Fator1)))
    media=tapply(response,paste(Fator1,Fator2), mean, na.rm=TRUE)[unique(paste(f1,f2))]
    desvio=tapply(response,paste(Fator1,Fator2), sd, na.rm=TRUE)[unique(paste(f1,f2))]
    f1=factor(f1,levels = unique(f1))
    f2=factor(f2,levels = unique(f2))

    graph=data.frame(f1=f1,
                     f2=f2,
                     media,
                     desvio,
                     letra,letra1,
                     numero=format(media,digits = dec))
    numero=graph$numero
    letras=paste(graph$letra, graph$letra1, sep="")
    matriz=data.frame(t(matrix(paste(format(graph$media,digits = dec),letras),
                               ncol = length(levels(Fator1)))))
    rownames(matriz)=levels(Fator1)
    colnames(matriz)=levels(Fator2)
    print(matriz)
    message(black("\nAverages followed by the same lowercase letter in the column and \nuppercase in the row do not differ by the", mcomp, "(p<",alpha.t,")\n"))


    #Checar o Fator3
    if(as.numeric(anavap[7,5])>alpha.f && as.numeric(anavap[8,5])>alpha.f) {
      if(pvalor[3]<=alpha.f) {
        cat(green(bold("\n------------------------------------------\n")))
        cat(green(italic('Analyzing the simple effects of the factor ',fac.names[3])))
        cat(green(bold("\n------------------------------------------\n")))
        cat(fac.names[i])
        if(mcomp=="tukey"){letra=TUKEY(response,fatores[,i],GL[3],qmres[3],alpha.t)}
        if(mcomp=="lsd"){letra=LSD(response,fatores[,i],GL[3],qmres[3],alpha.t)}
        if(mcomp=="duncan"){letra=duncan(response,fatores[,i],GL[3],qmres[3],alpha.t)}
        letra1 <- letra$groups; colnames(letra1)=c("resp","groups")
        print(letra1)
        cat(green(bold("\n-----------------------------------------------------------------")))
      }
    }
  }

  #####################################################################################################
  #Interacao Fator1*Fator3       + fator2
  #####################################################################################################
  # corrigir para variancia complexa
  qmresf1f3=(qmres[1]+(nv3-1)*qmres[3])/nv3

  # gl de Satterthwaite
  (nf1f3=((qmres[1]+(nv3-1)*qmres[3])^2)/
      ((qmres[1]^2)/GL[1]+(((nv3-1)*qmres[3])^2)/GL[3]))
  nf1f3=round(nf1f3)
  if(as.numeric(anavap[9,5])>alpha.f &&
     as.numeric(anavap[7,5])<=alpha.f){
    cat(green(bold("\n------------------------------------------\n")))
    cat(green(bold("\nInteraction",paste(fac.names[1],'*',fac.names[3],sep='')," significant: unfolding the interaction\n")))
    cat(green(bold("\n------------------------------------------\n")))

    #################################################################
    #### desdobramento f3 dentro de f1
    #################################################################
    cat(green(bold("\n------------------------------------------\n")))
    cat("Analyzing ", fac.names[3], " inside of the level of ",fac.names[1])
    cat(green(bold("\n------------------------------------------\n")))

    mod=aov(response~Fator1/Fator3+Fator1*Fator2+
              Fator1:Fator2:Fator3+
              Error(bloco/Fator1/Fator2))
    summary(mod)
    l3<-vector('list',nv1)
    names(l3)<-names(summary(Fator1))
    v<-numeric(0)
    for(j in 1:nv1) {
      for(i in 0:(nv3-2)) v<-cbind(v,i*nv1+j)
      l3[[j]]<-v
      v<-numeric(0)
    }
    desd1.tab<-summary(mod,split=list('Fator1:Fator3'=l3))
    desdf3f1=data.frame(desd1.tab$`Error: Within`[[1]])
    nlines=nrow(desdf3f1)
    desdf3f1=desdf3f1[-(nlines-1),]
    colnames(desdf3f1)=c("Df","Sum sq","Mean Sq", "F value","Pr(>F)")
    print(as.matrix(desdf3f1),na.print="")

    ###############################################
    #### desdobramento f1 dentro de f3
    ###############################################
    cat(green(bold("\n------------------------------------------\n")))
    cat("Analyzing ", fac.names[1], " inside of the level of ",fac.names[3])
    cat(green(bold("\n------------------------------------------\n")))

    mod=aov(response~Fator3/Fator1+Fator2*Fator3+Fator1:Fator2:Fator3+
              Error(bloco/Fator2))
    summary(mod)
    l1<-vector('list',nv3)
    names(l1)<-names(summary(Fator3))
    v<-numeric(0)
    for(j in 1:nv3) {
      for(i in 0:(nv1-2)) v<-cbind(v,i*nv3+j)
      l1[[j]]<-v
      v<-numeric(0)
    }
    desd1.tab<-summary(mod,split=list('Fator3:Fator1'=l1))
    desd1.tab
    desd=data.frame(desd1.tab$`Error: Within`[[1]])
    nlinhas=nrow(desd)
    desd=desd[-c(1,2,nlinhas-2,nlinhas-1,nlinhas),]
    qmresf1f3=(qmres[1]+(nv3-1)*qmres[3])/nv3
    nf1f3=((qmres[1]+(nv3-1)*qmres[3])^2)/
      ((qmres[1]^2)/GL[1]+(((nv3-1)*qmres[3])^2)/GL[3])
    nf1f3=round(nf1f3)
    desd$F.value=desd$Mean.Sq/qmresf1f3
    nline=nrow(desd)
    for(i in 1:nline){
      desd$Pr..F.[i]=1-pf(desd$F.value[i],desd$Df[i],nf1f3)}
    f1f3=data.frame(desd1.tab$`Error: Within`[[1]])[2,]
    desdf1f3=rbind(f1f3,desd,c(nf1f3,qmresf1f3/nf1f3,qmresf1f3,NA,NA))
    nline1=nrow(desdf1f3)
    rownames(desdf1f3)[nline1]="Residuals combined"
    colnames(desdf1f3)=c("Df","Sum sq","Mean Sq", "F value","Pr(>F)")
    print(as.matrix(desdf1f3),na.print = "")


    cat(green(bold("\n------------------------------------------\n")))
    cat(green(bold("Final table")))
    cat(green(bold("\n------------------------------------------\n")))


    if(mcomp=="tukey"){
    tukeygrafico=c()
    ordem=c()
    for (i in 1:nv3) {
      trati=fatores[, 1][Fator3 == lf3[i]]
      trati=factor(trati,levels = unique(trati))
      respi=response[Fator3 == lf3[i]]
      tukey=TUKEY(respi,trati,nf1f3,qmresf1f3,alpha.t)
      tukeygrafico[[i]]=tukey$groups[levels(trati),2]
      ordem[[i]]=rownames(tukey$groups[levels(trati),])}
    letra=unlist(tukeygrafico)
    datag=data.frame(letra,ordem=unlist(ordem))
    datag$ordem=factor(datag$ordem,levels = unique(datag$ordem))
    datag=datag[order(datag$ordem),]
    letra=datag$letra

    tukeygrafico1=c()
    for (i in 1:nv1) {
      trati=fatores[, 3][Fator1 == lf1[i]]
      trati=factor(trati,levels = unique(trati))
      respi=response[Fator1 == lf1[i]]
      tukey=TUKEY(respi,trati,GL[3],qmres[3],alpha.t)
      tukeygrafico1[[i]]=tukey$groups[levels(trati),2]}
    letra1=unlist(tukeygrafico1)
    letra1=toupper(letra1)}
    if(mcomp=="lsd"){
      lsdgrafico=c()
      ordem=c()
      for (i in 1:nv3) {
        trati=fatores[, 1][Fator3 == lf3[i]]
        trati=factor(trati,levels = unique(trati))
        respi=response[Fator3 == lf3[i]]
        lsd=LSD(respi,trati,nf1f3,qmresf1f3,alpha.t)
        lsdgrafico[[i]]=lsd$groups[levels(trati),2]
        ordem[[i]]=rownames(lsd$groups[levels(trati),])}
      letra=unlist(lsdgrafico)
      datag=data.frame(letra,ordem=unlist(ordem))
      datag$ordem=factor(datag$ordem,levels = unique(datag$ordem))
      datag=datag[order(datag$ordem),]
      letra=datag$letra

      lsdgrafico1=c()
      for (i in 1:nv1) {
        trati=fatores[, 3][Fator1 == lf1[i]]
        trati=factor(trati,levels = unique(trati))
        respi=response[Fator1 == lf1[i]]
        lsd=LSD(respi,trati,GL[3],qmres[3],alpha.t)
        lsdgrafico1[[i]]=lsd$groups[levels(trati),2]}
      letra1=unlist(lsdgrafico1)
      letra1=toupper(letra1)}
    if(mcomp=="duncan"){
      duncangrafico=c()
      ordem=c()
      for (i in 1:nv3) {
        trati=fatores[, 1][Fator3 == lf3[i]]
        trati=factor(trati,levels = unique(trati))
        respi=response[Fator3 == lf3[i]]
        duncan=duncan(respi,trati,nf1f3,qmresf1f3,alpha.t)
        duncangrafico[[i]]=duncan$groups[levels(trati),2]
        ordem[[i]]=rownames(duncan$groups[levels(trati),])}
      letra=unlist(duncangrafico)
      datag=data.frame(letra,ordem=unlist(ordem))
      datag$ordem=factor(datag$ordem,levels = unique(datag$ordem))
      datag=datag[order(datag$ordem),]
      letra=datag$letra

      duncangrafico1=c()
      for (i in 1:nv1) {
        trati=fatores[, 3][Fator1 == lf1[i]]
        trati=factor(trati,levels = unique(trati))
        respi=response[Fator1 == lf1[i]]
        duncan=duncan(respi,trati,GL[3],qmres[3],alpha.t)
        duncangrafico1[[i]]=duncan$groups[levels(trati),2]}
      letra1=unlist(duncangrafico1)
      letra1=toupper(letra1)}
    if(mcomp=="sk"){
      skgrafico=c()
      ordem=c()
      for (i in 1:nv3) {
        trati=fatores[, 1][Fator3 == lf3[i]]
        trati=factor(trati,levels = unique(trati))
        respi=response[Fator3 == lf3[i]]
        nrep=table(trati)[1]
        medias=sort(tapply(respi,trati,mean),decreasing = TRUE)
        sk=scottknott(means = medias,
                      df1 = nf1f3,
                      nrep = nrep,
                      QME = qmresf1f3,
                      alpha = alpha.t)
        sk=data.frame(respi=medias,groups=sk)
        # sk=sk(respi,trati,nf1f3,qmresf1f3/nf1f3,alpha.t)
        skgrafico[[i]]=sk[levels(trati),2]
        ordem[[i]]=rownames(sk[levels(trati),])}
      letra=unlist(skgrafico)
      datag=data.frame(letra,ordem=unlist(ordem))
      datag$ordem=factor(datag$ordem,levels = unique(datag$ordem))
      datag=datag[order(datag$ordem),]
      letra=datag$letra

      skgrafico1=c()
      for (i in 1:nv1) {
        trati=fatores[, 3][Fator1 == lf1[i]]
        trati=factor(trati,levels = unique(trati))
        respi=response[Fator1 == lf1[i]]
        nrep=table(trati)[1]
        medias=sort(tapply(respi,trati,mean),decreasing = TRUE)
        sk=scottknott(means = medias,
                      df1 = GL[3],
                      nrep = nrep,
                      QME = qmres[3],
                      alpha = alpha.t)
        sk=data.frame(respi=medias,groups=sk)
        # sk=sk(respi,trati,GL[3],qmres[3]/GL[3],alpha.t)
        skgrafico1[[i]]=sk[levels(trati),2]}
      letra1=unlist(skgrafico1)
      letra1=toupper(letra1)}
    f1=rep(levels(Fator1),e=length(levels(Fator3)))
    f3=rep(unique(as.character(Fator3)),length(levels(Fator1)))
    media=tapply(response,paste(Fator1,Fator3), mean, na.rm=TRUE)[unique(paste(f1,f3))]
    desvio=tapply(response,paste(Fator1,Fator3), sd, na.rm=TRUE)[unique(paste(f1,f3))]
    f1=factor(f1,levels = unique(f1))
    f3=factor(f3,levels = unique(f3))
    graph=data.frame(f1=f1, f3=f3,
                     media, desvio,
                     letra,letra1,
                     numero=format(media,digits = dec))
    numero=graph$numero
    letras=paste(graph$letra,graph$letra1,sep="")
    matriz=data.frame(t(matrix(paste(format(graph$media,digits = dec),letras),ncol = length(levels(Fator1)))))
    rownames(matriz)=levels(Fator1)
    colnames(matriz)=levels(Fator3)
    print(matriz)
    message(black("\nAverages followed by the same lowercase letter in the column and \nuppercase in the row do not differ by the", mcomp, "(p<",alpha.t,")\n"))

    #Checar o Fator2
    if(as.numeric(anavap[4,5])>alpha.f && as.numeric(anavap[6,5])>alpha.f) {
      i=2
      cat(green(bold("\n------------------------------------------\n")))
      cat(green(italic('Analyzing the simple effects of the factor ',fac.names[2])))
      cat(green(bold("\n------------------------------------------\n")))
      cat(fac.names[i])
      if(mcomp=="tukey"){letra=TUKEY(response,fatores[,i], GL[2], qmres[2],alpha.t)}
      if(mcomp=="lsd"){letra=LSD(response,fatores[,i], GL[2], qmres[2],alpha.t)}
      if(mcomp=="duncan"){letra=duncan(response,fatores[,i], GL[2], qmres[2],alpha.t)}
      letra1 <- letra$groups; colnames(letra1)=c("resp","groups")
      print(letra1)
      cat(green(bold("\n-----------------------------------------------------------------")))
    }
  }

  ######################################################################################################################
  #Interacao Fator2*Fator3     + fator1
  ######################################################################################################################
  # corrigir para variancia complexa
  qmresf2f3=(qmres[2]+(nv3-1)*qmres[3])/nv3

  # gl de Satterthwaite
  (nf2f3=((qmres[2]+(nv3-1)*qmres[3])^2)/
    ((qmres[2]^2)/GL[2]+(((nv3-1)*qmres[3])^2)/GL[3]))
  nf2f3=round(nf2f3)
  if(as.numeric(anavap[9,5])>alpha.f &&
     as.numeric(anavap[8,5])<=alpha.f){
    cat(green(bold("\n------------------------------------------\n")))
    cat(green(bold("Interaction",paste(fac.names[2],'*',fac.names[3],sep='')," significant: unfolding the interaction")))
    cat(green(bold("\n------------------------------------------\n")))

    cat(green(bold("\n------------------------------------------\n")))
    cat("Analyzing ", fac.names[2], ' inside of each level of ', fac.names[3])
    cat(green(bold("\n------------------------------------------\n")))

    #### desdobramento f3 dentro de f2
    mod=aov(response~Fator1*Fator2+Fator2/Fator3+Fator1:Fator2:Fator3+
              Error(bloco/Fator1/Fator2))
    summary(mod)
    l3<-vector('list',nv2)
    names(l3)<-names(summary(Fator2))
    v<-numeric(0)
    for(j in 1:nv2) {
      for(i in 0:(nv3-2)) v<-cbind(v,i*nv2+j)
      l3[[j]]<-v
      v<-numeric(0)
    }
    des1.tab<-summary(mod,split=list('Fator2:Fator3'=l3))
    desd=data.frame(des1.tab$`Error: Within`[[1]])
    nlinhas=nrow(desd)
    desdf3f2=desd[-c(nlinhas-1),]
    colnames(desdf3f2)=c("Df","Sum sq","Mean Sq", "F value","Pr(>F)")
    print(as.matrix(desdf3f2),na.print="")

    #### desdobramento f2 dentro de f3
    cat(green(bold("\n------------------------------------------\n")))
    cat("Analyzing ", fac.names[3], ' inside of each level of ', fac.names[2])
    cat(green(bold("\n------------------------------------------\n")))
    mod=aov(response~Fator1+Fator3/Fator2+Fator1:Fator2+Fator1:Fator2:Fator3+
              Error(bloco/Fator1))
    l2<-vector('list',nv3)
    names(l2)<-names(summary(Fator3))
    v<-numeric(0)
    for(j in 1:nv3) {
      for(i in 0:(nv2-2)) v<-cbind(v,i*nv3+j)
      l2[[j]]<-v
      v<-numeric(0)
    }
    des1.tab<-summary(mod,split=list('Fator3:Fator2'=l2))
    desd=data.frame(des1.tab$`Error: Within`[[1]])
    nlinhas=nrow(desd)
    desd=desd[-c(1,2,nlinhas-1,nlinhas-2,nlinhas),]
    qmresf2f3=(qmres[2]+(nv3-1)*qmres[3])/nv3
    nf2f3=((qmres[2]+(nv3-1)*qmres[3])^2)/
      ((qmres[2]^2)/GL[2]+(((nv3-1)*qmres[3])^2)/GL[3])
    nf2f3=round(nf2f3)
    desd$F.value=desd$Mean.Sq/qmresf2f3
    nline=nrow(desd)
    for(i in 1:nline){
      desd$Pr..F.[i]=1-pf(desd$F.value[i],desd$Df[i],nf2f3)}
    f3f2=data.frame(des1.tab$`Error: Within`[[1]])[2,]
    desdf2f3=rbind(f3f2,desd,c(nf2f3,qmresf2f3/nf2f3,qmresf2f3,NA,NA))
    nline1=nrow(desdf2f3)
    rownames(desdf2f3)[nline1]="Residuals combined"
    colnames(desdf2f3)=c("Df","Sum sq","Mean Sq", "F value","Pr(>F)")
    print(as.matrix(desdf2f3),na.print ="")

    cat(green(bold("\n------------------------------------------\n")))
    cat(green(bold("Final table")))
    cat(green(bold("\n------------------------------------------\n")))

    if(mcomp=="tukey"){
    tukeygrafico=c()
    ordem=c()
    for (i in 1:nv3) {
      trati=fatores[, 2][Fator3 == lf3[i]]
      trati=factor(trati,levels = unique(trati))
      respi=response[Fator3 == lf3[i]]
      tukey=TUKEY(respi,trati,nf2f3,qmresf2f3,alpha.t)
      tukeygrafico[[i]]=tukey$groups[levels(trati),2]
      ordem[[i]]=rownames(tukey$groups[levels(trati),])}
    letra=unlist(tukeygrafico)
    datag=data.frame(letra,ordem=unlist(ordem))
    datag$ordem=factor(datag$ordem,levels = unique(datag$ordem))
    datag=datag[order(datag$ordem),]
    letra=datag$letra
    tukeygrafico1=c()
    for (i in 1:nv2) {
      trati=fatores[, 3][Fator2 == lf2[i]]
      trati=factor(trati,levels = unique(trati))
      respi=response[Fator2 == lf2[i]]
      tukey=TUKEY(respi,trati,GL[3],qmres[3],alpha.t)
      tukeygrafico1[[i]]=tukey$groups[levels(trati),2]}
    letra1=unlist(tukeygrafico1)
    letra1=toupper(letra1)}
    if(mcomp=="lsd"){
      lsdgrafico=c()
      ordem=c()
      for (i in 1:nv3) {
        trati=fatores[, 2][Fator3 == lf3[i]]
        trati=factor(trati,levels = unique(trati))
        respi=response[Fator3 == lf3[i]]
        lsd=LSD(respi,trati,nf2f3,qmresf2f3,alpha.t)
        lsdgrafico[[i]]=lsd$groups[levels(trati),2]
        ordem[[i]]=rownames(lsd$groups[levels(trati),])}
      letra=unlist(lsdgrafico)
      datag=data.frame(letra,ordem=unlist(ordem))
      datag$ordem=factor(datag$ordem,levels = unique(datag$ordem))
      datag=datag[order(datag$ordem),]
      letra=datag$letra
      lsdgrafico1=c()
      for (i in 1:nv2) {
        trati=fatores[, 3][Fator2 == lf2[i]]
        trati=factor(trati,levels = unique(trati))
        respi=response[Fator2 == lf2[i]]
        lsd=LSD(respi,trati,GL[3],qmres[3],alpha.t)
        lsdgrafico1[[i]]=lsd$groups[levels(trati),2]}
      letra1=unlist(lsdgrafico1)
      letra1=toupper(letra1)}
    if(mcomp=="duncan"){
      duncangrafico=c()
      ordem=c()
      for (i in 1:nv3) {
        trati=fatores[, 2][Fator3 == lf3[i]]
        trati=factor(trati,levels = unique(trati))
        respi=response[Fator3 == lf3[i]]
        duncan=duncan(respi,trati,nf2f3,qmresf2f3,alpha.t)
        duncangrafico[[i]]=duncan$groups[levels(trati),2]
        ordem[[i]]=rownames(duncan$groups[levels(trati),])}
      letra=unlist(duncangrafico)
      datag=data.frame(letra,ordem=unlist(ordem))
      datag$ordem=factor(datag$ordem,levels = unique(datag$ordem))
      datag=datag[order(datag$ordem),]
      letra=datag$letra
      duncangrafico1=c()
      for (i in 1:nv2) {
        trati=fatores[, 3][Fator2 == lf2[i]]
        trati=factor(trati,levels = unique(trati))
        respi=response[Fator2 == lf2[i]]
        duncan=duncan(respi,trati,GL[3],qmres[3],alpha.t)
        duncangrafico1[[i]]=duncan$groups[levels(trati),2]}
      letra1=unlist(duncangrafico1)
      letra1=toupper(letra1)}
    if(mcomp=="sk"){
      skgrafico=c()
      ordem=c()
      for (i in 1:nv3) {
        trati=fatores[, 2][Fator3 == lf3[i]]
        trati=factor(trati,levels = unique(trati))
        respi=response[Fator3 == lf3[i]]
        nrep=table(trati)[1]
        medias=sort(tapply(respi,trati,mean),decreasing = TRUE)
        sk=scottknott(means = medias,
                      df1 = nf2f3,
                      nrep = nrep,
                      QME = qmresf2f3,
                      alpha = alpha.t)
        sk=data.frame(respi=medias,groups=sk)
        # sk=sk(respi,trati,nf2f3,qmresf2f3/nf2f3,alpha.t)
        skgrafico[[i]]=sk[levels(trati),2]
        ordem[[i]]=rownames(sk[levels(trati),])}
      letra=unlist(skgrafico)
      datag=data.frame(letra,ordem=unlist(ordem))
      datag$ordem=factor(datag$ordem,levels = unique(datag$ordem))
      datag=datag[order(datag$ordem),]
      letra=datag$letra
      skgrafico1=c()
      for (i in 1:nv2) {
        trati=fatores[, 3][Fator2 == lf2[i]]
        trati=factor(trati,levels = unique(trati))
        respi=response[Fator2 == lf2[i]]
        nrep=table(trati)[1]
        medias=sort(tapply(respi,trati,mean),decreasing = TRUE)
        sk=scottknott(means = medias,
                      df1 = GL[3],
                      nrep = nrep,
                      QME = qmres[3],
                      alpha = alpha.t)
        sk=data.frame(respi=medias,groups=sk)
        # sk=sk(respi,trati,GL[3],qmres[3]/GL[3],alpha.t)
        skgrafico1[[i]]=sk[levels(trati),2]}
      letra1=unlist(skgrafico1)
      letra1=toupper(letra1)}
    f2=rep(levels(Fator2),e=length(levels(Fator3)))
    f3=rep(unique(as.character(Fator3)),length(levels(Fator2)))
    media=tapply(response,paste(Fator2,Fator3), mean, na.rm=TRUE)[unique(paste(f2,f3))]
    desvio=tapply(response,paste(Fator2,Fator3), sd, na.rm=TRUE)[unique(paste(f2,f3))]
    f2=factor(f2,levels = unique(f2))
    f3=factor(f3,levels = unique(f3))
    graph=data.frame(f2=f2,
                     f3=f3,
                     media,
                     desvio,
                     letra,letra1,
                     numero=format(media,digits = dec))
    numero=graph$numero
    letras=paste(graph$letra,graph$letra1,sep="")
    matriz=data.frame(t(matrix(paste(format(graph$media,digits = dec),letras),ncol = length(levels(Fator2)))))
    rownames(matriz)=levels(Fator2)
    colnames(matriz)=levels(Fator3)
    print(matriz)
    message(black("\nAverages followed by the same lowercase letter in the column and \nuppercase in the row do not differ by the Tukey (p<",alpha.t,")\n"))
    #Checar o Fator1
    if(as.numeric(anavap[4,5])>alpha.f && as.numeric(anavap[7,5])>alpha.f) {
      i<-1
      if(pvalor[i]<=alpha.f) {
        cat(green(bold("\n------------------------------------------\n")))
        cat(green(italic('Analyzing the simple effects of the factor ',fac.names[1])))
        cat(green(bold("\n------------------------------------------\n")))
        cat(fac.names[i])
        if(mcomp=="tukey"){letra=TUKEY(response,fatores[,i],GL[1],qmres[1],alpha.t)}
        if(mcomp=="lsd"){letra=LSD(response,fatores[,i],GL[1],qmres[1],alpha.t)}
        if(mcomp=="duncan"){letra=duncan(response,fatores[,i],GL[1],qmres[1],alpha.t)}
        letra1 <- letra$groups; colnames(letra1)=c("resp","groups")
        print(letra1)
        cat(green(bold("\n------------------------------------------\n")))
      }
    }
  }

  #########################################################################################################################
  #Para interacao tripla significativa, desdobramento
  #########################################################################################################################
  qmresf2f3=(qmres[2]+(nv3-1)*qmres[3])/nv3
  (nf2f3=((qmres[2]+(nv3-1)*qmres[3])^2)/
      ((qmres[2]^2)/GL[2]+(((nv3-1)*qmres[3])^2)/GL[3]))
  nf2f3=round(nf2f3)

  glconj=c(nf1f2,nf1f3,nf2f3)
  qmconj=c(qmresf1f2,qmresf1f3,qmresf2f3)
  if(as.numeric(anavap[9,5])<=alpha.f){
    cat(green(bold("\n------------------------------------------\n")))
    cat(green(bold("Interaction",paste(fac.names[1],'*',fac.names[2],'*',fac.names[3],sep='')," significant: unfolding the interaction")))
    cat(green(bold("\n------------------------------------------\n")))


    cat(green(bold("\n------------------------------------------\n")))
    cat("Analyzing ", fac.names[2], ' inside of each level of ', fac.names[1], 'and',fac.names[3])
    cat(green(bold("\n------------------------------------------\n")))
    # testando f2
    # substituir qmres
    qmresf1f3=(qmres[1]+(nv3-1)*qmres[3])/nv3
    nf1f3=((qmres[1]+(nv3-1)*qmres[3])^2)/
      ((qmres[1]^2)/GL[1]+(((nv3-1)*qmres[3])^2)/GL[3])
    nf1f3=round(nf1f3)
    mod=aov(response~Fator3/Fator2/Fator1)
    summary(mod)
    l1<-vector('list',(nv2*nv3))
    nomes=expand.grid(names(summary(Fator2)),
                      names(summary(Fator3)))
    names(l1)<-paste(nomes$Var1,nomes$Var2)
    v<-numeric(0)
    for(j in 1:(nv3*nv2)) {
      for(i in 0:(nv1-2)) v<-cbind(v,i*nv3*nv2+j)
      l1[[j]]<-v
      v<-numeric(0)
    }
    dtf2a=data.frame(summary(mod,split=list('Fator3:Fator2:Fator1'=l1))[[1]])
    nl=nrow(dtf2a)
    dtf2=dtf2a[-c(1,2,3,nl),]
    nline=nrow(dtf2)
    for(i in 1:nline){
      dtf2$F.value[i]=dtf2$Mean.Sq[i]/qmresf1f3
      dtf2$Pr..F.[i]=1-pf(dtf2$F.value[i],dtf2$Df[i],nf1f3)
    }
    f11=dtf2a[3,]
    desd=rbind(f11,
               dtf2,
               c(nf2f3,qmresf1f3/nf1f3,qmresf1f3,NA,NA))
    nline1=nrow(desd)
    rownames(desd)[nline1]="Residuals combined"
    colnames(desd)=c("Df","Sum sq","Mean Sq", "F value","Pr(>F)")
    print(as.matrix(desd),na.print="")

    ii<-0
    for(i in 1:nv2) {
      for(j in 1:nv3) {
        ii<-ii+1
        cat("\n\n------------------------------------------")
        cat('\n',fac.names[1],' within the combination of levels ',lf2[i],' of  ',fac.names[2],' and ',lf3[j],' of  ',fac.names[3],"\n")
        cat("------------------------------------------\n")
        if(mcomp=="tukey"){
        tukey=TUKEY(response[fatores[,2]==lf2[i] & fatores[,3]==lf3[j]],
                       fatores[,1][Fator2==lf2[i] & Fator3==lf3[j]],
                       nf1f3,
                       qmresf1f3,
                       alpha.t)
        tukey=tukey$groups;colnames(tukey)=c("resp","letters")
        print(tukey)}
        if(mcomp=="lsd"){
          lsd=LSD(response[fatores[,2]==lf2[i] & fatores[,3]==lf3[j]],
                         fatores[,1][Fator2==lf2[i] & Fator3==lf3[j]],
                         nf1f3,
                         qmresf1f3,
                         alpha.t)
          lsd=lsd$groups;colnames(lsd)=c("resp","letters")
          print(lsd)}
        if(mcomp=="duncan"){
          duncan=duncan(response[fatores[,2]==lf2[i] & fatores[,3]==lf3[j]],
                         fatores[,1][Fator2==lf2[i] & Fator3==lf3[j]],
                         nf1f3,
                         qmresf1f3,
                         alpha.t)
          duncan=duncan$groups;colnames(duncan)=c("resp","letters")
          print(duncan)}
        if(mcomp=="sk"){
          respi=response[fatores[,2]==lf2[i] & fatores[,3]==lf3[j]]
          trati=fatores[,1][Fator2==lf2[i] & Fator3==lf3[j]]
          nrep=table(trati)[1]
          medias=sort(tapply(respi,trati,mean),decreasing = TRUE)
          sk=scottknott(means = medias,
                        df1 = nf1f3,
                        nrep = nrep,
                        QME = qmresf1f3,
                        alpha = alpha.t)
          sk=data.frame(respi=medias,groups=sk)
          # sk=sk(response[fatores[,2]==lf2[i] & fatores[,3]==lf3[j]],
          #               fatores[,1][Fator2==lf2[i] & Fator3==lf3[j]],
          #               nf1f3,
          #               qmresf1f3,
          #               alpha.t);colnames(sk)=c("resp","letters")
          print(sk)}

        }
    }

    cat('\n\n')

    cat(green(bold("\n------------------------------------------\n")))
    cat("Analyzing ", fac.names[1], ' inside of each level of ', fac.names[2], 'and',fac.names[3])
    cat(green(bold("\n------------------------------------------\n")))

    # desdobrando f1
    qmresf2f3=(qmres[2]+(nv3-1)*qmres[3])/nv3
    nf2f3=((qmres[2]+(nv3-1)*qmres[3])^2)/
      ((qmres[2]^2)/GL[2]+(((nv3-1)*qmres[3])^2)/GL[3])
    nf2f3=round(nf2f3)
    mod=aov(response~Fator3/Fator1/Fator2)
    summary(mod)
    l2<-vector('list',(nv1*nv3))
    nomes=expand.grid(names(summary(Fator1)),
                      names(summary(Fator3)))
    names(l2)<-paste(nomes$Var1,nomes$Var2)
    v<-numeric(0)
    for(j in 1:(nv3*nv1)) {
      for(i in 0:(nv2-2)) v<-cbind(v,i*nv1*nv3+j)
      l2[[j]]<-v
      v<-numeric(0)
    }
    dtf1a=data.frame(summary(mod,split=list('Fator3:Fator1:Fator2'=l2))[[1]])
    nl=nrow(dtf1a)
    dtf1=dtf1a[-c(1,2,3,nl),]
    nline=nrow(dtf1)
    for(i in 1:nline){
      dtf1$F.value[i]=dtf1$Mean.Sq[i]/qmresf2f3
      dtf1$Pr..F.[i]=1-pf(dtf1$F.value[i],dtf1$Df[i],nf2f3)
    }
    f11=dtf1a[3,]
    desd=rbind(f11,
               dtf1,
               c(nf2f3,qmresf2f3/nf2f3,qmresf2f3,NA,NA))
    nline1=nrow(desd)
    rownames(desd)[nline1]="Residuals combined"
    colnames(desd)=c("Df","Sum sq","Mean Sq", "F value","Pr(>F)")
    print(as.matrix(desd),na.print = "")

    ii<-0
    for(k in 1:nv1) {
      for(j in 1:nv3) {
        ii<-ii+1
        cat("\n\n------------------------------------------")
        cat('\n',fac.names[2],' within the combination of levels ',lf1[k],' of  ',fac.names[1],' and ',lf3[j],' of  ',fac.names[3],'\n')
        cat("------------------------------------------\n")
        if(mcomp=="tukey"){
        tukey=TUKEY(response[fatores[,1]==lf1[k] & fatores[,3]==lf3[j]],
                       fatores[,2][Fator1==lf1[k] & fatores[,3]==lf3[j]],
                       nf2f3,
                       qmresf2f3,
                       alpha.t)
        tukey=tukey$groups;colnames(tukey)=c("resp","letters")
        print(tukey)}
        if(mcomp=="lsd"){
          lsd=LSD(response[fatores[,1]==lf1[k] & fatores[,3]==lf3[j]],
                         fatores[,2][Fator1==lf1[k] & fatores[,3]==lf3[j]],
                         nf2f3,
                         qmresf2f3,
                         alpha.t)
          lsd=lsd$groups;colnames(lsd)=c("resp","letters")
          print(lsd)}
        if(mcomp=="duncan"){
          duncan=duncan(response[fatores[,1]==lf1[k] & fatores[,3]==lf3[j]],
                         fatores[,2][Fator1==lf1[k] & fatores[,3]==lf3[j]],
                         nf2f3,
                         qmresf2f3,
                         alpha.t)
          duncan=duncan$groups;colnames(duncan)=c("resp","letters")
          print(duncan)}
        if(mcomp=="sk"){
          respi=response[fatores[,1]==lf1[k] & fatores[,3]==lf3[j]]
          trati=fatores[,2][Fator1==lf1[k] & Fator3==lf3[j]]
          nrep=table(trati)[1]
          medias=sort(tapply(respi,trati,mean),decreasing = TRUE)
          sk=scottknott(means = medias,
                        df1 = nf2f3,
                        nrep = nrep,
                        QME = qmresf2f3,
                        alpha = alpha.t)
          sk=data.frame(respi=medias,groups=sk)

          # sk=sk(response[fatores[,1]==lf1[k] & fatores[,3]==lf3[j]],
          #               fatores[,2][Fator1==lf1[k] & fatores[,3]==lf3[j]],
          #               nf2f3,
          #               qmresf2f3/nf2f3,
          #               alpha.t);colnames(sk)=c("resp","letters")
          print(sk)}
      }
    }

    cat(green(bold("\n------------------------------------------\n")))
    cat("Analyzing ", fac.names[3], ' inside of each level of ', fac.names[1], 'and',fac.names[2])
    cat(green(bold("\n------------------------------------------\n")))
    mod=aov(response~Fator1/Fator2/Fator3+bloco+Error(bloco/fator1/fator2))
    summary(mod)
    l3<-vector('list',(nv2*nv1))
    nomes=expand.grid(names(summary(Fator1)),
                      names(summary(Fator2)))
    names(l3)<-paste(nomes$Var1,nomes$Var2)
    v<-numeric(0)
    for(j in 1:(nv1*nv2)) {
      for(i in 0:(nv3-2)) v<-cbind(v,i*nv2*nv1+j)
      l3[[j]]<-v
      v<-numeric(0)
    }
    dtf33=summary(mod,split=list('Fator1:Fator2:Fator3'=l3))
    dtf3=data.frame(dtf33$`Error: Within`[[1]])
    colnames(dtf3)=c("Df","Sum sq","Mean Sq", "F value","Pr(>F)")
    print(as.matrix(dtf3),na.print = "")

    ii<-0
    for(k in 1:nv1) {
      for(i in 1:nv2) {
        ii<-ii+1
        cat("\n\n------------------------------------------")
        cat('\n',fac.names[3],' within the combination of levels ',lf1[k],' of ',fac.names[1],' and ',lf2[i],' of  ',fac.names[2],'\n')
        cat("------------------------------------------\n")
        if(mcomp=="tukey"){
        tukey=TUKEY(response[fatores[,1]==lf1[k] & fatores[,2]==lf2[i]],
                       fatores[,3][fatores[,1]==lf1[k] & fatores[,2]==lf2[i]],
                       GL[3],
                       qmres[3],
                       alpha.t)
        tukey=tukey$groups;colnames(tukey)=c("resp","letters")
        print(tukey)}
        if(mcomp=="lsd"){
          lsd=LSD(response[fatores[,1]==lf1[k] & fatores[,2]==lf2[i]],
                         fatores[,3][fatores[,1]==lf1[k] & fatores[,2]==lf2[i]],
                         GL[3],
                         qmres[3],
                         alpha.t)
          lsd=lsd$groups;colnames(lsd)=c("resp","letters")
          print(lsd)}
        if(mcomp=="duncan"){
          duncan=duncan(response[fatores[,1]==lf1[k] & fatores[,2]==lf2[i]],
                         fatores[,3][fatores[,1]==lf1[k] & fatores[,2]==lf2[i]],
                         GL[3],
                         qmres[3],
                         alpha.t)
          duncan=duncan$groups;colnames(duncan)=c("resp","letters")
          print(duncan)}
        if(mcomp=="sk"){
        respi=response[fatores[,1]==lf1[k] & fatores[,2]==lf2[i]]
        trati=fatores[,3][fatores[,1]==lf1[k] & fatores[,2]==lf2[i]]
        nrep=table(trati)[1]
        medias=sort(tapply(respi,trati,mean),decreasing = TRUE)
        sk=scottknott(means = medias,
                      df1 = nf1f3,
                      nrep = nrep,
                      QME = qmresf1f3,
                      alpha = alpha.t)
        sk=data.frame(respi=medias,groups=sk)
          #
          # sk=sk(response[fatores[,1]==lf1[k] & fatores[,2]==lf2[i]],
          #               fatores[,3][fatores[,1]==lf1[k] & fatores[,2]==lf2[i]],
          #               GL[3],
          #               qmres[3]/GL[3],
          #               alpha.t)
          # colnames(sk)=c("resp","letters")
          print(sk)}

        }
    }
  }
}
