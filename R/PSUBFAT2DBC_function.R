#' Analysis: Plot subdivided into randomized blocks with a subplot in a double factorial scheme
#'
#' @description This function performs the analysis of a randomized block design in a split-plot with a subplot in a double factorial scheme.
#'
#' @param f1 Numeric or complex vector with plot levels
#' @param f2 Numeric or complex vector with splitplot levels
#' @param f3 Numeric or complex vector with splitsplitplot levels
#' @param block Numeric or complex vector with blocks
#' @param mcomp Multiple comparison test (Tukey (\emph{default}), LSD and Duncan)
#' @param resp Numeric vector with responses
#' @param alpha.f Level of significance of the F test (\emph{default} is 0.05)
#' @param alpha.t Significance level of the multiple comparison test (\emph{default} is 0.05)
#'
#' @export
#' @return Analysis of variance of fixed effects and multiple comparison test of Tukey, Scott-Knott, LSD or Duncan.
#'
#' @examples
#' f1=rep(c("PD","PDE","C"), e = 40);f1=factor(f1,unique(f1))
#' f2=rep(c(300,400), e = 20,3);f2=factor(f2,unique(f2))
#' f3=rep(c("c1", "c2", "c3", "c4"), e = 5,6);f3=factor(f3,unique(f3))
#' bloco=rep(paste("B",1:5),24); bloco=factor(bloco,unique(bloco))
#' set.seed(10)
#' resp=rnorm(120,50,5)
#' PSUBFAT2DBC(f1,f2,f3,bloco,resp,alpha.f = 0.5) # para soltar interação tripla
#' PSUBFAT2DBC(f1,f2,f3,bloco,resp,alpha.f = 0.4) # para soltar interação dupla

PSUBFAT2DBC=function(f1,
                     f2,
                     f3,
                     block,
                     resp,
                     alpha.f=0.05,
                     alpha.t=0.05,
                     mcomp="tukey"){
  requireNamespace("crayon")
  fac.names=c("F1","F2","F3")
  Fator1=f1
  Fator2=f2
  Fator3=f3
  f1<-factor(f1,levels=unique(f1))
  f2<-factor(f2,levels=unique(f2))
  f3<-factor(f3,levels=unique(f3))
  nv1=length(levels(f1))
  nv2 = length(levels(f2))
  nv3 = length(levels(f3))
  bloco=factor(block,levels = unique(block))
  nbl<-length(summary(bloco))
  j<-(length(resp))/(nv1*nv2*nv3)
  lf1<-levels(f1); lf2<-levels(f2); lf3<-levels(f3)

  mod=aov(resp~f1 * f2 * f3 + Error(bloco:f1)+bloco)
  a = summary(mod)
  anava = rbind(data.frame(a$`Error: bloco:f1`[[1]]),
                data.frame(a$`Error: Within`[[1]]))
  rownames(anava) = c("F1",
                      "block",
                      "Error A",
                      "F2",
                      "F3",
                      "F1 x F2",
                      "F1 x F3",
                      "F2 x F3",
                      "F1 x F2 x F3",
                      "Residuals")
  colnames(anava) = c("df", "SS", "MS", "F-value", "p")

  cat("\n==================================\n")
  cat("Analysis of variance")
  cat("\n==================================\n")
  print(as.matrix(anava), na.print = "")
  qmres = c(anava$MS[3], anava$MS[10])
  GLres = c(anava$df[3], anava$df[10])
  if (anava$p[9] > alpha.f) {
    if (anava$p[1] < alpha.f &
        anava$p[6] > alpha.f &
        anava$p[7] > alpha.f){
      if(mcomp=="tukey"){comp=TUKEY(resp, f1, DFerror = anava$df[3], MSerror = anava$MS[3])
      comp=comp$groups
      colnames(comp)=c("resp","groups")}
      if(mcomp=="lsd"){comp=LSD(resp, f1, DFerror = anava$df[3], MSerror = anava$MS[3])
      comp=comp$groups
      colnames(comp)=c("resp","groups")}
      if(mcomp=="duncan"){comp=duncan(resp, f1, DFerror = anava$df[3], MSerror = anava$MS[3])
      comp=comp$groups
      colnames(comp)=c("resp","groups")}
      if(mcomp=="sk"){
        nrep=table(f1)[1]
        medias=sort(tapply(resp,f1,mean),decreasing = TRUE)
        comp=scottknott(medias,alpha = alpha.t, df1 = anava$df[3], nrep = nrep,
                        QME = anava$MS[3])
        comp=data.frame(resp=medias,groups=comp)}
      cat("\n==================================\n")
      cat("F1")
      cat("\n==================================\n")
      print(comp)}

    if (anava$p[4] < alpha.f&
        anava$p[6] > alpha.f&
        anava$p[8] > alpha.f){
      if(mcomp=="tukey"){comp=TUKEY(resp, f2, DFerror = anava$df[10], MSerror = anava$MS[10])
      comp=comp$groups
      colnames(comp)=c("resp","groups")}
      if(mcomp=="lsd"){comp=LSD(resp, f2, DFerror = anava$df[10], MSerror = anava$MS[10])
      comp=comp$groups
      colnames(comp)=c("resp","groups")}
      if(mcomp=="duncan"){comp=duncan(resp, f2, DFerror = anava$df[10], MSerror = anava$MS[10])
      comp=comp$groups
      colnames(comp)=c("resp","groups")}
      if(mcomp=="sk"){
        nrep=table(f2)[1]
        medias=sort(tapply(resp,f2,mean),decreasing = TRUE)
        comp=scottknott(medias,alpha = alpha.t, df1 = anava$df[10],
                        nrep = nrep, QME = anava$MS[10])
        comp=data.frame(resp=medias,groups=comp)}
      cat("\n==================================\n")
      cat("F2")
      cat("\n==================================\n")
      print(comp)}
    if (anava$p[5] < alpha.f&
        anava$p[7] > alpha.f&
        anava$p[8] > alpha.f){
      if(mcomp=="tukey"){comp=TUKEY(resp, f3, DFerror = anava$df[10], MSerror = anava$MS[10])
      comp=comp$groups
      colnames(comp)=c("resp","groups")}
      if(mcomp=="lsd"){comp=LSD(resp, f3, DFerror = anava$df[10], MSerror = anava$MS[10])
      comp=comp$groups
      colnames(comp)=c("resp","groups")}
      if(mcomp=="duncan"){comp=duncan(resp, f3, DFerror = anava$df[10], MSerror = anava$MS[10])
      comp=comp$groups
      colnames(comp)=c("resp","groups")}
      if(mcomp=="sk"){
        nrep=table(f3)[1]
        medias=sort(tapply(resp,f3,mean),decreasing = TRUE)
        comp=scottknott(medias,alpha = alpha.t, df1 = anava$df[10],
                        nrep = nrep, QME = anava$MS[10])
        comp=data.frame(resp=medias,groups=comp)}
      cat("\n==================================\n")
      cat("F3")
      cat("\n==================================\n")
      print(comp)}
    cat("\n")
  }

  #=============================================
  # desdobramento de f1 x f2
  #=============================================
  if (anava$p[9] > alpha.f &
      anava$p[6] < alpha.f) {
    mod = aov(resp ~ f1 / f2 + f2:f3 + f1:f3 + f1:f2:f3 + Error(bloco:f1)+bloco)
    l2 <- vector('list', nv1)
    names(l2) <- names(summary(f1))
    v <- numeric(0)
    for (j in 1:nv1) {
      for (i in 0:(nv2 - 2))
        v <- cbind(v, i * nv1 + j)
      l2[[j]] <- v
      v <- numeric(0)
    }
    des1.tab <- summary(mod, split = list('f1:f2' = l2))
    desdf2f1 = data.frame(des1.tab$`Error: Within`[[1]])
    colnames(desdf2f1) = c("Df", "Sum sq", "Mean Sq", "F value", "Pr(>F)")
    nlin = nrow(desdf2f1)
    desdf2f1 = desdf2f1[-c(nlin - 1, nlin - 2, nlin - 3), ]
    cat(green(bold("\n-----------------------------------------------------\n")))
    cat("Analyzing ", fac.names[2], ' inside of each level of ', fac.names[1])
    cat(green(bold("\n-----------------------------------------------------\n")))
    print(as.matrix(desdf2f1), na.print = "")
    compf1f2=c()
    letterf1f2=c()
    for (i in 1:nv1) {
      trat1 = f2[f1 == levels(f1)[i]]
      resp1 = resp[f1 == levels(f1)[i]]
      nrep=table(trat1)[1]
      if(mcomp=="tukey"){comp = TUKEY(resp1, trat1, DFerror = anava$df[10], MSerror = anava$MS[10])
      comp=comp$groups
      colnames(comp)=c("resp","groups")}
      if(mcomp=="lsd"){comp = LSD(resp1, trat1, DFerror = anava$df[10], MSerror = anava$MS[10])
      comp=comp$groups
      colnames(comp)=c("resp","groups")}
      if(mcomp=="duncan"){comp = duncan(resp1, trat1, DFerror = anava$df[10], MSerror = anava$MS[10])
      comp=comp$groups
      colnames(comp)=c("resp","groups")}
      if(mcomp=="sk"){
        medias=sort(tapply(resp1,trat1,mean),decreasing = TRUE)
        comp = scottknott(medias,df1 = anava$df[10], QME = anava$MS[10],nrep = nrep)
        comp=data.frame(resp=medias,groups=comp)}
      comp=comp[unique(as.character(f2)),]
      compf1f2[[i]]=comp$resp
      letterf1f2[[i]]=comp$groups
      # cat("\n======================\n")
      # cat(levels(f1)[i])
      # cat("\n======================\n")
      # print(comp)
    }
    #===========================================
    mod = aov(resp ~ f2 / f1 + f2:f3 + f1:f3 + f1:f2:f3 + Error(bloco / f2))
    summary(mod)
    l1 <- vector('list', nv2)
    names(l1) <- names(summary(f2))
    v <- numeric(0)
    for (j in 1:nv2) {
      for (i in 0:(nv1 - 2))
        v <- cbind(v, i * nv2 + j)
      l1[[j]] <- v
      v <- numeric(0)
    }
    desd1.tab <- summary(mod, split = list('f2:f1' = l1))
    desd = data.frame(desd1.tab$`Error: Within`[[1]])
    nlinhas = nrow(desd)
    desd = desd[-c(1, nlinhas - 3, nlinhas - 2, nlinhas - 1, nlinhas), ]
    qmresf1f2 = (qmres[1] + (nv2 - 1) * qmres[2]) / nv2
    nf1f2 = ((qmres[1] + (nv2 - 1) * qmres[2]) ^ 2) /
      ((qmres[1] ^ 2) / GLres[1] + (((nv2 - 1) * qmres[2]) ^ 2) / GLres[2])
    nf1f2 = round(nf1f2)
    desd$F.value = desd$Mean.Sq / qmresf1f2
    nline = nrow(desd)
    for (i in 1:nline) {
      desd$Pr..F.[i] = 1 - pf(desd$F.value[i], desd$Df[i], nf1f2)
    }
    f1f2 = data.frame(desd1.tab$`Error: Within`[[1]])[1, ]
    desdf1f2 = rbind(f1f2, desd, c(nf1f2, qmresf1f2 / nf1f2, qmresf1f2, NA, NA))
    nline1 = nrow(desdf1f2)
    rownames(desdf1f2)[nline1] = "Residuals combined"
    colnames(desdf1f2) = c("Df", "Sum sq", "Mean Sq", "F value", "Pr(>F)")
    cat(green(bold("\n-----------------------------------------------------\n")))
    cat("Analyzing ", fac.names[1], ' inside of each level of ', fac.names[2])
    cat(green(bold("\n-----------------------------------------------------\n")))
    print(as.matrix(desdf1f2), na.print = "")
    compf2f1=c()
    letterf2f1=c()
    for (i in 1:nv2) {
      trat1 = f1[f2 == levels(f2)[i]]
      resp1 = resp[f2 == levels(f2)[i]]
      if(mcomp=="tukey"){comp = TUKEY(resp1, trat1, DFerror = desdf1f2[nrow(desdf1f2), 1],
                                      MSerror = desdf1f2[nrow(desdf1f2), 3])
      comp=comp$groups
      colnames(comp)=c("resp","groups")}
      if(mcomp=="lsd"){comp = LSD(resp1, trat1, DFerror = desdf1f2[nrow(desdf1f2), 1],
                                  MSerror = desdf1f2[nrow(desdf1f2), 3])
      comp=comp$groups
      colnames(comp)=c("resp","groups")}
      if(mcomp=="duncan"){comp = duncan(resp1, trat1, DFerror = desdf1f2[nrow(desdf1f2), 1],
                                        MSerror = desdf1f2[nrow(desdf1f2), 3])
      comp=comp$groups
      colnames(comp)=c("resp","groups")}
      if(mcomp=="sk"){
        medias=sort(tapply(resp1,trat1,mean),decreasing = TRUE)
        comp = scottknott(medias,df1 = desdf1f2[nrow(desdf1f2), 1],
                          QME = desdf1f2[nrow(desdf1f2), 3],nrep = nrep)
        comp=data.frame(resp=medias,groups=comp)}
      comp=comp[unique(as.character(f1)),]
      compf2f1[[i]]=comp$resp
      letterf2f1[[i]]=comp$groups

      # cat("\n======================\n")
      # cat(levels(f2)[i])
      # cat("\n======================\n")
      # print(comp)
    }
    final=paste(round(unlist(compf1f2),3),
                paste(unlist(letterf1f2),
                      toupper(c(t(matrix(unlist(letterf2f1),ncol=length(levels(f2)))))),sep = ""))
    final=data.frame(matrix(final,ncol=length(unique(f1))))
    colnames(final)=as.character(unique(f1))
    rownames(final)=as.character(unique(f2))

    # cat("\n======================\n")
    # cat("Multiple comparasion")
    # cat("\n======================\n")
    # print(final)
  }

  #=========================================================
  # desdobramento de f1 x f3
  #=========================================================
  if (anava$p[9] > alpha.f &
      anava$p[7] < alpha.f) {
    mod = aov(resp ~ f1 / f3 + f2:f3 + f1:f2 + f1:f2:f3 + Error(bloco / f1))
    l3 <- vector('list', nv1)
    names(l3) <- names(summary(f1))
    v <- numeric(0)
    for (j in 1:nv1) {
      for (i in 0:(nv3 - 2))
        v <- cbind(v, i * nv1 + j)
      l3[[j]] <- v
      v <- numeric(0)
    }
    des1.tab <- summary(mod, split = list('f1:f3' = l3))
    desdf3f1 = data.frame(des1.tab$`Error: Within`[[1]])
    colnames(desdf3f1) = c("Df", "Sum sq", "Mean Sq", "F value", "Pr(>F)")
    nlin = nrow(desdf3f1)
    desdf3f1 = desdf3f1[-c(nlin - 1, nlin - 2, nlin - 3), ]
    cat(green(bold("\n-----------------------------------------------------\n")))
    cat("Analyzing ", fac.names[3], ' inside of each level of ', fac.names[1])
    cat(green(bold("\n-----------------------------------------------------\n")))
    print(as.matrix(desdf3f1), na.print = "")

    compf1f3=c()
    letterf1f3=c()
    for (i in 1:nv1) {
      trat1 = f3[f1 == levels(f1)[i]]
      resp1 = resp[f1 == levels(f1)[i]]
      nrep=table(trat1)[1]
      if(mcomp=="tukey"){comp = TUKEY(resp1, trat1, DFerror = anava$df[10],
                                      MSerror = anava$MS[10])
      comp=comp$groups
      colnames(comp)=c("resp","groups")}
      if(mcomp=="lsd"){comp = LSD(resp1, trat1, DFerror = anava$df[10],
                                  MSerror = anava$MS[10])
      comp=comp$groups
      colnames(comp)=c("resp","groups")}
      if(mcomp=="duncan"){comp = duncan(resp1, trat1, DFerror = anava$df[10],
                                        MSerror = anava$MS[10])
      comp=comp$groups
      colnames(comp)=c("resp","groups")}
      if(mcomp=="sk"){
        medias=sort(tapply(resp1,trat1,mean),decreasing = TRUE)
        comp = scottknott(medias,df1 = anava$df[10],
                          QME = anava$MS[10],nrep = nrep)
        comp=data.frame(resp=medias,groups=comp)}
      comp=comp[unique(as.character(f3)),]
      compf1f3[[i]]=comp$resp
      letterf1f3[[i]]=comp$groups

      # cat("\n======================\n")
      # cat(levels(f1)[i])
      # cat("\n======================\n")
      # print(comp)
    }

    #===========================================
    mod = aov(resp ~ f3 / f1 + f1:f2 + f2:f3 + f1:f2:f3 + Error(bloco / f3))
    summary(mod)
    l1 <- vector('list', nv3)
    names(l1) <- names(summary(f3))
    v <- numeric(0)
    for (j in 1:nv3) {
      for (i in 0:(nv1 - 2))
        v <- cbind(v, i * nv3 + j)
      l1[[j]] <- v
      v <- numeric(0)
    }
    desd1.tab <- summary(mod, split = list('f3:f1' = l1))
    desd = data.frame(desd1.tab$`Error: Within`[[1]])
    nlinhas = nrow(desd)
    desd = desd[-c(1, nlinhas - 3, nlinhas - 2, nlinhas - 1, nlinhas), ]
    qmresf1f3 = (qmres[1] + (nv3 - 1) * qmres[2]) / nv3
    nf1f3 = ((qmres[1] + (nv3 - 1) * qmres[2]) ^ 2) /
      ((qmres[1] ^ 2) / GLres[1] + (((nv3 - 1) * qmres[2]) ^ 2) / GLres[2])
    nf1f3 = round(nf1f3)
    desd$F.value = desd$Mean.Sq / qmresf1f3
    nline = nrow(desd)
    for (i in 1:nline) {
      desd$Pr..F.[i] = 1 - pf(desd$F.value[i], desd$Df[i], nf1f3)
    }
    f1f3 = data.frame(desd1.tab$`Error: Within`[[1]])[1, ]
    desdf1f3 = rbind(f1f3, desd, c(nf1f3, qmresf1f3 / nf1f3, qmresf1f3, NA, NA))
    nline1 = nrow(desdf1f3)
    rownames(desdf1f3)[nline1] = "Residuals combined"
    colnames(desdf1f3) = c("Df", "Sum sq", "Mean Sq", "F value", "Pr(>F)")
    cat(green(bold("\n-----------------------------------------------------\n")))
    cat("Analyzing ", fac.names[1], ' inside of each level of ', fac.names[3])
    cat(green(bold("\n-----------------------------------------------------\n")))
    print(as.matrix(desdf1f3), na.print = "")
    compf3f1=c()
    letterf3f1=c()

    for (i in 1:nv3) {
      trat1 = f1[f3 == levels(f3)[i]]
      resp1 = resp[f3 == levels(f3)[i]]
      nrep=table(trat1)[1]
      if(mcomp=="tukey"){comp = TUKEY(resp1, trat1, DFerror = desdf1f3[nrow(desdf1f3), 1],
                                      MSerror = desdf1f3[nrow(desdf1f3), 3])
      comp=comp$groups
      colnames(comp)=c("resp","groups")}
      if(mcomp=="lsd"){comp = LSD(resp1, trat1, DFerror = desdf1f3[nrow(desdf1f3), 1],
                                  MSerror = desdf1f3[nrow(desdf1f3), 3])
      comp=comp$groups
      colnames(comp)=c("resp","groups")}
      if(mcomp=="duncan"){comp = duncan(resp1, trat1, DFerror = desdf1f3[nrow(desdf1f3), 1],
                                        MSerror = desdf1f3[nrow(desdf1f3), 3])
      comp=comp$groups
      colnames(comp)=c("resp","groups")}
      if(mcomp=="sk"){
        medias=sort(tapply(resp1,trat1,mean),decreasing = TRUE)
        comp = scottknott(medias,df1 = desdf1f3[nrow(desdf1f3), 1],
                          QME = desdf1f3[nrow(desdf1f3), 3],nrep = nrep)
        comp=data.frame(resp=medias,groups=comp)}
      comp=comp[unique(as.character(f1)),]
      compf3f1[[i]]=comp$resp
      letterf3f1[[i]]=comp$groups
      final=paste(round(unlist(compf1f3),3),
                  paste(unlist(letterf1f3),
                        toupper(c(t(matrix(unlist(letterf3f1),ncol=length(levels(f3)))))),sep = ""))
      final=data.frame(matrix(final,ncol=length(unique(f1))))
      colnames(final)=as.character(unique(f1))
      rownames(final)=as.character(unique(f3))
    }
    cat("\n======================\n")
    cat("Multiple comparasion")
    cat("\n======================\n")
    print(final)
  }

  # desdobramento de f2 x f3
  if (anava$p[9] > alpha.f &
      anava$p[8] < alpha.f) {
    mod = aov(resp ~ f2 / f3 + f1:f2 + f1:f3 + f1:f2:f3 + Error(bloco / f1))
    l3 <- vector('list', nv2)
    names(l3) <- names(summary(f2))
    v <- numeric(0)
    for (j in 1:nv2) {
      for (i in 0:(nv3 - 2))
        v <- cbind(v, i * nv2 + j)
      l3[[j]] <- v
      v <- numeric(0)
    }
    des1.tab <- summary(mod, split = list('f2:f3' = l3))
    desdf3f2 = data.frame(des1.tab$`Error: Within`[[1]])
    colnames(desdf3f2) = c("Df", "Sum sq", "Mean Sq", "F value", "Pr(>F)")
    nlin = nrow(desdf3f2)
    desdf3f2 = desdf3f2[-c(nlin - 1, nlin - 2, nlin - 3), ]
    cat(green(bold("\n-----------------------------------------------------\n")))
    cat("Analyzing ", fac.names[3], ' inside of each level of ', fac.names[2])
    cat(green(bold("\n-----------------------------------------------------\n")))
    print(as.matrix(desdf3f2), na.print = "")
    compf2f3=c()
    letterf2f3=c()
    for (i in 1:nv2) {
      trat1 = f3[f2 == levels(f2)[i]]
      resp1 = resp[f2 == levels(f2)[i]]
      nrep=table(trat1)[1]
      if(mcomp=="tukey"){comp = TUKEY(resp1, trat1, DFerror = anava$df[10],
                                      MSerror = anava$MS[10])
      comp=comp$groups
      colnames(comp)=c("resp","groups")}
      if(mcomp=="lsd"){comp = LSD(resp1, trat1, DFerror = anava$df[10],
                                  MSerror = anava$MS[10])
      comp=comp$groups
      colnames(comp)=c("resp","groups")}
      if(mcomp=="duncan"){comp = duncan(resp1, trat1, DFerror = anava$df[10],
                                        MSerror = anava$MS[10])
      comp=comp$groups
      colnames(comp)=c("resp","groups")}
      if(mcomp=="sk"){
        medias=sort(tapply(resp1,trat1,mean),decreasing = TRUE)
        comp = scottknott(medias,df1 = anava$df[10],
                          QME = anava$MS[10],nrep = nrep)
        comp=data.frame(resp=medias,groups=comp)}
      comp=comp[unique(as.character(f3)),]
      compf2f3[[i]]=comp$resp
      letterf2f3[[i]]=comp$groups
    }

    mod = aov(resp ~ f3 / f2 + f1:f2 + f1:f3 + f1:f2:f3 + Error(bloco / f1))
    l2 <- vector('list', nv3)
    names(l2) <- names(summary(f3))
    v <- numeric(0)
    for (j in 1:nv3) {
      for (i in 0:(nv2 - 2))
        v <- cbind(v, i * nv3 + j)
      l2[[j]] <- v
      v <- numeric(0)
    }
    des1.tab <- summary(mod, split = list('f3:f2' = l2))
    desdf2f3 = data.frame(des1.tab$`Error: Within`[[1]])
    colnames(desdf2f3) = c("Df", "Sum sq", "Mean Sq", "F value", "Pr(>F)")
    nlin = nrow(desdf2f3)
    desdf2f3 = desdf2f3[-c(nlin - 1, nlin - 2, nlin - 3), ]
    cat(green(bold("\n-----------------------------------------------------\n")))
    cat("Analyzing ", fac.names[2], ' inside of each level of ', fac.names[3])
    cat(green(bold("\n-----------------------------------------------------\n")))
    print(as.matrix(desdf2f3), na.print = "")
    compf3f2=c()
    letterf3f2=c()
    for (i in 1:nv3) {
      trat1 = f2[f3 == levels(f3)[i]]
      resp1 = resp[f3 == levels(f3)[i]]
      nrep=table(trat1)[1]
      if(mcomp=="tukey"){comp = TUKEY(resp1, trat1, DFerror = anava$df[10],
                                      MSerror = anava$MS[10])
      comp=comp$groups
      colnames(comp)=c("resp","groups")}
      if(mcomp=="lsd"){comp = LSD(resp1, trat1, DFerror = anava$df[10],
                                  MSerror = anava$MS[10])
      comp=comp$groups
      colnames(comp)=c("resp","groups")}
      if(mcomp=="duncan"){comp = duncan(resp1, trat1, DFerror = anava$df[10],
                                        MSerror = anava$MS[10])
      comp=comp$groups
      colnames(comp)=c("resp","groups")}
      if(mcomp=="sk"){
        medias=sort(tapply(resp1,trat1,mean),decreasing = TRUE)
        comp = scottknott(medias,df1 = anava$df[10],
                          QME = anava$MS[10],nrep = nrep)
        comp=data.frame(resp=medias,groups=comp)}
      comp=comp[unique(as.character(f2)),]
      compf3f2[[i]]=comp$resp
      letterf3f2[[i]]=comp$groups
      final=paste(round(unlist(compf2f3),3),
                  paste(unlist(letterf2f3),
                        toupper(c(t(matrix(unlist(letterf3f2),ncol=length(levels(f3)))))),sep = ""))
      final=data.frame(matrix(final,ncol=length(unique(f2))))
      colnames(final)=as.character(unique(f2))
      rownames(final)=as.character(unique(f3))
    }
    cat("\n======================\n")
    cat("Multiple comparasion")
    cat("\n======================\n")
    print(final)
  }

  #==================================================
  # desdobramento de f1 x f2 x f3
  #==================================================
  if (anava$p[9] < alpha.f) {
    # desdobramnto de f2
    m1=aov(resp~(f1*f3)/f2+Error(bloco/f1))
    summary(m1)
    pattern <- c(outer(levels(f1), levels(f3),
                       function(x,y) paste("f1",x,":f3",y,":",sep="")))
    des.tab <- sapply(pattern, simplify=FALSE,
                      grep, x=names(coef(m1$Within)[m1$Within$assign==4]))
    des1.tab <- summary(m1, split = list("f1:f3:f2" = des.tab))
    des1.tab=data.frame(des1.tab$`Error: Within`[[1]])
    nomes=expand.grid(levels(f1),levels(f3))
    nomes=paste(nomes$Var1,nomes$Var2)
    nomes=c("f3","f1:f3","f1:f3:f2",
            paste("   f1:f3:f2",nomes),"residuals")
    colnames(des1.tab) = c("Df", "Sum sq", "Mean Sq", "F value", "Pr(>F)")
    rownames(des1.tab)=nomes

    cat(green(bold("\n-----------------------------------------------------\n")))
    cat("Analyzing ", fac.names[2], ' inside of each level of ', fac.names[1], 'and',fac.names[3])
    cat(green(bold("\n-----------------------------------------------------\n")))
    print(as.matrix(des1.tab[-c(1,2),]), na.print = "")
    fatores=data.frame(f1,f2,f3)
    ii<-0
    for(k in 1:nv1) {
      for(i in 1:nv2) {
        ii<-ii+1
        cat("\n\n------------------------------------------")
        cat('\n',fac.names[2],' within the combination of levels ',lf1[k],' of ',fac.names[1],' and ',lf3[i],' of  ',fac.names[3],'\n')
        cat("------------------------------------------\n")
        respi=resp[fatores[,1]==lf1[k] & fatores[,3]==lf3[i]]
        trati=fatores[,2][Fator1==lf1[k] & fatores[,3]==lf3[i]]
        nlinhas=nrow(des1.tab)
        nrep=table(trati)[1]
        if(mcomp=="tukey"){comp = TUKEY(respi, trati,
                                        DFerror = des1.tab[nlinhas,1],
                                        MSerror = des1.tab[nlinhas,3])
        comp=comp$groups
        colnames(comp)=c("resp","groups")}
        if(mcomp=="lsd"){comp = LSD(respi, trati, DFerror = des1.tab[nlinhas,1],
                                    MSerror = des1.tab[nlinhas,3])
        comp=comp$groups
        colnames(comp)=c("resp","groups")}
        if(mcomp=="duncan"){comp = duncan(respi, trati,
                                          DFerror = des1.tab[nlinhas,1],
                                          MSerror = des1.tab[nlinhas,3])
        comp=comp$groups
        colnames(comp)=c("resp","groups")}
        if(mcomp=="sk"){
          medias=sort(tapply(respi,trati,mean),decreasing = TRUE)
          comp = scottknott(medias,df1 = des1.tab[nlinhas,1],
                            QME = des1.tab[nlinhas,3],nrep = nrep)
          comp=data.frame(resp=medias,groups=comp)}
        print(comp)}
    }

    ####################################################
    # dentro de f3, testa 1 dentro 2 - usar qm normal
    m1=aov(resp~(f1*f2)/f3)
    anova(m1)
    pattern <- c(outer(levels(f1), levels(f2),
                       function(x,y) paste("f1",x,":f2",y,":",sep="")))
    des.tab <- sapply(pattern, simplify=FALSE,
                      grep, x=names(coef(m1)[m1$assign==4]))
    des1.tab <- summary(m1, split = list("f1:f2:f3" = des.tab))
    des1.tab=data.frame(des1.tab[[1]])
    nomes=expand.grid(levels(f1),levels(f2))
    nomes=paste(nomes$Var1,nomes$Var2)
    nomes=c("f1","f2","f1:f2","f1:f2:f3",
            paste("   f1:f2:f3:",nomes),"residuals")
    colnames(des1.tab) = c("Df", "Sum sq", "Mean Sq", "F value", "Pr(>F)")
    rownames(des1.tab)=nomes
    cat(green(bold("\n-----------------------------------------------------\n")))
    cat("Analyzing ", fac.names[3], ' inside of each level of ', fac.names[1], 'and',fac.names[2])
    cat(green(bold("\n-----------------------------------------------------\n")))
    print(as.matrix(des1.tab[-c(1,2,3),]), na.print = "")
    ii<-0
    for(k in 1:nv1) {
      for(j in 1:nv2) {
        ii<-ii+1
        cat("\n\n------------------------------------------")
        cat('\n',fac.names[3],' within the combination of levels ',lf1[k],
            ' of  ',fac.names[1],' and ',lf2[j],' of  ',fac.names[2],'\n')
        cat("------------------------------------------\n")
        if(mcomp=="tukey"){
          trati=fatores[,3][fatores[,1]==lf1[k] & fatores[,2]==lf2[i]]
          respi=resp[fatores[,1]==lf1[k] & fatores[,2]==lf2[i]]
          nlinhas=nrow(des1.tab)
          if(mcomp=="tukey"){comp = TUKEY(respi, trati,
                                          DFerror = des1.tab[nlinhas,1],
                                          MSerror = des1.tab[nlinhas,3])
          comp=comp$groups
          colnames(comp)=c("resp","groups")}
          if(mcomp=="lsd"){comp = LSD(respi, trati, DFerror = des1.tab[nlinhas,1],
                                      MSerror = des1.tab[nlinhas,3])
          comp=comp$groups
          colnames(comp)=c("resp","groups")}
          if(mcomp=="duncan"){comp = duncan(respi, trati,
                                            DFerror = des1.tab[nlinhas,1],
                                            MSerror = des1.tab[nlinhas,3])
          comp=comp$groups
          colnames(comp)=c("resp","groups")}
          if(mcomp=="sk"){
            medias=sort(tapply(respi,trati,mean),decreasing = TRUE)
            comp = scottknott(medias,df1 = des1.tab[nlinhas,1],
                              QME = des1.tab[nlinhas,3],nrep = nrep)
            comp=data.frame(resp=medias,groups=comp)}
          print(comp)}
      }
    }
    #=======================
    # teste de f1, testa 2 e 3 - QM normal
    m1=aov(resp~(f2*f3)/f1)
    anova(m1)
    pattern <- c(outer(levels(f2), levels(f3),
                       function(x,y) paste("f2",x,":f3",y,":",sep="")))
    des.tab <- sapply(pattern, simplify=FALSE,
                      grep, x=names(coef(m1)[m1$assign==4]))
    des1.tab <- summary(m1, split = list("f2:f3:f1" = des.tab))
    nv23 = nv3 * nv2
    qmresf1f3 = (qmres[1] + (nv23 - 1) * qmres[2]) / nv23
    nf1f3 = ((qmres[1] + (nv23 - 1) * qmres[2]) ^ 2) /
      ((qmres[1] ^ 2) / GLres[1] + (((nv23 - 1) * qmres[2]) ^ 2) / GLres[2])
    nf1f3 = round(nf1f3)
    des1.tab=data.frame(des1.tab[[1]])
    nl = nrow(des1.tab)
    dtf2 = des1.tab[-c(1, 2, 3, 4,nl), ]
    nline = nrow(dtf2)
    for (i in 1:nline) {
      dtf2$F.value[i] = dtf2$Mean.Sq[i] / qmresf1f3
      dtf2$Pr..F.[i] = 1 - pf(dtf2$F.value[i], dtf2$Df[i], nf1f3)
    }
    f11 = dtf2[3, ]
    desd = rbind(f11,
                 dtf2,
                 c(nf1f3, qmresf1f3 * nf1f3, qmresf1f3, NA, NA))
    nline1 = nrow(desd)
    rownames(desd)[nline1] = "Residuals combined"
    colnames(desd) = c("Df", "Sum sq", "Mean Sq", "F value", "Pr(>F)")
    nomes=expand.grid(levels(f2),levels(f3))
    nomes=paste(nomes$Var1,nomes$Var2)
    nomes=c("f3:f2:f1",
            paste("   f3:f2:f1:",nomes),"Residuals combined")
    colnames(desd) = c("Df", "Sum sq", "Mean Sq", "F value", "Pr(>F)")
    rownames(desd)=nomes
    cat(green(bold("\n-----------------------------------------------------\n")))
    cat("Analyzing ", fac.names[1], ' inside of each level of ', fac.names[2], 'and',fac.names[3])
    cat(green(bold("\n-----------------------------------------------------\n")))
    print(as.matrix(desd), na.print = "")
    ii<-0
    for(i in 1:nv2) {
      for(j in 1:nv3) {
        ii<-ii+1
        cat("\n\n------------------------------------------")
        cat('\n',fac.names[1],' within the combination of levels ',lf2[i],' of  ',fac.names[2],' and ',lf3[j],' of  ',fac.names[3],"\n")
        cat("------------------------------------------\n")
        respi=resp[fatores[,2]==lf2[i] & fatores[,3]==lf3[j]]
        trati=fatores[,1][Fator2==lf2[i] & Fator3==lf3[j]]
        nlinhas=nrow(desd)
        if(mcomp=="tukey"){comp = TUKEY(respi, trati,
                                        DFerror = des1.tab[nlinhas,1],
                                        MSerror = des1.tab[nlinhas,3])
        comp=comp$groups
        colnames(comp)=c("resp","groups")}
        if(mcomp=="lsd"){comp = LSD(respi, trati, DFerror = des1.tab[nlinhas,1],
                                    MSerror = des1.tab[nlinhas,3])
        comp=comp$groups
        colnames(comp)=c("resp","groups")}
        if(mcomp=="duncan"){comp = duncan(respi, trati,
                                          DFerror = des1.tab[nlinhas,1],
                                          MSerror = des1.tab[nlinhas,3])
        comp=comp$groups
        colnames(comp)=c("resp","groups")}
        if(mcomp=="sk"){
          medias=sort(tapply(respi,trati,mean),decreasing = TRUE)
          comp = scottknott(medias,df1 = des1.tab[nlinhas,1],
                            QME = des1.tab[nlinhas,3],nrep = nrep)
          comp=data.frame(resp=medias,groups=comp)}
        print(comp)}
    }
  }
}


