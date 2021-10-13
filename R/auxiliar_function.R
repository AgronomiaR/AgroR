# funcoes auxiliares para as funcoes de analise
# algumas sao do pacote agricolae
# Mendiburu, F., and de Mendiburu, M. F. (2019). Package ‘agricolae’. R Package, Version, 1-2.

mean.stat <-
  function (y, x, stat = "mean")
  {k<-0
    numerico<- NULL
    if(is.null(ncol(x))){
      if(is.numeric(x)){ k<-1
      numerico[1]<-1}}
    else{
      ncolx<-ncol(x)
      for (i in 1:ncolx) {
        if(is.numeric(x[,i])){
          k<-k+1
          numerico[k]<-i
        }}}
    cx <- deparse(substitute(x))
    cy <- deparse(substitute(y))
    x <- data.frame(c1 = 1, x)
    y <- data.frame(v1 = 1, y)
    nx <- ncol(x)
    ny <- ncol(y)
    namex <- names(x)
    namey <- names(y)
    if (nx == 2)
      namex <- c("c1", cx)
    if (ny == 2)
      namey <- c("v1", cy)
    namexy <- c(namex, namey)
    for (i in 1:nx) {
      x[, i] <- as.character(x[, i])}
    z <- NULL
    for (i in 1:nx){z <- paste(z, x[, i], sep = "&")}
    w <- NULL
    for (i in 1:ny) {
      m <- tapply(y[, i], z, stat)
      m <- as.matrix(m)
      w <- cbind(w, m)}
    nw <- nrow(w)
    c <- rownames(w)
    v <- rep("", nw * nx)
    dim(v) <- c(nw, nx)
    for (i in 1:nw) {
      for (j in 1:nx) {
        v[i, j] <- strsplit(c[i], "&")[[1]][j + 1]}}
    rownames(w) <- NULL
    junto <- data.frame(v[, -1], w)
    junto <- junto[, -nx]
    names(junto) <- namexy[c(-1, -(nx + 1))]
    if(k==1 & nx==2) {
      junto[,numerico[1]]<-as.character(junto[,numerico[1]])
      junto[,numerico[1]]<-as.numeric(junto[,numerico[1]])
      junto<-junto[order(junto[,1]),]}
    if (k>0 & nx > 2) {
      for (i in 1:k){
        junto[,numerico[i]]<-as.character(junto[,numerico[i]])
        junto[,numerico[i]]<-as.numeric(junto[,numerico[i]])}
      junto<-junto[do.call("order", c(junto[,1:(nx-1)])),]}
    rownames(junto)<-1:(nrow(junto))
    return(junto)}

levenehomog <- function (y, ...) {
  UseMethod("levenehomog")}

levenehomog.default <- function (y, group, center=median, ...) {
  if (!is.numeric(y))
    stop(deparse(substitute(y)), " is not a numeric variable")
  if (!is.factor(group)){warning(deparse(substitute(group)), " coerced to factor.")
    group <- as.factor(group)}
  valid <- complete.cases(y, group)
  meds <- tapply(y[valid], group[valid], center, ...)
  resp <- abs(y - meds[group])
  table <- anova(lm(resp ~ group))[, c(1, 4, 5)]
  rownames(table)[2] <- " "
  dots <- deparse(substitute(...))
  attr(table, "heading") <- paste("Levene's Test  (center = ",
                                  deparse(substitute(center)),
                                  if(!(dots == "NULL")) paste(":", dots),  ")", sep="")
  table}


levenehomog.formula <- function(y, data, ...) {
  form <- y
  mf <- if (missing(data)) model.frame(form) else model.frame(form, data)
  if (any(sapply(2:dim(mf)[2], function(j) is.numeric(mf[[j]]))))
    stop("Levene's test is not appropriate with quantitative explanatory variables.")
  y <- mf[,1]
  if(dim(mf)[2]==2) group <- mf[,2]
  else {
    if (length(grep("\\+ | \\| | \\^ | \\:",form))>0) stop("Model must be completely crossed formula only.")
    group <- interaction(mf[,2:dim(mf)[2]])}
  levenehomog.default(y=y, group=group, ...)}

levenehomog.lm <- function(y, ...) {
  m <- model.frame(y)
  m$..y <- model.response(m)
  f <- formula(y)
  f[2] <- expression(..y)
  levenehomog.formula(f, data=m, ...)}

ordenacao=function (treatment, means, alpha, pvalue, console){
    n <- length(means)
    z <- data.frame(treatment, means)
    letras<-c(letters[1:26],LETTERS[1:26],1:9,
              c(".","+","-","*","/","#","$","%","&","^","[","]",":",
                "@",";","_","?","!","=","#",rep(" ",2000)))
    w <- z[order(z[, 2], decreasing = TRUE), ]
    M<-rep("",n)
    k<-1
    k1<-0
    j<-1
    i<-1
    cambio<-n
    cambio1<-0
    chequeo=0
    M[1]<-letras[k]
    q <- as.numeric(rownames(w)) #Check
    while(j<n) {
      chequeo<-chequeo+1
      if (chequeo > n) break
      for(i in j:n) {
        s<-pvalue[q[i],q[j]]>alpha
        if(s) {
          if(lastC(M[i]) != letras[k])M[i]<-paste(M[i],letras[k],sep="")
        }
        else {
          k<-k+1
          cambio<-i
          cambio1<-0
          ja<-j
          for(jj in cambio:n) M[jj]<-paste(M[jj],"",sep="") # El espacio
          M[cambio]<-paste(M[cambio],letras[k],sep="")
          for( v in ja:cambio) {
            if(pvalue[q[v],q[cambio]]<=alpha) {j<-j+1
            cambio1<-1
            }
            else break
          }
          break
        }
      }
      if (cambio1 ==0 )j<-j+1
    }
    w<-data.frame(w,stat=M)
    trt <- as.character(w$treatment)
    means <- as.numeric(w$means)
    output <- data.frame(means, groups=M)
    rownames(output)<-trt
    if(k>81)
      cat("\n",k,"groups are estimated.The number of groups exceeded the maximum of 81 labels. change to group=FALSE.\n")
    invisible(output)
  }
lastC <-
  function(x) {
    y<-sub(" +$", "",x)
    p1<-nchar(y)
    cc<-substr(y,p1,p1)
    return(cc)}
duncan <- function(y,
                   trt,
                   DFerror,
                   MSerror,
                   alpha = 0.05,
                   group = TRUE,
                   main = NULL,
                   console = FALSE)
{name.y <- paste(deparse(substitute(y)))
    name.t <- paste(deparse(substitute(trt)))
    if(is.null(main))main<-paste(name.y,"~", name.t)
    clase<-c("aov","lm")
    if("aov"%in%class(y) | "lm"%in%class(y)){
      if(is.null(main))main<-y$call
      A<-y$model
      DFerror<-df.residual(y)
      MSerror<-deviance(y)/DFerror
      y<-A[,1]
      ipch<-pmatch(trt,names(A))
      nipch<- length(ipch)
      for(i in 1:nipch){
        if (is.na(ipch[i]))
          return(if(console)cat("Name: ", trt, "\n", names(A)[-1], "\n"))}
      name.t<- names(A)[ipch][1]
      trt <- A[, ipch]
      if (nipch > 1){
        trt <- A[, ipch[1]]
        for(i in 2:nipch){
          name.t <- paste(name.t,names(A)[ipch][i],sep=":")
          trt <- paste(trt,A[,ipch[i]],sep=":")
        }}
      name.y <- names(A)[1]
    }
    junto <- subset(data.frame(y, trt), is.na(y) == FALSE)
    Mean<-mean(junto[,1])
    CV<-sqrt(MSerror)*100/Mean
    medians<-mean.stat(junto[,1],junto[,2],stat="median")
    for(i in c(1,5,2:4)) {
      x <- mean.stat(junto[,1],junto[,2],function(x)quantile(x)[i])
      medians<-cbind(medians,x[,2])
    }
    medians<-medians[,3:7]
    names(medians)<-c("Min","Max","Q25","Q50","Q75")
    means <- mean.stat(junto[,1],junto[,2],stat="mean") # change
    sds <-   mean.stat(junto[,1],junto[,2],stat="sd") #change
    nn <-   mean.stat(junto[,1],junto[,2],stat="length") # change
    means<-data.frame(means,std=sds[,2],r=nn[,2],medians)
    names(means)[1:2]<-c(name.t,name.y)
    ntr<-nrow(means)
    Tprob<-NULL
    k<-0
    for(i in 2:ntr){
      k<-k+1
      x <- suppressWarnings(warning(qtukey((1-alpha)^(i-1), i, DFerror)))
      if(x=="NaN")break
      else Tprob[k]<-x
    }
    if(k<(ntr-1)){
      for(i in k:(ntr-1)){
        f <- Vectorize(function(x)ptukey(x,i+1,DFerror)-(1-alpha)^i)
        Tprob[i]<-uniroot(f, c(0,100))$root
      }
    }
    Tprob<-as.numeric(Tprob)
    nr <- unique(nn[,2])
    if(console){
      cat("\nStudy:", main)
      cat("\n\nDuncan's new multiple range test\nfor",name.y,"\n")
      cat("\nMean Square Error: ",MSerror,"\n\n")
      cat(paste(name.t,",",sep="")," means\n\n")
      print(data.frame(row.names = means[,1], means[,2:6]))
    }
    if(length(nr) == 1 ) sdtdif <- sqrt(MSerror/nr)
    else {
      nr1 <-  1/mean(1/nn[,2])
      sdtdif <- sqrt(MSerror/nr1)
    }
    DUNCAN <- Tprob * sdtdif
    names(DUNCAN)<-2:ntr
    duncan<-data.frame(Table=Tprob,CriticalRange=DUNCAN)
    if ( group & length(nr) == 1 & console){
      cat("\nAlpha:",alpha,"; DF Error:",DFerror,"\n")
      cat("\nCritical Range\n")
      print(DUNCAN)}
    if ( group & length(nr) != 1 & console) cat("\nGroups according to probability of means differences and alpha level(",alpha,")\n")
    if ( length(nr) != 1) duncan<-NULL
    Omeans<-order(means[,2],decreasing = TRUE) #correccion 2019, 1 abril.
    Ordindex<-order(Omeans)
    comb <-utils::combn(ntr,2)
    nn<-ncol(comb)
    dif<-rep(0,nn)
    DIF<-dif
    LCL<-dif
    UCL<-dif
    pvalue<-dif
    odif<-dif
    sig<-NULL
    for (k in 1:nn) {
      i<-comb[1,k]
      j<-comb[2,k]
      dif[k]<-means[i,2]-means[j,2]
      DIF[k]<-abs(dif[k])
      nx<-abs(i-j)+1
      odif[k] <- abs(Ordindex[i]- Ordindex[j])+1
      pvalue[k]<- round(1-ptukey(DIF[k]/sdtdif,odif[k],DFerror)^(1/(odif[k]-1)),4)
      LCL[k] <- dif[k] - DUNCAN[odif[k]-1]
      UCL[k] <- dif[k] + DUNCAN[odif[k]-1]
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
      if(console){cat("\nComparison between treatments means\n\n")
        print(comparison)}
      groups=NULL
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
      groups <- ordenacao(means[, 1], means[, 2],alpha, Q,console)
      names(groups)[1]<-name.y
      if(console) {
        cat("\nMeans with the same letter are not significantly different.\n\n")
        print(groups)
      }
    }
    parameters<-data.frame(test="Duncan",name.t=name.t,ntr = ntr,alpha=alpha)
    statistics<-data.frame(MSerror=MSerror,Df=DFerror,Mean=Mean,CV=CV)
    rownames(parameters)<-" "
    rownames(statistics)<-" "
    rownames(means)<-means[,1]
    means<-means[,-1]
    output<-list(statistics=statistics,parameters=parameters, duncan=duncan,
                 means=means,comparison=comparison,groups=groups)
    class(output)<-"group"
    invisible(output)
  }

TUKEY <- function(y, trt, DFerror,
                      MSerror, alpha=0.05, group=TRUE,
                      main = NULL,unbalanced=FALSE,console=FALSE){
    name.y <- paste(deparse(substitute(y)))
    name.t <- paste(deparse(substitute(trt)))
    if(is.null(main))main<-paste(name.y,"~", name.t)
    clase<-c("aov","lm")
    if("aov"%in%class(y) | "lm"%in%class(y)){
      if(is.null(main))main<-y$call
      A<-y$model
      DFerror<-df.residual(y)
      MSerror<-deviance(y)/DFerror
      y<-A[,1]
      ipch<-pmatch(trt,names(A))
      nipch<- length(ipch)
      for(i in 1:nipch){
        if (is.na(ipch[i]))
          return(if(console)cat("Name: ", trt, "\n", names(A)[-1], "\n"))
      }
      name.t<- names(A)[ipch][1]
      trt <- A[, ipch]
      if (nipch > 1){
        trt <- A[, ipch[1]]
        for(i in 2:nipch){
          name.t <- paste(name.t,names(A)[ipch][i],sep=":")
          trt <- paste(trt,A[,ipch[i]],sep=":")
        }}
      name.y <- names(A)[1]
    }
    junto <- subset(data.frame(y, trt), is.na(y) == FALSE)
    Mean<-mean(junto[,1])
    CV<-sqrt(MSerror)*100/Mean
    medians<-mean.stat(junto[,1],junto[,2],stat="median")
    for(i in c(1,5,2:4)) {
      x <- mean.stat(junto[,1],junto[,2],function(x)quantile(x)[i])
      medians<-cbind(medians,x[,2])
    }
    medians<-medians[,3:7]
    names(medians)<-c("Min","Max","Q25","Q50","Q75")
    means <- mean.stat(junto[,1],junto[,2],stat="mean")
    sds <-   mean.stat(junto[,1],junto[,2],stat="sd")
    nn <-   mean.stat(junto[,1],junto[,2],stat="length")
    means<-data.frame(means,std=sds[,2],r=nn[,2],medians)
    names(means)[1:2]<-c(name.t,name.y)
    ntr<-nrow(means)
    Tprob <- qtukey(1-alpha,ntr, DFerror)
    nr<-unique(nn[, 2])
    nr1<-1/mean(1/nn[,2])
    if(console){
      cat("\nStudy:", main)
      cat("\n\nHSD Test for",name.y,"\n")
      cat("\nMean Square Error: ",MSerror,"\n\n")
      cat(paste(name.t,",",sep="")," means\n\n")
      print(data.frame(row.names = means[,1], means[,2:6]))
      cat("\nAlpha:",alpha,"; DF Error:",DFerror,"\n")
      cat("Critical Value of Studentized Range:", Tprob,"\n")
    }
    HSD <- Tprob * sqrt(MSerror/nr)
    statistics<-data.frame(MSerror=MSerror,Df=DFerror,Mean=Mean,CV=CV,MSD=HSD)
    if ( group & length(nr) == 1 & console) cat("\nMinimun Significant Difference:",HSD,"\n")
    if ( group & length(nr) != 1 & console) cat("\nGroups according to probability of means differences and alpha level(",alpha,")\n")
    if ( length(nr) != 1) statistics<-data.frame(MSerror=MSerror,Df=DFerror,Mean=Mean,CV=CV)
    comb <-utils::combn(ntr,2)
    nn<-ncol(comb)
    dif<-rep(0,nn)
    sig<-NULL
    LCL<-dif
    UCL<-dif
    pvalue<-rep(0,nn)
    for (k in 1:nn) {
      i<-comb[1,k]
      j<-comb[2,k]
      dif[k]<-means[i,2]-means[j,2]
      sdtdif<-sqrt(MSerror * 0.5*(1/means[i,4] + 1/means[j,4]))
      if(unbalanced)sdtdif<-sqrt(MSerror /nr1)
      pvalue[k]<- round(1-ptukey(abs(dif[k])/sdtdif,ntr,DFerror),4)
      LCL[k] <- dif[k] - Tprob*sdtdif
      UCL[k] <- dif[k] + Tprob*sdtdif
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
      if(console){cat("\nComparison between treatments means\n\n")
        print(comparison)}
      groups=NULL
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
      groups <- ordenacao(means[, 1], means[, 2],alpha, Q,console)
      names(groups)[1]<-name.y
      if(console) {
        cat("\nTreatments with the same letter are not significantly different.\n\n")
        print(groups)
      }
    }
    parameters<-data.frame(test="Tukey",name.t=name.t,ntr = ntr, StudentizedRange=Tprob,alpha=alpha)
    rownames(parameters)<-" "
    rownames(statistics)<-" "
    rownames(means)<-means[,1]
    means<-means[,-1]
    output<-list(statistics=statistics,parameters=parameters,
                 means=means,comparison=comparison,groups=groups)
    class(output)<-"group"
    invisible(output)
  }

LSD = function(y,
               trt,
               DFerror,
               MSerror,
               alpha = 0.05,
               p.adj = c("none",
                         "holm",
                         "hommel",
                         "hochberg",
                         "bonferroni",
                         "BH",
                         "BY",
                         "fdr"),
               group = TRUE,
               main = NULL,
               console = FALSE) {

  p.adj <- match.arg(p.adj)
  clase <- c("aov", "lm")
  name.y <- paste(deparse(substitute(y)))
  name.t <- paste(deparse(substitute(trt)))
  if(is.null(main))main<-paste(name.y,"~", name.t)
  if ("aov" %in% class(y) | "lm" %in% class(y)) {
    if(is.null(main))main<-y$call
    A <- y$model
    DFerror <- df.residual(y)
    MSerror <- deviance(y)/DFerror
    y <- A[, 1]
    ipch <- pmatch(trt, names(A))
    nipch<- length(ipch)
    for(i in 1:nipch){
      if (is.na(ipch[i]))
        return(if(console)cat("Name: ", trt, "\n", names(A)[-1], "\n"))}
    name.t<- names(A)[ipch][1]
    trt <- A[, ipch]
    if (nipch > 1){
      trt <- A[, ipch[1]]
      for(i in 2:nipch){
        name.t <- paste(name.t,names(A)[ipch][i],sep=":")
        trt <- paste(trt,A[,ipch[i]],sep=":")
      }}
    name.y <- names(A)[1]}
  junto <- subset(data.frame(y, trt), is.na(y) == FALSE)
  Mean<-mean(junto[,1])
  CV<-sqrt(MSerror)*100/Mean
  medians<-mean.stat(junto[,1],junto[,2],stat="median")
  for(i in c(1,5,2:4)) {
    x <- mean.stat(junto[,1],junto[,2],function(x)quantile(x)[i])
    medians<-cbind(medians,x[,2])}
  medians<-medians[,3:7]
  names(medians)<-c("Min","Max","Q25","Q50","Q75")
  means <- mean.stat(junto[, 1], junto[, 2], stat = "mean")
  sds <- mean.stat(junto[, 1], junto[, 2], stat = "sd")
  nn <- mean.stat(junto[, 1], junto[, 2], stat = "length")
  std.err <- sqrt(MSerror)/sqrt(nn[, 2])
  Tprob <- qt(1 - alpha/2, DFerror)
  LCL <- means[, 2] - Tprob * std.err
  UCL <- means[, 2] + Tprob * std.err
  means <- data.frame(means, std=sds[,2], r = nn[, 2], LCL, UCL,medians)
  names(means)[1:2] <- c(name.t, name.y)
  ntr <- nrow(means)
  nk <- choose(ntr, 2)
  if (p.adj != "none") {
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
  nr <- unique(nn[, 2])
  if(console){
    cat("\nStudy:", main)
    if(console)cat("\n\nLSD t Test for", name.y, "\n")
    if (p.adj != "none")cat("P value adjustment method:", p.adj, "\n")
    cat("\nMean Square Error: ", MSerror, "\n\n")
    cat(paste(name.t, ",", sep = ""), " means and individual (",
        (1 - alpha) * 100, "%) CI\n\n")
    print(data.frame(row.names = means[, 1], means[, 2:8]))
    cat("\nAlpha:", alpha, "; DF Error:", DFerror)
    cat("\nCritical Value of t:", Tprob, "\n")}
  statistics<-data.frame(MSerror=MSerror,Df=DFerror,Mean=Mean,CV=CV)
  if (length(nr) == 1)  LSD <- Tprob * sqrt(2 * MSerror/nr)
  if ( group & length(nr) == 1 & console) {
    if(p.adj=="none") cat("\nleast Significant Difference:",LSD,"\n")
    else cat("\nMinimum Significant Difference:",LSD,"\n")}
  if ( group & length(nr) != 1 & console)
    cat("\nGroups according to probability of means differences and alpha level(",alpha,")\n")

  if ( length(nr) == 1 & p.adj=="none") statistics<-data.frame(statistics, t.value=Tprob,LSD=LSD)
  if ( length(nr) == 1 & p.adj!="none") statistics<-data.frame(statistics, t.value=Tprob,MSD=LSD)
  LSD=" "
  comb <- utils::combn(ntr, 2)
  nn <- ncol(comb)
  dif <- rep(0, nn)
  pvalue <- dif
  sdtdif <- dif
  sig <- rep(" ", nn)
  for (k in 1:nn) {
    i <- comb[1, k]
    j <- comb[2, k]
    dif[k] <-means[i, 2] - means[j, 2]
    sdtdif[k] <- sqrt(MSerror * (1/means[i, 4] + 1/means[j,4]))
    pvalue[k] <- 2 * (1 - pt(abs(dif[k])/sdtdif[k], DFerror))}
  if (p.adj != "none")
    pvalue <- p.adjust(pvalue, p.adj)
  pvalue <- round(pvalue,4)
  for (k in 1:nn) {
    if (pvalue[k] <= 0.001)
      sig[k] <- "***"
    else if (pvalue[k] <= 0.01)
      sig[k] <- "**"
    else if (pvalue[k] <= 0.05)
      sig[k] <- "*"
    else if (pvalue[k] <= 0.1)
      sig[k] <- "."}
  tr.i <- means[comb[1, ], 1]
  tr.j <- means[comb[2, ], 1]
  LCL <- dif - Tprob * sdtdif
  UCL <- dif + Tprob * sdtdif
  comparison <- data.frame(difference = dif, pvalue = pvalue, "signif."=sig, LCL, UCL)
  if (p.adj !="bonferroni" & p.adj !="none"){
    comparison<-comparison[,1:3]
  }
  rownames(comparison) <- paste(tr.i, tr.j, sep = " - ")
  if (!group) {
    if(console){
      cat("\nComparison between treatments means\n\n")
      print(comparison)}
    groups <- NULL}
  if (group){
    comparison=NULL
    Q<-matrix(1,ncol=ntr,nrow=ntr)
    p<-pvalue
    k<-0
    for(i in 1:(ntr-1)){
      for(j in (i+1):ntr){
        k<-k+1
        Q[i,j]<-p[k]
        Q[j,i]<-p[k]}}
    groups <- ordenacao(means[, 1], means[, 2],alpha, Q,console)
    names(groups)[1]<-name.y
    if(console) {
      cat("\nTreatments with the same letter are not significantly different.\n\n")
      print(groups)}
  }
  parameters<-data.frame(test="Fisher-LSD",p.ajusted=p.adj,name.t=name.t,ntr = ntr,alpha=alpha)
  rownames(parameters)<-" "
  rownames(statistics)<-" "
  rownames(means)<-means[,1]
  means<-means[,-1]
  output<-list(statistics=statistics,parameters=parameters,
               means=means,comparison=comparison,groups=groups)
  class(output)<-"group"
  invisible(output)
}

sk<-function(y,
             trt,
             DFerror,
             SSerror,
             alpha = 0.05,
             group = TRUE,
             main = NULL){
  sk <- function(medias,s2,dfr,prob){
    bo <- 0
    si2 <- s2
    defr <- dfr
    parou <- 1
    np <- length(medias) - 1
    for (i in 1:np){
      g1 <- medias[1:i]
      g2 <- medias[(i+1):length(medias)]
      B0 <- sum(g1)^2/length(g1) + sum(g2)^2/length(g2) - (sum(g1) + sum(g2))^2/length(c(g1,g2))
      if (B0 > bo)
      {bo <- B0
      parou <- i}
    }

    g1 <- medias[1:parou]
    g2 <- medias[(parou+1):length(medias)]
    teste <- c(g1,g2)
    sigm2 <- (sum(teste^2) - sum(teste)^2/length(teste) + defr*si2)/(length(teste) + defr)
    lamb <- pi*bo/(2*sigm2*(pi-2))
    v0 <- length(teste)/(pi-2)
    p <- pchisq(lamb,v0,lower.tail = FALSE)
    if (p < prob) {
      for (i in 1:length(g1)){
        cat(names(g1[i]),"\n",file="sk_groups",append=TRUE)}
      cat("*","\n",file="sk_groups",append=TRUE)}
    if (length(g1)>1){sk(g1,s2,dfr,prob)}
    if (length(g2)>1){sk(g2,s2,dfr,prob)}
  }
  trt=factor(trt,unique(trt))
  trt1=trt
  levels(trt)=paste("T",1:length(levels(trt)),sep = "")
  medias <- sort(tapply(y,trt,mean),decreasing=TRUE)
  dfr <- DFerror

  rep <- tapply(y,trt,length)
  s0 <- MSerror <-SSerror/DFerror
  s2 <- s0/rep[1]
  prob <- alpha
  sk(medias,s2,dfr,prob)
  f <- names(medias)
  names(medias) <- 1:length(medias)
  resultado <- data.frame("r"=0,"f"=f,"m"=medias)
  if (file.exists("sk_groups") == FALSE) {stop} else{
    xx <- read.table("sk_groups")
    file.remove("sk_groups")
    x <- xx[[1]]
    x <- as.vector(x)
    z <- 1

    for (j in 1:length(x)){
      if (x[j] == "*")	{z <- z+1}
      for (i in 1:length(resultado$f)){
        if (resultado$f[i]==x[j]){
          resultado$r[i] <- z;}
      }
    }

  }
  letras<-letters
  if(length(resultado$r)>26) {
    l<-floor(length(resultado$r)/26)
    for(i in 1:l) letras<-c(letras,paste(letters,i,sep=''))
  }
  res <- 1
  for (i in 1:(length(resultado$r)-1))
  {
    if (resultado$r[i] != resultado$r[i+1]){
      resultado$r[i] <- letras[res]
      res <- res+1
      if (i == (length(resultado$r)-1)){
        resultado$r[i+1] <- letras[res]}
    }
    else{
      resultado$r[i] <- letras[res]
      if (i == (length(resultado$r)-1)){
        resultado$r[i+1] <- letras[res]
      }
    }
  }
  names(resultado) <- c("groups","Tratamentos","Means")
  resultado1=resultado[,c(3,1)]
  rownames(resultado1)=resultado$Tratamentos
  final=list(resultado1)[[1]]
  final=final[as.character(unique(trt)),]
  rownames(final)=as.character(unique(trt1))
  final
}
