#' Utils: Experimental sketch
#'
#' @description Experimental sketching function
#' @author Gabriel Danilo Shimizu, \email{shimizu@uel.br}
#' @author Leandro Simoes Azeredo Goncalves
#' @author Rodrigo Yudi Palhaci Marubayashi
#' @name sketch
#' @param trat Vector with factor A levels
#' @param trat1 Vector with levels of factor B (Set to NULL if not factorial or psub)
#' @param trat2 Vector with levels of factor C (Set to NULL if not factorial)
#' @param r Number of repetitions
#' @param design Experimental design ("dic", "dbc", "dql","psubdic","psubdbc","fat2dic","fat2dbc")
#' @param pos Repeat position (line or column)
#' @keywords croqui
#' @keywords experimental
#' @return Returns an experimental sketch according to the specified design.
#' @note The sketches have only a rectangular shape, and the blocks (in the case of randomized blocks) can be in line or in a column.
#' @references
#' Mendiburu, F., & de Mendiburu, M. F. (2019). Package ‘agricolae’. R Package, Version, 1-2.
#' @export
#' @examples
#' Trat=paste("Tr",1:6)
#'
#' #=============================
#' # Completely randomized design
#' #=============================
#' sketch(Trat,r=3)
#' sketch(Trat,r=3,pos="column")
#'
#' #=============================
#' # Randomized block design
#' #=============================
#' sketch(Trat, r=3, design="dbc")
#' sketch(Trat, r=3, design="dbc",pos="column")
#'
#' #=============================
#' # Completely randomized experiments in double factorial
#' #=============================
#' sketch(trat=c("A","B"),
#'        trat1=c("A","B","C"),
#'        design = "fat2dic",
#'        r=3)
#'
#' sketch(trat=c("A","B"),
#'        trat1=c("A","B","C"),
#'        design = "fat2dic",
#'        r=3,
#'        pos="column")

sketch=function(trat,
                trat1=NULL,
                trat2=NULL,
                r,
                design="dic",
                pos="line"){
  requireNamespace("ggplot2")
  design.crd <-
    function(trt,r,serie=2,seed=0,kinds="Super-Duper",randomization=TRUE)
    {
      number<-0
      if(serie>0) number<-10^serie
      junto<-data.frame(trt,r)
      junto<-junto[order(junto[,1]),]
      TR<-as.character(junto[,1])
      r<-as.numeric(junto[,2])
      y <- rep(TR[1], r[1])
      tr <- length(TR)
      if (seed == 0) {
        genera<-runif(1)
        seed <-.Random.seed[3]
      }
      set.seed(seed,kinds)
      parameters<-list(design="crd",trt=trt,r=r,serie=serie,seed=seed,kinds=kinds,randomization)
      for (i in 2:tr) y <- c(y, rep(TR[i], r[i]))
      trat<-y
      if(randomization)trat <- sample(y, length(y), replace = FALSE)
      plots <- number+1:length(trat)
      dca<-data.frame(plots, trat)
      dca[,1]<-as.numeric(dca[,1])
      xx<-dca[order(dca[,2],dca[,1]),]
      r1<-seq(1,r[1])
      for (i in 2:length(r)) {
        r1<-c(r1,seq(1,r[i]))
      }
      yy<-data.frame(plots=xx[,1],r=r1,xx[,2])
      book<-yy[order(yy[,1]),]
      rownames(book)<-rownames(yy)
      names(book)[3]<-c(paste(deparse(substitute(trt))))
      outdesign<-list(parameters=parameters,book=book)
      return(outdesign)
    }

  design.rcbd <-
    function (trt, r,serie=2,seed=0,kinds="Super-Duper",first=TRUE,continue=FALSE,randomization=TRUE )
    {
      number<-10
      if(serie>0) number<-10^serie
      ntr <- length(trt)
      if (seed == 0) {
        genera<-runif(1)
        seed <-.Random.seed[3]
      }
      set.seed(seed,kinds)
      parameters<-list(design="rcbd",trt=trt,r=r,serie=serie,seed=seed,kinds=kinds,randomization)
      mtr <-trt
      if(randomization)mtr <- sample(trt, ntr, replace = FALSE)
      block <- c(rep(1, ntr))
      for (y in 2:r) {
        block <- c(block, rep(y, ntr))
        if(randomization)mtr <- c(mtr, sample(trt, ntr, replace = FALSE))
      }
      if(randomization){
        if(!first) mtr[1:ntr]<-trt
      }
      plots <- block*number+(1:ntr)
      book <- data.frame(plots, block = as.factor(block), trt = as.factor(mtr))
      names(book)[3] <- c(paste(deparse(substitute(trt))))
      names(book)[3]<-c(paste(deparse(substitute(trt))))
      if(continue){
        start0<-10^serie
        if(serie==0) start0<-0
        book$plots<-start0+1:nrow(book)
      }
      outdesign<-list(parameters=parameters,sketch=matrix(book[,3], byrow = TRUE, ncol = ntr),book=book)
      return(outdesign)
    }

  design.lsd <-
    function (trt,serie=2,seed=0,kinds="Super-Duper",first=TRUE,randomization=TRUE)
    {
      number<-10
      if(serie>0) number<-10^serie
      r <- length(trt)
      if (seed == 0) {
        genera<-runif(1)
        seed <-.Random.seed[3]
      }
      set.seed(seed,kinds)
      parameters<-list(design="lsd",trt=trt,r=r,serie=serie,seed=seed,kinds=kinds,randomization)
      a <- 1:(r * r)
      dim(a) <- c(r, r)
      for (i in 1:r) {
        for (j in 1:r) {
          k <- i + j - 1
          if (k > r)
            k <- i + j - r - 1
          a[i, j] <- k
        }
      }
      m<-2:r
      if(randomization)m<-sample(2:r,r-1)
      a<-a[,c(1,m)]
      if(randomization){
        if (first) {
          m<-sample(1:r,r)
          a<-a[m,]
        }}
      trat<-trt[a]
      columna <- rep(gl(r, 1), r)
      fila <- gl(r, r)
      fila <- as.character(fila)
      fila <- as.numeric(fila)
      plots <- fila*number+(1:r)
      book <- data.frame(plots, row = as.factor(fila), col = as.factor(columna),
                         trat = as.factor(trat))
      names(book)[4] <- c(paste(deparse(substitute(trt))))
      outdesign<-list(parameters=parameters,sketch=matrix(book[,4], byrow = TRUE, ncol = r),book=book)
      return(outdesign)
    }

  design.split <-
    function (trt1, trt2,r=NULL, design=c("rcbd","crd","lsd"),serie = 2, seed = 0, kinds = "Super-Duper",
              first=TRUE,randomization=TRUE )
    {
      n1<-length(trt1)
      n2<-length(trt2)
      if (seed == 0) {
        genera<-runif(1)
        seed <-.Random.seed[3]
      }
      set.seed(seed,kinds)
      design <- match.arg(design)
      number<-10^serie +1
      if (design == "crd") {
        plan<-design.crd(trt1,r,serie, seed, kinds,randomization)
        k<-3
      }
      if (design == "rcbd"){
        plan<-design.rcbd(trt1,r,serie, seed, kinds, first,randomization)
        k<-3
      }
      if (design == "lsd") {
        plan<-design.lsd(trt1,serie, seed, kinds, first,randomization)
        r<-n1
        k<-4
      }
      book<-plan$book
      parameters<-plan$parameters
      names(parameters)[2]<-"trt1"
      parameters$applied<-parameters$design
      parameters$design<-"split"
      parameters$trt2<-trt2
      j<-0
      B<-list()
      for(i in c(1,7,2,8,3:6)){
        j<-j+1
        B[[j]]<-parameters[[i]]
        names(B)[j]<-names(parameters)[i]
      }
      nplot<-nrow(book)
      d<-NULL
      if(randomization){
        for(i in 1:nplot)d<-rbind(d,sample(trt2,n2))
      }
      else{
        d<-rbind(d,trt2[1:n2])
      }
      aa<-data.frame(book,trt2=d[,1])
      for(j in 2:n2) aa<-rbind(aa,data.frame(book,trt2=d[,j]))
      aa<-aa[order(aa[,1]),]
      splots<-rep(gl(n2,1),nplot)
      book <- data.frame(plots=aa[,1],splots,aa[,-1])
      rownames(book)<-1:(nrow(book))
      names(book)[k+1] <- c(paste(deparse(substitute(trt1))))
      names(book)[k+2] <- c(paste(deparse(substitute(trt2))))
      outdesign<-list(parameters=B,book=book)
      return(outdesign)
    }

  design.ab <-
    function(trt, r=NULL,serie=2,design=c("rcbd","crd","lsd"),seed=0,kinds="Super-Duper",
             first=TRUE,randomization=TRUE ){
      design <- match.arg(design)
      if( design=="rcbd" | design=="crd") posicion <- 3
      else posicion <- 4
      serie<-serie; seed<-seed; kinds<-kinds; first<-first;
      ntr<-length(trt)
      fact<-NULL
      tr0<-1:trt[1]
      k<-0
      a<-trt[1];b<-trt[2]
      for(i in 1:a){
        for(j in 1:b){
          k<-k+1
          fact[k]<-paste(tr0[i],j)
        }
      }

      if(ntr >2) {
        for(m in 3:ntr){
          k<-0
          tr0<-fact
          fact<-NULL
          a<-a*b
          b<-trt[m]
          for(i in 1:a){
            for(j in 1:b){
              k<-k+1
              fact[k]<-paste(tr0[i],j)
            }
          }
        }
      }
      if(design=="rcbd")plan<-design.rcbd(trt=fact, r, serie, seed, kinds, first,randomization )
      if(design=="crd")plan<-design.crd(trt=fact, r, serie, seed, kinds,randomization)
      if(design=="lsd")plan<-design.lsd(trt=fact, serie, seed, kinds, first,randomization )
      parameters<-plan$parameters
      parameters$applied<-parameters$design
      parameters$design<-"factorial"
      plan<-plan$book
      trt<-as.character(plan[,posicion])
      nplan<-nrow(plan)
      A<-rep(" ",nplan*ntr)
      dim(A)<-c(nplan,ntr)
      colnames(A)<-LETTERS[1:ntr]

      for(i in 1:nplan) {
        A[i,]<-unlist(strsplit(trt[i], " "))
      }
      A<-as.data.frame(A)
      book<-data.frame(plan[,1:(posicion-1)],A)
      outdesign<-list(parameters=parameters,book=book)
      return(outdesign)
    }

  #=================
  if(design=="dic"){sort=design.crd(trat,r,serie=0)
  data=sort$book
  data$x=rep(1:length(unique(data$trat)),r)
  data$x=factor(data$x,unique(data$x))
  data$y=rep(1:r,e=length(unique(data$trat)))
  data$y=factor(data$y,unique(data$y))
  x=data$x
  y=data$y
  if(pos=="line"){graph=ggplot(data,aes(x=x,y=y,fill=trat))+
    geom_tile(color="black")+labs(x="",y="",fill="Treatments")+
    geom_text(aes(label=trat))+theme_bw()}
  if(pos=="column"){graph=ggplot(data,aes(y=x,x=y,fill=trat))+
    geom_tile(color="black")+labs(x="",y="",fill="Treatments")+
    geom_text(aes(label=trat))+theme_bw()}}

  #=================
  if(design=="dbc"){sort=design.rcbd(trat,r,serie=0)
  data=sort$book
  data$x=rep(1:length(unique(data$trat)),r)
  data$x=factor(data$x,unique(data$x))
  x=data$x
  block=data$block
  if(pos=="line"){graph=ggplot(data,aes(x=x,y=block,fill=trat))+
    geom_tile(color="black")+labs(y="Block",x="",fill="Treatments")+
    geom_text(aes(label=trat))+theme_bw()}
  if(pos=="column"){graph=ggplot(data,aes(y=x,x=block,fill=trat))+
    geom_tile(color="black")+labs(y="",x="Block",fill="Treatments")+
    geom_text(aes(label=trat))+theme_bw()}}

  #=================
  if(design=="dql"){sort=design.lsd(trat,r,serie=0)
  data=sort$book
  graph=ggplot(data,aes(x=row,y=col,fill=trat))+
    geom_tile(color="black")+labs(x="Row",y="Column",fill="Treatments")+
    geom_text(aes(label=trat))+theme_bw()}

  #=================
  if(design=="psubdic"){sort=design.split(trat,trat1,r,design = "crd",serie=0)
  data=sort$book
  data$x=rep(1:length(unique(paste(data$trat,data$trat1))),r)
  data$x=factor(data$x,unique(data$x))
  data$y=rep(1:r,e=length(unique(paste(data$trat,data$trat1))))
  data$y=factor(data$y,unique(data$y))
  x=data$x
  y=data$y

  if(pos=="column"){graph=ggplot(data,aes(x=y,y=x,fill=paste(trat,trat1)))+
    geom_tile(color="black")+labs(x="",y="",fill="Treatments")+
    geom_text(aes(label=paste(trat,trat1)))+theme_bw()}
  if(pos=="line"){graph=ggplot(data,aes(x=x,y=y,fill=paste(trat,trat1)))+
    geom_tile(color="black")+labs(x="",y="",fill="Treatments")+
    geom_text(aes(label=paste(trat,trat1)))+theme_bw()}}

  #=================
  if(design=="psubdbc"){sort=design.split(trat,trat1,r,design = "rcbd",serie=0)
  data=sort$book
  data$x=rep(1:length(unique(paste(data$trat,data$trat1))),r)
  data$x=factor(data$x,unique(data$x))
  x=data$x
  block=data$block
  if(pos=="column"){graph=ggplot(data,aes(y=x,x=block,fill=paste(trat,trat1)))+
    geom_tile(color="black")+labs(x="Block",y="",fill="Treatments")+
    geom_text(aes(label=paste(trat,trat1)))+theme_bw()}
  if(pos=="line"){graph=ggplot(data,aes(y=block,x=x,fill=paste(trat,trat1)))+
    geom_tile(color="black")+labs(x="",y="Block",fill="Treatments")+
    geom_text(aes(label=paste(trat,trat1)))+theme_bw()}}

  #=================
  if(design=="fat2dic"){sort=design.ab(c(length(trat),length(trat1)),r,design = "crd",serie=0)
  sort$book$A=as.factor(sort$book$A)
  sort$book$B=as.factor(sort$book$B)
  levels(sort$book$A)=trat
  levels(sort$book$B)=trat1
  sort$book$trat=paste(sort$book$A,sort$book$B)
  data=sort$book
  data$x=rep(1:length(unique(paste(data$trat,data$trat1))),r)
  data$x=factor(data$x,unique(data$x))
  data$y=rep(1:r,e=length(unique(paste(data$trat,data$trat1))))
  data$y=factor(data$y,unique(data$y))
  A=data$A
  B=data$B
  x=data$x
  y=data$y
  if(pos=="column"){graph=ggplot(data,aes(x=y,y=x,fill=paste(A,B)))+
    geom_tile(color="black")+labs(x="",y="",fill="Treatments")+
    geom_text(aes(label=paste(A,B)))+theme_bw()}
  if(pos=="line"){graph=ggplot(data,aes(x=x,y=y,fill=paste(A,B)))+
    geom_tile(color="black")+labs(x="",y="",fill="Treatments")+
    geom_text(aes(label=paste(A,B)))+theme_bw()}}

  #=================
  if(design=="fat2dbc"){sort=design.ab(c(length(trat),length(trat1)),r,design = "rcbd",serie=0)
  sort$book$A=as.factor(sort$book$A)
  sort$book$B=as.factor(sort$book$B)
  levels(sort$book$A)=trat
  levels(sort$book$B)=trat1
  sort$book$trat=paste(sort$book$A,sort$book$B)
  data=sort$book
  data$x=rep(1:length(unique(paste(data$trat,data$trat1))),r)
  data$x=factor(data$x,unique(data$x))
  A=data$A
  B=data$B
  x=data$x
  block=data$block
  if(pos=="column"){graph=ggplot(data,aes(y=x,x=block,fill=paste(A,B)))+
    geom_tile(color="black")+labs(x="Block",y="",fill="Treatments")+
    geom_text(aes(label=paste(A,B)))+theme_bw()}
  if(pos=="line"){graph=ggplot(data,aes(y=block,x=x,fill=paste(A,B)))+
    geom_tile(color="black")+labs(x="",y="Block",fill="Treatments")+
    geom_text(aes(label=paste(A,B)))+theme_bw()}}

  if(design=="fat3dic"){
    trat=expand.grid(trat,trat1,trat2)
    tr=paste(trat$Var1,trat$Var2,trat$Var3)
    trats=rep(tr,r)
    sorteio=sample(trats)
    x=rep(1:(length(sorteio)/r),r)
    y=rep(1:r,e=(length(sorteio)/r))
    data=data.frame(x,y,sorteio)
    data$x=factor(data$x,unique(data$x))
    data$y=factor(data$y,unique(data$y))

    if(pos=="line"){graph=ggplot(data,aes(x=y,y=x,fill=sorteio))+
      geom_tile(color="black")+labs(x="",y="",fill="Treatments")+
      geom_text(aes(label=sorteio))+theme_bw()}
    if(pos=="column"){graph=ggplot(data,aes(x=x,y=y,fill=sorteio))+
      geom_tile(color="black")+labs(x="",y="",fill="Treatments")+
      geom_text(aes(label=sorteio))+theme_bw()}
    print(graph)}

  if(design=="fat3dbc"){
    trat=expand.grid(trat,trat1,trat2)
    tr=paste(trat$Var1,trat$Var2,trat$Var3)

    sorteio=matrix(NA,ncol=length(tr),nrow=r)
    for(i in 1:r){
      sorteio[i,]=sample(tr)
    }
    sorteio=as.vector(sorteio)
    x=rep(1:(length(sorteio)/r),e=r)
    y=rep(1:r,(length(sorteio)/r))
    data=data.frame(x,y,sorteio)
    data$x=factor(data$x,unique(data$x))
    data$y=factor(data$y,unique(data$y))
    if(pos=="line"){graph=ggplot(data,aes(x=y,y=x,fill=sorteio))+
      geom_tile(color="black")+labs(x="block",y="",fill="Treatments")+
      geom_text(aes(label=sorteio))+theme_bw()}
    if(pos=="column"){graph=ggplot(data,aes(x=x,y=y,fill=sorteio))+
      geom_tile(color="black")+labs(x="",y="block",fill="Treatments")+
      geom_text(aes(label=sorteio))+theme_bw()}
  }
  #=================
  print(graph)
}

