#' Utils: Summary of Analysis of Variance and Test of Means
#' @description Summarizes the output of the analysis of variance and the multiple comparisons test for completely randomized (DIC), randomized block (DBC) and Latin square (DQL) designs.
#' @author Gabriel Danilo Shimizu
#' @param analysis List with the analysis outputs of the DIC, DBC, DQL, FAT2DIC, FAT2DBC, PSUBDIC and PSUBDBC functions
#' @param design Type of experimental project (DIC, DBC, DQL, FAT2DIC, FAT2DBC, PSUBDIC or PSUBDBC)
#' @param round Number of decimal places
#' @param divisor Add divider between columns
#' @param inf Analysis of variance information (can be "p", "f", "QM" or "SQ")
#' @note Adding table divider can help to build tables in microsoft word. Copy console output, paste into MS Word, Insert, Table, Convert text to table, Separated text into:, Other: |.
#' @note The column names in the final output are imported from the ylab argument within each function.
#' @note This function is only for declared qualitative factors. In the case of a quantitative factor and the other qualitative in projects with two factors, this function will not work.
#' @note Triple factorials and split-split-plot do not work in this function.
#'
#' @export
#' @examples
#'
#'
#' library(AgroR)
#'
#' #=====================================
#' # DIC
#' #=====================================
#' data(pomegranate)
#' attach(pomegranate)
#' a=DIC(trat, WL, geom = "point", ylab = "WL")
#' b=DIC(trat, SS, geom = "point", ylab="SS")
#' c=DIC(trat, AT, geom = "point", ylab = "AT")
#' summarise_anova(analysis = list(a,b,c), divisor = TRUE)
#'
#' #=====================================
#' # DBC
#' #=====================================
#' data(soybean)
#' attach(soybean)
#' a=DBC(cult,bloc,prod,ylab = "Yield")
#' summarise_anova(list(a),design = "DBC")
#'
#' #=====================================
#' # FAT2DIC
#' #=====================================
#' data(corn)
#' attach(corn)
#' a=FAT2DIC(A, B, Resp, quali=c(TRUE, TRUE))
#' summarise_anova(list(a),design="FAT2DIC")

summarise_anova=function(analysis,
                         inf="p",
                         design="DIC",
                         round=3,
                         divisor=TRUE){
  if(design=="DIC"){
    nlinhas=length(analysis[[1]][[1]]$plot$dadosm$groups)
    infor=data.frame(matrix(ncol=length(analysis),nrow = nlinhas))
    trats=analysis[[1]][[1]]$plot$dadosm$trats
    variable=1:length(analysis)
    for(i in 1:length(analysis)){
      tests=analysis[[i]][[1]]$plot$test

      if(tests=="parametric"){
      variable[i]=analysis[[i]][[1]]$plot$ylab
      letra=analysis[[i]][[1]]$plot$dadosm$groups
      transf=analysis[[i]][[1]]$plot$transf
      if(transf==1){media=round(analysis[[i]][[1]]$plot$dadosm$resp,round)}
      if(transf!=1){media=round(analysis[[i]][[1]]$plot$dadosm$respO,round)}
      infor[,i]=paste(media,letra)}
      if(tests=="noparametric"){
        variable[i]=analysis[[i]][[1]]$plot$ylab
        letra=analysis[[i]][[1]]$plot$dadosm$groups
        media=round(analysis[[i]][[1]]$plot$dadosm$media,round)
        infor[,i]=paste(media,letra)}
    }
    names(infor)=variable
    rownames(infor)=trats

    cvs=data.frame(matrix(ncol=length(analysis),nrow = 1))
    variable=1:length(analysis)
    for(i in 1:length(analysis)){
      tests=analysis[[i]][[1]]$plot$test
      if(tests=="parametric"){
      variable[i]=analysis[[i]][[1]]$plot$ylab
      cvs[,i]=round(sqrt(analysis[[i]][[1]]$plot$a$`Mean Sq`[2])/
                      mean(analysis[[i]][[1]]$plot$resp)*100,round)}
      if(tests=="noparametric"){
        variable[i]=analysis[[i]][[1]]$plot$ylab
        cvs[,i]=" "}}
    rownames(cvs)="CV(%)"
    names(cvs)=variable

    transf=data.frame(matrix(ncol=length(analysis),nrow = 1))
    variable=1:length(analysis)
    for(i in 1:length(analysis)){
      tests=analysis[[i]][[1]]$plot$test
      variable[i]=analysis[[i]][[1]]$plot$ylab

      if(tests=="parametric"){
      tran=analysis[[i]][[1]]$plot$transf
      if(tran==1){tran=1}else{tran=tran}
      if(tran==0){tran="log"}
      if(tran==0.5){tran="sqrt(x)"}
      if(tran==-0.5){tran="1/sqrt(x)"}
      if(tran==-1){tran="1/x"}
      transf[,i]=ifelse(tran==1,"No transf",tran)}

      if(tests=="noparametric"){
        transf[,i]=""}
      }
    rownames(transf)="Transformation"
    names(transf)=variable

    n=5
    nc=1
    infor1=data.frame(matrix(ncol=length(analysis),nrow = 1))
    variable=1:length(analysis)
    for(i in 1:length(analysis)){
      tests=analysis[[i]][[1]]$plot$test
      if(tests=="parametric"){
      variable[i]=analysis[[i]][[1]]$plot$ylab
      pvalor=round(analysis[[i]][[1]]$plot$anava[1,n],round)
      infor1[,i]=ifelse(pvalor<0.001,"p<0.001",pvalor)}
      if(tests=="noparametric"){
        variable[i]=analysis[[i]][[1]]$plot$ylab
        pvalor=round(analysis[[i]][[1]]$plot$krusk$statistics[3][[1]],round)
        infor1[,i]=ifelse(pvalor<0.001,"p<0.001",pvalor)}}

    names(infor1)=variable
    rownames(infor1)="p-value"

    n=4
    nc=1
    infor2=data.frame(matrix(ncol=length(analysis),nrow = 1))
    variable=1:length(analysis)
    for(i in 1:length(analysis)){
      tests=analysis[[i]][[1]]$plot$test
      if(tests=="parametric"){
        variable[i]=analysis[[i]][[1]]$plot$ylab
      infor2[,i]=round(analysis[[i]][[1]]$plot$anava[1,n],round)}
      if(tests=="noparametric"){
        variable[i]=analysis[[i]][[1]]$plot$ylab
        infor2[,i]=paste(round(analysis[[i]][[1]]$plot$krusk$statistics[1][[1]],round),
                         "(Chisq)")}}
    names(infor2)=variable
    rownames(infor2)="F"

    n=3
    nc=2
    infor3=data.frame(matrix(ncol=length(analysis),nrow = 2))
    variable=1:length(analysis)
    for(i in 1:length(analysis)){
      tests=analysis[[i]][[1]]$plot$test
      if(tests=="parametric"){
      variable[i]=analysis[[i]][[1]]$plot$ylab
      infor3[,i]=round(analysis[[i]][[1]]$plot$anava[1:nc,n],round)}
      if(tests=="noparametric"){
        variable[i]=analysis[[i]][[1]]$plot$ylab
        infor3[,i]=""}
      }
    names(infor3)=variable
    rownames(infor3)=c("QM_tr","QM_r")

    n=2
    nc=2
    infor4=data.frame(matrix(ncol=length(analysis),nrow = 2))
    variable=1:length(analysis)
    for(i in 1:length(analysis)){
      tests=analysis[[i]][[1]]$plot$test
      if(tests=="parametric"){
      variable[i]=analysis[[i]][[1]]$plot$ylab
      infor4[,i]=round(analysis[[i]][[1]]$plot$anava[1:nc,n],round)}
      if(tests=="noparametric"){
        variable[i]=analysis[[i]][[1]]$plot$ylab
        infor4[,i]=""}}
    names(infor4)=variable
    rownames(infor4)=c("SQ_tr","SQ_r")
    if(inf=="p"){juntos=rbind(infor,cvs,infor1,transf)}
    if(inf=="f"){juntos=rbind(infor,cvs,infor2,transf)}
    if(inf=="QM"){juntos=rbind(infor,cvs,infor3,transf)}
    if(inf=="SQ"){juntos=rbind(infor,cvs,infor4,transf)}
    if(inf=="all"){juntos=rbind(infor,cvs,infor1,infor2,infor3,infor4,transf)}
  if(divisor==TRUE){
  nl=nrow(juntos)
  nc=ncol(juntos)
  market=data.frame(matrix(rep("|",nl*nc),ncol=nc,nrow = nl))
  juntosnovo=cbind("|"=rep("|",nl),juntos,market)
  for(i in 1:nc){
    ordem=matrix(1:(nc*2),nrow=2)
    nomes=colnames(juntos)
    juntosnovo[,c(ordem[,i]+1)]=cbind(juntos[,i],market[,i])
    colnames(juntosnovo)[1]="|"
    colnames(juntosnovo)[c(ordem[,i]+1)]=c(nomes[i],"|")}
  juntos=juntosnovo}
  }




  if(design=="DBC"){
    nlinhas=length(analysis[[1]][[1]]$plot$dadosm$groups)
    infor=data.frame(matrix(ncol=length(analysis),nrow = nlinhas))
    trats=analysis[[1]][[1]]$plot$dadosm$trats
    variable=1:length(analysis)
    for(i in 1:length(analysis)){
      variable[i]=analysis[[i]][[1]]$plot$ylab
      tests=analysis[[i]][[1]]$plot$test

      if(tests=="parametric"){
      letra=analysis[[i]][[1]]$plot$dadosm$groups
      transf=analysis[[i]][[1]]$plot$transf
      if(transf==1){media=round(analysis[[i]][[1]]$plot$dadosm$resp,round)}
      if(transf!=1){media=round(analysis[[i]][[1]]$plot$dadosm$respO,round)}
      infor[,i]=paste(media,letra)}
      if(tests=="noparametric"){
        letra=analysis[[i]][[1]]$plot$dadosm$groups
        media=round(analysis[[i]][[1]]$plot$dadosm$media,round)
        infor[,i]=paste(media,letra)}
    }
    names(infor)=variable
    rownames(infor)=trats
    cvs=data.frame(matrix(ncol=length(analysis),nrow = 1))
    variable=1:length(analysis)
    for(i in 1:length(analysis)){
      variable[i]=analysis[[i]][[1]]$plot$ylab
      tests=analysis[[i]][[1]]$plot$test

      if(tests=="parametric"){
      cvs[,i]=round(sqrt(analysis[[i]][[1]]$plot$a$`Mean Sq`[3])/
                      mean(analysis[[i]][[1]]$plot$resp)*100,round)}
      if(tests=="noparametric"){
        cvs[,i]="-"}}
    rownames(cvs)="CV(%)"
    names(cvs)=variable

    transf=data.frame(matrix(ncol=length(analysis),nrow = 1))
    variable=1:length(analysis)
    for(i in 1:length(analysis)){
      variable[i]=analysis[[i]][[1]]$plot$ylab
      tests=analysis[[i]][[1]]$plot$test

      if(tests=="parametric"){
      tran=analysis[[i]][[1]]$plot$transf
      if(tran==1){tran=1}else{tran=tran}
      if(tran==0){tran="log"}
      if(tran==0.5){tran="sqrt(x)"}
      if(tran==-0.5){tran="1/sqrt(x)"}
      if(tran==-1){tran="1/x"}
      transf[,i]=ifelse(tran==1,"No transf",tran)}
      if(tests=="noparametric"){transf[,i]="-"}}
    rownames(transf)="Transformation"
    names(transf)=variable

    n=5
    nc=2
    infor1=data.frame(matrix(ncol=length(analysis),nrow = 2))
    variable=1:length(analysis)
    for(i in 1:length(analysis)){
      variable[i]=analysis[[i]][[1]]$plot$ylab
      tests=analysis[[i]][[1]]$plot$test

      if(tests=="parametric"){
        pvalor=round(analysis[[i]][[1]]$plot$anava[1:2,n],round)
      infor1[,i]=ifelse(pvalor<0.001,"p<0.001",pvalor)}
      if(tests=="noparametric"){
        pvalor=round(analysis[[i]][[1]]$plot$fried$statistics[6],round)
        infor1[,i]=c(ifelse(pvalor<0.001,"p<0.001",pvalor[[1]]),
                     "-")}
      }
    names(infor1)=variable
    rownames(infor1)=c("p_tr","p_bl")


    n=4
    nc=2
    infor2=data.frame(matrix(ncol=length(analysis),nrow = 2))
    variable=1:length(analysis)
    for(i in 1:length(analysis)){
      variable[i]=analysis[[i]][[1]]$plot$ylab
      tests=analysis[[i]][[1]]$plot$test

      if(tests=="parametric"){
        infor2[,i]=round(analysis[[i]][[1]]$plot$anava[1:2,n],round)}
      if(tests=="noparametric"){
        infor2[,i]=c(round(analysis[[i]][[1]]$plot$fried$statistics[4][[1]],round),
                     "-")}}
    names(infor2)=variable
    rownames(infor2)=c("F_tr","F_bl")

    n=3
    nc=3
    infor3=data.frame(matrix(ncol=length(analysis),nrow = 3))
    variable=1:length(analysis)
    for(i in 1:length(analysis)){
      tests=analysis[[i]][[1]]$plot$test

      if(tests=="parametric"){
        variable[i]=analysis[[i]][[1]]$plot$ylab
      infor3[,i]=round(analysis[[i]][[1]]$plot$anava[1:nc,n],round)}
      if(tests=="noparametric"){
        variable[i]=analysis[[i]][[1]]$plot$ylab
        infor3[,i]=c("-","-","-")}
      }
    names(infor3)=variable
    rownames(infor3)=c("QM_tr","QM_bl","QM_r")

    n=2
    nc=3
    infor4=data.frame(matrix(ncol=length(analysis),nrow = 3))
    variable=1:length(analysis)
    for(i in 1:length(analysis)){
      variable[i]=analysis[[i]][[1]]$plot$ylab
      tests=analysis[[i]][[1]]$plot$test
      if(tests=="parametric"){
        infor4[,i]=round(analysis[[i]][[1]]$plot$anava[1:nc,n],round)}
      if(tests=="noparametric"){
        infor4[,i]=c("-","-","-")}}
    names(infor4)=variable
    rownames(infor4)=c("SQ_tr","SQ_bl","SQ_r")
    if(inf=="p"){juntos=rbind(infor,cvs,infor1,transf)}
    if(inf=="f"){juntos=rbind(infor,cvs,infor2,transf)}
    if(inf=="QM"){juntos=rbind(infor,cvs,infor3,transf)}
    if(inf=="SQ"){juntos=rbind(infor,cvs,infor4,transf)}
    if(inf=="all"){juntos=rbind(infor,cvs,infor1,infor2,infor3,infor4,transf)}

    if(divisor==TRUE){
      nl=nrow(juntos)
      nc=ncol(juntos)
      market=data.frame(matrix(rep("|",nl*nc),ncol=nc,nrow = nl))
      juntosnovo=cbind("|"=rep("|",nl),juntos,market)
      for(i in 1:nc){
        ordem=matrix(1:(nc*2),nrow=2)
        nomes=colnames(juntos)
        juntosnovo[,c(ordem[,i]+1)]=cbind(juntos[,i],market[,i])
        colnames(juntosnovo)[1]="|"
        colnames(juntosnovo)[c(ordem[,i]+1)]=c(nomes[i],"|")}
      juntos=juntosnovo}
  }
  if(design=="DQL"){
    nlinhas=length(analysis[[1]][[1]]$plot$dadosm$groups)
    infor=data.frame(matrix(ncol=length(analysis),nrow = nlinhas))
    trats=analysis[[1]][[1]]$plot$dadosm$trats
    variable=1:length(analysis)
    for(i in 1:length(analysis)){
      variable[i]=analysis[[i]][[1]]$plot$ylab
      letra=analysis[[i]][[1]]$plot$dadosm$groups
      transf=analysis[[i]][[1]]$plot$transf
      if(transf==1){media=round(analysis[[i]][[1]]$plot$dadosm$resp,round)}
      if(transf!=1){media=round(analysis[[i]][[1]]$plot$dadosm$respO,round)}
      infor[,i]=paste(media,letra)
    }
    names(infor)=variable
    rownames(infor)=trats
    cvs=data.frame(matrix(ncol=length(analysis),nrow = 1))
    variable=1:length(analysis)
    for(i in 1:length(analysis)){
      variable[i]=analysis[[i]][[1]]$plot$ylab
      cvs[,i]=round(sqrt(analysis[[i]][[1]]$plot$a$`Mean Sq`[4])/
                      mean(analysis[[i]][[1]]$plot$resp)*100,round)
    }
    rownames(cvs)="CV(%)"
    names(cvs)=variable

    transf=data.frame(matrix(ncol=length(analysis),nrow = 1))
    variable=1:length(analysis)
    for(i in 1:length(analysis)){
      variable[i]=analysis[[i]][[1]]$plot$ylab
      tran=analysis[[i]][[1]]$plot$transf
      if(tran==1){tran=1}else{tran=tran}
      if(tran==0){tran="log"}
      if(tran==0.5){tran="sqrt(x)"}
      if(tran==-0.5){tran="1/sqrt(x)"}
      if(tran==-1){tran="1/x"}
      transf[,i]=ifelse(tran==1,"No transf",tran)}
    rownames(transf)="Transformation"
    names(transf)=variable

    n=5
    nc=3
    infor1=data.frame(matrix(ncol=length(analysis),nrow = 3))
    variable=1:length(analysis)
    for(i in 1:length(analysis)){
      variable[i]=analysis[[i]][[1]]$plot$ylab
      pvalor=round(analysis[[i]][[1]]$plot$anava[1:nc,n],round)
      infor1[,i]=ifelse(pvalor<0.001,"p<0.001",pvalor)}
    names(infor1)=variable
    rownames(infor1)=c("p_tr","p_l","p_c")

    n=4
    nc=3
    infor2=data.frame(matrix(ncol=length(analysis),nrow = 3))
    variable=1:length(analysis)
    for(i in 1:length(analysis)){
      variable[i]=analysis[[i]][[1]]$plot$ylab
      infor2[,i]=round(analysis[[i]][[1]]$plot$anava[1:nc,n],round)}
    names(infor2)=variable
    rownames(infor2)=c("F_tr","F_l","F_c")
    n=3
    nc=4
    infor3=data.frame(matrix(ncol=length(analysis),nrow = 4))
    variable=1:length(analysis)
    for(i in 1:length(analysis)){
      variable[i]=analysis[[i]][[1]]$plot$ylab
      infor3[,i]=round(analysis[[i]][[1]]$plot$anava[1:nc,n],round)}
    names(infor3)=variable
    rownames(infor3)=c("QM_tr","QM_l","QM_c","QM_r")
    n=2
    nc=4
    infor4=data.frame(matrix(ncol=length(analysis),nrow = 4))
    variable=1:length(analysis)
    for(i in 1:length(analysis)){
      variable[i]=analysis[[i]][[1]]$plot$ylab
      infor4[,i]=round(analysis[[i]][[1]]$plot$anava[1:nc,n],round)}
    names(infor4)=variable
    rownames(infor4)=c("SQ_tr","SQ_l","SQ_c","SQ_r")
    if(inf=="p"){juntos=rbind(infor,cvs,infor1,transf)}
    if(inf=="f"){juntos=rbind(infor,cvs,infor2,transf)}
    if(inf=="QM"){juntos=rbind(infor,cvs,infor3,transf)}
    if(inf=="SQ"){juntos=rbind(infor,cvs,infor4,transf)}
    if(inf=="all"){juntos=rbind(infor,cvs,infor1,infor2,infor3,infor4,transf)}

    if(divisor==TRUE){
      nl=nrow(juntos)
      nc=ncol(juntos)
      market=data.frame(matrix(rep("|",nl*nc),ncol=nc,nrow = nl))
      juntosnovo=cbind("|"=rep("|",nl),juntos,market)
      for(i in 1:nc){
        ordem=matrix(1:(nc*2),nrow=2)
        nomes=colnames(juntos)
        juntosnovo[,c(ordem[,i]+1)]=cbind(juntos[,i],market[,i])
        colnames(juntosnovo)[1]="|"
        colnames(juntosnovo)[c(ordem[,i]+1)]=c(nomes[i],"|")}
      juntos=juntosnovo}
  }
  if(design=="FAT2DIC"){
    cvs=data.frame(matrix(ncol=length(analysis),nrow = 1))
    variable=1:length(analysis)
    for(i in 1:length(analysis)){
      variable[i]=analysis[[i]][[1]]$plot$ylab
      cvs[,i]=round(sqrt(analysis[[i]][[1]]$plot$a$`Mean Sq`[4])/
                      mean(analysis[[i]][[1]]$plot$resp)*100,round)}
    rownames(cvs)="CV(%)"
    names(cvs)=variable

    transf=data.frame(matrix(ncol=length(analysis),nrow = 1))
    variable=1:length(analysis)
    for(i in 1:length(analysis)){
      variable[i]=analysis[[i]][[1]]$plot$ylab
      tran=analysis[[i]][[1]]$plot$transf
      if(tran==1){tran=1}else{tran=tran}
      if(tran==0){tran="log"}
      if(tran==0.5){tran="sqrt(x)"}
      if(tran==-0.5){tran="1/sqrt(x)"}
      if(tran==-1){tran="1/x"}
      transf[,i]=ifelse(tran==1,"No transf",tran)}
    rownames(transf)="Transformation"
    names(transf)=variable

    n=5
    nc=3
    infor1=data.frame(matrix(ncol=length(analysis),nrow = 3))
    variable=1:length(analysis)
    for(i in 1:length(analysis)){
      variable[i]=analysis[[i]][[1]]$plot$ylab
      pvalor=round(analysis[[i]][[1]]$plot$anava[1:nc,n],round)
      infor1[,i]=ifelse(pvalor<0.001,"p<0.001",pvalor)}
    names(infor1)=variable
    rownames(infor1)=c("p_F1","p_F2","p_F1xF2")
    n=4
    nc=3
    infor2=data.frame(matrix(ncol=length(analysis),nrow = 3))
    variable=1:length(analysis)
    for(i in 1:length(analysis)){
      variable[i]=analysis[[i]][[1]]$plot$ylab
      infor2[,i]=round(analysis[[i]][[1]]$plot$anava[1:nc,n],round)}
    names(infor2)=variable
    rownames(infor2)=c("f_F1","f_F2","f_F1xF2")

    n=3
    nc=4
    infor3=data.frame(matrix(ncol=length(analysis),nrow = 4))
    variable=1:length(analysis)
    for(i in 1:length(analysis)){
      variable[i]=analysis[[i]][[1]]$plot$ylab
      infor3[,i]=round(analysis[[i]][[1]]$plot$anava[1:nc,n],round)}
    names(infor3)=variable
    rownames(infor3)=c("QM_F1","QM_F2","QM_F1xF2","QM_r")

    n=2
    nc=4
    infor4=data.frame(matrix(ncol=length(analysis),nrow = 4))
    variable=1:length(analysis)
    for(i in 1:length(analysis)){
      variable[i]=analysis[[i]][[1]]$plot$ylab
      infor4[,i]=round(analysis[[i]][[1]]$plot$anava[1:nc,n],round)}
    names(infor4)=variable
    rownames(infor4)=c("SQ_F1","SQ_F2","SQ_F1xF2","SQ_r")
    variable=1:length(analysis)
    nf1=analysis[[1]][[1]]$plot$anava[1,1]+1
    nf2=analysis[[1]][[1]]$plot$anava[2,1]+1
    nf1f2=nf1*nf2
    f1mean=data.frame(matrix(rep("-",(length(analysis)*nf1)),ncol=length(analysis),nrow = nf1))
    f2mean=data.frame(matrix(rep("-",(length(analysis)*nf2)),ncol=length(analysis),nrow = nf2))
    f1f2mean=data.frame(matrix(rep("-",(length(analysis)*nf1f2)),ncol=length(analysis),nrow = nf1f2))
    rownames(f1mean)=unique(analysis[[1]][[1]]$plot$f1)
    rownames(f2mean)=unique(analysis[[1]][[1]]$plot$f2)
    nomes=expand.grid(unique(analysis[[1]][[1]]$plot$f2),unique(analysis[[1]][[1]]$plot$f1))
    rownames(f1f2mean)=paste(nomes$Var2,nomes$Var1)
    for(i in 1:length(analysis)){
      sigF=analysis[[i]][[1]]$plot$alpha.f
      variable[i]=analysis[[i]][[1]]$plot$ylab
      pvalue=round(analysis[[i]][[1]]$plot$anava[1:3,5],round)
      if(pvalue[3]<sigF){f1f2mean[,i]=data.frame(analysis[[i]][[1]]$plot$graph$numero)
      rownames(f1f2mean)=rownames(analysis[[i]][[1]]$plot$graph)}else{
        if(pvalue[1]<sigF){f1mean[,i]=data.frame(analysis[[i]][[2]]$data$letra)
        rownames(f1mean)=rownames(analysis[[i]][[2]]$data)}
        if(pvalue[2]<sigF){f2mean[,i]=data.frame(analysis[[i]][[3]]$data$letra)
        rownames(f2mean)=rownames(analysis[[i]][[3]]$data)}}
    }
    names(f1mean)=variable
    names(f2mean)=variable
    names(f1f2mean)=variable
    # juntos=rbind(f1mean,f2mean,f1f2mean,cvs,infor1,infor2,infor3,infor4,transf)
    if(inf=="p"){juntos=rbind(f1mean,f2mean,f1f2mean,cvs,infor1,transf)}
    if(inf=="f"){juntos=rbind(f1mean,f2mean,f1f2mean,cvs,infor2,transf)}
    if(inf=="QM"){juntos=rbind(f1mean,f2mean,f1f2mean,cvs,infor3,transf)}
    if(inf=="SQ"){juntos=rbind(f1mean,f2mean,f1f2mean,cvs,infor4,transf)}
    if(inf=="all"){juntos=rbind(f1mean,f2mean,f1f2mean,cvs,infor1,infor2,infor3,infor4,transf)}

    if(divisor==TRUE){
      nl=nrow(juntos)
      nc=ncol(juntos)
      market=data.frame(matrix(rep("|",nl*nc),ncol=nc,nrow = nl))
      juntosnovo=cbind("|"=rep("|",nl),juntos,market)
      for(i in 1:nc){
        ordem=matrix(1:(nc*2),nrow=2)
        nomes=colnames(juntos)
        juntosnovo[,c(ordem[,i]+1)]=cbind(juntos[,i],market[,i])
        colnames(juntosnovo)[1]="|"
        colnames(juntosnovo)[c(ordem[,i]+1)]=c(nomes[i],"|")}
      juntos=juntosnovo}
  }

  if(design=="FAT2DBC"){
    cvs=data.frame(matrix(ncol=length(analysis),nrow = 1))
    variable=1:length(analysis)
    for(i in 1:length(analysis)){
      variable[i]=analysis[[i]][[1]]$plot$ylab
      cvs[,i]=round(sqrt(analysis[[i]][[1]]$plot$a$`Mean Sq`[5])/
                      mean(analysis[[i]][[1]]$plot$resp)*100,round)}
    rownames(cvs)="CV(%)"
    names(cvs)=variable
    transf=data.frame(matrix(ncol=length(analysis),nrow = 1))
    variable=1:length(analysis)
    for(i in 1:length(analysis)){
      variable[i]=analysis[[i]][[1]]$plot$ylab
      tran=analysis[[i]][[1]]$plot$transf
      if(tran==1){tran=1}else{tran=tran}
      if(tran==0){tran="log"}
      if(tran==0.5){tran="sqrt(x)"}
      if(tran==-0.5){tran="1/sqrt(x)"}
      if(tran==-1){tran="1/x"}
      transf[,i]=ifelse(tran==1,"No transf",tran)}
    rownames(transf)="Transformation"
    names(transf)=variable

    n=5
    nc=4
    infor1=data.frame(matrix(ncol=length(analysis),nrow = 4))
    variable=1:length(analysis)
    for(i in 1:length(analysis)){
      variable[i]=analysis[[i]][[1]]$plot$ylab
      pvalor=round(analysis[[i]][[1]]$plot$anava[1:nc,n],round)
      infor1[,i]=ifelse(pvalor<0.001,"p<0.001",pvalor)}
    names(infor1)=variable
    rownames(infor1)=c("p_F1","p_F2","p_bl","p_F1xF2")

    n=4
    nc=4
    infor2=data.frame(matrix(ncol=length(analysis),nrow = 4))
    variable=1:length(analysis)
    for(i in 1:length(analysis)){
      variable[i]=analysis[[i]][[1]]$plot$ylab
      infor2[,i]=round(analysis[[i]][[1]]$plot$anava[1:nc,n],round)}
    names(infor2)=variable
    rownames(infor2)=c("f_F1","f_F2","p_bl","f_F1xF2")

    n=3
    nc=5
    infor3=data.frame(matrix(ncol=length(analysis),nrow = 5))
    variable=1:length(analysis)
    for(i in 1:length(analysis)){
      variable[i]=analysis[[i]][[1]]$plot$ylab
      infor3[,i]=round(analysis[[i]][[1]]$plot$anava[1:nc,n],round)}
    names(infor3)=variable
    rownames(infor3)=c("QM_F1","QM_F2","p_bl","QM_F1xF2","QM_r")

    n=2
    nc=5
    infor4=data.frame(matrix(ncol=length(analysis),nrow = 5))
    variable=1:length(analysis)
    for(i in 1:length(analysis)){
      variable[i]=analysis[[i]][[1]]$plot$ylab
      infor4[,i]=round(analysis[[i]][[1]]$plot$anava[1:nc,n],round)}
    names(infor4)=variable
    rownames(infor4)=c("SQ_F1","SQ_F2","p_bl","SQ_F1xF2","SQ_r")

    variable=1:length(analysis)
    nf1=analysis[[1]][[1]]$plot$anava[1,1]+1
    nf2=analysis[[1]][[1]]$plot$anava[2,1]+1
    nf1f2=nf1*nf2
    f1mean=data.frame(matrix(rep("-",(length(analysis)*nf1)),ncol=length(analysis),nrow = nf1))
    f2mean=data.frame(matrix(rep("-",(length(analysis)*nf2)),ncol=length(analysis),nrow = nf2))
    f1f2mean=data.frame(matrix(rep("-",(length(analysis)*nf1f2)),ncol=length(analysis),nrow = nf1f2))
    rownames(f1mean)=unique(analysis[[1]][[1]]$plot$f1)
    rownames(f2mean)=unique(analysis[[1]][[1]]$plot$f2)
    nomes=expand.grid(unique(analysis[[1]][[1]]$plot$f2),unique(analysis[[1]][[1]]$plot$f1))
    rownames(f1f2mean)=paste(nomes$Var2,nomes$Var1)
    for(i in 1:length(analysis)){
      sigF=analysis[[i]][[1]]$plot$alpha.f
      variable[i]=analysis[[i]][[1]]$plot$ylab
      pvalue=round(analysis[[i]][[1]]$plot$anava[1:4,5],round)
      if(pvalue[4]<sigF){f1f2mean[,i]=data.frame(analysis[[i]][[1]]$plot$graph$numero)
      rownames(f1f2mean)=rownames(analysis[[i]][[1]]$plot$graph)}else{
        if(pvalue[1]<sigF){f1mean[,i]=data.frame(analysis[[i]][[2]]$data$letra)
        rownames(f1mean)=rownames(analysis[[i]][[2]]$data)}
        if(pvalue[2]<sigF){f2mean[,i]=data.frame(analysis[[i]][[3]]$data$letra)
        rownames(f2mean)=rownames(analysis[[i]][[3]]$data)}}
    }
    names(f1mean)=variable
    names(f2mean)=variable
    names(f1f2mean)=variable
    if(inf=="p"){juntos=rbind(f1mean,f2mean,f1f2mean,cvs,infor1,transf)}
    if(inf=="f"){juntos=rbind(f1mean,f2mean,f1f2mean,cvs,infor2,transf)}
    if(inf=="QM"){juntos=rbind(f1mean,f2mean,f1f2mean,cvs,infor3,transf)}
    if(inf=="SQ"){juntos=rbind(f1mean,f2mean,f1f2mean,cvs,infor4,transf)}
    if(inf=="all"){juntos=rbind(f1mean,f2mean,f1f2mean,cvs,infor1,infor2,infor3,infor4,transf)}


    if(divisor==TRUE){
      nl=nrow(juntos)
      nc=ncol(juntos)
      market=data.frame(matrix(rep("|",nl*nc),ncol=nc,nrow = nl))
      juntosnovo=cbind("|"=rep("|",nl),juntos,market)
      for(i in 1:nc){
        ordem=matrix(1:(nc*2),nrow=2)
        nomes=colnames(juntos)
        juntosnovo[,c(ordem[,i]+1)]=cbind(juntos[,i],market[,i])
        colnames(juntosnovo)[1]="|"
        colnames(juntosnovo)[c(ordem[,i]+1)]=c(nomes[i],"|")}
      juntos=juntosnovo}
  }
  if(design=="PSUBDIC"){
    cvs=data.frame(matrix(ncol=length(analysis),nrow = 2))
    variable=1:length(analysis)
    for(i in 1:length(analysis)){
      variable[i]=analysis[[i]][[1]]$plot$ylab
      cvs[1,i]=round(sqrt(analysis[[i]][[1]]$plot$anova1$`Mean Sq`[2])/
                       mean(analysis[[i]][[1]]$plot$resp)*100,round)
      cvs[2,i]=round(sqrt(analysis[[i]][[1]]$plot$anova1$`Mean Sq`[5])/
                       mean(analysis[[i]][[1]]$plot$resp)*100,round)}
    rownames(cvs)=c("CV1 (%)","CV2 (%)")
    names(cvs)=variable

    transf=data.frame(matrix(ncol=length(analysis),nrow = 1))
    variable=1:length(analysis)
    for(i in 1:length(analysis)){
      variable[i]=analysis[[i]][[1]]$plot$ylab
      tran=analysis[[i]][[1]]$plot$transf
      if(tran==1){tran=1}else{tran=tran}
      if(tran==0){tran="log"}
      if(tran==0.5){tran="sqrt(x)"}
      if(tran==-0.5){tran="1/sqrt(x)"}
      if(tran==-1){tran="1/x"}
      transf[,i]=ifelse(tran==1,"No transf",tran)}
    rownames(transf)="Transformation"
    names(transf)=variable

    # p-value
    n=5
    infor1=data.frame(matrix(ncol=length(analysis),nrow = 3))
    variable=1:length(analysis)
    for(i in 1:length(analysis)){
      variable[i]=analysis[[i]][[1]]$plot$ylab
      pvalor=round(analysis[[i]][[1]]$plot$anova1[c(1,3,4),n],round)
      infor1[,i]=ifelse(pvalor<0.001,"p<0.001",pvalor)}
    names(infor1)=variable
    rownames(infor1)=c("p_F1","p_F2","p_F1xF2")
    n=4
    infor2=data.frame(matrix(ncol=length(analysis),nrow = 3))
    variable=1:length(analysis)
    for(i in 1:length(analysis)){
      variable[i]=analysis[[i]][[1]]$plot$ylab
      infor2[,i]=round(analysis[[i]][[1]]$plot$anova1[c(1,3,4),n],round)}
    names(infor2)=variable
    rownames(infor2)=c("f_F1","f_F2","f_F1xF2")

    n=3
    infor3=data.frame(matrix(ncol=length(analysis),nrow = 5))
    variable=1:length(analysis)
    for(i in 1:length(analysis)){
      variable[i]=analysis[[i]][[1]]$plot$ylab
      infor3[,i]=round(analysis[[i]][[1]]$plot$anova1[c(1:5),n],round)}
    names(infor3)=variable
    rownames(infor3)=c("QM_F1","QM_error1","QM_F2","QM_F1xF2","QM_error2")

    n=2
    infor4=data.frame(matrix(ncol=length(analysis),nrow = 5))
    variable=1:length(analysis)
    for(i in 1:length(analysis)){
      variable[i]=analysis[[i]][[1]]$plot$ylab
      infor4[,i]=round(analysis[[i]][[1]]$plot$anova1[1:5,n],round)}
    names(infor4)=variable
    rownames(infor4)=c("SQ_F1","SQ_error1","SQ_F2","SQ_F1xF2","SQ_error2")

    variable=1:length(analysis)
    nf1=analysis[[1]][[1]]$plot$anova1[1,1]+1
    nf2=analysis[[1]][[1]]$plot$anova1[3,1]+1
    nf1f2=nf1*nf2
    f1mean=data.frame(matrix(rep("-",(length(analysis)*nf1)),ncol=length(analysis),nrow = nf1))
    f2mean=data.frame(matrix(rep("-",(length(analysis)*nf2)),ncol=length(analysis),nrow = nf2))
    f1f2mean=data.frame(matrix(rep("-",(length(analysis)*nf1f2)),ncol=length(analysis),nrow = nf1f2))
    rownames(f1mean)=unique(analysis[[1]][[1]]$plot$f1)
    rownames(f2mean)=unique(analysis[[1]][[1]]$plot$f2)
    nomes=expand.grid(unique(analysis[[1]][[1]]$plot$f2),unique(analysis[[1]][[1]]$plot$f1))
    rownames(f1f2mean)=paste(nomes$Var2,nomes$Var1)
    for(i in 1:length(analysis)){
      sigF=analysis[[i]][[1]]$plot$alpha.f
      variable[i]=analysis[[i]][[1]]$plot$ylab
      pvalue=round(analysis[[i]][[1]]$plot$anova1[c(1,3,4),5],round)
      if(pvalue[3]<0.05){f1f2mean[,i]=data.frame(analysis[[i]][[1]]$plot$graph$numero)
      rownames(f1f2mean)=rownames(analysis[[i]][[1]]$plot$graph)}else{
        if(pvalue[1]<0.05){f1mean[,i]=data.frame(analysis[[i]][[2]]$data$letra)
        rownames(f1mean)=rownames(analysis[[i]][[2]]$data)}
        if(pvalue[2]<0.05){f2mean[,i]=data.frame(analysis[[i]][[3]]$data$letra)
        rownames(f2mean)=rownames(analysis[[i]][[3]]$data)}}
    }
    names(f1mean)=variable
    names(f2mean)=variable
    names(f1f2mean)=variable
    if(inf=="p"){juntos=rbind(f1mean,f2mean,f1f2mean,cvs,infor1,transf)}
    if(inf=="f"){juntos=rbind(f1mean,f2mean,f1f2mean,cvs,infor2,transf)}
    if(inf=="QM"){juntos=rbind(f1mean,f2mean,f1f2mean,cvs,infor3,transf)}
    if(inf=="SQ"){juntos=rbind(f1mean,f2mean,f1f2mean,cvs,infor4,transf)}
    if(inf=="all"){juntos=rbind(f1mean,f2mean,f1f2mean,cvs,infor1,infor2,infor3,infor4,transf)}

    if(divisor==TRUE){
      nl=nrow(juntos)
      nc=ncol(juntos)
      market=data.frame(matrix(rep("|",nl*nc),ncol=nc,nrow = nl))
      juntosnovo=cbind("|"=rep("|",nl),juntos,market)
      for(i in 1:nc){
        ordem=matrix(1:(nc*2),nrow=2)
        nomes=colnames(juntos)
        juntosnovo[,c(ordem[,i]+1)]=cbind(juntos[,i],market[,i])
        colnames(juntosnovo)[1]="|"
        colnames(juntosnovo)[c(ordem[,i]+1)]=c(nomes[i],"|")}
      juntos=juntosnovo}
  }
  if(design=="PSUBDBC"){
    cvs=data.frame(matrix(ncol=length(analysis),nrow = 2))
    variable=1:length(analysis)
    for(i in 1:length(analysis)){
      variable[i]=analysis[[i]][[1]]$plot$ylab
      cvs[1,i]=round(sqrt(analysis[[i]][[1]]$plot$anova1$`Mean Sq`[3])/
                       mean(analysis[[i]][[1]]$plot$resp)*100,round)
      cvs[2,i]=round(sqrt(analysis[[i]][[1]]$plot$anova1$`Mean Sq`[6])/
                       mean(analysis[[i]][[1]]$plot$resp)*100,round)}
    rownames(cvs)=c("CV1 (%)","CV2 (%)")
    names(cvs)=variable

    transf=data.frame(matrix(ncol=length(analysis),nrow = 1))
    variable=1:length(analysis)
    for(i in 1:length(analysis)){
      variable[i]=analysis[[i]][[1]]$plot$ylab
      tran=analysis[[i]][[1]]$plot$transf
      if(tran==1){tran=1}else{tran=tran}
      if(tran==0){tran="log"}
      if(tran==0.5){tran="sqrt(x)"}
      if(tran==-0.5){tran="1/sqrt(x)"}
      if(tran==-1){tran="1/x"}
      transf[,i]=ifelse(tran==1,"No transf",tran)}
    rownames(transf)="Transformation"
    names(transf)=variable

    # p-value
    n=5
    infor1=data.frame(matrix(ncol=length(analysis),nrow = 4))
    variable=1:length(analysis)
    for(i in 1:length(analysis)){
      variable[i]=analysis[[i]][[1]]$plot$ylab
      pvalor=round(analysis[[i]][[1]]$plot$anova1[c(1,2,4,5),n],round)
      infor1[,i]=ifelse(pvalor<0.001,"p<0.001",pvalor)}
    names(infor1)=variable
    rownames(infor1)=c("p_F1","p_bl","p_F2","p_F1xF2")
    n=4
    infor2=data.frame(matrix(ncol=length(analysis),nrow = 4))
    variable=1:length(analysis)
    for(i in 1:length(analysis)){
      variable[i]=analysis[[i]][[1]]$plot$ylab
      infor2[,i]=round(analysis[[i]][[1]]$plot$anova1[c(1,2,4,5),n],round)}
    names(infor2)=variable
    rownames(infor2)=c("f_F1","f_bl","f_F2","f_F1xF2")

    n=3
    infor3=data.frame(matrix(ncol=length(analysis),nrow = 6))
    variable=1:length(analysis)
    for(i in 1:length(analysis)){
      variable[i]=analysis[[i]][[1]]$plot$ylab
      infor3[,i]=round(analysis[[i]][[1]]$plot$anova1[c(1:6),n],round)}
    names(infor3)=variable
    rownames(infor3)=c("QM_F1","QM_bl","QM_error1","QM_F2","QM_F1xF2","QM_error2")

    n=2
    infor4=data.frame(matrix(ncol=length(analysis),nrow = 6))
    variable=1:length(analysis)
    for(i in 1:length(analysis)){
      variable[i]=analysis[[i]][[1]]$plot$ylab
      infor4[,i]=round(analysis[[i]][[1]]$plot$anova1[1:6,n],round)}
    names(infor4)=variable
    rownames(infor4)=c("SQ_F1","SQ_bl","SQ_error1","SQ_F2","SQ_F1xF2","SQ_error2")

    variable=1:length(analysis)
    nf1=analysis[[1]][[1]]$plot$anova1[1,1]+1
    nf2=analysis[[1]][[1]]$plot$anova1[4,1]+1
    nf1f2=nf1*nf2
    f1mean=data.frame(matrix(rep("-",(length(analysis)*nf1)),ncol=length(analysis),nrow = nf1))
    f2mean=data.frame(matrix(rep("-",(length(analysis)*nf2)),ncol=length(analysis),nrow = nf2))
    f1f2mean=data.frame(matrix(rep("-",(length(analysis)*nf1f2)),ncol=length(analysis),nrow = nf1f2))
    rownames(f1mean)=unique(analysis[[1]][[1]]$plot$f1)
    rownames(f2mean)=unique(analysis[[1]][[1]]$plot$f2)
    nomes=expand.grid(unique(analysis[[1]][[1]]$plot$f2),unique(analysis[[1]][[1]]$plot$f1))
    rownames(f1f2mean)=paste(nomes$Var2,nomes$Var1)
    for(i in 1:length(analysis)){
      sigF=analysis[[i]][[1]]$plot$alpha.f
      variable[i]=analysis[[i]][[1]]$plot$ylab
      pvalue=round(analysis[[i]][[1]]$plot$anova1[c(1,4,5),5],round)
      if(pvalue[3]<0.05){f1f2mean[,i]=data.frame(analysis[[i]][[1]]$plot$graph$numero)
      rownames(f1f2mean)=rownames(analysis[[i]][[1]]$plot$graph)}else{
        if(pvalue[1]<0.05){f1mean[,i]=data.frame(analysis[[i]][[2]]$data$letra)
        rownames(f1mean)=rownames(analysis[[i]][[2]]$data)}
        if(pvalue[2]<0.05){f2mean[,i]=data.frame(analysis[[i]][[3]]$data$letra)
        rownames(f2mean)=rownames(analysis[[i]][[3]]$data)}}
    }
    names(f1mean)=variable
    names(f2mean)=variable
    names(f1f2mean)=variable
    if(inf=="p"){juntos=rbind(f1mean,f2mean,f1f2mean,cvs,infor1,transf)}
    if(inf=="f"){juntos=rbind(f1mean,f2mean,f1f2mean,cvs,infor2,transf)}
    if(inf=="QM"){juntos=rbind(f1mean,f2mean,f1f2mean,cvs,infor3,transf)}
    if(inf=="SQ"){juntos=rbind(f1mean,f2mean,f1f2mean,cvs,infor4,transf)}
    if(inf=="all"){juntos=rbind(f1mean,f2mean,f1f2mean,cvs,infor1,infor2,infor3,infor4,transf)}



    if(divisor==TRUE){
      nl=nrow(juntos)
      nc=ncol(juntos)
      market=data.frame(matrix(rep("|",nl*nc),ncol=nc,nrow = nl))
      juntosnovo=cbind("|"=rep("|",nl),juntos,market)
      for(i in 1:nc){
        ordem=matrix(1:(nc*2),nrow=2)
        nomes=colnames(juntos)
        juntosnovo[,c(ordem[,i]+1)]=cbind(juntos[,i],market[,i])
        colnames(juntosnovo)[1]="|"
        colnames(juntosnovo)[c(ordem[,i]+1)]=c(nomes[i],"|")}
      juntos=juntosnovo}
  }

  if(design=="FAT2DBC.ad"){
    cvs=data.frame(matrix(ncol=length(analysis),nrow = 1))
    variable=1:length(analysis)
    for(i in 1:length(analysis)){
      variable[i]=analysis[[i]][[1]]$plot$ylab
      cvs[,i]=round(sqrt(analysis[[i]][[1]]$plot$anava$`Mean Sq`[6])/
                      mean(c(analysis[[i]][[1]]$plot$resp,
                             analysis[[i]][[1]]$plot$respAd))*100,round)}
    rownames(cvs)="CV(%)"
    names(cvs)=variable

    respAd=data.frame(matrix(ncol=length(analysis),nrow = 2))
    variable=1:length(analysis)
    for(i in 1:length(analysis)){
      variable[i]=analysis[[i]][[1]]$plot$ylab
      respAd[1,i]=round(mean(analysis[[i]][[1]]$plot$response),round)
      respAd[2,i]=round(mean(analysis[[i]][[1]]$plot$responseAd),round)}
    rownames(respAd)=c("Factorial","RespAd")
    names(respAd)=variable

    transf=data.frame(matrix(ncol=length(analysis),nrow = 1))
    variable=1:length(analysis)
    for(i in 1:length(analysis)){
      variable[i]=analysis[[i]][[1]]$plot$ylab
      tran=analysis[[i]][[1]]$plot$transf
      if(tran==1){tran=1}else{tran=tran}
      if(tran==0){tran="log"}
      if(tran==0.5){tran="sqrt(x)"}
      if(tran==-0.5){tran="1/sqrt(x)"}
      if(tran==-1){tran="1/x"}
      transf[,i]=ifelse(tran==1,"No transf",tran)}
    rownames(transf)="Transformation"
    names(transf)=variable

    n=5
    nc=5
    infor1=data.frame(matrix(ncol=length(analysis),nrow = 5))
    variable=1:length(analysis)
    for(i in 1:length(analysis)){
      variable[i]=analysis[[i]][[1]]$plot$ylab
      pvalor=round(analysis[[i]][[1]]$plot$anava[1:nc,n],round)
      infor1[,i]=ifelse(pvalor<0.001,"p<0.001",pvalor)}
    names(infor1)=variable
    rownames(infor1)=c("p_F1","p_F2","p_bl","p_F1xF2","p_Fat x Ad")

    n=4
    nc=5
    infor2=data.frame(matrix(ncol=length(analysis),nrow = 5))
    variable=1:length(analysis)
    for(i in 1:length(analysis)){
      variable[i]=analysis[[i]][[1]]$plot$ylab
      infor2[,i]=round(analysis[[i]][[1]]$plot$anava[1:nc,n],round)}
    names(infor2)=variable
    rownames(infor2)=c("f_F1","f_F2","p_bl","f_F1xF2","f_Fat x Ad")

    n=3
    nc=6
    infor3=data.frame(matrix(ncol=length(analysis),nrow = 6))
    variable=1:length(analysis)
    for(i in 1:length(analysis)){
      variable[i]=analysis[[i]][[1]]$plot$ylab
      infor3[,i]=round(analysis[[i]][[1]]$plot$anava[1:nc,n],round)}
    names(infor3)=variable
    rownames(infor3)=c("QM_F1","QM_F2","p_bl","QM_F1xF2","QM_Fat x Ad","QM_r")

    n=2
    nc=6
    infor4=data.frame(matrix(ncol=length(analysis),nrow = 6))
    variable=1:length(analysis)
    for(i in 1:length(analysis)){
      variable[i]=analysis[[i]][[1]]$plot$ylab
      infor4[,i]=round(analysis[[i]][[1]]$plot$anava[1:nc,n],round)}
    names(infor4)=variable
    rownames(infor4)=c("SQ_F1","SQ_F2","p_bl","SQ_F1xF2","SQ_Fat x Ad","SQ_r")

    variable=1:length(analysis)
    nf1=analysis[[1]][[1]]$plot$anava[1,1]+1
    nf2=analysis[[1]][[1]]$plot$anava[2,1]+1
    nf1f2=nf1*nf2
    f1mean=data.frame(matrix(rep("-",(length(analysis)*nf1)),ncol=length(analysis),nrow = nf1))
    f2mean=data.frame(matrix(rep("-",(length(analysis)*nf2)),ncol=length(analysis),nrow = nf2))
    f1f2mean=data.frame(matrix(rep("-",(length(analysis)*nf1f2)),ncol=length(analysis),nrow = nf1f2))
    rownames(f1mean)=unique(analysis[[1]][[1]]$plot$f1)
    rownames(f2mean)=unique(analysis[[1]][[1]]$plot$f2)
    nomes=expand.grid(unique(analysis[[1]][[1]]$plot$f2),unique(analysis[[1]][[1]]$plot$f1))
    rownames(f1f2mean)=paste(nomes$Var2,nomes$Var1)
    for(i in 1:length(analysis)){
      sigF=analysis[[i]][[1]]$plot$alpha.f
      variable[i]=analysis[[i]][[1]]$plot$ylab
      pvalue=round(analysis[[i]][[1]]$plot$anava[1:4,5],round)
      if(pvalue[4]<sigF){f1f2mean[,i]=data.frame(analysis[[i]][[1]]$plot$graph$numero)
      rownames(f1f2mean)=rownames(analysis[[i]][[1]]$plot$graph)}else{
        if(pvalue[1]<sigF){f1mean[,i]=data.frame(analysis[[i]][[2]]$data$letra)
        rownames(f1mean)=rownames(analysis[[i]][[2]]$data)}
        if(pvalue[2]<sigF){f2mean[,i]=data.frame(analysis[[i]][[3]]$data$letra)
        rownames(f2mean)=rownames(analysis[[i]][[3]]$data)}}
    }
    names(f1mean)=variable
    names(f2mean)=variable
    names(f1f2mean)=variable
    if(inf=="p"){juntos=rbind(f1mean,f2mean,f1f2mean,cvs,infor1,respAd,transf)}
    if(inf=="f"){juntos=rbind(f1mean,f2mean,f1f2mean,cvs,infor2,respAd,transf)}
    if(inf=="QM"){juntos=rbind(f1mean,f2mean,f1f2mean,cvs,infor3,respAd,transf)}
    if(inf=="SQ"){juntos=rbind(f1mean,f2mean,f1f2mean,cvs,infor4,respAd,transf)}
    if(inf=="all"){juntos=rbind(f1mean,f2mean,f1f2mean,cvs,infor1,infor2,infor3,infor4,respAd,transf)}


    if(divisor==TRUE){
      nl=nrow(juntos)
      nc=ncol(juntos)
      market=data.frame(matrix(rep("|",nl*nc),ncol=nc,nrow = nl))
      juntosnovo=cbind("|"=rep("|",nl),juntos,market)
      for(i in 1:nc){
        ordem=matrix(1:(nc*2),nrow=2)
        nomes=colnames(juntos)
        juntosnovo[,c(ordem[,i]+1)]=cbind(juntos[,i],market[,i])
        colnames(juntosnovo)[1]="|"
        colnames(juntosnovo)[c(ordem[,i]+1)]=c(nomes[i],"|")}
      juntos=juntosnovo}
  }

  if(design=="FAT2DIC.ad"){
    cvs=data.frame(matrix(ncol=length(analysis),nrow = 1))
    variable=1:length(analysis)
    for(i in 1:length(analysis)){
      variable[i]=analysis[[i]][[1]]$plot$ylab
      cvs[,i]=round(sqrt(analysis[[i]][[1]]$plot$anava$`Mean Sq`[5])/
                      mean(c(analysis[[i]][[1]]$plot$resp,
                             analysis[[i]][[1]]$plot$respAd))*100,round)}
    rownames(cvs)="CV(%)"
    names(cvs)=variable

    respAd=data.frame(matrix(ncol=length(analysis),nrow = 2))
    variable=1:length(analysis)
    for(i in 1:length(analysis)){
      variable[i]=analysis[[i]][[1]]$plot$ylab
      respAd[1,i]=round(mean(analysis[[i]][[1]]$plot$response),round)
      respAd[2,i]=round(mean(analysis[[i]][[1]]$plot$responseAd),round)}
    rownames(respAd)=c("Factorial","RespAd")
    names(respAd)=variable

    transf=data.frame(matrix(ncol=length(analysis),nrow = 1))
    variable=1:length(analysis)
    for(i in 1:length(analysis)){
      variable[i]=analysis[[i]][[1]]$plot$ylab
      tran=analysis[[i]][[1]]$plot$transf
      if(tran==1){tran=1}else{tran=tran}
      if(tran==0){tran="log"}
      if(tran==0.5){tran="sqrt(x)"}
      if(tran==-0.5){tran="1/sqrt(x)"}
      if(tran==-1){tran="1/x"}
      transf[,i]=ifelse(tran==1,"No transf",tran)}
    rownames(transf)="Transformation"
    names(transf)=variable

    n=5
    nc=4
    infor1=data.frame(matrix(ncol=length(analysis),nrow = 4))
    variable=1:length(analysis)
    for(i in 1:length(analysis)){
      variable[i]=analysis[[i]][[1]]$plot$ylab
      pvalor=round(analysis[[i]][[1]]$plot$anava[1:4,n],round)
      infor1[,i]=ifelse(pvalor<0.001,"p<0.001",pvalor)}
    names(infor1)=variable
    rownames(infor1)=c("p_F1","p_F2","p_F1xF2","p_Fat x Ad")

    n=4
    nc=4
    infor2=data.frame(matrix(ncol=length(analysis),nrow = 4))
    variable=1:length(analysis)
    for(i in 1:length(analysis)){
      variable[i]=analysis[[i]][[1]]$plot$ylab
      infor2[,i]=round(analysis[[i]][[1]]$plot$anava[1:nc,n],round)}
    names(infor2)=variable
    rownames(infor2)=c("f_F1","f_F2","f_F1xF2", "f_Fat x Ad")

    n=3
    nc=5
    infor3=data.frame(matrix(ncol=length(analysis),nrow = 5))
    variable=1:length(analysis)
    for(i in 1:length(analysis)){
      variable[i]=analysis[[i]][[1]]$plot$ylab
      infor3[,i]=round(analysis[[i]][[1]]$plot$anava[1:nc,n],round)}
    names(infor3)=variable
    rownames(infor3)=c("QM_F1","QM_F2","QM_F1xF2","QM_Fat x Ad","QM_r")

    n=2
    nc=5
    infor4=data.frame(matrix(ncol=length(analysis),nrow = 5))
    variable=1:length(analysis)
    for(i in 1:length(analysis)){
      variable[i]=analysis[[i]][[1]]$plot$ylab
      infor4[,i]=round(analysis[[i]][[1]]$plot$anava[1:nc,n],round)}
    names(infor4)=variable
    rownames(infor4)=c("SQ_F1","SQ_F2","SQ_F1xF2","SQ_Fat x Ad","SQ_r")
    variable=1:length(analysis)
    nf1=analysis[[1]][[1]]$plot$anava[1,1]+1
    nf2=analysis[[1]][[1]]$plot$anava[2,1]+1
    nf1f2=nf1*nf2
    f1mean=data.frame(matrix(rep("-",(length(analysis)*nf1)),ncol=length(analysis),nrow = nf1))
    f2mean=data.frame(matrix(rep("-",(length(analysis)*nf2)),ncol=length(analysis),nrow = nf2))
    f1f2mean=data.frame(matrix(rep("-",(length(analysis)*nf1f2)),ncol=length(analysis),nrow = nf1f2))
    rownames(f1mean)=unique(analysis[[1]][[1]]$plot$f1)
    rownames(f2mean)=unique(analysis[[1]][[1]]$plot$f2)
    nomes=expand.grid(unique(analysis[[1]][[1]]$plot$f2),unique(analysis[[1]][[1]]$plot$f1))
    rownames(f1f2mean)=paste(nomes$Var2,nomes$Var1)
    for(i in 1:length(analysis)){
      sigF=analysis[[i]][[1]]$plot$alpha.f
      variable[i]=analysis[[i]][[1]]$plot$ylab
      pvalue=round(analysis[[i]][[1]]$plot$anava[1:3,5],round)
      if(pvalue[3]<sigF){f1f2mean[,i]=data.frame(analysis[[i]][[1]]$plot$graph$numero)
      rownames(f1f2mean)=rownames(analysis[[i]][[1]]$plot$graph)}else{
        if(pvalue[1]<sigF){f1mean[,i]=data.frame(analysis[[i]][[2]]$data$letra)
        rownames(f1mean)=rownames(analysis[[i]][[2]]$data)}
        if(pvalue[2]<sigF){f2mean[,i]=data.frame(analysis[[i]][[3]]$data$letra)
        rownames(f2mean)=rownames(analysis[[i]][[3]]$data)}}
    }
    names(f1mean)=variable
    names(f2mean)=variable
    names(f1f2mean)=variable
    if(inf=="p"){juntos=rbind(f1mean,f2mean,f1f2mean,cvs,infor1,respAd,transf)}
    if(inf=="f"){juntos=rbind(f1mean,f2mean,f1f2mean,cvs,infor2,respAd,transf)}
    if(inf=="QM"){juntos=rbind(f1mean,f2mean,f1f2mean,cvs,infor3,respAd,transf)}
    if(inf=="SQ"){juntos=rbind(f1mean,f2mean,f1f2mean,cvs,infor4,respAd,transf)}
    if(inf=="all"){juntos=rbind(f1mean,f2mean,f1f2mean,cvs,infor1,infor2,infor3,infor4,respAd,transf)}

    if(divisor==TRUE){
      nl=nrow(juntos)
      nc=ncol(juntos)
      market=data.frame(matrix(rep("|",nl*nc),ncol=nc,nrow = nl))
      juntosnovo=cbind("|"=rep("|",nl),juntos,market)
      for(i in 1:nc){
        ordem=matrix(1:(nc*2),nrow=2)
        nomes=colnames(juntos)
        juntosnovo[,c(ordem[,i]+1)]=cbind(juntos[,i],market[,i])
        colnames(juntosnovo)[1]="|"
        colnames(juntosnovo)[c(ordem[,i]+1)]=c(nomes[i],"|")}
      juntos=juntosnovo}
  }
  if(design=="dunnett"){
    nlinhas=length(analysis[[1]]$plot$teste[,1])
    trats=rownames(analysis[[1]]$plot$data)
    infor=data.frame(matrix(ncol=length(analysis),nrow = nlinhas))
    variable=1:length(analysis)
    for(i in 1:length(analysis)){
      variable[i]=analysis[[i]]$plot$label
      infor[i]=analysis[[i]]$plot$teste[6]}
    names(infor)=variable
    rownames(infor)=trats
    juntos=infor}

  # print(juntos)
  list(juntos)[[1]]
}
