#' Analysis: Dunnett test
#' @export
#' @description The function performs the Dunnett test
#' @param trat Numerical or complex vector with treatments
#' @param resp Numerical vector containing the response of the experiment.
#' @param control Treatment considered control (write identical to the name in the vector)
#' @param model Experimental design (DIC, DBC or DQL)
#' @param block Numerical or complex vector with blocks
#' @param line Numerical or complex vector with lines
#' @param column Numerical or complex vector with columns
#' @param alpha.t Significance level (\emph{default} is 0.05)
#' @param label Variable label
#' @note Do not use the "-" symbol or space in treatment names
#' @return I return the Dunnett test for experiments in a completely randomized design, randomized blocks or Latin square.
#' @importFrom multcomp glht
#' @importFrom multcomp mcp
#' @examples
#'
#' #====================================================
#' # complete randomized design
#' #====================================================
#' data("pomegranate")
#' with(pomegranate,dunnett(trat=trat,resp=WL,control="T1"))
#'
#' #====================================================
#' # randomized block design in factorial double
#' #====================================================
#' library(AgroR)
#' data(cloro)
#' attach(cloro)
#' respAd=c(268, 322, 275, 350, 320)
#' a=FAT2DBC.ad(f1, f2, bloco, resp, respAd,
#'              ylab="Number of nodules",
#'              legend = "Stages",mcomp="sk")
#' data=rbind(data.frame(trat=paste(f1,f2,sep = ""),bloco=bloco,resp=resp),
#'            data.frame(trat=c("Test","Test","Test","Test","Test"),
#'                       bloco=unique(bloco),resp=respAd))
#' with(data,dunnett(trat = trat,
#'                   resp = resp,
#'                   control = "Test",
#'                   block=bloco,model = "DBC"))


dunnett=function(trat, resp, control, model="DIC",
                 block=NA, column=NA, line=NA, alpha.t=0.05,
                 label="Response"){
  trat1=as.factor(trat)
  trat=as.factor(trat)
  levels(trat1)=paste("T",1:length(levels(trat1)),sep = "")
  controle=as.character(trat1[trat==control][1])
  if(model=="DIC"){mod=aov(resp~trat1)}
  if(model=="DBC"){
    block=as.factor(block)
    mod=aov(resp~trat1+block)}
  if(model=="DQL"){
    column=as.factor(column)
    lines=as.factor(lines)
    mod=aov(resp~trat1+column+line)}
  requireNamespace("multcomp")
  dados=data.frame(trat1,resp)
  contras=unique(trat1)[!unique(trat1)==controle]
  a=confint(glht(mod,
               linfct = mcp(trat1=paste(contras,"-",
                                       controle,
                                       "==0",sep=""))),
          level = 1-alpha.t,)
  a=summary(a)
  teste=cbind(a$confint,
        round(a$test$tstat,4),
        round(a$test$pvalues,4))
  nomes=rownames(teste)
  nomes1=as.factor(t(matrix(unlist(strsplit(nomes," - ")),nrow=2))[,1])
  levels(nomes1)=levels(trat)[!levels(trat)==control]
  rownames(teste)=paste(control," - ",nomes1)
  teste=data.frame(teste)
  colnames(teste)=c("Estimate","IC-lwr","IC-upr","t value","p-value")
  teste$sig=ifelse(teste$`p-value`>0.05,"ns",
                   ifelse(teste$`p-value`<0.01,"**","*"))
  print(teste)
  data=data.frame(teste)
  `IC-lwr`=data$IC.lwr
  `IC-upr`=data$IC.upr
  sig=data$sig
  Estimate=data$Estimate
  graph=ggplot(data,aes(y=rownames(data),x=Estimate))+
    geom_errorbar(aes(xmin=`IC-lwr`,xmax=`IC-upr`),width=0.2,size=1)+
    geom_point(shape=21,size=5,color="black",fill="gray")+
    theme_classic()+
    labs(y="")+
    geom_vline(xintercept = 0,lty=2,size=1)+
    geom_label(aes(label=paste(round(Estimate,3),
                               sig)),fill="lightyellow",
               vjust=-0.5)+
    theme(axis.text = element_text(size=12))
  plot(graph)
  }
