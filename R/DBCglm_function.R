#' Analysis: Randomized block design by glm
#' @description Statistical analysis of experiments conducted in a randomized block design using a generalized linear model. It performs the deviance analysis and the effect is tested by a chi-square test. Multiple comparisons are adjusted by Tukey.
#' @author Gabriel Danilo Shimizu, \email{shimizu@uel.br}
#' @author Leandro Simoes Azeredo Goncalves
#' @author Rodrigo Yudi Palhaci Marubayashi
#' @param trat Numerical or complex vector with treatments
#' @param block Numerical or complex vector with blocks
#' @param response Numerical vector containing the response of the experiment. Use cbind(resp, n-resp) for binomial or quasibinomial family.
#' @param glm.family distribution family considered (\emph{default} is binomial)
#' @param quali Defines whether the factor is quantitative or qualitative (\emph{default} is qualitative)
#' @param alpha.f Level of significance of the F test (\emph{default} is 0.05)
#' @param alpha.t Significance level of the multiple comparison test (\emph{default} is 0.05)
#' @param geom Graph type (columns, boxes or segments)
#' @param theme ggplot2 theme (\emph{default} is theme_classic())
#' @param sup Number of units above the standard deviation or average bar on the graph
#' @param ylab Variable response name (Accepts the \emph{expression}() function)
#' @param xlab Treatments name (Accepts the \emph{expression}() function)
#' @param textsize Font size
#' @param labelsize Label size
#' @param fill Defines chart color (to generate different colors for different treatments, define fill = "trat")
#' @param angle x-axis scale text rotation
#' @param family Font family
#' @param dec Number of cells
#' @param addmean Plot the average value on the graph (\emph{default} is TRUE)
#' @param errorbar Plot the standard deviation bar on the graph (In the case of a segment and column graph) - \emph{default} is TRUE
#' @param posi Legend position
#' @param point Defines whether to plot mean ("mean"), mean with standard deviation ("mean_sd" - \emph{default}) or mean with standard error (\emph{default} - "mean_se").
#' @param angle.label label angle
#' @export
#' @examples
#' data("aristolochia")
#' attach(aristolochia)
#' # Assuming the same aristolochia data set, but considering randomized blocks
#' bloco=rep(paste("B",1:16),5)
#' resp=resp/2
#' DBC.glm(trat,bloco, cbind(resp,50-resp), glm.family="binomial")

DBC.glm=function(trat,
                 block,
                 response,
                 glm.family="binomial",
                 quali=TRUE,
                 alpha.f=0.05,
                 alpha.t=0.05,
                 geom="bar",
                 theme=theme_classic(),
                 sup=NA,
                 ylab="Response",
                 xlab="",
                 fill="lightblue",
                 angle=0,
                 family="sans",
                 textsize=12,
                 labelsize=5,
                 dec=3,
                 addmean=TRUE,
                 errorbar=TRUE,
                 posi="top",
                 point="mean_sd",
                 angle.label=0){
  if(angle.label==0){hjust=0.5}else{hjust=0}
  requireNamespace("emmeans")
  # requireNamespace("stringr")
  requireNamespace("multcomp")
  requireNamespace("crayon")
  requireNamespace("ggplot2")
  trat=as.factor(trat)
  block=as.factor(block)
  resp=response
  if(glm.family=="binomial" | glm.family=="quasibinomial"){
    if(glm.family=="binomial"){a = glm(resp ~ trat+block,family = "binomial")}
    if(glm.family=="quasibinomial"){a = glm(resp ~ trat+block,family = "quasibinomial")}
    anava1 = anova(a, test="Chisq")
    anava2 = summary(a)
    anava=rbind("Null deviance" = anava1$`Resid. Dev`[1],
                "Df Null deviance" = round(anava1$`Resid. Df`[1]),
                "-----"=NA,
                "Treatment effects"=NA,
                "Residual deviance" = anava1$`Resid. Dev`[2],
                "Df residual deviance" = round(anava1$`Resid. Df`[2],0),
                "p-value(Chisq)" = anava1$`Pr(>Chi)`[2],
                "-----"=NA,
                "Block effects"=NA,
                "Residual deviance" = anava1$`Resid. Dev`[3],
                "Df residual deviance" = round(anava1$`Resid. Df`[3],0),
                "p-value(Chisq)" = anava1$`Pr(>Chi)`[3],
                "-----"=NA,
                AIC = anava2$aic)
    colnames(anava)=""

    if(is.na(sup==TRUE)){sup=0.1*mean(a$fitted.values*100)}
    cat(green(bold("\n\n-----------------------------------------------------------------\n")))
    cat(green(bold("Analysis of deviance")))
    cat(green(bold("\n-----------------------------------------------------------------\n")))
    print(as.matrix(round(anava,3)),na.print = "",quote=F)
    cat("\n\n")
    message(if (anava1$`Pr(>Chi)`[2]<alpha.f){
      black("As the calculated p-value, it is less than the 5% significance level.The hypothesis H0 of equality of means is rejected. Therefore, at least two treatments differ")}
      else {"As the calculated p-value is greater than the 5% significance level, H0 is not rejected"})
    cat(green(bold("\n\n-----------------------------------------------------------------\n")))
    if(quali==TRUE){cat(green(bold("Multiple Comparison Test")))}else{cat(green(bold("Regression")))}
    cat(green(bold("\n-----------------------------------------------------------------\n")))
    if(quali==TRUE){
      letra <- cld(regrid(emmeans(a, "trat", alpha=alpha.t)),
                   Letters=letters,reversed=TRUE,adjusted="tukey")
      rownames(letra)=letra$trat
      letra=letra[unique(as.character(trat)),]
      out=letra[,-4]
      out[,c(2,3,4,5)]=round(out[,c(2,3,4,5)],2)
      print(out)

      prob=letra$prob*100
      superior=letra$asymp.UCL*100
      inferior=letra$asymp.LCL*100
      # grupo=str_trim(letra$.group)
      grupo=gsub(" ","",letra$.group)
      trat=letra$trat

      desvio=letra$asymp.UCL-letra$prob
      media=prob
      trats=trat
      dadosm=data.frame(trat,prob,superior,inferior,grupo,desvio,
                        media,trats)
      if(addmean==TRUE){dadosm$letra=paste(format(round(prob,3),digits = dec),
                                           grupo)}
      if(addmean==FALSE){dadosm$letra=grupo}
      letra=dadosm$letra
      if(geom=="bar"){
        grafico=ggplot(dadosm,aes(x=trat,y=prob))
        if(fill=="trat"){grafico=grafico+
          geom_col(aes(fill=trat),color=1)}else{grafico=grafico+
            geom_col(aes(fill=trat),fill=fill,color=1)}
        if(errorbar==TRUE){grafico=grafico+
          geom_text(aes(y=superior+sup,label=letra),
                    family=family,size=labelsize,
                    angle=angle.label,
                    hjust=hjust)}
        if(errorbar==FALSE){grafico=grafico+geom_text(aes(y=prob+sup,label=letra),family=family,
                                                      angle=angle.label, size=labelsize,hjust=hjust)}
        if(errorbar==TRUE){grafico=grafico+
          geom_errorbar(data=dadosm,aes(ymin=inferior, ymax=superior,color=1),
                        color="black",width=0.3)}}
      if(geom=="point"){grafico=ggplot(dadosm,aes(x=trat,
                                                  y=prob))
      if(errorbar==TRUE){grafico=grafico+
        geom_text(aes(y=superior+sup, label=letra),
                  family=family,angle=angle.label,size=labelsize, hjust=hjust)}
      if(errorbar==FALSE){grafico=grafico+
        geom_text(aes(y=prob+sup,label=letra),size=labelsize,family=family,angle=angle.label, hjust=hjust)}
      if(errorbar==TRUE){grafico=grafico+
        geom_errorbar(data=dadosm,
                      aes(ymin=inferior,
                          ymax=superior,color=1),
                      color="black",width=0.3)}
      if(fill=="trat"){grafico=grafico+
        geom_point(aes(color=trat),size=5)}
      else{grafico=grafico+
        geom_point(aes(color=trat),
                   color="black",
                   fill=fill,shape=21,size=5)}}
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
      grafico=as.list(grafico)
      print(grafico)
    }
  }
  if(glm.family=="poisson" | glm.family=="quasipoisson"){
    if(glm.family=="poisson"){a = glm(resp ~ trat+block,family = "poisson")}
    if(glm.family=="quasipoisson"){a = glm(resp ~ trat+block,family = "quasipoisson")}
    anava1 = anova(a, test="Chisq")
    anava2 = summary(a)
    anava=rbind("Null deviance" = anava1$`Resid. Dev`[1],
                "Df Null deviance" = round(anava1$`Resid. Df`[1]),
                "-----"=NA,
                "Treatment effects"=NA,
                "Residual deviance" = anava1$`Resid. Dev`[2],
                "Df residual deviance" = round(anava1$`Resid. Df`[2]),
                "p-value(Chisq)" = anava1$`Pr(>Chi)`[2],
                "-----"=NA,
                "Block effects"=NA,
                "Residual deviance" = anava1$`Resid. Dev`[3],
                "Df residual deviance" = round(anava1$`Resid. Df`[3]),
                "p-value(Chisq)" = anava1$`Pr(>Chi)`[3],
                "-----"=NA,
                AIC = anava2$aic)
    colnames(anava)=""
    if(is.na(sup==TRUE)){sup=0.1*mean(a$fitted.values)}
    cat(green(bold("\n\n-----------------------------------------------------------------\n")))
    cat(green(bold("Analysis of deviance")))
    cat(green(bold("\n-----------------------------------------------------------------\n")))
    print(as.matrix(round(anava,3)),na.print = "",quote=F)
    cat("\n\n")
    message(if (anava1$`Pr(>Chi)`[2]<alpha.f){
      black("As the calculated p-value, it is less than the 5% significance level.The hypothesis H0 of equality of means is rejected. Therefore, at least two treatments differ")}
      else {"As the calculated p-value is greater than the 5% significance level, H0 is not rejected"})
    cat(green(bold("\n\n-----------------------------------------------------------------\n")))
    if(quali==TRUE){cat(green(bold("Multiple Comparison Test")))}else{cat(green(bold("Regression")))}
    cat(green(bold("\n-----------------------------------------------------------------\n")))
    if(quali==TRUE){
      letra <- cld(regrid(emmeans(a, "trat", alpha=alpha.t,adjusted="tukey")),
                   Letters=letters,reversed=TRUE)
      rownames(letra)=letra$trat
      letra=letra[unique(as.character(trat)),]
      print(letra[,-4])
      rate=letra$rate
      superior=letra$asymp.UCL
      inferior=letra$asymp.LCL
      # grupo=str_trim(letra$.group)
      grupo=gsub(" ","",letra$.group)
      trat=letra$trat
      desvio=letra$asymp.UCL-letra$rate
      media=rate
      trats=trat
      dadosm=data.frame(trat,rate,superior,inferior,
                        grupo,desvio,
                        media,trats)
      if(addmean==TRUE){dadosm$letra=paste(format(rate,digits = dec),
                                           grupo)}
      if(addmean==FALSE){dadosm$letra=grupo}
      letra=dadosm$letra
      if(geom=="bar"){
        grafico=ggplot(dadosm,aes(x=trat,y=rate))
        if(fill=="trat"){grafico=grafico+
          geom_col(aes(fill=trat),color=1)}else{grafico=grafico+
            geom_col(aes(fill=trat),fill=fill,color=1)}
        if(errorbar==TRUE){grafico=grafico+
          geom_text(aes(y=superior+sup,label=letra),
                    family=family,size=labelsize,
                    angle=angle.label,
                    hjust=hjust)}
        if(errorbar==FALSE){grafico=grafico+geom_text(aes(y=rate+sup,label=letra),family=family,
                                                      angle=angle.label,size=labelsize, hjust=hjust)}
        if(errorbar==TRUE){grafico=grafico+
          geom_errorbar(data=dadosm,aes(ymin=inferior, ymax=superior,color=1),
                        color="black",width=0.3)}}
      if(geom=="point"){grafico=ggplot(dadosm,aes(x=trat,
                                                  y=rate))
      if(errorbar==TRUE){grafico=grafico+
        geom_text(aes(y=superior+sup, label=letra),
                  family=family,size=labelsize,angle=angle.label, hjust=hjust)}
      if(errorbar==FALSE){grafico=grafico+
        geom_text(aes(y=rate+sup,label=letra),family=family,size=labelsize,angle=angle.label, hjust=hjust)}
      if(errorbar==TRUE){grafico=grafico+
        geom_errorbar(data=dadosm,
                      aes(ymin=inferior,
                          ymax=superior,color=1),
                      color="black",width=0.3)}
      if(fill=="trat"){grafico=grafico+
        geom_point(aes(color=trat),size=5)}
      else{grafico=grafico+
        geom_point(aes(color=trat),
                   color="black",
                   fill=fill,shape=21,size=5)}}
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
      grafico=as.list(grafico)
      print(grafico)
    }
  }

}
