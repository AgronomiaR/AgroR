#' Analysis: Analysis of Variance of Aligned Rank Transformed Data in FAT2DIC
#' @author Gabriel Danilo Shimizu, \email{shimizu@uel.br}
#' @author Leandro Simoes Azeredo Goncalves
#' @author Rodrigo Yudi Palhaci Marubayashi
#' @description Apply the aligned rank transform to a factorial model (with optional grouping terms). Usually done in preparation for a nonparametric analyses of variance on models with numeric or ordinal responses, which can be done by following up with anova.art.
#' @param f1 Numeric or complex vector with factor 1 levels
#' @param f2 Numeric or complex vector with factor 2 levels
#' @param response Numerical vector containing the response of the experiment.
#' @param ylab Variable response name (Accepts the \emph{expression}() function)
#' @param xlab Treatments name (Accepts the \emph{expression}() function)
#' @param theme ggplot2 theme (\emph{default} is theme_classic())
#' @param legend.title Legend title name
#' @param decreasing Letter order (\emph{default} is TRUE)
#' @param sup Number of units above the standard deviation or average bar on the graph
#' @return The function returns the Anova of aligned ranks, the multiple comparison test and the interaction graph.
#' @references
#'
#' Wobbrock, J. O., Findlater, L., Gergle, D., Higgins, J. J. (2011, May). The aligned rank transform for nonparametric factorial analyses using only anova procedures. In Proceedings of the SIGCHI conference on human factors in computing systems (pp. 143-146).
#'
#' Kay, M., Wobbrock, J. O. (2020). Package ‘ARTool’.
#'
#' @seealso \link{FAT2DIC}
#' @import ARTool
#' @export
#' @examples
#' data(cloro)
#' with(cloro, FAT2DIC.art(f1,f2,resp))

FAT2DIC.art=function(f1,
                     f2,
                     response,
                     decreasing=TRUE,
                     sup=NA,
                     xlab=" ",
                     ylab="Sum of posts",
                     legend.title="Factor",
                     theme=theme_classic()){
  requireNamespace("ARTool")
  requireNamespace("ggplot2")
  requireNamespace("crayon")
  requireNamespace("emmeans")
  requireNamespace("multcomp")
  requireNamespace("stringr")
  fator1=f1
  fator2=f2
  fator1=as.factor(fator1)
  fator2=as.factor(fator2)
  data=data.frame(fator1,fator2,response)
  mod=art(response~fator1+fator2+fator1:fator2,data=data)
  anova(mod, response="aligned")
  summary(mod)
  cat(green(bold("\n-----------------------------------------------------------------\n")))
  cat(green(italic("Analysis of Variance of Aligned Rank Transformed Data")))
  cat(green(bold("\n-----------------------------------------------------------------\n")))
  an=data.frame(anova(mod));colnames(an)=c("FV","Df","Df.res","SQ","SQres","Fvalue","p-value")
  print(an)
  a=anova(mod)
  if(a$`Pr(>F)`[3]>0.05 && a$`Pr(>F)`[1]<0.05){mod1=artlm(mod,"fator1")
  cat(green(bold("\n\n-----------------------------------------------------------------\n")))
  a1=cld(emmeans(mod1,~fator1),Letters=letters,reversed=decreasing)
  fator1=a1$fator1
  emmean=a1$emmean
  if(is.na(sup==TRUE)){sup=0.1*mean(emmean)}
  lower.CL=a1$lower.CL
  upper.CL=a1$upper.CL
  .group=a1$.group
  g1=ggplot(a1,aes(x=fator1,
                   y=emmean))+
    geom_col(fill="lightblue",
             color="black")+
    geom_errorbar(aes(ymin=lower.CL,
                      ymax=upper.CL),
                  width=0.2,size=0.8)+
    geom_label(aes(label=.group,
                   y=upper.CL+sup),
               fill="lightyellow",
               show.legend = FALSE)+
    theme+
    theme(axis.text = element_text(size=12,color="black"))+
    labs(y=ylab,x=xlab)
  print(a1)
  print(g1)
  }
  if(a$`Pr(>F)`[3]>0.05 && a$`Pr(>F)`[2]<0.05){mod1=artlm(mod,"fator2")
  cat(green(bold("\n\n-----------------------------------------------------------------\n")))
  a2=cld(emmeans(mod1,~fator2),Letters=letters,reversed=decreasing)
  fator2=a2$fator2
  emmean=a2$emmean
  if(is.na(sup==TRUE)){sup=0.1*mean(emmean)}
  lower.CL=a2$lower.CL
  upper.CL=a2$upper.CL
  .group=a2$.group
  g2=ggplot(a2,aes(x=fator2,y=emmean))+
    geom_col(fill="lightblue",
             color="black")+
    geom_errorbar(aes(ymin=lower.CL,
                      ymax=upper.CL),
                  width=0.2,size=0.8)+
    geom_label(aes(label=.group,
                   y=upper.CL+sup),
               fill="lightyellow",
               show.legend = FALSE)+
    theme+
    theme(axis.text = element_text(size=12,color="black"))+
    labs(y=ylab,x=xlab)
  print(a2)
  print(g2)
  }
  if(a$`Pr(>F)`[3]<0.05){mod1=artlm(mod,"fator1:fator2")
  cat(green(bold("\n\n-----------------------------------------------------------------\n")))
  a=cld(emmeans(mod1,~fator1|fator2),Letters=LETTERS,reversed=decreasing)
  a=data.frame(a)
  print(a)
  cat(green(bold("\n\n-----------------------------------------------------------------\n")))
  b=cld(emmeans(mod1,~fator2|fator1),Letters=letters,reversed=decreasing)
  b=data.frame(b)
  print(b)
  a$trat=paste(a$fator1,a$fator2)
  b$trat=paste(b$fator1,b$fator2)
  tabela=merge(a,b,by.x="trat",by.y="trat")
  tabela=tabela[,c(2,3,4,7,8,9,17)]
  colnames(tabela)=c("fator1","fator2","emmean","lower","upper","group1","group2")
  tabela
  tabela$group1=str_trim(tabela$group1)
  tabela$group2=str_trim(tabela$group2)
  fator2=tabela$fator2
  fator1=tabela$fator1
  emmean=tabela$emmean
  if(is.na(sup==TRUE)){sup=0.1*mean(emmean)}
  lower=tabela$lower
  upper=tabela$upper
  group1=tabela$group1
  group2=tabela$group2
  grafico=ggplot(tabela,
                   aes(x=fator2,
                       y=emmean))+
    geom_col(aes(fill=fator1,
                 color=fator1,
                 group=fator1),
             position = position_dodge(width = 0.9),
             color="black")+
    geom_errorbar(aes(ymin=lower,
                      ymax=upper,
                      group=fator1),
                  position = position_dodge(width = 0.9),
                  width=0.2,size=0.8)+
    geom_label(aes(label=paste(group1,
                               group2,
                               sep=""),
                   y=upper+sup,
                   group=fator1),
               position = position_dodge(width = 0.8),
               fill="lightyellow",show.legend = FALSE)+
    theme+
    theme(axis.text = element_text(size=12,color="black"))+
    labs(y=ylab,x=xlab,fill=legend.title)
  print(grafico)}
}

