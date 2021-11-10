#' Analysis: Test for two samples
#' @author Gabriel Danilo Shimizu, \email{shimizu@uel.br}
#' @author Leandro Simoes Azeredo Goncalves
#' @author Rodrigo Yudi Palhaci Marubayashi
#' @description Test for two samples (paired and unpaired t test, paired and unpaired Wilcoxon test)
#' @param trat Categorical vector with the two treatments
#' @param resp Numeric vector with the response
#' @param test Test used (t for test t or w for Wilcoxon test)
#' @param alternative A character string specifying the alternative hypothesis, must be one of "two.sided" (default), "greater" or "less". You can specify just the initial letter.
#' @param paired A logical indicating whether you want a paired t-test.
#' @param var.equal A logical variable indicating whether to treat the two variances as being equal. If TRUE then the pooled variance is used to estimate the variance otherwise the Welch (or Satterthwaite) approximation to the degrees of freedom is used.
#' @param conf.level Confidence level of the interval.
#' @param theme ggplot2 theme (\emph{default} is theme_classic())
#' @param ylab Variable response name (Accepts the \emph{expression}() function)
#' @param xlab Treatments name (Accepts the \emph{expression}() function)
#' @param pointsize Point size
#' @param fill fill box
#' @details Alternative = "greater" is the alternative that x has a larger mean than y. For the one-sample case: that the mean is positive.
#' @details If paired is TRUE then both x and y must be specified and they must be the same length. Missing values are silently removed (in pairs if paired is TRUE). If var.equal is TRUE then the pooled estimate of the variance is used. By default, if var.equal is FALSE then the variance is estimated separately for both groups and the Welch modification to the degrees of freedom is used.
#' @details If the input data are effectively constant (compared to the larger of the two means) an error is generated.
#' @return Returns the test for two samples (paired or unpaired t test, paired or unpaired Wilcoxon test)
#' @export
#' @examples
#' resp=rnorm(100,100,5)
#' trat=rep(c("A","B"),e=50)
#' test_two(trat,resp)
#' test_two(trat,resp,paired = TRUE)

test_two=function(trat,
                resp,
                paired=FALSE,
                test="t",
                alternative = c("two.sided", "less", "greater"),
                conf.level=0.95,
                theme=theme_classic(),
                ylab="Response",
                xlab="",
                var.equal=FALSE,
                pointsize=2,
                fill="white"){
  if(test=="t"){teste=t.test(resp~trat,
                paired=paired,
                alternative=alternative,
                conf.level=conf.level,
                var.equal=var.equal)}
  if(test=="w"){teste=wilcox.test(resp~trat,
                              paired=paired,
                              alternative=alternative,
                              conf.level=conf.level,
                              var.equal=var.equal,
                              exact=FALSE)}
  dados=data.frame(resp,trat)
  media=data.frame(media=tapply(resp,trat,mean, na.rm=TRUE))
  media$trat=rownames(media)
  pvalor=sprintf("italic(\"p-value\")==%0.3f",teste$p.value)
  requireNamespace("ggplot2")
  grafico=ggplot(dados,aes(y=resp,x=trat))+
  geom_boxplot(size=0.8,outlier.colour = "white", fill=fill)+
  geom_jitter(width=0.1,alpha=0.2,size=pointsize)+
  stat_boxplot(geom="errorbar", size=0.8,width=0.2)+
  geom_label(data=media,aes(y=media,
                            x=trat,
                            label=round(media,2)),fill="lightyellow")+
  theme+ylab(ylab)+xlab(xlab)+
  theme(axis.text = element_text(size=12,color="black"))
  if(teste$p.value<0.001){grafico=grafico+annotate(geom="text",
           x=1.5,y=min(resp),hjust=0.5,vjust=0,
           label="italic(\"p-value <0.001\")",parse=TRUE)}
  if(teste$p.value>0.001){grafico=grafico+
    annotate(geom="text",x=1.5,y=min(resp),hjust=0.5,vjust=0,
             label=pvalor,parse=TRUE)}
  print(teste)
  print(grafico)
  graficos=list(grafico)[[1]]
  }
