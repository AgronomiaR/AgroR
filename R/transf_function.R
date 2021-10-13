#' Utils: Data transformation (Box-Cox, 1964)
#'
#' @description Estimates the lambda value for data transformation
#' @author Gabriel Danilo Shimizu, \email{shimizu@uel.br}
#' @author Leandro Simoes Azeredo Goncalves
#' @author Rodrigo Yudi Palhaci Marubayashi
#' @param response Numerical vector containing the response of the experiment.
#' @param f1 Numeric or complex vector with factor 1 levels
#' @param f2 Numeric or complex vector with factor 2 levels
#' @param f3 Numeric or complex vector with factor 3 levels
#' @param block Numerical or complex vector with blocks
#' @param line Numerical or complex vector with lines
#' @param column Numerical or complex vector with columns
#' @keywords Transformation
#' @keywords Experimental
#' @export
#' @return Returns the value of lambda and/or data transformation approximation, according to Box-Cox (1964)
#' @references
#'
#' Box, G. E., Cox, D. R. (1964). An analysis of transformations. Journal of the Royal Statistical Society: Series B (Methodological), 26(2), 211-243.
#' @examples
#' data("pomegranate")
#' with(pomegranate, transf(WL,trat))

transf=function(response,
                f1,
                f2=NA,
                f3=NA,
                block=NA,
                line=NA,
                column=NA){
  # DIC simples
  requireNamespace("MASS")

  if(is.na(f2)==TRUE && is.na(f3)==TRUE && is.na(block)==TRUE &&
     is.na(line)==TRUE && is.na(column)==TRUE){
    vero=MASS::boxcox(response~f1)}

  # DBC simples
  if(is.na(f2)==TRUE && is.na(f3)==TRUE && is.na(block)==FALSE &&
     is.na(line)==TRUE && is.na(column)==TRUE){
    vero=MASS::boxcox(response~f1+block)}

  # DQL
  if(is.na(f2)==TRUE && is.na(f3)==TRUE && is.na(block)==TRUE &&
     is.na(line)==FALSE && is.na(column)==FALSE){
    vero=MASS::boxcox(response~f1+column+line)}

  # fat2.dic
  if(is.na(f2)==FALSE && is.na(f3)==TRUE && is.na(block)==TRUE &&
     is.na(line)==TRUE && is.na(column)==TRUE){
    vero=MASS::boxcox(response~f1*f2)}

  #fat2dbc
  if(is.na(f2)==FALSE && is.na(f3)==TRUE && is.na(block)==FALSE &&
     is.na(line)==TRUE && is.na(column)==TRUE){
    vero=MASS::boxcox(response~f1*f2+block)}

  # fat3dic
  if(is.na(f2)==FALSE && is.na(f3)==FALSE && is.na(block)==TRUE &&
     is.na(line)==TRUE && is.na(column)==TRUE){
    vero=MASS::boxcox(response~f1*f2*f3)}

  # fat3dic
  if(is.na(f2)==FALSE && is.na(f3)==FALSE && is.na(block)==FALSE &&
     is.na(line)==TRUE && is.na(column)==TRUE){
    vero=MASS::boxcox(response~f1*f2*f3+block)}

  maxvero = vero$x[which.max(vero$y)]
  cat("\n------------------------------------------------\n")
  cat("Box-Cox Transformation (1964)\n")
  cat("\nLambda=")
  cat(lambda=maxvero)
  cat("\n------------------------------------------------\n")
  cat("Suggestion:\n")
  cat(if(round(maxvero,2)>0.75 & round(maxvero,2)<1.25){"Do not transform data"}
      else{if(round(maxvero,2)>0.25 & round(maxvero,2)<=0.75){"square root (sqrt(Y))"}
        else{if(round(maxvero,2)>-0.75 & round(maxvero,2)<=-0.25){"1/sqrt(Y)"}
          else{if(round(maxvero,2)>-1.25 & round(maxvero,2)<=-0.75){"1/Y"}
            else{if(round(maxvero,2)>=-0.25 & round(maxvero,2)<0.25){"log(Y)"}}}}}
  )
  cat("\n\n")
  cat("ou:")
  cat("Yt=(Y^lambda-1)/lambda")
}
