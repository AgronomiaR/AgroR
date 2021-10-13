#' Dataset: Pig development and production
#'
#' @description An experiment whose objective was to study the effect of castration
#' age on the development and production of pigs, evaluating the weight
#' of the piglets. Four treatments were studied: A - castration at
#' 56 days of age; B - castration at 7 days of age; C - castration at
#' 36 days of age; D - whole (not castrated); E - castration at 21 days
#' of age. The Latin square design was used in order to control the
#' variation between litters (lines) and the variation in the initial
#' weight of the piglets (columns), with the experimental portion
#' consisting of a piglet.
#'
#' @docType data
#'
#' @usage data(porco)
#'
#' @keywords datasets
#' @format data.frame containing data set
#'   \describe{
#'   \item{\code{trat}}{Categorical vector with treatments}
#'   \item{\code{linhas}}{Categorical vector with lines}
#'   \item{\code{colunas}}{Categorical vector with columns}
#'   \item{\code{resp}}{Numeric vector}
#'   }
#' @seealso \link{cloro}, \link{enxofre}, \link{laranja}, \link{mirtilo}, \link{pomegranate}, \link{sensorial}, \link{simulate1}, \link{simulate2}, \link{simulate3}, \link{tomate}, \link{weather}, \link{phao}, \link{passiflora}, \link{aristolochia}
#' @keywords datasets
#'
#' @examples
#' data(porco)
"porco"
