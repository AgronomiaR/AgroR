#' Dataset: Cutting blueberry data
#'
#' An experiment was carried out in order to evaluate the rooting
#' (resp1) of blueberry cuttings as a function of the cutting size
#' (Treatment Colume). This experiment was repeated three times
#' (Location column) and a randomized block design with four
#' replications was adopted.
#'
#' @docType data
#'
#' @usage data(mirtilo)
#'
#' @format data.frame containing data set
#'   \describe{
#'   \item{\code{trat}}{Categorical vector with treatments}
#'   \item{\code{exp}}{Categorical vector with experiment}
#'   \item{\code{bloco}}{Categorical vector with block}
#'   \item{\code{resp}}{Numeric vector}
#'   }
#' @keywords datasets
#' @seealso \link{cloro}, \link{enxofre}, \link{laranja}, \link{pomegranate}, \link{porco}, \link{sensorial}, \link{simulate1}, \link{simulate2}, \link{simulate3}, \link{tomate}, \link{weather}
#' @examples
#' data(mirtilo)
#' attach(mirtilo)
"mirtilo"
