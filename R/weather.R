#' Dataset: Weather data
#'
#' @docType data
#'
#' @description Climatic data from 01 November 2019 to 30 June 2020 in the municipality of Londrina-PR, Brazil. Data from the Instituto de Desenvolvimento Rural do Parana (IDR-PR)
#'
#' @usage data(weather)
#'
#' @format data.frame containing data set
#'   \describe{
#'   \item{\code{Data}}{POSIXct vector with dates}
#'   \item{\code{tempo}}{Numeric vector with time}
#'   \item{\code{Tmax}}{Numeric vector with maximum temperature}
#'   \item{\code{Tmed}}{Numeric vector with mean temperature}
#'   \item{\code{Tmin}}{Numeric vector with minimum temperature}
#'   \item{\code{UR}}{Numeric vector with relative humidity}
#'   }
#' @keywords datasets
#' @seealso \link{cloro}, \link{enxofre}, \link{laranja}, \link{mirtilo}, \link{pomegranate}, \link{porco}, \link{sensorial}, \link{simulate1}, \link{simulate2}, \link{simulate3}, \link{tomate}, \link{aristolochia}, \link{phao}, \link{passiflora}
#' @examples
#' data(weather)
"weather"
