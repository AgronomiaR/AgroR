% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plotcor_function.R
\name{plot_cor}
\alias{plot_cor}
\title{Graph: Plot correlation}
\usage{
plot_cor(
  x,
  y,
  method = "pearson",
  ylab = "Dependent",
  xlab = "Independent",
  theme = theme_classic(),
  pointsize = 5,
  shape = 21,
  fill = "gray",
  color = "black",
  axis.size = 12,
  ic = TRUE,
  title = NA,
  family = "sans"
)
}
\arguments{
\item{x}{Numeric vector with independent variable}

\item{y}{Numeric vector with dependent variable}

\item{method}{Method correlation (\emph{default} is Pearson)}

\item{ylab}{Variable response name (Accepts the \emph{expression}() function)}

\item{xlab}{Treatments name (Accepts the \emph{expression}() function)}

\item{theme}{ggplot2 theme (\emph{default} is theme_classic())}

\item{pointsize}{Point size}

\item{shape}{shape format}

\item{fill}{Fill point}

\item{color}{Color point}

\item{axis.size}{Axis text size}

\item{ic}{add interval of confidence}

\item{title}{title}

\item{family}{Font family}
}
\value{
The function returns a graph for correlation
}
\description{
Correlation analysis function (Pearson or Spearman)
}
\examples{
data("pomegranate")
with(pomegranate, plot_cor(WL, SS, xlab="WL", ylab="SS"))
}
\author{
Gabriel Danilo Shimizu, \email{shimizu@uel.br}

Leandro Simoes Azeredo Goncalves

Rodrigo Yudi Palhaci Marubayashi
}
