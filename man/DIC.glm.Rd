% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/DICglm_function.R
\name{DIC.glm}
\alias{DIC.glm}
\title{Analysis: Completely randomized design by glm}
\usage{
DIC.glm(
  trat,
  response,
  glm.family = "binomial",
  quali = TRUE,
  alpha.f = 0.05,
  alpha.t = 0.05,
  geom = "bar",
  theme = theme_classic(),
  sup = NA,
  ylab = "Response",
  xlab = "",
  fill = "lightblue",
  angle = 0,
  family = "sans",
  textsize = 12,
  labelsize = 5,
  dec = 3,
  addmean = TRUE,
  errorbar = TRUE,
  posi = "top",
  point = "mean_sd",
  angle.label = 0
)
}
\arguments{
\item{trat}{Numerical or complex vector with treatments}

\item{response}{Numerical vector containing the response of the experiment. Use cbind(resp, n-resp) for binomial or quasibinomial family.}

\item{glm.family}{distribution family considered (\emph{default} is binomial)}

\item{quali}{Defines whether the factor is quantitative or qualitative (\emph{default} is qualitative)}

\item{alpha.f}{Level of significance of the F test (\emph{default} is 0.05)}

\item{alpha.t}{Significance level of the multiple comparison test (\emph{default} is 0.05)}

\item{geom}{Graph type (columns, boxes or segments)}

\item{theme}{ggplot2 theme (\emph{default} is theme_classic())}

\item{sup}{Number of units above the standard deviation or average bar on the graph}

\item{ylab}{Variable response name (Accepts the \emph{expression}() function)}

\item{xlab}{Treatments name (Accepts the \emph{expression}() function)}

\item{fill}{Defines chart color (to generate different colors for different treatments, define fill = "trat")}

\item{angle}{x-axis scale text rotation}

\item{family}{Font family}

\item{textsize}{Font size}

\item{labelsize}{Label size}

\item{dec}{Number of cells}

\item{addmean}{Plot the average value on the graph (\emph{default} is TRUE)}

\item{errorbar}{Plot the standard deviation bar on the graph (In the case of a segment and column graph) - \emph{default} is TRUE}

\item{posi}{Legend position}

\item{point}{Defines whether to plot mean ("mean"), mean with standard deviation ("mean_sd" - \emph{default}) or mean with standard error (\emph{default} - "mean_se").}

\item{angle.label}{label angle}
}
\description{
Statistical analysis of experiments conducted in a completely randomized design using a generalized linear model. It performs the deviance analysis and the effect is tested by a chi-square test. Multiple comparisons are adjusted by Tukey.
}
\examples{
data("aristolochia")
attach(aristolochia)
#=============================
# Use the DIC function
#=============================
DIC(trat, resp)

#=============================
# Use the DIC function noparametric
#=============================
DIC(trat, resp, test="noparametric")

#=============================
# Use the DIC.glm function
#=============================

resp=resp/4 # total germinated seeds

# the value 25 is the total of seeds in the repetition
DIC.glm(trat, cbind(resp,25-resp), glm.family="binomial")
}
\author{
Gabriel Danilo Shimizu, \email{shimizu@uel.br}

Leandro Simoes Azeredo Goncalves

Rodrigo Yudi Palhaci Marubayashi
}
