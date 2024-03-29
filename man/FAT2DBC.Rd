% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/FAT2DBC_function.R
\name{FAT2DBC}
\alias{FAT2DBC}
\title{Analysis: DBC experiments in double factorial}
\usage{
FAT2DBC(
  f1,
  f2,
  block,
  response,
  norm = "sw",
  homog = "bt",
  alpha.f = 0.05,
  alpha.t = 0.05,
  quali = c(TRUE, TRUE),
  mcomp = "tukey",
  grau = NA,
  transf = 1,
  constant = 0,
  geom = "bar",
  theme = theme_classic(),
  ylab = "Response",
  xlab = "",
  legend = "Legend",
  fill = "lightblue",
  angle = 0,
  textsize = 12,
  dec = 3,
  family = "sans",
  point = "mean_sd",
  addmean = TRUE,
  errorbar = TRUE,
  CV = TRUE,
  sup = NA,
  color = "rainbow",
  posi = "right",
  ylim = NA,
  angle.label = 0
)
}
\arguments{
\item{f1}{Numeric or complex vector with factor 1 levels}

\item{f2}{Numeric or complex vector with factor 2 levels}

\item{block}{Numerical or complex vector with blocks}

\item{response}{Numerical vector containing the response of the experiment.}

\item{norm}{Error normality test (\emph{default} is Shapiro-Wilk)}

\item{homog}{Homogeneity test of variances (\emph{default} is Bartlett)}

\item{alpha.f}{Level of significance of the F test (\emph{default} is 0.05)}

\item{alpha.t}{Significance level of the multiple comparison test (\emph{default} is 0.05)}

\item{quali}{Defines whether the factor is quantitative or qualitative (\emph{qualitative})}

\item{mcomp}{Multiple comparison test (Tukey (\emph{default}), LSD, Scott-Knott and Duncan)}

\item{grau}{Degree of polynomial in case of quantitative factor (\emph{default} is 1)}

\item{transf}{Applies data transformation (default is 1; for log consider 0)}

\item{constant}{Add a constant for transformation (enter value)}

\item{geom}{Graph type (columns or segments (For simple effect only))}

\item{theme}{ggplot2 theme (\emph{default} is theme_classic())}

\item{ylab}{Variable response name (Accepts the \emph{expression}() function)}

\item{xlab}{Treatments name (Accepts the \emph{expression}() function)}

\item{legend}{Legend title name}

\item{fill}{Defines chart color (to generate different colors for different treatments, define fill = "trat")}

\item{angle}{x-axis scale text rotation}

\item{textsize}{font size}

\item{dec}{number of cells}

\item{family}{font family}

\item{point}{if quali=FALSE, defines whether to plot all points ("all"), mean ("mean"), standard deviation ("mean_sd" - \emph{default}) or mean with standard error (\emph{default} - "mean_se").}

\item{addmean}{Plot the average value on the graph (\emph{default} is TRUE)}

\item{errorbar}{Plot the standard deviation bar on the graph (In the case of a segment and column graph) - \emph{default} is TRUE}

\item{CV}{Plotting the coefficient of variation and p-value of Anova (\emph{default} is TRUE)}

\item{sup}{Number of units above the standard deviation or average bar on the graph}

\item{color}{Column chart color (\emph{default} is "rainbow")}

\item{posi}{Legend position}

\item{ylim}{y-axis scale}

\item{angle.label}{label angle}
}
\value{
The table of analysis of variance, the test of normality of errors (Shapiro-Wilk, Lilliefors, Anderson-Darling, Cramer-von Mises, Pearson and Shapiro-Francia), the test of homogeneity of variances (Bartlett or Levene), the test of independence of Durbin-Watson errors, the test of multiple comparisons (Tukey, LSD, Scott-Knott or Duncan) or adjustment of regression models up to grade 3 polynomial, in the case of quantitative treatments. The column chart for qualitative treatments is also returned.
}
\description{
Analysis of an experiment conducted in a randomized block design in a double factorial scheme using analysis of variance of fixed effects.
}
\note{
The ordering of the graph is according to the sequence in which the factor levels are arranged in the data sheet. The bars of the column and segment graphs are standard deviation.

The function does not perform multiple regression in the case of two quantitative factors.

In the final output when transformation (transf argument) is different from 1, the columns resp and respo in the mean test are returned, indicating transformed and non-transformed mean, respectively.
}
\examples{

#================================================
# Example cloro
#================================================
library(AgroR)
data(cloro)
attach(cloro)
FAT2DBC(f1, f2, bloco, resp, ylab="Number of nodules", legend = "Stages")
FAT2DBC(f1, f2, bloco, resp, mcomp="sk", ylab="Number of nodules", legend = "Stages")
#================================================
# Example covercrops
#================================================
library(AgroR)
data(covercrops)
attach(covercrops)
FAT2DBC(A, B, Bloco, Resp, ylab=expression("Yield"~(Kg~"100 m"^2)),
legend = "Cover crops")
FAT2DBC(A, B, Bloco, Resp, mcomp="sk", ylab=expression("Yield"~(Kg~"100 m"^2)),
legend = "Cover crops")
}
\references{
Principles and procedures of statistics a biometrical approach Steel, Torry and Dickey. Third Edition 1997

Multiple comparisons theory and methods. Departament of statistics the Ohio State University. USA, 1996. Jason C. Hsu. Chapman Hall/CRC.

Practical Nonparametrics Statistics. W.J. Conover, 1999

Ramalho M.A.P., Ferreira D.F., Oliveira A.C. 2000. Experimentacao em Genetica e Melhoramento de Plantas. Editora UFLA.

Scott R.J., Knott M. 1974. A cluster analysis method for grouping mans in the analysis of variance. Biometrics, 30, 507-512.

Mendiburu, F., and de Mendiburu, M. F. (2019). Package ‘agricolae’. R Package, Version, 1-2.
}
\seealso{
\link{FAT2DBC.ad}
}
\author{
Gabriel Danilo Shimizu, \email{shimizu@uel.br}

Leandro Simoes Azeredo Goncalves

Rodrigo Yudi Palhaci Marubayashi
}
\keyword{DBC}
\keyword{Factorial}
