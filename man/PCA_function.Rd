% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pca_function.R
\name{PCA_function}
\alias{PCA_function}
\title{Analysis: Principal components analysis}
\usage{
PCA_function(
  data,
  scale = TRUE,
  text = TRUE,
  pointsize = 5,
  textsize = 12,
  labelsize = 4,
  linesize = 0.6,
  repel = TRUE,
  ylab = NA,
  xlab = NA,
  groups = NA,
  sc = 1,
  font.family = "sans",
  theme = theme_bw(),
  label.legend = "Cluster",
  type.graph = "biplot"
)
}
\arguments{
\item{data}{Data.frame with data set. Line name must indicate the treatment}

\item{scale}{Performs data standardization (\emph{default} is TRUE)}

\item{text}{Add label (\emph{default} is TRUE)}

\item{pointsize}{Point size (\emph{default} is 5)}

\item{textsize}{Text size (\emph{default} is 12)}

\item{labelsize}{Label size (\emph{default} is 4)}

\item{linesize}{Line size (\emph{default} is 0.8)}

\item{repel}{Avoid text overlay (\emph{default} is TRUE)}

\item{ylab}{Names y-axis}

\item{xlab}{Names x-axis}

\item{groups}{Define grouping}

\item{sc}{Secondary axis scale ratio (\emph{default} is 1)}

\item{font.family}{Font family (\emph{default} is sans)}

\item{theme}{Theme ggplot2 (\emph{default} is theme_bw())}

\item{label.legend}{Legend title (when group is not NA)}

\item{type.graph}{Type of chart (\emph{default} is biplot)}
}
\value{
The eigenvalues and eigenvectors, the explanation percentages of each principal component, the correlations between the vectors with the principal components, as well as graphs are returned.
}
\description{
This function performs principal component analysis.
}
\details{
The type.graph argument defines the graph that will be returned,
in the case of "biplot" the biplot graph is returned with the
first two main components and with eigenvalues and eigenvectors.
In the case of "scores" only the treatment scores are returned,
while for "cor" the correlations are returned. For "corPCA" a
correlation between the vectors with the components is returned.
}
\examples{
data(pomegranate)
medias=tabledesc(pomegranate)
PCA_function(medias)
}
\author{
Gabriel Danilo Shimizu
}
