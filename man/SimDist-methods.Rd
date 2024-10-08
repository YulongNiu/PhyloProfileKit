% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/AllGenerics.R, R/SimDistBatch.R
\docType{methods}
\name{SimDist}
\alias{SimDist}
\alias{SimDist,PPIdx-method}
\alias{SimDist,PPTreeIdx-method}
\title{Batch process of similarity and distance}
\usage{
SimDist(x, method, ..., n = 1)

\S4method{SimDist}{PPIdx}(x, method, ..., n = 1)

\S4method{SimDist}{PPTreeIdx}(x, method, ..., n = 1)
}
\arguments{
\item{x}{A \code{PPIdx}/\code{PPTreeIdx} object.}

\item{method}{A character string.
\itemize{
  \item \code{"SimCor"}: Person's correlation coefficient.
  \item \code{"SimJaccard"}: Jaccard similarity.
  \item \code{"SimMI"}: Mutual information.
  \item \code{"DistHamming"}: Hamming distance.
  \item \code{"DistEuclidean"}: Euclidean distance.
}

\code{"SimCor"} and \code{"DistHamming"} are not recommended for binning profiles.}

\item{...}{Additional parameters passed to \code{method}.
\itemize{
  \item \code{bin}: Integer. The number of breaks in the \code{"MI"} method for continuous profiles.
}}

\item{n}{The number of CPUs or processors.}
}
\value{
A \code{PPResult} object.
}
\description{
Similarity or distance
}
\details{
Similarity and distance of paired profiles. If the input is a \code{PPTreeIdx} object, the paired profile is collapsed according to the phylogenetic tree.
}
\examples{
require('magrittr')
require('ape')

tree <- system.file('extdata', 'bioinfoTree.nex', package = "PhyloProfileKit") \%>\% read.nexus
ppPath <- system.file('extdata', 'bioinfoProfile.csv', package = "PhyloProfileKit")

sceP <- ppPath \%>\% read.csv(row.names = 1) \%>\% as.matrix \%>\% PP
scePI <- PPIdx(sceP, 1:6, 1:6)
scePTI <- sceP \%>\% PPTree(tree) \%>\% PPTreeIdx(1:6, 1:6)

## Mutual information
SimDist(scePI, 'SimMI', n = 2)
SimDist(scePTI, 'SimMI', n = 2)

}
\author{
Yulong Niu \email{yulong.niu@hotmail.com}
}
