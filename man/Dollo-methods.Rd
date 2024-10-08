% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/AllGenerics.R, R/DolloBatch.R
\docType{methods}
\name{Dollo}
\alias{Dollo}
\alias{Dollo,PPTreeIdx-method}
\title{Batch process of Dollo's parsimony distance}
\usage{
Dollo(x, ..., n = 1)

\S4method{Dollo}{PPTreeIdx}(x, ..., n = 1)
}
\arguments{
\item{x}{\code{PPTreeIdx} object.}

\item{...}{Additional parameters.}

\item{n}{The number of CPUs or processors.}
}
\value{
A \code{PPResult} object.
}
\description{
Dollo's parsimony distance
}
\details{
Dollo's parsimony distance of paired profiles.
}
\examples{
require('magrittr')
require('ape')

tree <- system.file('extdata', 'bioinfoTree.nex', package = "PhyloProfileKit") \%>\% read.nexus
ppPath <- system.file('extdata', 'bioinfoProfile.csv', package = "PhyloProfileKit")

sceP <- ppPath \%>\% read.csv(row.names = 1) \%>\% as.matrix \%>\% PP
scePTI <- sceP \%>\% PPTree(tree) \%>\% PPTreeIdx(1:6, 1:6)
Dollo(scePTI, n = 2)

}
\references{
\href{https://www.ncbi.nlm.nih.gov/pubmed/?term=17535793}{Dollo's parsimony description}
}
\author{
Yulong Niu \email{yulong.niu@hotmail.com}
}
