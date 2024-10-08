% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/PPIdx.R
\name{PPIdx}
\alias{PPIdx}
\title{Constructor of \code{PPIdx}}
\usage{
PPIdx(p, x, ...)
}
\arguments{
\item{p}{A \code{PP} object.}

\item{x}{A matrix with two columns or a numeric vector. Proteins not in \code{p} (profile) are removed. The numeric vector indicates the indices of proteins.
\itemize{
  \item a \code{character matrix}: Each row containing the protein names is a linkage.
  \item a \code{numeric matrix} or \code{big.matrix}: Each row containing the protein indices is a linkage.
  \item a \code{numeric vector}: Each element is the index of a interested protein.
}}

\item{...}{Additional parameters if \code{x} is a numeric vector.
\itemize{
  \item \code{y}: Another numeric vector used to generate paired linkages with \code{x}. Each element of \code{x} should be in \code{y}.
  \item \code{self}: Whether include self pairs, and default set is \code{FALSE}.
  \item \code{bidirect}: Whether to include two directions, and default set is \code{FALSE}.
}}
}
\value{
A \code{PPIdx} object.
}
\description{
The constructor the \code{PPIdx} class
}
\details{
Construct a \code{PPIdx} object.
}
\examples{
require('magrittr')

ppBinning <- sample(0:1, 10 * 20, replace = TRUE) \%>\% matrix(ncol = 20) \%>\% PP

## pre-built linkages
linkM <- sample(1:30, 20 * 2, replace = TRUE) \%>\% paste0('protein', .) \%>\% matrix(ncol = 2)
PPIdx(ppBinning, linkM)

## within top 3 proteins
PPIdx(ppBinning, 1:3, 1:3)
## top 3 proteins with whole profiles
PPIdx(ppBinning, 1:3, 1:nrow(ppBinning))
## with self linkages
PPIdx(ppBinning, 1:3, 1:3, self = TRUE)
## with bidirectional linkages
PPIdx(ppBinning, 1:3, 1:3, self = TRUE, bidirect = TRUE)

}
\seealso{
PPTreeIdx
}
\author{
Yulong Niu \email{yulong.niu@hotmail.com}
}
