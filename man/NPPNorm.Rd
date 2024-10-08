% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppExports.R
\name{NPPNorm}
\alias{NPPNorm}
\title{z-score normalization of phylogenetic profile}
\usage{
NPPNorm(rawBitM, bitCutoff, bitReset, minConserve)
}
\arguments{
\item{rawBitM}{Raw bit score matrix.}

\item{bitCutoff}{Minimum value of the bit score.}

\item{bitReset}{Reset the bit score for ones lower than the `bitCutoff`.}

\item{minConserve}{Minimum number of homologous. The proteins with homologous less than this value are discarded.}
}
\value{
NPPNorm(): NPP normalized bit score matrix.
}
\description{
NPP normalization
}
\details{
Algorithm:

Step1: rawBitM < hitCutoff to hitReset;

Step2: filter genes without enough homologys

Step3: in each row (species), log2(x/max(x));

Step4: z-score for each column.
}
\references{
\href{http://www.nature.com/nature/journal/v493/n7434/extref/nature11779-s1.pdf}{NPP description}
}
\author{
Yulong Niu \email{yulong.niu@hotmail.com}
}
\keyword{internal}
