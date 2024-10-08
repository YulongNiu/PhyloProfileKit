% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppExports.R
\name{eachMI}
\alias{eachMI}
\alias{Info}
\alias{HistTwo}
\alias{FindInter}
\alias{FindInterSingle}
\alias{gInter}
\alias{CountRepeat}
\title{Utilities for MI}
\usage{
eachMI(p1, p2, p3, n)

Info(v, n)

HistTwo(x, y, bin)

FindInter(x, interval)

FindInterSingle(value, interval)

gInter(x, bin)

CountRepeat(x)
}
\arguments{
\item{p1, }{p2, p3: Counts of variables in cells.}

\item{n}{Total variables.}

\item{v}{Histogram of counts.}

\item{x, }{y Numeric vector.}

\item{bin}{Integer.}

\item{value}{Number.}

\item{internal}{Interval numeric vector.}
}
\description{
\code{eachMI()}: Info for a cell.
\code{Info()}: Entropy.
\code{HistTwo()}: Joint counts of two vectors.
\code{FindInter()}: Interval indices of a vector.
\code{FindInterSingle()}: Interval index of a value.
\code{gInter()}: Generate an interval vector, inspired from the \code{cut()} function of the \code{base} package.
\code{CountRepeat()}: Repeat counts of a vector.
}
\author{
Yulong Niu \email{yulong.niu@hotmail.com}
}
\keyword{internal}
