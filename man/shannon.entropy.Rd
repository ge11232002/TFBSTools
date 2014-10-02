\name{shannon.entropy}
\alias{shannon.entropy}

\title{
  Calculate the Shannon entropy
}
\description{
  This function calculates the Shannon entropy for a discrete random variable
  with finite n values sample.
}
\usage{
shannon.entropy(p)
}

\arguments{
  \item{p}{
  A \code{numeric} object with non-negative values.
}
}
\details{
  The entropy is calculated by H(x) = -sum_i^n(P(x_i)log_b(P(x_i))).
}
\value{
  A \code{numeric} value of entropy is returned.
}

\author{
  Ge Tan
}

\examples{
  x <- c(1, 1, 1, 1)
  shannon.entropy(x)
  x <- c(1, 0, 0, 0)
  shannon.entropy(x)
}
