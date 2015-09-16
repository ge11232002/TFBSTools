\name{getPosProb}
\alias{getPosProb}
\alias{getPosProb,TFFMDetail-method}
\alias{getPosProb,TFFMFirst-method}
\title{
  Get the emission probabilities of nucleotides
}
\description{
  Get the emission probabilities of ACGT at each position of TFFM.
}
\usage{
  getPosProb(tffm)
}

\arguments{
  \item{tffm}{
    A \code{\link{TFFMFirst}} object or a \code{\link{TFFMDetail}} object.
  }
}
\details{
  This function calculates the probabilities of emitting nucleotides ACGT
  at each position of TFFM.
}
\value{
  A matrix of \code{numeric} with dimensions of 4 * \code{ncol(tffm)}.
}

\author{
  Ge Tan
}

\seealso{
  \code{\link{getEmissionProb}}
}

\examples{
  xmlFirst <- file.path(system.file("extdata", package="TFBSTools"),
                        "tffm_first_order.xml")
  tffmFirst <- readXMLTFFM(xmlFirst, type="First")
  getPosProb(tffmFirst)
  
  xmlDetail <- file.path(system.file("extdata", package="TFBSTools"),
                         "tffm_detailed.xml")
  tffmDetail <- readXMLTFFM(xmlDetail, type="Detail")
  getPosProb(tffmDetail)
}

\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
