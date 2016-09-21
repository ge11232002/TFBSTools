\name{getEmissionProb}
\alias{getEmissionProb}
\alias{getEmissionProb,TFFMDetail-method}
\alias{getEmissionProb,TFFMFirst-method}
\title{
  Get the emission distribution parameters.
}
\description{
  This function accesses the emission distribution parameters of the TFFM.
}
\usage{
  getEmissionProb(tffm)
}

\arguments{
  \item{tffm}{
    A \code{\link{TFFMFirst}} object or a \code{\link{TFFMDetail}} object.
  }
}
\details{
  This function accesses the emission distribution parameters for each position 
  of the TFFM.
  It returns the probability of emitting certain nucleotide based on the
  nucleotide on the previous site.
}
\value{
  A matrix of \code{numeric} with dimensions of 16 * \code{ncol(tffm)}.
}

\author{
  Ge Tan
}

\seealso{
  \code{\link{getPosProb}}
}
\examples{
  xmlFirst <- file.path(system.file("extdata", package="TFBSTools"),
                        "tffm_first_order.xml")
  tffmFirst <- readXMLTFFM(xmlFirst, type="First")
  getEmissionProb(tffmFirst)
  
  xmlDetail <- file.path(system.file("extdata", package="TFBSTools"),
                         "tffm_detailed.xml")
  tffmDetail <- readXMLTFFM(xmlDetail, type="Detail")
  getEmissionProb(tffmDetail)
}
