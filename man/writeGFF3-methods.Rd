\name{writeGFF3-methods}
\docType{methods}
\alias{writeGFF3}
\alias{writeGFF2}
\alias{writeGFF3-methods}
\alias{writeGFF2-methods}
\alias{writeGFF3,SitePairSet-method}
\alias{writeGFF2,SitePairSet-method}
\alias{writeGFF3,SitePairSetList-method}
\alias{writeGFF2,SitePairSetList-method}
\alias{writeGFF3,SiteSet-method}
\alias{writeGFF2,SiteSet-method}
\alias{writeGFF3,SiteSetList-method}
\alias{writeGFF2,SiteSetList-method}

\title{ \code{writeGFF3}, \code{writeGFF2} functions}
\description{
  write the \code{SiteSet}, \code{SitePairSet}, \code{SiteSetList},
  \code{SitePairSetList} into the GFF3 or GFF2 format.
}

\usage{
  writeGFF3(x, scoreType=c("absolute", "relative"))
  writeGFF2(x, scoreType=c("absolute", "relative"))
}

\arguments{
  \item{x}{
    A \code{SiteSet}, \code{SitePairSet}, \code{SiteSetList},
    or \code{SitePairSetList} object.
  }
  \item{scoreType}{
    The score column can have absolute value or relative value.
  }
}

\author{
  Ge Tan
}

\keyword{methods}
\keyword{ ~~ other possible keyword(s) ~~ }
