\name{XMatrixList}
\Rdversion{1.1}
\docType{class}
\alias{XMatrixList}
\alias{XMatrixList,list-method}
\alias{PFMatrixList}
\alias{PFMatrixList-class}
\alias{PWMatrixList}
\alias{PWMatrixList-class}
\alias{ICMatrixList}
\alias{ICMatrixList-class}
\alias{bg,XMatrixList-method}
\alias{ID,XMatrixList-method}
\alias{Matrix,XMatrixList-method}
\alias{matrixClass,XMatrixList-method}
\alias{matrixType,XMatrixList-method}
\alias{name,XMatrixList-method}
\alias{strand,XMatrixList-method}
\alias{tags,XMatrixList-method}

\title{Class \code{"XMatrixList"}}
\description{
  The XMatrixList virtual class is a container for storing a collection of XMatrix objects.
  No object can be constructed directly from this virtual and it has three subclasses: PFMatrixList, PWMatrixList and ICMatrixList.
  Basically it is a SimpleList and is designed for manipulating the set of XMatrix objects as a whole.
}

\usage{
## Constructors:
PFMatrixList(..., use.names=TRUE)
PWMatrixList(..., use.names=TRUE)
ICMatrixList(..., use.names=TRUE)

## Accessor-like methods:
\S4method{ID}{XMatrixList}(x)
\S4method{name}{XMatrixList}(x)

}

\arguments{
  \item{...}{
    The \code{XMatrix} objects are supplied in ....
  }
  \item{use.names}{
    A logical value. When TRUE, the names of the \code{XMatrix} will be kept.
  }
  \item{x}{
    A \code{XMatrixList} object.
  }
}

\author{
Ge Tan
}


\seealso{
  \code{\link{XMatrix}},
}

\examples{
  \dontrun{
    PFMatrixList(pfm1=pfm, pfm2=pfm, use.names=TRUE)
  }
}
\keyword{classes}

