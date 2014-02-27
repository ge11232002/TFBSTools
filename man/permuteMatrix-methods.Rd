\name{permuteMatrix-methods}
\docType{methods}
\alias{permuteMatrix}
\alias{permuteMatrix-methods}
\alias{permuteMatrix,matrix-method}
\alias{permuteMatrix,PFMatrix-method}
\alias{permuteMatrix,PFMatrixList-method}
\title{permuteMatrix method}
\description{
  This method simply shuffles the columns in matrices. 
  This can either be done by just shuffling columns within each selected matrix, 
  or by shuffling columns almong all selected matrices.
}

\usage{
  permuteMatrix(x, type="intra")
}

\arguments{
  \item{x}{
    A \code{matrix} which meets the PFM standard, \code{PFMatrix} object, or 
    \code{PFMatrixList} obejct.
  }
  \item{type}{
    The type of shuffling. 
    It can be "intra" or "inter", which shuffle within each matrix, or
    between all the matrix.
  }
}

\value{
  A object with shuffled matrix.
}

\author{
  Ge Tan
}

\examples{
    library(JASPAR2014)
    pfmSubject = getMatrixByID(JASPAR2014, ID="MA0055")
    pfmQuery = getMatrixByID(JASPAR2014, ID="MA0048")  
    opts = list()
    opts[["class"]] = "Ig-fold"
    pfmList = getMatrixSet(JASPAR2014, opts)
    foo = permuteMatrix(pfmQuery)
    foo1 = permuteMatrix(pfmList, type="intra")
    foo2 = permuteMatrix(pfmList, type="inter")
}

\keyword{methods}
