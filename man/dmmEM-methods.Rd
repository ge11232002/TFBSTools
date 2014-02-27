\name{dmmEM-methods}
\docType{methods}
\alias{dmmEM}
\alias{dmmEM-methods}
\alias{dmmEM,ANY-method}
\alias{dmmEM,matrix-method}
\alias{dmmEM,PFMatrixList-method}
\title{dmmEM method}
\description{
  Train Dirichlet mixture models parameters for matrix.
}

\usage{
  dmmEM(x, K, alpha0=NULL, pmix=NULL)
}

\arguments{
  \item{x}{
    x can be a \code{matrix}, \code{PFMatrixList} or 
    \code{JASPAR2014} to be trained.
  }
  \item{K}{
    The number of components in the mixture model.
  }
  \item{alpha0}{
    The initial Dirichlet mixture parameters matrix with dimensions A * K.
    A is the row number of x. It can be missing.
  }
  \item{pmix}{
    The mixing proportions with length K and the sum should be 1.
    It can be missing.
  }
}

\value{
  A list of trainned alpha0, pmix and likelihood during the training.
}

\section{Methods}{
\describe{

\item{\code{signature(x = "ANY")}}{
}

\item{\code{signature(x = "matrix")}}{
}

\item{\code{signature(x = "PFMatrixList")}}{
}
}}

\seealso{
  \code{\link{rPWMDmm}}
}
\author{
  Ge Tan
}

\examples{
  \donttest{
    data(MA0003.2)
    data(MA0004.1)
    pfmList <- PFMatrixList(pfm1=MA0003.2, pfm2=MA0004.1, use.names=TRUE)
    dmmParameters <- dmmEM(pfmList, 6)
  }
}

\keyword{methods}
\keyword{ ~~ other possible keyword(s) ~~ }
