\name{dmmEM-methods}
\docType{methods}
\alias{dmmEM}
\alias{dmmEM-methods}
\alias{dmmEM,ANY-method}
\alias{dmmEM,matrix-method}
\alias{dmmEM,PFMatrixList-method}
\title{dmmEM method}
\description{
  This function trains the Dirichlet multinomial mixture models parameters for
  a set of profile matrices.
}

\usage{
  dmmEM(x, K=6, alg=c("C", "R"))
}

\arguments{
  \item{x}{
    x can be a \code{matrix}, \code{PFMatrixList} or 
    \code{JASPAR2014} to be trained.
  }
  \item{K}{
    The maximal number of components to test in the mixture model when 
    alg is "C". 
    Then an optimal number of components between 1 and K will be chosen 
    based on the fitness of the model.
    
    The fixed number of components to use when alg is "R".
    The default is 6.
  }
  \item{alg}{
    The algorithm to use.
    "C" uses the implementation from \code{DirichletMultinomial} package
    which has more advanced feature and performance.
    "R" uses our own implemention in R.
  }
}

\details{
  When using the implementation from \code{DirichletMultinomial} package,
  the final number of components can be 1:K.
  An internal selection will be made based on the maximum likelihood.

  When using the implementation of R, the number of component is fixed to K.
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
    dmmParameters <- dmmEM(pfmList, K=6, alg="C")
  }
}

\keyword{methods}

