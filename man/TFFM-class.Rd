\name{TFFM}
\alias{TFFM}
\alias{TFFMFirst}
\alias{TFFMDetail}
\alias{ncol,TFFMDetail-method}
\alias{ncol,TFFMFirst-method}

\title{
  The TFFM class
}
\description{
  The TFFM is a virtual class. 
  Two classes are derived from this class: TFFMFirst and TFFMDetail.
  
  TFFMFirst class stands for the first-order TFFMs and 
  TFFMDetail stands for the more detailed and descriptive TFFMs.
}

\usage{
## constructors:
  TFFMFirst(ID="Unknown", name="Unknown", matrixClass="Unknown",
            strand="+", bg=c(A=0.25, C=0.25, G=0.25, T=0.25),
            tags=list(), profileMatrix=matrix(),
            type=character(), emission=list(),
            transition=matrix())
  TFFMDetail(ID="Unknown", name="Unknown", matrixClass="Unknown",
             strand="+", bg=c(A=0.25, C=0.25, G=0.25, T=0.25),
             tags=list(), profileMatrix=matrix(),
             type=character(), emission=list(),
             transition=matrix())

## Accessor-like mothods:
  \S4method{ncol}{TFFMFirst}(x)
  \S4method{ncol}{TFFMDetail}(x)
}

\arguments{
  \item{ID,name,matrixClass,strand,bg,tags,profileMatrix}{
    See \code{\link{XMatrix}}
  }
  \item{type}{
    The type of TFFM.
  }
  \item{emission}{
    The emission distribution parameters.
  }
  \item{transition}{
    The transition probability matrix.
  }
  \item{x}{
    a \code{\link{TFFMFirst}} object or a \code{\link{TFFMDetail}} object.
  }
}

\references{
  Mathelier, A., and Wasserman, W.W. (2013). The next generation of transcription factor binding site prediction. PLoS Comput. Biol. 9, e1003214.

\url{http://cisreg.cmmt.ubc.ca/TFFM/doc/#}
}

\author{
  Ge Tan
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
