\name{PFMSimilarity-methods}
\docType{methods}
\alias{PFMSimilarity}
\alias{PFMSimilarity-methods}
\alias{PFMSimilarity,matrix,character-method}
\alias{PFMSimilarity,matrix,matrix-method}
\alias{PFMSimilarity,matrix,PFMatrix-method}
\alias{PFMSimilarity,PFMatrix,character-method}
\alias{PFMSimilarity,PFMatrix,matrix-method}
\alias{PFMSimilarity,PFMatrix,PFMatrix-method}
\alias{PFMSimilarity,PFMatrixList,character-method}
\alias{PFMSimilarity,PFMatrixList,matrix-method}
\alias{PFMSimilarity,PFMatrixList,PFMatrix-method}
\title{PFMSimilarity method}
\description{
  Given a \code{PFMatrix} or a normal matrix, 
  align it with another set of \code{PFMatrix} to assess the similarity.
}

\usage{
  PFMSimilarity(pfmSubject, pfmQuery, openPenalty=3, extPenalty=0.01)
}

\arguments{
  \item{pfmSubject}{
    A \code{matrix}, \code{PFMatrix} or \code{PFMatrixList} object,
    which is compared with query matrix.
  }
  \item{pfmQuery}{
    A \code{matrix}, \code{PFMatrix} or IUPAC \code{character} object.
  }
  \item{openPenalty}{
    The gap open penalty used in the modified Needleman-Wunsch algorithm.
    By default, it is 3.
  }
  \item{extPenalty}{
    The gap extension penalty used in the modified Needleman-Wunsch algorithm.
    By default, it is 0.01.
  }
}

\value{
  For each \code{pfmSubject}, an absolute score and a relative percentage score is returned. 
  The maximum absolute score is 2*the width of the smaller matrix in the comparison pair.
}

\references{
  Sandelin, A., H glund, A., Lenhard, B., & Wasserman, W. W. (2003). Integrated analysis of yeast regulatory sequences for biologically linked clusters of genes. Functional & Integrative Genomics, 3(3), 125-134. doi:10.1007/s10142-003-0086-6

}

\author{
  Ge Tan
}

\examples{
    library(JASPAR2014)
    ## one to one comparison
    pfmSubject = getMatrixByID(JASPAR2014, ID="MA0055")
    pfmQuery = getMatrixByID(JASPAR2014, ID="MA0048")
    score = PFMSimilarity(pfmSubject, pfmQuery)
    ## one to several comparsion
    opts = list()
    opts[["class"]] = "Ig-fold"
    pfmList = getMatrixSet(JASPAR2014, opts)
    scores = PFMSimilarity(pfmList, pfmQuery)
}

\keyword{methods}
\keyword{ ~~ other possible keyword(s) ~~ }
