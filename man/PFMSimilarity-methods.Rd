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
  For each \code{pfmSubject}, an absolute score and 
  a relative percentage score is returned. 
  The maximum absolute score is 2*the width of the smaller matrix 
  in the comparison pair.
}

\references{
  Sandelin, A., H glund, A., Lenhard, B., & Wasserman, W. W. (2003). Integrated analysis of yeast regulatory sequences for biologically linked clusters of genes. Functional & Integrative Genomics, 3(3), 125-134. doi:10.1007/s10142-003-0086-6
}

\author{
  Ge Tan
}

\examples{
  library(Biostrings)
  library(JASPAR2016)
  ## Example matrix from JASPAR database
  profileMatrix <- matrix(as.integer(
    c(13, 13,  3,  1, 54,  1,  1,  1,  0,  3,  2,  5,
      13, 39,  5, 53,  0,  1, 50,  1,  0, 37,  0, 17,
      17,  2, 37,  0,  0, 52,  3,  0, 53,  8, 37, 12,
      11,  0,  9,  0,  0,  0,  0, 52,  1,  6, 15, 20)),
    nrow=4, byrow=TRUE, dimnames=list(DNA_BASES))
  pfmQuery <- PFMatrix(profileMatrix=profileMatrix)
  pfmSubjects <- getMatrixSet(JASPAR2016,
                              opts=list(ID=c("MA0500", "MA0499", "MA0521",
                                             "MA0697", "MA0048", "MA0751",
                                             "MA0832")))
  PFMSimilarity(pfmSubjects, pfmQuery)
}

\keyword{methods}