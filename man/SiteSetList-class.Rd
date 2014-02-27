\name{SiteSetList}
\Rdversion{1.1}
\docType{class}
\alias{SiteSetList}
\alias{SiteSetList-class}
\alias{relScore,SiteSetList-method}
\alias{writeGFF2,SiteSetList-method}
\alias{writeGFF3,SiteSetList-method}
\alias{pvalues,SiteSetList-method}

\title{Class \code{"SiteSetList"}}
\description{
%%  ~~ A concise (1-5 lines) description of what the class is. ~~
  The SiteSetList class is a container for storing a collection of SiteSet objects.
  Basically it is a SimpleList and is designed for manipulating the set of SiteSet objects as a whole.
}

\usage{
## Constructors:
  SiteSetList(..., use.names=TRUE)
}

\arguments{
  \item{...}{
    The SiteSet objects are supplied in .... 
    A list of SiteSet objects is also acceptable.
  }
  \item{use.names}{
    A logical value. When TRUE, the names of the \code{SiteSet} will be kept.
  }
}

\section{Methods}{
  \describe{
    \item{writeGFF2}{\code{signature(x = "SiteSetList")}(x, scoreType=c("absolute", "relative")): 
      Returns a standard GFF2 format data frame.}
    \item{writeGFF3}{\code{signature(x = "SiteSetList")}(x, scoreType=c("absolute", "relative")): 
      Returns a standard GFF3 format data frame.}
    \item{pvalues}{\code{signature(x = "SiteSetList")}(x):
      Calculates the empirical p-values for the scores.
	 }
  }
}
\author{
  Ge Tan
}

\seealso{
  \code{\linkS4class{SiteSet}},
  \code{\link{searchSeq}},
  \code{\link{searchAln}}
}

\examples{
    data(MA0003.2)
    data(MA0004.1)
    pwmList = PWMatrixList(MA0003.2=toPWM(MA0003.2), MA0004.1=toPWM(MA0004.1))
    sitesetList = searchSeq(pwmList, "GAATTCTCTCTTGTTGTAGTCTCTTGACAAAATG", 
                            min.score="50\%")
    writeGFF3(sitesetList, scoreType="absolute")
    pvalues(sitesetList)

}
\keyword{classes}
