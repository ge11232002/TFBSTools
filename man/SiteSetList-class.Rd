\name{SiteSetList}
\Rdversion{1.1}
\docType{class}
\alias{SiteSetList}
\alias{SiteSetList-class}
\alias{relScore,SiteSetList-method}
\alias{writeGFF2,SiteSetList-method}
\alias{writeGFF3,SiteSetList-method}

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
    \item{writeGFF2}{\code{signature(x = "SiteSetList")}: Returns a standard GFF2 format data frame.}
    \item{writeGFF3}{\code{signature(x = "SiteSetList")}: Returns a standard GFF3 format data frame.}
	 }
}
\author{
  Ge Tan
}

\seealso{
%%  ~~objects to See Also as \code{\link{~~fun~~}}, ~~~
%%  ~~or \code{\linkS4class{CLASSNAME}} for links to other classes ~~~
  \code{\linkS4class{SiteSet}},
  \code{\link{searchSeq}},
  \code{\link{searchAln}}
}

\examples{
  \dontrun{
    sitesetList = SiteSetList(siteset1, siteset2)
    sitesetList = searchSeq(pwmList, "GAATTCTCTCTTGTTGTAGTCTCTTGACAAAATG", 
                            min.score="60%")
  }
}
\keyword{classes}
