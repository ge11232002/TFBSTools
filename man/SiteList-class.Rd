\name{SiteList}
\Rdversion{1.1}
\docType{class}
\alias{SiteList}
\alias{SiteList-class}
\alias{writeGFF2,SiteList-method}
\alias{writeGFF3,SiteList-method}

\title{Class \code{"SiteList"}}
\description{
%%  ~~ A concise (1-5 lines) description of what the class is. ~~
  The SiteList class is a container for storing a collection of Site objects.
  Basically it is a SimpleList and is designed for manipulating the set of Site objects as a whole.
}

\usage{
## Constructors:
  SiteList(..., use.names=TRUE)
}

\arguments{
  \item{...}{
    The Site objects are supplied in .... 
    A list of Site objects is also acceptable.
  }
  \item{use.names}{
    A logical value. When TRUE, the names of the \code{Site} will be kept.
  }
}

\section{Methods}{
  \describe{
    \item{writeGFF2}{\code{signature(x = "SiteList")}: Returns a standard GFF2 format data frame.}
    \item{writeGFF3}{\code{signature(x = "SiteList")}: Returns a standard GFF3 format data frame.}
	 }
}
\author{
  Ge Tan
}

\seealso{
%%  ~~objects to See Also as \code{\link{~~fun~~}}, ~~~
%%  ~~or \code{\linkS4class{CLASSNAME}} for links to other classes ~~~
  \code{\linkS4class{Site}},
  \code{\link{searchSeq}},
  \code{\link{searchAln}}
}

\examples{
  \dontrun{
    siteList = SiteList(site, site)
    siteList = searchSeq(pwmList, "GAATTCTCTCTTGTTGTAGTCTCTTGACAAAATG", 
                        min.score="60%")
  }
}
\keyword{classes}
