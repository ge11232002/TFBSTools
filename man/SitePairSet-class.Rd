\name{SitePairSet}
\Rdversion{1.1}
\docType{class}
\alias{SitePairSet}
\alias{SitePairSet-class}
\alias{length,SitePairSet-method}
\alias{as.data.frame,SitePairSet-method}
\alias{site1,SitePairSet-method}
\alias{site2,SitePairSet-method}

\title{Class \code{"SitePairSet"}}
\description{
%%  ~~ A concise (1-5 lines) description of what the class is. ~~
  The SitePairSet object is a container for storing two SiteSet objects.
  Usually it is used to hold the results returned by \link{searchAln}.
}

\usage{
  ## Constructor
  SitePairSet(siteset1, siteset2)
}
\arguments{
  \item{siteset1, siteset2}{
    Each SiteSet object is from one sequence in the pairwise alignment.
  }
}

\section{Methods}{
  \describe{
%    \item{length}{\code{signature(x = "SitePair")}: ... }
    \item{siteset1}{\code{signature(x = "SitePairSet")}: Gets the first SiteSet object.}
    \item{siteset2}{\code{signature(x = "SitePairSet")}: Gets the second SiteSet object.}
  }
}

\value{
  A \code{\link{SitePairSet}} object.
}

\author{
  Ge Tan
}

\seealso{
%%  ~~objects to See Also as \code{\link{~~fun~~}}, ~~~
%%  ~~or \code{\linkS4class{CLASSNAME}} for links to other classes ~~~
  \code{\link{SitePairSet}},
  \code{\link{searchAln}}
}

\keyword{classes}
