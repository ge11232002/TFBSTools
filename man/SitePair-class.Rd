\name{SitePair}
\Rdversion{1.1}
\docType{class}
\alias{SitePair}
\alias{SitePair-class}
\alias{length,SitePair-method}
\alias{site1,SitePair-method}
\alias{site2,SitePair-method}
\alias{writeGFF2,SitePair-method}
\alias{writeGFF3,SitePair-method}

\title{Class \code{"SitePair"}}
\description{
%%  ~~ A concise (1-5 lines) description of what the class is. ~~
  The SitePair object is a container for storing two Site objects.
  Usually it is used to hold the results returned by \link{searchAln}.
}

\usage{
  ## Constructor
  SitePair(site1, site2)
}
\arguments{
  \item{site1, site2}{
    Each Site object is from one sequence in the pairwise alignment.
  }
}

\section{Methods}{
  \describe{
%    \item{length}{\code{signature(x = "SitePair")}: ... }
    \item{site1}{\code{signature(x = "SitePair")}: Gets the first Site object.}
    \item{site2}{\code{signature(x = "SitePair")}: Gets the second Site object.}
    \item{writeGFF2}{\code{signature(x = "SitePair")}: Returns a standard GFF2 format data frame.}
    \item{writeGFF3}{\code{signature(x = "SitePair")}: Returns a standard GFF3 format data frame.}
	 }
}

\author{
%%  ~~who you are~~
  Ge Tan
}

\seealso{
%%  ~~objects to See Also as \code{\link{~~fun~~}}, ~~~
%%  ~~or \code{\linkS4class{CLASSNAME}} for links to other classes ~~~
  \code{\link{SitePair}},
  \code{\link{searchAln}}
}

\keyword{classes}
