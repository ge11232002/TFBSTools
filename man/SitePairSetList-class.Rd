\name{SitePairSetList-class}
\Rdversion{1.1}
\docType{class}
\alias{SitePairSetList}
\alias{SitePairSetList-class}

\title{Class \code{"SitePairSetList"}}
\description{
%%  ~~ A concise (1-5 lines) description of what the class is. ~~
  The SitePairSetList class is a container for storing a collection of SitePairSet objects.
  Basically it is a SimpleList and is designed for manipulating the set of SitePairSet objects as a whole.
}

\usage{
## Constructors:
  SitePairSetList(..., use.names=TRUE)
}

\arguments{
  \item{...}{
    The SitePairSet objects are supplied in ....
    A list of SitePairSet objects is also acceptable.
  }
  \item{use.names}{
    A logical value. When TRUE, the names of the \code{SitePairSet} will be kept.
  }
}

\author{
Ge Tan
}

\seealso{
%%  ~~objects to See Also as \code{\link{~~fun~~}}, ~~~
%%  ~~or \code{\linkS4class{CLASSNAME}} for links to other classes ~~~
\code{\linkS4class{SitePairSet}},
}
\keyword{classes}
