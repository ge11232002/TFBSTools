\name{readXMLTFFM}
\alias{readXMLTFFM}

\title{
  readXMLTFFM
}
\description{
  Read the ouput xml files from Puython module "TFFM" into R.
}
\usage{
  readXMLTFFM(fn, type=c("First", "Detail"))
}

\arguments{
  \item{fn}{
    The path of xml file.
  }
  \item{type}{
    The type of xml file. It can be one of the two types of xml files,
    "First" or "Detail".
  }
}


\value{
  A \code{\link{TFFMFirst}} object or a \code{\link{TFFMDetail}} object is returned.
}

\author{
  Ge Tan
}

\seealso{
  \code{\link{TFFM}}
}
\examples{
  xmlFirst <- file.path(system.file("extdata", package="TFBSTools"),
                        "tffm_first_order.xml")
  tffmFirst <- readXMLTFFM(xmlFirst, type="First")
}

\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
