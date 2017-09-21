\name{JASPAR2018-class}
\docType{class}
\alias{JASPAR2018-class}
\alias{JASPAR2018}

\title{JASPAR2018 object}
\description{
  The JASPAR2018 object class is a thin class for storing the path of
  JASPAR2018 style SQLite file.
}

\section{Slots}{
  \describe{
    \item{\code{db}:}{Object of class \code{"character"}:
      a character string of the path of SQLite file.}
  }
}

\author{
  Ge Tan
}

\examples{
  \dontrun{
    library(JASPAR2018)
    JASPAR2018

    library(TFBSTools)

    opts <- list()
    opts[["species"]] <- 9606
    opts[["type"]] <- "SELEX"
    opts[["all_versions"]] <- TRUE
    PFMatrixList <- getMatrixSet(JASPAR2018, opts)

    opts2 <- list()
    opts2[["type"]] <- "SELEX"
    PFMatrixList2 <- getMatrixSet(JASPAR201, opts2)
  }
}
\keyword{classes}
