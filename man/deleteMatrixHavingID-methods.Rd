\name{deleteMatrixHavingID}
\docType{methods}
\alias{deleteMatrixHavingID}
\alias{deleteMatrixHavingID,character-method}
\alias{deleteMatrixHavingID,JASPAR2014-method}
\alias{deleteMatrixHavingID,JASPAR2016-method}
\alias{deleteMatrixHavingID,SQLiteConnection-method}
\alias{storeMatrix}
\alias{storeMatrix-methods}
\alias{storeMatrix,character,PFMatrix-method}
\alias{storeMatrix,character,PFMatrixList-method}
\alias{storeMatrix,JASPAR2014,PFMatrix-method}
\alias{storeMatrix,JASPAR2016,PFMatrix-method}
\alias{storeMatrix,JASPAR2014,PFMatrixList-method}
\alias{storeMatrix,JASPAR2016,PFMatrixList-method}
\alias{storeMatrix,SQLiteConnection,PFMatrix-method}
\alias{storeMatrix,SQLiteConnection,PFMatrixList-method}
\alias{initializeJASPARDB}
\alias{initializeJASPARDB,character-method}
\alias{initializeJASPARDB,JASPAR2014-method}
\alias{initializeJASPARDB,JASPAR2016-method}
\alias{initializeJASPARDB,SQLiteConnection-method}

\title{JASPAR database operations}
\description{
  The functions to initialize, store matrix or delete matrix in JASPAR database.
}

\usage{
  \S4method{deleteMatrixHavingID}{character}(x, IDs)
  \S4method{deleteMatrixHavingID}{SQLiteConnection}(x, IDs)
  \S4method{deleteMatrixHavingID}{JASPAR2014}(x, IDs)
  \S4method{storeMatrix}{character,PFMatrixList}(x, pfmList)
  \S4method{storeMatrix}{SQLiteConnection,PFMatrixList}(x, pfmList)
  \S4method{storeMatrix}{JASPAR2014,PFMatrixList}(x, pfmList)
  \S4method{storeMatrix}{character,PFMatrix}(x, pfmList)
  \S4method{storeMatrix}{SQLiteConnection,PFMatrix}(x, pfmList)
  \S4method{storeMatrix}{JASPAR2014,PFMatrix}(x, pfmList)
  \S4method{initializeJASPARDB}{SQLiteConnection}(x, version=c("2014", "2016"))
  \S4method{initializeJASPARDB}{character}(x, version=c("2014", "2016"))
  \S4method{initializeJASPARDB}{JASPAR2014}(x, version)
  \S4method{initializeJASPARDB}{JASPAR2016}(x, version)
}

\arguments{
  \item{x}{
    A character vector of length 1 for the path of JASPAR SQLite file, 
    or a \code{SQLiteConnection} object.
  }
  \item{IDs}{
    JASPAR stable IDs.
  }
  \item{pfmList}{
    The PFMatrixList object, or pfm object.
  }
  \item{version}{
    Which version of JASPAR to create. So far, it supports 2014 and 2016.
  }
}

\value{
  If the operation works, a "success" will be returned.
}

\examples{
  initializeJASPARDB("jaspar.sqlite", version="2014")
  data("MA0043")
  storeMatrix("jaspar.sqlite", MA0043)
  deleteMatrixHavingID("jaspar.sqlite","MA0043.1")
  file.remove("jaspar.sqlite")
}

\keyword{methods}
\keyword{ ~~ other possible keyword(s) ~~ }
