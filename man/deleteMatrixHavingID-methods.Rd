\name{deleteMatrixHavingID}
\docType{methods}
\alias{deleteMatrixHavingID}
\alias{deleteMatrixHavingID,character-method}
\alias{deleteMatrixHavingID,JASPAR2014-method}
\alias{deleteMatrixHavingID,SQLiteConnection-method}
\alias{storeMatrix}
\alias{storeMatrix-methods}
\alias{storeMatrix,character,PFMatrix-method}
\alias{storeMatrix,character,PFMatrixList-method}
\alias{storeMatrix,JASPAR2014,PFMatrix-method}
\alias{storeMatrix,JASPAR2014,PFMatrixList-method}
\alias{storeMatrix,SQLiteConnection,PFMatrix-method}
\alias{storeMatrix,SQLiteConnection,PFMatrixList-method}
\alias{initializeJASPARDB}
\alias{initializeJASPARDB,character-method}
\alias{initializeJASPARDB,JASPAR2014-method}
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
  \S4method{initializeJASPARDB}{SQLiteConnection}(x)
  \S4method{initializeJASPARDB}{character}(x)
  \S4method{initializeJASPARDB}{JASPAR2014}(x)
}

\arguments{
  \item{x}{
    a character vector of length 1 for the path of JASPAR SQLite file, or a \code{SQLiteConnection} object.
  }
  \item{IDs}{
    JASPAR stable IDs.
  }
  \item{pfmList}{
    The PFMatrixList object, or pfm object.
  }
}

\examples{
  \dontrun{
    library(JASPAR2014)
    deleteMatrixHavingID(JASPAR2014, "MA0003")
    storeMatrix(JASPAR2014, pfm)
    initializeJASPARDB("jaspar.sqlite")
  }
}

\keyword{methods}
\keyword{ ~~ other possible keyword(s) ~~ }
