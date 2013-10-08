\name{getMatrixByID}
\docType{methods}
\alias{getMatrixByID}
\alias{getMatrixByID,character-method}
\alias{getMatrixByID,JASPAR2014-method}
\alias{getMatrixByID,SQLiteConnection-method}
\alias{getMatrixByName}
\alias{getMatrixByName,character-method}
\alias{getMatrixByName,JASPAR2014-method}
\alias{getMatrixByName,SQLiteConnection-method}

\title{Basic JASPAR database search functionis \code{getMatrixByID}, 
  \code{getMatrixByName}}
\description{
  This method fetches matrix data under the given ID or name from the database and returns a XMatrix object.
}

\usage{
  \S4method{getMatrixByID}{character}(x, ID)
  \S4method{getMatrixByID}{SQLiteConnection}(x, ID)
  \S4method{getMatrixByID}{JASPAR2014}(x, ID)
  \S4method{getMatrixByName}{character}(x, name)
  \S4method{getMatrixByName}{SQLiteConnection}(x, name)
  \S4method{getMatrixByName}{JASPAR2014}(x, name)
}

\arguments{
  \item{x}{
    a character vector of length 1 for the path of JASPAR SQLite file, a \code{SQLiteConnection} object, or a \code{JASPAR2014} object.
  }
  \item{ID}{
    a character vector of length 1 for the JASPAR stable ID. 
    See more details below.
  } 
  \item{name}{
    a character vector of length 1 for the JASPAR profile name.
  }
}

\details{
  For getMatrixByID, 
  ID is a string which refers to the stable JASPAR ID 
  (usually something like "MA0001") with or without version numbers. 
  "MA0001" will give the latest version on MA0001, while "MA0001.2" will give the second version, if existing. 

  For getMatrixByName, according to the current JASPAR data model, 
  name is not necessarily a unique identifier. 
  Also, names change over time. 
  In the case where there are several matrices with the same name in the database, 
  the function fetches the first one and prints a warning. 
  You've been warned. 
  Some matrices have multiple versions. 
  The function will return the latest version. 
  For specific versions, use getMatrixByID(ID.version)
}

\value{
  A XMatrix object is returned. 
  The exact type of the object depending on the \code{type} argument.
  NA is returned if matrix with the given ID/name is not found.
}

\author{
  Ge Tan
}

\seealso{
  \code{\link{getMatrixSet}}
}

\examples{
  \dontrun{
    db = "/Users/gtan/CSC/JASPAR/JASPAR2014/inst/extdata/JASPAR2014.sqlite"
    pfm = getMatrixByID(db, ID="MA0003")
    pfm2 = getMatrixByName(db, name="TFAP2A")

    library(JASPAR2014)
    pfm3 = getMatrixByID(JASPAR2014, ID="MA0003")
  }
}

\keyword{methods}
\keyword{ ~~ other possible keyword(s) ~~ }
