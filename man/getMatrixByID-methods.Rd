\name{getMatrixByID}
\docType{methods}
\alias{getMatrixByID}
\alias{getMatrixByID,character-method}
\alias{getMatrixByID,JASPAR2014-method}
\alias{getMatrixByID,JASPAR2016-method}
\alias{getMatrixByID,SQLiteConnection-method}
\alias{getMatrixByName}
\alias{getMatrixByName,character-method}
\alias{getMatrixByName,JASPAR2014-method}
\alias{getMatrixByName,JASPAR2016-method}
\alias{getMatrixByName,SQLiteConnection-method}

\title{Basic JASPAR database search functionis \code{getMatrixByID}, 
  \code{getMatrixByName}}
\description{
  This method fetches matrix data under the given ID or name from 
  the database and returns a XMatrix object.
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
    \code{character}(1) for the path of JASPAR SQLite file, 
    a \code{SQLiteConnection} object, a \code{JASPAR2014}, 
    or a \code{JASPAR2016}object.
  }
  \item{ID}{
    \code{character}() of JASPAR stable ID(s).
    See more details below.
  } 
  \item{name}{
    \code{character}() of JASPAR stable name(s).
  }
}

\details{
  For getMatrixByID, 
  ID is a string which refers to the stable JASPAR ID 
  (usually something like "MA0001") with or without version numbers. 
  "MA0001" will give the latest version on MA0001, 
  while "MA0001.2" will give the second version, if existing. 

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
  A PFMMatrix object is returned when input ID or name is length 1.
  Otherwise, \code{PFMMatrixList} is returned.
}

\author{
  Ge Tan
}

\seealso{
  \code{\link{getMatrixSet}}
}

\examples{
    library(JASPAR2014)
    db <- file.path(system.file("extdata", package="JASPAR2014"), 
                    "JASPAR2014.sqlite")
                    
    ## character and ID
    pfm <- getMatrixByID(db, ID="MA0003")
    
    ## character and IDs
    pfmList <- getMatrixByID(db, ID=c("MA0003", "MA0004"))
    
    ## character and name
    pfm <- getMatrixByName(db, name="TFAP2A")
    
    ##
    ## character and name
    pfmList <- getMatrixByName(db, name=c("TFAP2A", "Arnt"))
    
    ## JASPAR2014 and ID
    pfm <- getMatrixByID(JASPAR2014, ID="MA0003")
}
