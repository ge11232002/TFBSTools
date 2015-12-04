\name{makeFlatFileDir}
\alias{makeFlatFileDir}

\title{
  Generate "FlatFileDir" directory 
}
\description{
  On JASPAR web service, "FlatFileDir" includes all the *.pfm and a 
  matrix_list.txt file
}
\usage{
  makeFlatFileDir(JASPAR)
}

\arguments{
  \item{JASPAR}{
    A JASPAR object. Now it can be JASPAR2014 or JASPAR2016.
}
}
\details{
  The matrix_list.txt file contains each pfm per line.
  Each line has the ID, total information content, name, class and tags of
  one pfm.
}
\value{
  The generated files are under "FlatFileDir" directory.
}

\author{
  Ge Tan
}

\examples{
  \donttest{
    library(JASPAR2014)
    makeFlatFileDir(JASPAR2014)
  }
}

