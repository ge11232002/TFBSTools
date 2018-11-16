\name{readJASPARMatrix}
\alias{readJASPARMatrix}

\title{
  Read JASPAR format matrix
}
\description{
  Read a JASPAR format matrix file with \sQuote{individual} matrix or
  \sQuote{all} matrices in one file.
}
\usage{
  readJASPARMatrix(fn, matrixClass=c("PFM", "PWM", "PWMProb"))
}
\arguments{
  \item{fn}{
    \code{character}(1): The filename of JASPAR format matrix file.
  }
  \item{matrixClass}{
    \code{character}(1): the matrices in the files are \sQuote{PFM}, \sQuote{PWM} or
    \sQuote{PFM} in probability.
  }
}
\details{
  An example of \sQuote{individual} format matrix file is available at 
  \url{http://jaspar.genereg.net/html/DOWNLOAD/JASPAR_CORE/pfm/individual/MA0001.1.pfm}
  
  An exmaple of \sQuote{all} format matrix file is available at 
  \url{http://jaspar.genereg.net/html/DOWNLOAD/JASPAR_CORE/pfm/nonredundant/pfm_all.txt}
}
\value{
  A \code{PFMatrixList} or \code{PWMatrixList} object is returned, depending on the
  matrix class.
}
\author{
  Ge Tan
}
\seealso{
  \code{\link{makeFlatFileDir}}
}
\examples{
  fn <- file.path(system.file("extdata", package="TFBSTools"),
                  "MA0001.1.pfm")
  readJASPARMatrix(fn, matrixClass="PFM")
  fn <- file.path(system.file("extdata", package="TFBSTools"),
                  "pfm_all_example.txt")
  readJASPARMatrix(fn, matrixClass="PFM")
}