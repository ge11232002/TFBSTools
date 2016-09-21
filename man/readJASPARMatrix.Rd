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
  readJASPARMatrix(fn, type = c("individual", "all"))
}
\arguments{
  \item{fn}{
    \code{character}(1): The filename of JASPAR format matrix file.
  }
  \item{type}{
    \code{character}(1): It can be \sQuote{individual} or \sQuote{all}.
  }
}
\details{
  An example of \sQuote{individual} format matrix file is available at 
  \url{http://jaspar.genereg.net/html/DOWNLOAD/JASPAR_CORE/pfm/individual/MA0001.1.pfm}
  
  An exmaple of \sQuote{all} format matrix file is available at 
  \url{http://jaspar.genereg.net/html/DOWNLOAD/JASPAR_CORE/pfm/nonredundant/pfm_all.txt}
}
\value{
  When \sQuote{individual}, a \code{PFMatrix} object is returned.
  When \sQuote{all}, a \code{PFMatrixList} object is returned.
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
  readJASPARMatrix(fn, type="individual")
  fn <- file.path(system.file("extdata", package="TFBSTools"),
                  "pfm_all_example.txt")
  readJASPARMatrix(fn, type="all")
}