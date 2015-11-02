\name{IUPAC2Matrix}
\alias{IUPAC2Matrix}
\title{
  IUPAC2Matrix
}
\description{
  Convert a IUPAC string into a Postion Weight Matirx
}
\usage{
  IUPAC2Matrix(x)
}

\arguments{
  \item{x}{
    The IUPAC string.
}
}

\details{
  The mapping between IUPAC Extended Genetic Alphabet and the DNA bases letters
  are from \code{IUPAC_CODE_MAP} in Biostrings package.
}

\value{
  A matrix with position weight.
}

\author{
  Ge Tan
}

\examples{
  x <- "RMGNV"
  IUPAC2Matrix(x)

}
