\name{calConservation-methods}
\docType{methods}
\alias{calConservation}
\alias{calConservation-methods}
\alias{calConservation,character,character-method}
\alias{calConservation,character,missing-method}
\alias{calConservation,DNAString,DNAString-method}
\alias{calConservation,DNAStringSet,missing-method}
\title{calConservation method}
\description{
  Calculate the conservation score for a pairwise alignment given a smooth window size.
}

\usage{
  calConservation(aln1, aln2, windowSize=51L, which="1")
}

\arguments{
  \item{aln1}{
    A \code{DNAString} object , a \code{DNAStringSet} or a \code{character} object, 
    which contains the pairwise alignments.
    When the last two objects have a length of 2, the argument \code{aln2} can be missing.
  }
  \item{aln2}{
    A \code{DNAString}, a \code{character} object or missing.
  }
  \item{windowSize}{
    The size of the sliding window (in nucleotides) for calculating local conservation in the alignment. This should be an odd value.
  }
  \item{which}{
    The conservation profile of Which sequence in the alignments is computed. It can be "1" or "2".
  }
}

\value{
  A \code{numeric} vector with the same length of alignment is returned.
}

\author{
  Ge Tan
}

\seealso{
  \code{\link{searchAln}}
}

\keyword{methods}
