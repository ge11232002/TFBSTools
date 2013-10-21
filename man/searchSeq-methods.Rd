\name{searchSeq}
\docType{methods}
\alias{searchSeq}
\alias{searchSeq-methods}
\alias{searchSeq,PWMatrix-method}
\alias{searchSeq,PWMatrixList-method}
\title{searchSeq method}
\description{
  It scans a nucleotide sequence with the pattern represented by a PWMatrix
  and identifies putative transcription factor binding sites.
}

\usage{
  searchSeq(x, subject, seqname="Unknown", strand="*", min.score="80\%")
}

\arguments{
  \item{x}{
    x can be a \link{PWMatrix} object or a \link{PWMatrixList} object.
  }
  \item{subject}{
    A \code{DNAString}, \code{XStringViews} or \code{MaskedDNAString} object which will be scanned.
  }
  \item{seqname}{
    This is sequence name of the target sequence.
  }
  \item{strand}{
    When searching the sequence, 
    we can search the positive strand or negative strand.
    While strand is "*", it will search both strands and
    return the results based on the positvie strand coordinate.
  }
  \item{min.score}{
    The minimum score for the hit. Can be given an character string in th format
    of "80%" or as a single absolute value.
  }
}

\value{
  A Site object is returned when x is a \code{PWMatrix} object, while 
  a SiteList object is returned when x is a \code{PWMatrixList}.
}

\author{
  Ge Tan
}

\seealso{
  \code{\link{searchAln}},
  \link[Biostrings]{matchPWM}
}

\examples{
  \dontrun{
    site = searchSeq(pwm, "GAATTCTCTCTTGTTGTAGTCTCTTGACAAAATG", min.score="60%")
    siteList = searchSeq(pwmList, "GAATTCTCTCTTGTTGTAGTCTCTTGACAAAATG", 
                         min.score="60%", strand="*", type="any")
  }
}

\keyword{methods}
