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
    A \code{DNAString}, \code{XStringViews} or \code{MaskedDNAString} object 
    which will be scanned.
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
    The minimum score for the hit. Can be given an character string 
    in the format of "80\%" or as a single absolute value.
    When it is percentage value, it means the percentage of 
    the maximal possible from the PWM.
  }
}

\value{
  A Site object is returned when x is a \code{PWMatrix} object, while 
  a SiteList object is returned when x is a \code{PWMatrixList}.
}

\references{
Wasserman, W. W., & Sandelin, A. (2004). Applied bioinformatics for the identification of regulatory elements. Nature Publishing Group, 5(4), 276-287. doi:10.1038/nrg1315
}

\author{
  Ge Tan
}

\seealso{
  \code{\link{searchAln}},
  \link[Biostrings]{matchPWM}
}

\examples{
  data(MA0003.2)
  pwm <- toPWM(MA0003.2)
  pwmList <- PWMatrixList(pwm1=pwm, pwm2=pwm)
  siteset <- searchSeq(pwm, "GAATTCTCTCTTGTTGTAGTCTCTTGACAAAATG", 
                       min.score="60\%")
  sitesetList <- searchSeq(pwmList, "GAATTCTCTCTTGTTGTAGTCTCTTGACAAAATG", 
                           min.score="60\%", strand="*")
}

\keyword{methods}
