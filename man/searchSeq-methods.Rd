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
    A \code{DNAStringSet}, \code{DNAString}, \code{XStringViews} or 
    \code{MaskedDNAString} object 
    which will be scanned.
  }
  \item{seqname}{
    This is sequence name of the target sequence.
    If subject is a \code{DNAStringSet}, the names of the \code{DNAStringSet}
    object will be used.
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
  A Site object is returned when x is a \code{PWMatrix} object.
  A SiteList object is returned when x is a \code{PWMatrixList} or
  subject is a \code{DNAStringSet}.
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
  seq1 <- "GAATTCTCTCTTGTTGTAGTCTCTTGACAAAATG"
  seq2 <- "GTTTCACTACTTCCTTTCGGGTAAGTAAATATATAAATATATAAAAATATAATTTTCATC"
  siteset <- searchSeq(pwm, seq1, min.score="60\%")
  as(siteset, "data.frame")
  as(siteset, "DataFrame")
  as(siteset, "GRanges")
  writeGFF3(siteset)
  sitesetList <- searchSeq(pwmList, seq1, min.score="60\%", strand="*")
  writeGFF2(sitesetList)
  as(sitesetList, "data.frame")
  as(sitesetList, "DataFrame")
  as(sitesetList, "GRanges")

  library(Biostrings)
  seqs <- DNAStringSet(c(seq1, seq2))
  sitesetList <- searchSeq(pwm, seqs, min.score="60\%")

}

\keyword{methods}
