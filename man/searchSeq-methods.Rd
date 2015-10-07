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
    \code{MaskedDNAString} object that will be scanned.
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
    in the format of "80\%" or as a single absolute value between 0 and 1.
    When it is percentage value, it represents the quantile 
    between the minimal and the maximal possible value from the PWM.
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
  data(MA0004.1)
  pwm1 <- toPWM(MA0003.2)
  pwm2 <- toPWM(MA0004.1)
  pwmList <- PWMatrixList(pwm1=pwm1, pwm2=pwm2)
  seq1 <- "GAATTCTCTCTTGTTGTAGCATTGCCTCAGGGCACACGTGCAAAATG"
  seq2 <- "GTTTCACCATTGCCTCAGGGCATAAATATATAAAAAAATATAATTTTCATC"
  
  # PWMatrix, character
  ## Only scan the positive strand of the input sequence
  siteset <- searchSeq(pwm1, seq1, seqname="seq1", strand="+", min.score="80\%")
  siteset <- searchSeq(pwm1, seq1, seqname="seq1", strand="+", min.score=0.8)
  ## Only scan the negative strand of the input sequence
  siteset <- searchSeq(pwm1, seq1, seqname="seq1", strand="-", min.score="80\%")
  ## Scan both strands of the input sequences
  siteset <- searchSeq(pwm1, seq1, seqname="seq1", strand="*", min.score="80\%")
  ## Convert the SiteSet object into other R objects
  as(siteset, "data.frame")
  as(siteset, "DataFrame")
  as(siteset, "GRanges")
  writeGFF3(siteset)
  writeGFF2(siteset)
  
  # PWMatrixList, character
  sitesetList <- searchSeq(pwmList, seq1, seqname="seq1", strand="*", 
                           min.score="80\%")
  ## Convert the SiteSteList object into other R objects
  as(sitesetList, "data.frame")
  as(sitesetList, "DataFrame")
  as(sitesetList, "GRanges")
  writeGFF3(sitesetList)
  writeGFF2(sitesetList)

  # PWMatrix, DNAStringSet
  library(Biostrings)
  seqs <- DNAStringSet(c(seq1=seq1, seq2=seq2))
  sitesetList <- searchSeq(pwm1, seqs, min.score="80\%")

  # PWMatrixList, DNAStringSet
  sitesetList <- searchSeq(pwmList, seqs, min.score="80\%")
}

\keyword{methods}
