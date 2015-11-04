\name{sampleRanges}
\alias{sampleRanges}

\title{
  sampleRanges
}
\description{
  Sample ranges with same widths of input rannges
  from a set of subject ranges.
}
\usage{
  sampleRanges(inputGRanges, subjectGRanges, ignore.strand=TRUE)
}

\arguments{
  \item{inputGRanges}{
    The input \code{GRanges}.
    }
  \item{subjectGRanges}{
    The subject \code{GRanges}.
  }
  \item{ignore.strand}{
    When set to \code{TRUE}, the strand information is ignoreed during the 
    sampling. Otherwise, the input ranges on positvie strand will only sample 
    from subject ranges on positvie strand.
  }
}


\value{
  A \code{GRanges} object with the same length and widths of 
  \code{inputGRanges}.
}

\author{
  Ge Tan
}

\examples{
  library(GenomicRanges)
  inputGRanges <- GRanges(seqnames=c("chr1", "chr2"), 
                         range=IRanges(start=c(2L, 10L), end=c(6L, 15L)),
                          strand=c("+", "-"))

  subjectGRanges <- GRanges(
                            seqnames=c("chr1", "chr1", "chr1", "chr1", 
                                       "chr2", "chr2"),
                            ranges=IRanges(start=c(20L, 20L, 30L, 30L, 7L, 25L),
                                           end=c(50L, 50L, 32L, 32L,9L, 55L)),
                            strand=c("+","-", "+", "-", "+","-"))
  set.seed(16)
  sampleRanges(inputGRanges, subjectGRanges, ignore.strand=TRUE)
  sampleRanges(inputGRanges, subjectGRanges, ignore.strand=FALSE)

}
