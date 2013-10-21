\name{PairwiseAlignmentTFBS}
\Rdversion{1.1}
\docType{class}
\alias{PairwiseAlignmentTFBS}
\alias{alignments}
\alias{alignments,PairwiseAlignmentTFBS-method}
\alias{alnlength,PairwiseAlignmentTFBS-method}
\alias{conservation1,PairwiseAlignmentTFBS-method}
%\alias{searchAln,PWMatrix,PairwiseAlignmentTFBS,missing-method}
%\alias{searchAln,PWMatrixList,PairwiseAlignmentTFBS,missing-method}
\alias{seqlength,PairwiseAlignmentTFBS-method}
\alias{seqname,PairwiseAlignmentTFBS-method}

\title{Class \code{"PairwiseAlignmentTFBS"}}
\description{
  The PairwiseAlignmentTFBS object is a container for pairwise alignment
  and its corresponding information.
}

\usage{
  ## Constructor
  PairwiseAlignmentTFBS(pattern, subject, type="global", 
                        substitutionMatrix=NULL, gapOpening=0,
                        gapExtension=-1, seqname1="Unknown", seqname2="Unknown",
                        windowSize=51L, cutoff=0.7)
}

\section{Slots}{
  \describe{
    \item{\code{alignments}:}{Object of class \code{"PairwiseAlignments"}: A PairwiseAlignments object representing the pairwise alignment sequence.}
    \item{\code{seqname1}:}{Object of class \code{"character"}: Sequence name of first sequence}
    \item{\code{seqname2}:}{Object of class \code{"character"}: Sequence name of second sequence}
    \item{\code{conservation1}:}{Object of class \code{"numeric"}: Conservation score based on the first sequence in the pairwise alignment.}
    \item{\code{windowSize}:}{Object of class \code{"integer"}: The window size used to compute the conservation profile.}
    \item{\code{cutoff}:}{Object of class \code{"numeric"}: The conservation cutoff (between 0 and 1) for including the region in the results of the pattern search.}
    \item{\code{seq1length}:}{Object of class \code{"integer"}: The first sequence's length}
    \item{\code{seq2length}:}{Object of class \code{"integer"}: The second sequence's length}
  }
}
\author{
  Ge Tan
}

\seealso{
%%  ~~objects to See Also as \code{\link{~~fun~~}}, ~~~
%%  ~~or \code{\linkS4class{CLASSNAME}} for links to other classes ~~~
  \code{\link{searchAln}}
}
\examples{
  \dontrun{
    aln = PairwiseAlignmentTFBS(
            pattern="ACTTCACCAGCTCCCTGGCGGTAAGTTGATC---AAAGG---AAACGCAAAGTTTTCAAG",
            subject="GTTTCACTACTTCCTTTCGGGTAAGTAAATATATAAATATATAAAAATATAATTTTCATC",            
            windowSize=13L) 
  }
}
\keyword{classes}
