\name{searchAln}
\docType{methods}
\alias{searchAln}
\alias{searchAln,PWMatrix,character,character-method}
\alias{searchAln,PWMatrix,character,missing-method}
\alias{searchAln,PWMatrix,DNAString,DNAString-method}
\alias{searchAln,PWMatrix,DNAStringSet,missing-method}
\alias{searchAln,PWMatrix,Axt,missing-method}
%\alias{searchAln,PWMatrix,PairwiseAlignmentTFBS,missing-method}
\alias{searchAln,PWMatrixList,character,character-method}
\alias{searchAln,PWMatrixList,character,missing-method}
\alias{searchAln,PWMatrixList,DNAString,DNAString-method}
\alias{searchAln,PWMatrixList,DNAStringSet,missing-method}
%\alias{searchAln,PWMatrixList,PairwiseAlignmentTFBS,missing-method}
\title{searchAln method}
\description{
  Scans a pairwise alignment of nucleotide sequences 
  with the pattern represented by the PWMatrix. 
  It reports only those hits that are overlapped in the alignment of 
  of the two sequences and exceed a specified threshold score in both, 
  AND are found in regions of the alignment above the specified
  conservation cutoff value.
}

\usage{
  searchAln(pwm, aln1, aln2, seqname1="Unknown1", seqname2="Unknown2",
            min.score="80\%", windowSize=51L, 
            cutoff=0.7, strand="*", type="any", conservation=NULL,
            mc.cores=1L)
}

\arguments{
  \item{pwm}{
    A \link{PWMatrix} object or a \link{PWMatrixList} object.
  }
  \item{aln1}{
    A \code{DNAString}, \code{character}, \code{DNAStringSet} or \code{Axt}
    object can be used to represent the pairwise alignment.
    When the last two objects are used and have a length of 2, 
    the argument \code{aln2} can be missing.
  }
  \item{aln2}{
    A \code{DNAString}, \code{character}.
    It can be missing when aln1 is \code{DNAStringSet} or \code{Axt} object.
  }
  \item{seqname1,seqname2}{
    A \code{chracter} object for the name of sequence. 
    "Unknown1" and "Unknown2" are used by default.
    These two arguments are ignored when aln1 is \code{Axt}, 
    or the seqnames are available from aln1.
  }
  \item{min.score}{
    The minimum score for the hit. Can be given an character string 
    in the format of "80\%" or as a single absolute value.
    When it is percentage value, it means the percentage of the 
    maximal possible from the PWM.
  }
  \item{windowSize}{
    The size of the sliding window (in nucleotides) for calculating 
    local conservation in the alignment. 
    This should be an odd value.
  }
  \item{cutoff}{
    The conservation cutoff can be from 0 (0\% identity) to 1 (100\% identity).
    The regions which have lower conservation than the cutoff 
    will be discarded from the results of the pattern searching.
    The conservation is calculated by comparing the alignments 
    within the \code{windowSize}:
    1 for match and 0 for mismatch and gap.
  }
  \item{strand}{
    When searching the alignment, we can search the positive strand or 
    negative strand. 
    While strand is "*", it will search both strands and 
    return the results based on the positvie strand coordinate.
  }
  \item{type}{
    This argument can be "any" or "all". 
    When it is "any", one motif will be kept 
    if the maximal conservation value of  the motif is larger than the cutoff.
    When it is "all", one motif will be kept 
    if the minimal conservation value of the motif is larger than the cutoff.
  }
  \item{conservation}{
    A vector of conservation profile. 
    If not supplied, the conservation profile will be computed internally 
    on the fly.
  }
  \item{mc.cores}{
    The number of cpu threads to use when searching \code{Axt}. 
    1L is assigned by default.
  }
}

\value{
  A \code{SitePairSet} object is returned when pwm is a \code{PWMatrix}, 
  while a \code{SitePairSetList} is returned when pwm is a \code{PWMatrixList}.
}

\details{
  In brief, given a pairwise alignment of two sequences, 
  first of all, we remove the gaps ("-", "-", ".").
  Then we scan both ungapped sequences with the pwm and 
  return the hits that above \code{min.score}.
  Since we only want to keep the conserved hits, 
  we choose the pair of motifs that overlap most in the alignment.
  Finally, the pair of motifs have to be conserved above 
  the threshold \code{cutoff}.
  
  In the returned \code{SitePairSet}, the coordinates of start, end are based
  on the ungapped sequences, instead of the original alignment.
  This is due to we are more concerned about the actual location of motif
  in the genome rather than in the alignment.
}

\author{
  Ge Tan
}

\seealso{
  \code{\link{searchSeq}}
}

\examples{
  data(MA0003.2)
  data(MA0004.1)
  pwm1 <- toPWM(MA0003.2)
  pwm2 <- toPWM(MA0004.1)
  pwmList <- PWMatrixList(pwm1=pwm1, pwm2=pwm2)
  # Two character objects
  aln1 <- "ACCACATTGCCTCAGGGCAGGTAAGTTGATC---AAAGG---AAACGCAAAGTTTTCAAG"
  aln2 <- "GTTTCACTACATTGCTTCAGGGCAGTAAATATATAAATATATAAAAATATAATTTTCATC"
  aln <- c(aln1=aln1, aln2=aln2)
  library(Biostrings)
  alnDNAStringSet <- DNAStringSet(c(aln1=aln1, aln2=aln2))
  
  # PWMatrix, character, character
  ## Only scan the positive strand of the alignments
  sitePairSet <- searchAln(pwm1, aln1, aln2, seqname1="aln1", seqname2="aln2",
                           min.score="70\%", cutoff=0.5, 
                           strand="+", type="any")
  ## Only scan the negative strand of the alignments
  sitePairSet <- searchAln(pwm1, aln1, aln2, seqname1="aln1", seqname2="aln2",
                           min.score="70\%", cutoff=0.5, 
                           strand="-", type="any")
  ## Scan the both strands of the alignments                         
  sitePairSet <- searchAln(pwm1, aln1, aln2, seqname1="aln1", seqname2="aln2",
                           min.score="70\%", cutoff=0.5, 
                           strand="*", type="any")
  ## Convert the SitePairSet object into other R objects
  as(sitePairSet, "data.frame")
  as.data.frame(sitePairSet)
  as(sitePairSet, "DataFrame")
  as(sitePairSet, "GRanges")
  writeGFF3(sitePairSet)
  writeGFF2(sitePairSet)
  
  # PWMatrix, character, missing
  sitePairSet <- searchAln(pwm1, aln,
                           min.score="70\%", cutoff=0.5, 
                           strand="*", type="any")
  
  # PWMatrix, DNAString, DNAString
  sitePairSet <- searchAln(pwm1, DNAString(aln1), DNAString(aln2), 
                           seqname1="aln1", seqname2="aln2",
                           min.score="70\%", cutoff=0.5, 
                           strand="*", type="any")

  # PWMatrix, DNAStringSet, missing
  sitePairSet <- searchAln(pwm1, alnDNAStringSet,
                           min.score="70\%", cutoff=0.5, 
                           strand="*", type="any")
                           
  # PWMatrixList, character, character
  sitePairSetList <- searchAln(pwmList, aln1, aln2, 
                           seqname1="aln1", seqname2="aln2",
                           min.score="70\%", cutoff=0.5, 
                           strand="*", type="any")
  ## elementLenths of each pwm hits
   elementNROWS(sitePairSetList)
                           
  ## output
  writeGFF2(sitePairSetList)
  writeGFF3(sitePairSetList)
  as(sitePairSetList, "DataFrame")
  as(sitePairSetList, "data.frame")
  as.data.frame(sitePairSetList)
  as(sitePairSetList, "GRanges")
  
  # PWMatrix, Axt, missing
  library(CNEr)
  axtFilesHg19DanRer7 <- file.path(system.file("extdata", package="TFBSTools"), 
                                   "hg19.danRer7.net.axt")
  axtHg19DanRer7 <- readAxt(axtFilesHg19DanRer7)
  sitePairSetList <-  searchAln(pwm1, axtHg19DanRer7, min.score="80\%", 
                            windowSize=51L, cutoff=0.7, strand="*", 
                            type="any", conservation=NULL, mc.cores=1)
  ## We may want to coordinates of motif in the genome
  GRangesTFBS <- toGRangesList(sitePairSetList, axtHg19DanRer7)
}

\keyword{methods}

