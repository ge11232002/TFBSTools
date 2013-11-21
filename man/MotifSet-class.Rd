\name{MotifSet}
\Rdversion{1.1}
\docType{class}
\alias{MotifSet}
\alias{MotifSet-class}
\alias{[,MotifSet-method}
\alias{[,MotifSet,ANY,ANY-method}
\alias{length,MotifSet-method}
\alias{sitesSeq}
\alias{sitesSeq,MotifSet-method}

\title{Class \code{"MotifSet"}}
\description{
%%  ~~ A concise (1-5 lines) description of what the class is. ~~
  This MotifSet object is a container for storing the generated 
  motifs from Motif identification softwares, such as MEME.
}

\usage{
  ## Constructor
  MotifSet(motifList=GRangesList(), motifEvalues=numeric(), subjectSeqs=DNAStringSet())
}

\arguments{
  \item{motifList}{
    A GRangesList. Each GRanges store the starts, ends, strand, seqnames and scores information of one motif sites sequences.
  }
  \item{motifEvalues}{
    A numeric vector of the E values generated from MEME for each motif.
  }
  \item{subjectSeqs}{
    A DNAStringSet object. It stores the original sequences which are scanned by the software.
  }
}

\section{Methods}{
  \describe{
    \item{[}{\code{signature(x = "MotifSet")}: Getter}
    \item{length}{\code{signature(x = "MotifSet")}: Returns the number of motifs.}
    \item{sitesSeq}{\code{signature(x = "MotifSet")}(x, n=10, type="none"): 
      Gets the sites sequences. 
      \code{n} is the number of bases to include from flanking region. 
      \code{type} controls "all", "left", "right" or "none" flanking sequences are included.}
	 }
}
\author{
  Ge Tan
}

\seealso{
%%  ~~objects to See Also as \code{\link{~~fun~~}}, ~~~
%%  ~~or \code{\linkS4class{CLASSNAME}} for links to other classes ~~~
  \code{\link{runMEME}}
}
\examples{
  \dontrun{
    motifSet = runMEME("/Users/gtan/src/meme_4.9.1/tests/crp0.s",
                    binary="/usr/local/Cellar/meme/4.9.0-p4/bin/meme",
                    arguments="-nmotifs 3")
    sitesSeq(motifSet, type="all")
  }  
}
\keyword{classes}
