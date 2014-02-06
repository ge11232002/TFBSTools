\name{SiteSet}
\Rdversion{1.1}
\docType{class}
\alias{SiteSet}
\alias{SiteSet-class}
\alias{[,SiteSet-method}
\alias{[,SiteSet,ANY,ANY-method}
\alias{c-methods}
\alias{c,SiteSet-method}
\alias{length,SiteSet-method}
\alias{pattern,SiteSet-method}
\alias{primary,SiteSet-method}
\alias{relScore}
\alias{relScore,SiteSet-method}
\alias{score,SiteSet-method}
\alias{seqname,SiteSet-method}
\alias{show,SiteSet-method}
\alias{sitesource,SiteSet-method}
\alias{strand,SiteSet-method}
\alias{views}
\alias{views,SiteSet-method}
\alias{writeGFF2}
\alias{writeGFF2,SiteSet-method}
\alias{writeGFF3}
\alias{writeGFF3,SiteSet-method}

\title{Class \code{"SiteSet"}}
\description{
%%  ~~ A concise (1-5 lines) description of what the class is. ~~
  The SiteSet object is a container for storing a set of putative 
  transcription factor binding sites on a nucleotide sequence
  (start, end, strand, score, pattern as a PWMatrix, etc.)
}

\usage{
## Constructors:
SiteSet(views, score, strand="*", seqname="Unknown", sitesource="TFBS", 
        primary="TF binding site", pattern)
}

\arguments{
    \item{views}{Object of class \code{"XStringViews"}: 
      It holds the start, end and the nucleotide sequence information 
      of the transcription factor binding sites.}
    \item{score}{Object of class \code{"numeric"}: A vector of PWM score for each putative binding site.}
    \item{strand}{Object of class \code{"character"}: 
      The binding site is from the positive ("+"), negative ("-") or unknown("*") strand.}
    \item{seqname}{Object of class \code{"character"}: the seqname of the sequence which contains these binding sites.}
    \item{sitesource}{Object of class \code{"character"}: Currently it is set to "TFBS"}
    \item{primary}{Object of class \code{"character"}: Currently it is set to "TF binding site"}
    \item{pattern}{Object of class \code{"PWMatrix"}: The PWMatrix object which is used to search the binding sites.}
}

\section{Methods}{
  \describe{
    \item{[}{\code{signature(x = "SiteSet")}: Getter function.}
    \item{length}{\code{signature(x = "SiteSet")}: The number of binding sites in this SiteSet.}
    \item{pattern}{\code{signature(x = "SiteSet")}: Returns the PWMatrix used.}
    \item{relScore}{\code{signature(x = "SiteSet")}: 
      Gets relative score (between 0.0 to 1.0) with respect of the score 
      range of the associated pattern (PWMatrix).}
    \item{score}{\code{signature(x = "SiteSet")}: Returns the score of each site.}
    \item{seqname}{\code{signature(x = "SiteSet")}: Returns the sequence name of the sequence which contains these sites.}
    \item{strand}{\code{signature(x = "SiteSet")}: Returns the strand information.}
    \item{views}{\code{signature(x = "SiteSet")}: Returns the views object.}
    \item{writeGFF2}{\code{signature(x = "SiteSet")}(x, scoreType=c("absolute", "relative")): 
    Returns a standard GFF2 format data frame.
    The score column can have absolute value or relative value.}
    \item{writeGFF3}{\code{signature(x = "SiteSet")}(x, scoreType=c("absolute", "relative")): 
    Returns a standard GFF3 format data frame.
    The score column can have absolute value or relative value.}
	 }
}

\author{
Ge Tan
}

\seealso{
  \code{\link{searchSeq}},
  \code{\link{searchAln}},
  \code{\linkS4class{PWMatrix}},
  \code{\linkS4class{SiteSetList}},
  \code{\linkS4class{SitePairSet}}
}


\examples{
    data(MA0003.2)
    pwm = toPWM(MA0003.2)
    siteset = searchSeq(pwm, "GAATTCTCTCTTGTTGTAGTCTCTTGACAAAATG", min.score="60\%")
    writeGFF3(siteset)
    relScore(siteset)
}

\keyword{classes}
