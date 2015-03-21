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
\alias{as.data.frame,SiteSet-method}
\alias{score,SiteSet-method}
\alias{seqname,SiteSet-method}
\alias{show,SiteSet-method}
\alias{sitesource,SiteSet-method}
\alias{strand,SiteSet-method}
\alias{views}
\alias{views,SiteSet-method}
\alias{start,SiteSet-method}
\alias{end,SiteSet-method}
\alias{pvalues}
\alias{pvalues,SiteSet-method}

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
    \item{score}{Object of class \code{"numeric"}: 
      A vector of PWM score for each putative binding site based on the PWM matrix..}
    \item{strand}{Object of class \code{"character"}: 
      The binding site is from the positive ("+"), negative ("-") or 
      unknown("*") strand.}
    \item{seqname}{Object of class \code{"character"}: 
      The seqname of the sequence which contains these binding sites.}
    \item{sitesource}{Object of class \code{"character"}: 
      Currently it is set to "TFBS"}
    \item{primary}{Object of class \code{"character"}: 
      Currently it is set to "TF binding site"}
    \item{pattern}{Object of class \code{"PWMatrix"}: 
      The PWMatrix object which is used to search the binding sites.}
}

\section{Methods}{
  \describe{
    \item{[}{\code{signature(x = "SiteSet")}: Getter function.}
    \item{length}{\code{signature(x = "SiteSet")}: 
      The number of binding sites in this SiteSet.}
    \item{pattern}{\code{signature(x = "SiteSet")}: Returns the PWMatrix used.}
    \item{relScore}{\code{signature(x = "SiteSet")}: 
      Gets relative score (between 0.0 to 1.0) with respect of the score 
      range of the associated pattern (PWMatrix).}
    \item{score}{\code{signature(x = "SiteSet")}: 
      Returns the score of each site.}
    \item{seqname}{\code{signature(x = "SiteSet")}: 
      Returns the sequence name of the sequence which contains these sites.}
    \item{strand}{\code{signature(x = "SiteSet")}: 
      Returns the strand information.}
    \item{views}{\code{signature(x = "SiteSet")}: Returns the views object.}
    \item{start}{\code{signature(x = "SiteSet")}: 
      Returns the start coordinates.}
    \item{end}{\code{signature(x = "SiteSet")}: Returns the end coordinates.}
   \item{pvalues}{\code{signature(x = "SiteSet")}
     (x, type=c("TFMPvalue", "sampling")):
     Calculates the empirical p-values for the scores with two methods:
     the exact method from TFMPaluve package or 
     implementation of sampling in this package.
     The background probability for sampling is based on 
     the PWM matrix in the \code{SiteSet} object.
   }
  }
}

\details{
  The score retuned in SiteSet is the absolute score of each putative TFBS scanned by the corresponding PWM.
  The way of calculating the score is shown on the refernce, Page 281.
}

\references{
Wasserman, W. W., & Sandelin, A. (2004). Applied bioinformatics for the identification of regulatory elements. Nature Publishing Group, 5(4), 276-287. doi:10.1038/nrg1315
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
    pwm <- toPWM(MA0003.2)
    siteset <- searchSeq(pwm, "GAATTCTCTCTTGTTGTAGTCTCTTGACAAAATG", 
                         min.score="60\%")
    writeGFF3(siteset, scoreType="absolute")
    as(siteset, "data.frame")
    as(siteset, "DataFrame")
    as(siteset, "GRanges")
    
    relScore(siteset)
    pvalues(siteset, type="TFMPvalue")
    pvalues(siteset, type="sampling")
}

\keyword{classes}
