\name{Site}
\Rdversion{1.1}
\docType{class}
\alias{Site}
\alias{Site-class}
\alias{[,Site-method}
\alias{c-methods}
\alias{c,Site-method}
\alias{length,Site-method}
\alias{pattern,Site-method}
\alias{primary,Site-method}
\alias{relScore}
\alias{relScore,Site-method}
\alias{score,Site-method}
\alias{seqname,Site-method}
\alias{show,Site-method}
\alias{sitesource,Site-method}
\alias{strand,Site-method}
\alias{views}
\alias{views,Site-method}
\alias{writeGFF2}
\alias{writeGFF2,Site-method}
\alias{writeGFF3}
\alias{writeGFF3,Site-method}

\title{Class \code{"Site"}}
\description{
%%  ~~ A concise (1-5 lines) description of what the class is. ~~
  The Site object is a container for storing a set of putative 
  transcription factor binding sites on a nucleotide sequence
  (start, end, strand, score, pattern as a PWMatrix, etc.)
}

\usage{
## Constructors:
Site(views, score, strand="*", seqname="Unknown", sitesource="TFBS", 
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
    \item{[}{\code{signature(x = "Site")}: Getter function.}
    \item{length}{\code{signature(x = "Site")}: The number of binding sites in this Site.}
    \item{pattern}{\code{signature(x = "Site")}: Returns the PWMatrix used.}
    \item{relScore}{\code{signature(x = "Site")}: 
      Gets relative score (between 0.0 to 1.0) with respect of the score 
      range of the associated pattern (PWMatrix).}
    \item{score}{\code{signature(x = "Site")}: Returns the score of each site.}
    \item{seqname}{\code{signature(x = "Site")}: Returns the sequence name of the sequence which contains these sites.}
    \item{strand}{\code{signature(x = "Site")}: Returns the strand information.}
    \item{views}{\code{signature(x = "Site")}: Returns the views object.}
    \item{writeGFF2}{\code{signature(x = "Site")}: Returns a standard GFF2 format data frame.}
    \item{writeGFF3}{\code{signature(x = "Site")}: Returns a standard GFF3 format data frame.}
	 }
}
\author{
Ge Tan
}

\seealso{
%%  ~~objects to See Also as \code{\link{~~fun~~}}, ~~~
%%  ~~or \code{\linkS4class{CLASSNAME}} for links to other classes ~~~
  \code{\link{searchSeq}},
  \code{\link{searchAln}},
  \code{\linkS4class{PWMatrix}},
  \code{\linkS4class{SiteList}},
  \code{\linkS4class{SitePair}}
}


\examples{
  \dontrun{
    site = searchSeq(pwm, "GAATTCTCTCTTGTTGTAGTCTCTTGACAAAATG", min.score="60%")
    writeGFF3(pwm)
    relScore(site)
  }
}
\keyword{classes}
