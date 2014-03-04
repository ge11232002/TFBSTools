\name{toGRangesList-methods}
\docType{methods}
\alias{toGRangesList}
\alias{toGRangesList-methods}
\alias{toGRangesList,SitePairSetList,Axt-method}
\title{toGRangesList function}
\description{
  Get the genomic coordinates from SitePairSetList.
}
\section{Methods}{
\describe{
\item{\code{signature(x = "SitePairSetList", axt = "Axt")}}{
  Convert the relative coordinates to absolute coordinates.
}
}}

\value{
  A list of two \code{GRanges} objects are returned, 
  one for the target sequences
  and another for query sequences.

  In the \code{GRanges}, strand is taken from the \code{Axt} object.
  In the meta-data columns, PWM matrix ID, the strand of matrix and 
  match score are also returned.
}

\author{
  Ge Tan
}

\examples{
  data(MA0003.2)
  pwm <- toPWM(MA0003.2)
  library(CNEr)
  axtFilesHg19DanRer7 <- file.path(system.file("extdata", package="CNEr"),
                                   "hg19.danRer7.net.axt")
  axtHg19DanRer7 <- readAxt(axtFilesHg19DanRer7)
  sitePairSet <-  searchAln(pwm, axtHg19DanRer7, min.score="80\%",
                            windowSize=51L, cutoff=0.7, strand="*",
                            type="any", conservation=NULL, mc.cores=2)
  toGRangesList(sitePairSet, axtHg19DanRer7)
}


\keyword{methods}
\keyword{ ~~ other possible keyword(s) ~~ }
