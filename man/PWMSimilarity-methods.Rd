\name{PWMSimilarity-methods}
\docType{methods}
\alias{PWMSimilarity}
\alias{PWMSimilarity-methods}
\alias{PWMSimilarity,matrix,matrix-method}
\alias{PWMSimilarity,matrix,PWMatrix-method}
\alias{PWMSimilarity,PWMatrix,matrix-method}
\alias{PWMSimilarity,PWMatrix,PWMatrix-method}
\alias{PWMSimilarity,PWMatrixList,matrix-method}
\alias{PWMSimilarity,PWMatrixList,PWMatrix-method}
\alias{PWMSimilarity,PWMatrixList,PWMatrixList-method}
\title{PWMSimilarity method}
\description{
  This function measures the similarity of two PWM matrix in three measurements:
  "normalised Euclidian distance", "Pearson correlation" and "Kullback Leibler divergence".
}

\usage{
  PWMSimilarity(pwmSubject, pwmQuery, method=c("Euclidian", "Pearson", "KL"))
}

\arguments{
  \item{pwmSubject}{
    A \code{matrix} or \code{PWMatrix} or \code{PWMatrixList} object.
  }
  \item{pwmQuery}{
    A \code{matrix} or \code{PWMatrix} object.
  }
  \item{method}{
    The method can be "Euclidian", "Pearson", "KL".
  }
}

\details{
  When pwmSubject and pwmQuery have different number of columns,
  the smaller PWM will be shifted from the start position of larger PWM 
  and compare all the possible alignments.
  Only the smallest distance, divergence or largest correlation will be reported.
  
  %Given two PWMs, P^1 and P^2, where l is the length. 
  %P_{i,b} is the values in column i with base b.
  %The normalised Euclidian distance is computed in
  %D(a,b) = {1 \over {\sqrt{2}l}} \cdot \sum_{i=1}^{l} \sqrt{\sum_{b \in {\{A,C,G,T\}}} (P_{i,b}^1-P_{i,b}^2)^2}
  %.
  %This distance is between 0 (perfect identity) and 1 (complete dis-similarity).
  
  %The pearson correlation coefficient is computed in
  %r(P^1, P^2) = {1 \over l} \cdot \sum_{i=1}^l {\sum_{b \in \{A,C,G,T\}} (P_{i,b}^1 - 0.25)(P_{i,b}^2-0.25) \over \sqrt{\sum_{b \in \{A,C,G,T\}} (P_{i,b}^1 - 0.25)^2 \cdot \sum_{b \in \{A,C,G,T\}} (P_{i,b}^2 - 0.25)^2}}.
  
  %The Kullback-Leibler divergence is computed in 
  %KL(P^1, P^2) = {1 \over {2l}} \cdot \sum_{i=1}^l \sum_{b \in \{A,C,G,T\}} (P_{i,b}^1\log{ P_{i,b}^1 \over P_{i,b}^2}+ P_{i,b}^2\log{P_{i,b}^2 \over {P_{i,b}^1}}).
}

\references{
  Linhart, C., Halperin, Y., & Shamir, R. (2008). Transcription factor and microRNA motif discovery: The Amadeus platform and a compendium of metazoan target sets. Genome Research, 18(7), 1180-1189. doi:10.1101/gr.076117.108
}

\section{Methods}{
\describe{

\item{\code{signature(pwmSubject = "matrix", pwmQuery = "matrix")}}{
%%  ~~describe this method here~~
}

\item{\code{signature(pwmSubject = "matrix", pwmQuery = "PWMatrix")}}{
%%  ~~describe this method here~~
}

\item{\code{signature(pwmSubject = "PWMatrix", pwmQuery = "matrix")}}{
%%  ~~describe this method here~~
}

\item{\code{signature(pwmSubject = "PWMatrix", pwmQuery = "PWMatrix")}}{
%%  ~~describe this method here~~
}

\item{\code{signature(pwmSubject = "PWMatrixList", pwmQuery = "matrix")}}{
%%  ~~describe this method here~~
}

\item{\code{signature(pwmSubject = "PWMatrixList", pwmQuery = "PWMatrix")}}{
%%  ~~describe this method here~~
}

\item{\code{signature(pwmSubject = "PWMatrixList", pwmQuery = "PWMatrixList")}}{
%%  ~~describe this method here~~
}
}}

\seealso{
  \code{\link{PFMSimilarity}}
}

\examples{
  data(MA0003.2)
  data(MA0004.1)
  pwm1 = toPWM(MA0003.2, type="prob")
  pwm2 = toPWM(MA0004.1, type="prob")
  PWMSimilarity(pwm1, pwm2, method="Euclidian")
}


\keyword{methods}
\keyword{ ~~ other possible keyword(s) ~~ }
