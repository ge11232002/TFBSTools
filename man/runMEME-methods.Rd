\name{runMEME}
\docType{methods}
\alias{runMEME}
\alias{runMEME,character-method}
\alias{runMEME,DNAStringSet-method}
\title{Wrapper function for MEME} 
\description{
  This function builds position frequency matrices using 
  an external program \strong{MEME} written by \emph{Bailey} and \emph{Elkan}.
}

\usage{
  \S4method{runMEME}{character}(x, binary="meme", seqtype="DNA", 
                                arguments=list(), tmpdir=tempdir())
  \S4method{runMEME}{DNAStringSet}(x, binary="meme", seqtype="DNA", 
                                   arguments=list(), tmpdir=tempdir())
}

\arguments{
  \item{x}{
    A \code{character}(1) vector of the path of fasta file 
    or a \code{XStringSet}.
  }
  \item{binary}{
    \code{character}(1): the path of MEME binary. 
    By default, we assume the meme is in the \emph{PATH}.
  }
  \item{seqtype}{
    The sequence type. "AA" and "DNA" are allowed.
  }
  \item{arguments}{
    A \code{list}: the addtional arguments for meme.
    This list takes the parameter of MEME as names of the elements,
    and the values of the parameters as the elements.
    For examples, arguments=list("-nmotifs"=3).

  }
  \item{tmpdir}{
    A \code{character}(1) vector to change the defult R's temp directory.
  }
}

\value{
  A \code{MotifSet} object is returned.
}


\references{
  Bailey, T. L., Boden, M., Buske, F. A., Frith, M., Grant, C. E., Clementi, L., et al. (2009). MEME SUITE: tools for motif discovery and searching. Nucleic acids research, 37(Web Server issue), W202-8. doi:10.1093/nar/gkp335
  
  \url{http://meme-suite.org/}
}

\author{
  Ge Tan
}

\note{
  This wrapper is tested on \dQuote{MEME} 4.10.1.
}
\seealso{
  \code{\linkS4class{MotifSet}}
}

\examples{
  \dontrun{
    motifSet <- runMEME(file.path(system.file("extdata", package="TFBSTools"), 
                        "crp0.s"), 
                    binary="/usr/local/Cellar/meme/4.10.1/bin/meme", 
                    arguments=list("-nmotifs"=3))
    ## Get the site sequences                
    sitesSeq(motifSet, type="all")
    sitesSeq(motifSet, type="none")
    
    ## Get the consensu matrix, then it can be used as a PFMatrix
    consensusMatrix(motifSet)
  }
}

\keyword{methods}
