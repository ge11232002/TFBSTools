\name{runMEME}
\docType{methods}
\alias{runMEME}
\alias{runMEME,character-method}
\alias{runMEME,DNAStringSet-method}
\title{Wrapper function for MEME} 
\description{
  This function builds position frequency matrices using an external program MEME written by Bailey and Elkan.
}

\usage{
  \S4method{runMEME}{character}(x, binary="meme", seqtype="DNA", arguments="", tmpdir=tempdir())
  \S4method{runMEME}{DNAStringSet}(x, binary="meme", seqtype="DNA", arguments="", tmpdir=tempdir())
}

\arguments{
  \item{x}{
    A chracter vector of the path of fasta file or a \code{XStringSet}.
  }
  \item{binary}{
    The path of MEME binary. By default, we assume the meme is in the PATH.
  }
  \item{seqtype}{
    The sequence type. "AA" and "DNA" are allowed.
  }
  \item{arguments}{
    A character vector. The addtional arguments for meme.
  }
  \item{tmpdir}{
    A character vector to change the defult tmp directory.
  }
}

\value{
  A MotifSet object is returned.
}


\references{
  Bailey, T. L., Boden, M., Buske, F. A., Frith, M., Grant, C. E., Clementi, L., et al. (2009). MEME SUITE: tools for motif discovery and searching. Nucleic acids research, 37(Web Server issue), W202-8. doi:10.1093/nar/gkp335
}

\author{
  Ge Tan
}

\seealso{
  \code{\linkS4class{MotifSet}}
}

\examples{
  \dontrun{
    motifSet = runMEME("/Users/gtan/src/meme_4.9.1/tests/crp0.s", 
                    binary="/usr/local/Cellar/meme/4.9.0-p4/bin/meme", 
                    arguments="-nmotifs 3")
  }
}

\keyword{methods}
