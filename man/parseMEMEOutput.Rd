\name{parseMEMEOutput}
\alias{parseMEMEOutput}

\title{
  parseMEMEOutput
}
\description{
  Parse the output file from \dQuote{MEME}.
}
\usage{
  parseMEMEOutput(x)
}
\arguments{
  \item{x}{
  \code{character}(1): the filename of the \dQuote{MEME} output.
  }
}
\value{
  A list of motifs and evalues is returned.
}
\author{
  Ge Tan
}
\seealso{
  \code{\link{runMEME}}
}
\examples{
  memeOutput <- file.path(system.file("extdata", package="TFBSTools"),
                          "meme.output")
  parseMEMEOutput(memeOutput)
}
