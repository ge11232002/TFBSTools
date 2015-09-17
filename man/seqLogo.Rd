\name{seqLogo}
\docType{methods}
\alias{seqLogo}
\alias{seqLogo,ICMatrix-method}
\alias{seqLogo,TFFM-method}
\title{Plot a sequence logo}
\description{
  This function takes a \code{ICMatrix} or \code{TFFM} object 
  and plot the sequence logo.
}

\usage{
  seqLogo(x, ic.scale = TRUE, xaxis = TRUE, yaxis = TRUE,
           xfontsize = 15, yfontsize = 15)
}

\arguments{
  \item{x}{
    x is a valid \code{\link{ICMatrix}} object or \code{\link{TFFM}} object.
  }
  \item{ic.scale}{
    A logical value. If TRUE, the total height of one column is proportional to the information content at that position.
    Otherwise, all the columns will have the same height.
    
    Ignored for \code{\link{TFFM}} object.
  }
  \item{xaxis}{
    A logical value. If TRUE, the x-axis will be plotted.
    
    Ignored for \code{\link{TFFM}} object.
  }
  \item{yaxis}{
    A logical value. If TRUE, the y-axis will be plotted.
    
    Ignored for \code{\link{TFFM}} object.
  }
  \item{xfontsize}{
    A numeric value. The font size for x-axis.
  }
  \item{yfontsize}{
    A numeric value. The font size for y-axis.
  }
}

\details{
  A sequence logo is a graphical representation of the matrix model, 
  based on the information content of each position.
  The information content ranges from 0 (no base preference) to 2 (only 1 base used).
  If \code{ic.scale} is TRUE, 
  the height of the logo at certain site is proportinal to 
  the information content value. 
  And each stacked base (A, C, G, T)'s height is also proportional to 
  the information content of each base at that position, 
  and sorted based on the character size.
  
  For a \code{\link{TFFM}} object, a novel graphical representation is used for 
  capturing the dinucleotide dependencies on the TFFM.
  For the upper part of the sequence logo, 
  we represent the nucleotide probabilities 
  at position p for each possible nucleotide at position p-1.
  Hence, each column represents a position within a TFBS and 
  each row the nucleotide probabilities found at that position. 
  Each row assumes a specific nucleotide has been emitted 
  by the previous hidden state.
  The intersection between a column corresponding to position p and 
  row corresponding to nucleotide n gives the probabilities of 
  getting each nucleotide at position p if n has been seen at position p-1.
  The opacity to represent the sequence logo is proportional to 
  the probablity of 
  possible row to be used by the TFFM.
}

\value{
  No return value.
}


\references{
  T D Schneider, R. M. S. (1990). Sequence logos: a new way to display consensus sequences. Nucleic acids research, 18(20), 6097.

  Mathelier, A., and Wasserman, W.W. (2013). The next generation of transcription factor binding site prediction. PLoS Comput. Biol. 9, e1003214.
}

\author{
  Ge Tan
}

\note{
  This function is based on the function \code{seqLogo} from the Bioconductor
  package \code{seqLogo}, especially for the plotting code of TFFM.
}

\seealso{
  \code{\link{toICM}},
  \code{\linkS4class{ICMatrix}},
}

\examples{
  ## ICMatrix
  data(MA0003.2)
  icm = toICM(MA0003.2)
  seqLogo(icm, ic.scale = TRUE)
  
  ## TFFM
  xmlFirst <- file.path(system.file("extdata", package="TFBSTools"),
                        "tffm_first_order.xml")
  tffmFirst <- readXMLTFFM(xmlFirst, type="First")
  seqLogo(tffmFirst)
}

\keyword{methods}
\keyword{ ~~ other possible keyword(s) ~~ }
