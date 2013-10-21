\name{seqLogo}
\docType{methods}
\alias{seqLogo}
\alias{seqLogo,ICMatrix-method}
\title{Plot a sequence logo for a ICMatrix object }
\description{
  This function takes a ICMatrix object and plot the sequence logo.
}

\usage{
  seqLogo(x, ic.scale = TRUE, xaxis = TRUE, yaxis = TRUE,
           xfontsize = 15, yfontsize = 15)
}

\arguments{
  \item{x}{
    x is a valid ICMatrix class.
  }
  \item{ic.scale}{
    A logical value. If TRUE, the total height of one column is proportional to the information content at that position.
    Otherwise, all the columns will have the same height.
  }
  \item{xaxis}{
    A logical value. If TRUE, the x-axis will be plotted.
  }
  \item{yaxis}{
    A logical value. If TRUE, the y-axis will be plotted.
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
  If \code{ic.scale} is TRUE, the height of the logo at certain site is proportinal to the information content value. 
  And each stacked base (A, C, G, T)'s height is also proportional to the information content of each base at that position, and sorted based on the character size.
}

\section{Methods}{
\describe{

\item{\code{signature(x = "ICMatrix")}}{
%%  ~~describe this method here~~
}
}}

\references{
T D Schneider, R. M. S. (1990). Sequence logos: a new way to display consensus sequences. Nucleic acids research, 18(20), 6097.
}

\author{
  Ge Tan
}

\seealso{
  \code{\link{toICM}},
  \code{\linkS4class{ICMatrix}},
}

\examples{
  \dontrun{
    seqLogo(icm, ic.scale = TRUE)
  }  
}

\keyword{methods}
\keyword{ ~~ other possible keyword(s) ~~ }
