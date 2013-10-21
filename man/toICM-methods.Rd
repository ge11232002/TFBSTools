\name{toICM}
\docType{methods}
\alias{toICM}
\alias{toICM,character-method}
\alias{toICM,DNAStringSet-method}
\alias{toICM,matrix-method}
\alias{toICM,PFMatrix-method}
\title{toICM method}
\description{
  Converts a raw frequency matrix (PFMatrix) to a information content matrix (ICMatrix).
  It takes the bases background frequencies, pseudocounts and schneider as parameters.
}

\usage{
toICM(x, pseudocounts=0.8, schneider=FALSE, bg=c(A=0.25, C=0.25, G=0.25, T=0.25))
}

\arguments{
  \item{x}{
    For \code{toPWM}, a \link{PFMatrix},
    rectangular \link{DNAStringSet} object ("rectangular" means that all elements have the same number of characters) with no IUPAC ambiguity letters,
    a rectangular \link{character} vector or
    a \link{matrix} with rownames containing at least A, C, G and T.
  }
  \item{pseudocounts}{
    A default value 0.8 is used.
  }
  \item{schneider}{
    This logical parameter controls whether a Schneider correction will be done.
    See more details below.
  }
  \item{bg}{
    bg is a vector of background frequencies of four bases with names containing A, C, G, T. When toPWM is applied to a \code{PFMatrix}, if bg is not specified, it will use the bg information contained in \code{PFMatrix}.
  }
}
\details{
  The information content matrix has a column sum between 0 (no base preference) and 2 (only 1 base used). Usually this information is used to plot sequence log.
  
  The information content at each position is computed
  \deqn{D = \log_2(nrow(pfm)) + colSums(postProbs * \log_2(postProbs)}{}
  \deqn{icm = postProbs * D}

  If a Schneider correction will be done if requested. 
  Please see the reference below for more comprehensive explanation.
}

\value{
  A \code{ICMatrix} object which contains the background frequency, pseudocounts and Schneider correction used.
}

\references{
Schneider, T. D., Stormo, G. D., Gold, L., & Ehrenfeucht, A. (1986). Information content of binding sites on nucleotide sequences. Journal of molecular biology, 188(3), 415-431.
}

\author{
  Ge Tan
}

\seealso{
  \code{\link{toPWM}},
  \code{\linkS4class{XMatrix}},
  \code{\link{seqLogo}}
}

\examples{
  ## Constructor a PFMatrix
  pfm = PFMatrix(ID="MA0004.1", name="Arnt", matrixClass="Zipper-Type", strand="+",
        bg=c(A=0.25, C=0.25, G=0.25, T=0.25),
        tags=list(family="Helix-Loop-Helix", species="10090", tax_group="vertebrates",
        medline="7592839", type="SELEX", ACC="P53762", pazar_tf_id="TF0000003",
        TFBSshape_ID="11", TFencyclopedia_ID="580"),
        matrix=matrix(c(4L,  19L, 0L,  0L,  0L,  0L,
                16L, 0L,  20L, 0L,  0L,  0L,
                0L,  1L,  0L,  20L, 0L,  20L,
                0L,  0L,  0L,  0L,  20L, 0L),
                byrow=TRUE, nrow=4, dimnames=list(c("A", "C", "G", "T")))
        )
   ## Convert it into a PWMatrix 
  icm = toICM(pfm, pseudocounts=0.8, schneider=TRUE)
}

\keyword{methods}
