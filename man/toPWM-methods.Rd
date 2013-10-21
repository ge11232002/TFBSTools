\name{toPWM}
\docType{methods}
\alias{toPWM}
\alias{toPWM,character-method}
\alias{toPWM,DNAStringSet-method}
\alias{toPWM,matrix-method}
\alias{toPWM,PFMatrix-method}

\title{toPWM method}
\description{
  Converts a raw frequency matrix (PFMatrix) to a position weight matrix (PWMatrix).
  It takes the type, bases background frequencies, pseudocounts as parameters.
}

\usage{
toPWM(x, type="log2probratio", pseudocounts=0.8, bg=c(A=0.25, C=0.25, G=0.25, T=0.25))
}

\arguments{
  \item{x}{
    For \code{toPWM}, a \link{PFMatrix}, 
    rectangular \link{DNAStringSet} object ("rectangular" means that all elements have the same number of characters) with no IUPAC ambiguity letters,
    a rectangular \link{character} vector or 
    a \link{matrix} with rownames containing at least A, C, G and T.
  }
  \item{type}{
    The type of PWM generated, should be one of "log2probratio" or "prob".
  }
  \item{pseudocounts}{
    pseudocounts is a numeric non-negative vector, which means you can specify different pseudocounts for each site. 
    The values will be recycled if shorter than the length of sites. 
    0.8 is recommended. See the reference below for more details.
  }
  \item{bg}{
    bg is a vector of background frequencies of four bases with names containing A, C, G, T. When toPWM is applied to a \code{PFMatrix}, if bg is not specified, it will use the bg information contained in \code{PFMatrix}. 
  }
}

\details{
  The raw position frequency matrix (PFM) is usually converted into 
  a position weight matrix (PWM), also known as position specific scoring matrix (PSSM).
  The PWM provides the probability of each base at certain position and used for scanning the genomic sequences.
  The implementation here is slightly different from \code{PWM} in \code{Biostrings} package by choosing the pseudocounts. 
  Pseudocounts is necessary for correcting the small number of counts 
  or eliminating the zero values before log transformation.

  \deqn{postProbs = {{PFM + bg * pseudocounts} \over {ncol(PFM) + sum(bg) * pseudocounts}}}{}
  \deqn{priorProbs = bg / sum(bg)}{}
  \deqn{PWM_{log2probratio} = \log_2{postProbs \over priorProbs}}{}
  \deqn{PWM_{prob} = postProbs}
}

\value{
  A \code{PWMatrix} object which contains the background frequency and pseudocounts used.
}

\references{
Wasserman, W. W., & Sandelin, A. (2004). Applied bioinformatics for the identification of regulatory elements. Nature Publishing Group, 5(4), 276-287. doi:10.1038/nrg1315

Nishida, K., Frith, M. C., & Nakai, K. (2009). Pseudocounts for transcription factor binding sites. Nucleic acids research, 37(3), 939-944. doi:10.1093/nar/gkn1019
}

\author{
Ge Tan
}

\seealso{
\code{\link{toICM}},
\code{\linkS4class{XMatrix}},
\link[Biostrings]{PWM}
}

\examples{
  ## Constructe a PFMatrix
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
  pwm = toPWM(pfm, pseudocounts=0.8)
}

\keyword{methods}
