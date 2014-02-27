\name{rPWMDmm-methods}
\docType{methods}
\alias{rPWMDmm}
\alias{rPWMDmm-methods}
\alias{rPWMDmm,matrix-method}
\alias{rPWMDmm,PFMatrix-method}
\alias{rPWMDmm,PFMatrixList-method}
\title{rPWMDmm method}
\description{
  This function sample matrice from trainned Dirichlet mixture model based on selected matrices.
}

\usage{
  rPWMDmm(x, alpha0, pmix, N=1, W=6)
}
\arguments{
  \item{x}{
    x can be a \code{matrix}, \code{PFMatrixList}. 
    The count matrix on which the sampling is based.
  }
  \item{alpha0}{
    The trained Dirichlet mixture parameters.
  }
  \item{pmix}{
    The trained mixing proportions.
  }
  \item{N}{
    The number of matrice to sample.
  }
  \item{W}{
    The desired width of matrice from the sampling.
  }
}

\value{
  A list of matrices from the sampling.
}

\details{
This feature enables the users to generate random Position Frequency Matrices 
(PFMs) from selected profiles.

We assume that each column in the profile is independent and 
described by a mixture of Dirichlet multinomials 
in which the letters are drawn from a multinomial and 
the multinomial parameters are drawn from a mixture of Dirichlets. 
Within this model each column has its own set of multinomial parameters 
but the higher level parameters -- 
those of the mixture prior is assumed to be common to all Jaspar matrices. 
We can therefore use a maximum likelihood approach to learn these 
from the observed column counts of all Jaspar matrices. 
The maximum likelihood approach automatically ensures that 
matrices receive a weight relative to the number of counts it contains.

Drawing samples from the prior distribution will generate PWMs
with the same statistical properties as the Jaspar matrices as a whole. 
PWMs with statistical properties like those of the selected profiles 
can be obtained by drawing from a posterior distribution 
which is proportional to the prior times a multinomial likelihood term 
with counts taken from one of the columns of the selected profiles.

Each 4-dimensional column is sampled by the following three-step procedure: 
1. draw the mixture component according to the distribution of 
mixing proportions, 
2. draw an input column randomly from the concatenated selected profiles 
and 3. draw the probability vector over nucleotides from a 
4-dimensional Dirichlet distribution. 
The parameter vector alpha of the Dirichlet is equal to the sum of the count 
(of the drawn input)  and the parameters of the Dirichlet prior 
(of the drawn component). 

Draws from a Dirichlet can be obtained in the following way from 
Gamma distributed samples:
(X1,X2,X3,X4) = (Y1/V,Y2/V,Y3/V,Y4/V) ~ Dir(a1,a2,a3,a4)
where V = sum(Yi) ~ Gamma(shape = sum(ai), scale = 1).
}

\section{Methods}{
\describe{

\item{\code{signature(x = "PFMatrix")}}{
%%  ~~describe this method here~~
}

\item{\code{signature(x = "matrix")}}{
%%  ~~describe this method here~~
}

\item{\code{signature(x = "PFMatrixList")}}{
%%  ~~describe this method here~~
}
}}

\references{
L. Devroye, "Non-Uniform Random Variate Generation", Springer-Verlag, 1986

Kimura, T., Tokuda, T., Nakada, Y., Nokajima, T., Matsumoto, T., & Doucet, A. (2011). Expectation-maximization algorithms for inference in Dirichlet processes mixture. Pattern Analysis and Applications, 16(1), 55-67. doi:10.1007/s10044-011-0256-4
}

\note{
  This code is based on the Matlab code original written by Ole Winther, binf.ku.dk, June 2006.
}

\seealso{
  \code{\link{dmmEM}}
}

\author{
  Ge Tan
}

\examples{
  \donttest{
    data(MA0003.2)
    data(MA0004.1)
    pfmList <- PFMatrixList(pfm1=MA0003.2, pfm2=MA0004.1, use.names=TRUE)
    dmmParameters <- dmmEM(pfmList, 6)
    rPWMDmm(MA0003.2, dmmParameters$alpha0, dmmParameters$pmix, N=1, W=6)
  }
}

\keyword{methods}
\keyword{ ~~ other possible keyword(s) ~~ }
