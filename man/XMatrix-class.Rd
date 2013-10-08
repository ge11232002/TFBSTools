\name{XMatrix}
\Rdversion{1.1}
\docType{class}
\alias{XMatrix}
\alias{XMatrix-class}
\alias{PFMatrix}
\alias{PFMatrix-class}
\alias{bg,XMatrix-method}
\alias{bg<-,XMatrix-method}
\alias{ID,XMatrix-method}
\alias{ID<-,XMatrix-method}
\alias{length,XMatrix-method}
\alias{Matrix,XMatrix-method}
\alias{Matrix<-,XMatrix-method}
\alias{matrixClass,XMatrix-method}
\alias{matrixClass<-,XMatrix-method}
\alias{name,XMatrix-method}
\alias{name<-,XMatrix-method}
\alias{reverseComplement,XMatrix-method}
\alias{strand,XMatrix-method}
\alias{strand<-,XMatrix-method}
\alias{tags}
\alias{tags,XMatrix-method}
\alias{matrixType-methods}
\alias{matrixType,ICMatrix-method}
\alias{matrixType,PFMatrix-method}
\alias{matrixType,PWMatrix-method}
\alias{[-methods}
\alias{[,Site-method}
\alias{length-methods}
\alias{length,Site-method}
\alias{length,SitePair-method}

\title{\code{"XMatrix"} objects}
\description{
XMatrix is a virtual class. No objects can be created from it directly. Three classes are derived from this class: PFMatrix, PWMatrix and ICMatrix.

PFMatrix is a class whose instances are objects representing raw position frequency matrices (PFMs).
}

\usage{
## Constructors:
PFMatrix(ID="Unknown", name="Unknown", matrixClass="Unknown", 
         strand="*", bg=c(A=0.25, C=0.25, G=0.25, T=0.25), 
         tags=list(), matrix=matrix())
## Accessor-like methods:
\S4method{ID}{XMatrix}(x)
\S4method{bg}{XMatrix}(x)

## ... and more (see Methods)
}
\arguments{
    \item{ID}{Object of class \code{"character"}: a unique identifier for each matrix.}
    \item{name}{Object of class \code{"character"}: The name of the transcription factor. 
      In JASPAR, as far as possible, the name is based on the standardized Entrez gene symbols. 
      In the case the model describes a transcription factor hetero-dimer, two names are concatenated, such as RXR-VDR. 
      In a few cases, different splice forms of the same gene have different binding specificity: 
      in this case the splice form information is added to the name, based on the relevant literature.}
    \item{matrixClass}{Object of class \code{"character"}: Structural class of the transcription factor, based on the TFCaT system}
    \item{strand}{Object of class \code{"character"}: which strand is the binding sites sequences from.}
    \item{bg}{Object of class \code{"numeric"}: Background frequencies of the four bases. By default, it is equally 0.25.}
    \item{tags}{Object of class \code{"list"}: Some tags information about this model. 
      Tags include:
      
      (1) "family": Structural sub-class of the transcription factor, based on the TFCaT system.
      
      (2) "species": The species source for the sequences, in NCBI tax IDs.
     
      (3) "tax_group": Group of species, currently consisting of 4 larger groups: vertebrate, insect, plant, chordate.
     
      (4) "medline": a ID to the relevant publication reporting the sites used in the mode building.
     
      (5) "type": Methodology used for matrix construction.

      (6) "ACC": A representative protein accession number in Genbank for the transcription factor. 
      Human takes precedence if several exists. 
     
      (6) "pazar_tf_id": a ID to PAZAR database.
     
      (7) "TFBSshape_ID": a ID to TFBSshape database.
      
      (8) "TFencyclopedia_ID": a ID to the Transcription Factor Encyclopedia.
     
      (9) "comment": For some matrices, a curator comment is added.
    }
    \item{matrix}{Object of class \code{"matrix"}: This is the matrix information.}
    \item{x}{Object of class \code{XMatrix}.}
 %   \item{pseudocounts}{Object of class \code{"numeric"}: This is the pseudocounts used when computing ICM or PWM from PFM. 
 %     By default, a threshold of 0.8 is used based on the previous research (doi:10.1093/nar/gkn1019).}
 %   \item{schneider}{Obejct of class \code{"logical"}: this logical value indicates whether the schneider correction is used during the conversion from PFM to ICM.}
}

\section{Methods}{
  \describe{
%    \item{bg}{\code{signature(x = "XMatrix")}: ... }
%    \item{bg<-}{\code{signature(x = "XMatrix")}: ... }
    \item{ID}{
      \code{signature(x = "XMatrix")}: Gets the ID slot.
    }
%    \item{ID<-}{\code{signature(x = "XMatrix")}: ... }
    \item{length}{
      \code{signature(x = "XMatrix")}: Gets the pattern length in nucleotides (i.e. number of columns in the matrix).
    }
    \item{reverseComplement}{
      \code{signature(x = "XMatrix")}: Generates the reverse complement matrix object. Note than the strand is XMatrix will also be changed to the opossite strand.
    }
    \item{as.matrix}{
      \code{signature(x = "XMatrix")}: Returns the matrix in the XMatrix class.
    }
%    \item{total_ic}{
%      \code{signature(x = "ICMatrix")}: Returns the information content vector.
%    }
%    \item{Matrix}{\code{signature(x = "XMatrix")}: ... }
%    \item{Matrix<-}{\code{signature(x = "XMatrix")}: ... }
%    \item{matrixClass}{\code{signature(x = "XMatrix")}: ... }
%    \item{matrixClass<-}{\code{signature(x = "XMatrix")}: ... }
%    \item{name}{\code{signature(x = "XMatrix")}: ... }
%    \item{name<-}{\code{signature(x = "XMatrix")}: ... }
%    \item{revcom}{\code{signature(x = "XMatrix")}: ... }
%    \item{show}{\code{signature(object = "XMatrix")}: ... }
%    \item{strand}{\code{signature(x = "XMatrix")}: ... }
%    \item{strand<-}{\code{signature(x = "XMatrix")}: ... }
%    \item{tags}{\code{signature(x = "XMatrix")}: ... }
%    \item{update}{\code{signature(object = "XMatrix")}: ... }
	 }
}
\author{
Ge Tan
}

\examples{
  ## ----------------------------------------------------------------
  ## Constructor
  ## Note that there is no XMatrix() constructor, 
  ## but an XMatrix family of constructors: PFMatrix(), PWMatrix(), ICMatrix()
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
  
  ## Coersion
  as.matrix(pfm)
  as(pfm, "matrix")

  ## Methods
  reverseComplement(pfm)
  length(pfm)
}
\keyword{classes}
