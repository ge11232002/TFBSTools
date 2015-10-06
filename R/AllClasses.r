### -----------------------------------------------------------------
### JASPAR classes for each release
### Not exported!
### This is ugly hack. 
### These classes have been defined in each release of JASPAR data package.
### However, import of class from the data package will require 
### all the data packages.

setClass("JASPAR2014",
         slots=c(
                 db="character"
                )
         )

setClass("JASPAR2016",
         slots=c(
                 db="character"
                 )
         )

### -----------------------------------------------------------------
### The position frequency matrix (PWM) class
###
setClass("PFMatrix",
         slots=c(ID="character",
                 name="character",
                 matrixClass="character",
                 strand="character",
                 bg="numeric",
                 tags="list",
                 profileMatrix="matrix")
         )

setClass("PWMatrix", contains="PFMatrix",
         slots=c(pseudocounts="numeric")
         )

setClass("ICMatrix", contains="PWMatrix",
         slots=c(
                 schneider="logical"
                 )
         )

setClassUnion("XMatrix", c("PFMatrix", "ICMatrix", "PWMatrix"))


setValidity("XMatrix",
            function(object){
              ## Check the length of arguments
              if(!isConstant(c(length(object@ID), length(object@name),
                               length(object@matrixClass), 
                               length(object@strand), 1L))){
                return("The lengths of ID, name, matrixClass, strand
                       must be length 1")
              }
              ## check the strand can only be "+", "-", "*"
              if(!object@strand %in% c("+", "-", "*")){
                return("The strand can only be '+', '-' or '*'")
              }
              ## Check the bg
              if(length(object@bg) != 0L){
                if (!is.numeric(object@bg)){
                  return("'bg' must be a numeric vector")
                }
                if (length(object@bg) != length(DNA_BASES) ||
                    !identical(names(object@bg), DNA_BASES)){
                  return(("'bg' elements must be named 
                          A, C, G and T"))
                }
                if (any(is.na(object@bg)) || any(object@bg < 0)){
                  return(("'bg' contains NAs and/or negative values"))
                }
              }
              ## Check the profileMatrix
              if(!identical(dim(object@profileMatrix), c(1L,1L))){
                if (!is.matrix(object@profileMatrix) || 
                    !is.numeric(object@profileMatrix)){
                  return("The profileMatrix must be a numeric matrix")
                }
                if (!identical(rownames(object@profileMatrix), DNA_BASES)){
                  return("The profileMatrix must be the 4 DNA bases 
                         ('DNA_BASES')")
                }
                if (any(is.na(object@profileMatrix))){
                  return("The profileMatrix contains NAs")
                }
              return(TRUE)
              }
            }
            )


### ----------------------------------------------------------------------
### The XMatrix constructor
### Exported!!
PFMatrix = function(ID="Unknown", name="Unknown", matrixClass="Unknown",
                    strand="+", bg=c(A=0.25, C=0.25, G=0.25, T=0.25), 
                    tags=list(), profileMatrix=matrix()){
  new("PFMatrix", ID=ID, name=name, matrixClass=matrixClass, 
      strand=strand, bg=bg,
      tags=tags,
      profileMatrix=profileMatrix)
}

ICMatrix = function(ID="Unknown", name="Unknown", matrixClass="Unknown",
                    strand="+", bg=c(A=0.25, C=0.25, G=0.25, T=0.25), 
                    tags=list(), profileMatrix=matrix(),
                    pseudocounts=numeric(), schneider=logical()){
  new("ICMatrix", ID=ID, name=name, matrixClass=matrixClass, 
      strand=strand, bg=bg,
      tags=tags,
      profileMatrix=profileMatrix, pseudocounts=pseudocounts, 
      schneider=schneider)
}

PWMatrix = function(ID="Unknown", name="Unknown", matrixClass="Unknown",
                    strand="+", bg=c(A=0.25, C=0.25, G=0.25, T=0.25), 
                    tags=list(), profileMatrix=matrix(),
                    pseudocounts=numeric()){
  new("PWMatrix", ID=ID, name=name, matrixClass=matrixClass, 
      strand=strand, bg=bg,
      tags=tags,
      profileMatrix=profileMatrix, pseudocounts=pseudocounts)
}

### -----------------------------------------------------------------
### The set of PWM class
###
setClass("PFMatrixList",
         contains="SimpleList",
         representation(
                        ),
         prototype(
                   elementType="PFMatrix"
                   )
         )

setClass("PWMatrixList",
         contains="SimpleList",
         representation(
                        ),
         prototype(
                   elementType="PWMatrix"
                   )
         )

setClass("ICMatrixList",
         contains="PWMatrixList",
         representation(
                        ),
         prototype(
                   elementType="ICMatrix"
                   )
         )

setClassUnion("XMatrixList", c("PFMatrixList", "ICMatrixList", "PWMatrixList"))

### -----------------------------------------------------------------------
### XMatrixList() constructor.
### Exported!!

setMethod("XMatrixList", "list",
          function(x, use.names=TRUE, type, matrixClass, ...){
            ok = sapply(x, class) == matrixClass
            if(!all(ok))
              stop(type,"() only accepts ", matrixClass, " objects!")
            if(!use.names)
              names(x) = NULL
            S4Vectors:::new_SimpleList_from_list(type, x)
          }
          )

PFMatrixList = function(..., use.names=TRUE){
  listData = list(...)
  XMatrixList(listData, use.names=use.names, type="PFMatrixList", 
              matrixClass="PFMatrix")
}

PWMatrixList = function(..., use.names=TRUE){
  listData = list(...)
  XMatrixList(listData, use.names=use.names, type="PWMatrixList",
              matrixClass="PWMatrix")
}

ICMatrixList = function(..., use.names=TRUE){
  listData = list(...)
  XMatrixList(listData, use.names=use.names, type="ICMatrixList",
              matrixClass="ICMatrix")
}

### ---------------------------------------------------------------------
### SiteSet object: a nucleotide sequence feature object representing 
### (possibly putative) transcription factor binding site. 
### Different from TFBS perl module, here one 
### Site object contains multiple sites.
###

setClass("SiteSet",
         slots=c(
                 views="XStringViews", 
                 score="numeric",    # vector
                 strand="character",  ## make it Rle() later.
                 seqname="character", # length 1
                 sitesource="character", # length 1
                 primary="character",  # length 1
                 pattern="PWMatrix"   # length 1
                 )
         )

### -------------------------------------------------------------------
### The SiteSet constructor
### Exportted!
SiteSet = function(views=Views(subject=DNAString("")), 
                   score=numeric(), strand="*",
                   seqname="Unknown",
                   sitesource="TFBS", primary="TF binding site",
                   pattern=PWMatrix()){
  new("SiteSet", views=views, seqname=seqname, score=score, strand=strand,
      sitesource=sitesource, primary=primary, pattern=pattern)
}


### --------------------------------------------------------------
### SiteSetList objects
###

setClass("SiteSetList",
         contains="SimpleList",
         representation(
                        ),
         prototype(
                   elementType="SiteSet"
                   )
         )

### -------------------------------------------------------------
### SiteSetList() constructor
###
SiteSetList = function(..., use.names=TRUE){
  listData = list(...)
  ok = sapply(listData, is, "SiteSet")
  if(!all(ok))
    stop("SiteSetList() only accepts SiteSet objects!")
  if(!use.names)
    names(listData) = NULL
  S4Vectors:::new_SimpleList_from_list("SiteSetList", listData)
}


### -------------------------------------------------------------------
### SitePairSet object: a nucleotide sequence feature object 
### representing (possibly putative) transcription factor binding site 
### from A alignment
### 
setClass("SitePairSet",
         slots=c(siteset1="SiteSet",
                 siteset2="SiteSet"
                 )
         )

### -----------------------------------------------------------------
### The SitePairSet constructor
###
SitePairSet = function(siteset1, siteset2){
  new("SitePairSet", siteset1=siteset1, siteset2=siteset2)
}


### ----------------------------------------------------------------
### SitePairSetList obejct: holds the list of SitePairSet.
###

setClass("SitePairSetList",
         contains="SimpleList",
         representation(
                        ),
         prototype(
                   elementType="SitePairSet"
                   )
         )

### ------------------------------------------------------------
### SitePairSetList constructor
###
SitePairSetList = function(..., use.names=TRUE){
  listData = list(...)
  #if(is(listData[[1]], "list")) # This is pretty ugly. better solution?
  #  listData = listData[[1]]
  ok = sapply(listData, is, "SitePairSet")
  if(!all(ok))
    stop("SitePairSetList() only accepts SitePairSet objects!")
  if(!use.names)
    names(listData) = NULL
  S4Vectors:::new_SimpleList_from_list("SitePairSetList", listData)
}

### ----------------------------------------------------------------
### The PairwiseAlignmentTFBS object
### Let's define our PairwiseAlignmentTFBS object, 
### based on the PairwiseAlignments object from Biostrings. 
### For simplicity, make it as a slot.
#setClass("PairwiseAlignmentTFBS",
#         slots=c(alignments="PairwiseAlignments",
#                 seqname1="character",
#                 seqname2="character",
#                 conservation1="numeric",
#                 #conservation2="numeric", 
#                 # because conservation2 is never used.
#                 windowSize="integer",
#                 cutoff="numeric",
#                 seq1length="integer",
#                 seq2length="integer"
#                 )
#         )

### ----------------------------------------------------------------
### The constructor
### Not exported!
#PairwiseAlignmentTFBS = function(pattern, subject, type="global",
#                                 substitutionMatrix=NULL, gapOpening=0,
#                                 gapExtension=-1,
#                                 seqname1="Unknown", seqname2="Unknown",
#                                 windowSize=51L, cutoff=0.7){
#  alignments = PairwiseAlignments(pattern, subject, type=type,
#                                  substitutionMatrix=substitutionMatrix,
#                                  gapOpening=gapOpening, gapExtension=gapExtension)
#  conservation1 = calConservation(as.character(pattern(alignments)), 
#                                  as.character(subject(alignments)), 
#                                  windowSize=windowSize)
#  seq1length = nchar(gsub("(-|_|\\.)", "", as.character(pattern(alignments))))
#  seq2length = nchar(gsub("(-|_|\\.)", "", as.character(subject(alignments))))
#  new("PairwiseAlignmentTFBS", alignments=alignments, seqname1=seqname1,
#      seqname2=seqname2, conservation1=conservation1, windowSize=windowSize,
#      cutoff=cutoff, seq1length=seq1length, seq2length=seq2length)
#}

### ---------------------------------------------------------
### Motif and MotifSet class
### Exported!
setClass("Motif",
         slots=c(
                motif="GRanges",
                motifEvalue="numeric",
                subjectSeqs="DNAStringSet"
                )
         )

setClass("MotifSet",
         slots=c(
                 motifList="GRangesList",
                 motifEvalues="numeric",
                 subjectSeqs="DNAStringSet"
                 )
         )

### ---------------------------------------------------------
### The constructor function
### Exported!
Motif = function(motif=GRanges(), motifEvalue=numeric(), subjectSeqs=DNAStringSet()){
  new("Motif", motif=motif, motifEvalue=motifEvalue, subjectSeqs=subjectSeqs)
}


MotifSet = function(motifList=GRangesList(), motifEvalues=numeric(), subjectSeqs=DNAStringSet()){
  new("MotifSet", motifList=motifList, motifEvalues=motifEvalues,
      subjectSeqs=subjectSeqs)
}



### ----------------------------------------------------------------
### The MEME object which holds the result of a MEME run
### Not Exported!
setClass("MEME",
         slots=c(
                 version="character",
                 command="character",
                 motifs="MotifSet"
                 )
         )

### -----------------------------------------------------------------
### The MEME object constructor
### Not Exported!
MEME = function(version=character(), alphabet=c("A", "C", "G", "T"),
                command=character(), motifs){
  new("MEME", version=version, alphabet=alphabet, command=command, motifs=motifs)
}


### -----------------------------------------------------------------
### The transcription factor flexible model (TFFM) class
### Exported!
setClass("TFFMFirst", contains="PFMatrix",
         slots=c(type="character",
                 emission="list",
                 transition="matrix")
         )

setClass("TFFMDetail", contains="PFMatrix",
         slots=c(type="character",
                 emission="list",
                 transition="matrix"))

setClassUnion("TFFM", c("TFFMDetail", "TFFMFirst"))

### ----------------------------------------------------------------------
### The TFFM constructor
### Exported!!
TFFMFirst <- function(ID="Unknown", name="Unknown", matrixClass="Unknown",
                 strand="+", bg=c(A=0.25, C=0.25, G=0.25, T=0.25),
                 tags=list(), profileMatrix=matrix(), 
                 type=character(), emission=list(),
                 transition=matrix()){
  new("TFFMFirst", ID=ID, name=name, matrixClass=matrixClass,
      strand=strand, bg=bg,
      tags=tags,
      profileMatrix=profileMatrix,
      type=type,
      emission=emission, transition=transition)
}

TFFMDetail <- function(ID="Unknown", name="Unknown", matrixClass="Unknown",
                       strand="+", bg=c(A=0.25, C=0.25, G=0.25, T=0.25),
                       tags=list(), profileMatrix=matrix(),
                       type=character(), emission=list(),
                       transition=matrix()){
  new("TFFMDetail", ID=ID, name=name, matrixClass=matrixClass,
      strand=strand, bg=bg,
      tags=tags,
      profileMatrix=profileMatrix,
      type=type,
      emission=emission, transition=transition)
}



