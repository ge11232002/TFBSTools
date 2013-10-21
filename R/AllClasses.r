## To define it again, then do not need to import JASPAR2014, as JASPAR2014 is not in Bioconductor yet.
## clean it later.
setClass("JASPAR2014",
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
                 matrix="matrix")
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


### ----------------------------------------------------------------------
### The XMatrix constructor
###
PFMatrix = function(ID="Unknown", name="Unknown", matrixClass="Unknown",
                    strand="*", bg=c(A=0.25, C=0.25, G=0.25, T=0.25), tags=list(), matrix=matrix()){
  new("PFMatrix", ID=ID, name=name, matrixClass=matrixClass, strand=strand, bg=bg,
      tags=tags,
      matrix=matrix)
}

ICMatrix = function(ID="Unknown", name="Unknown", matrixClass="Unknown",
                    strand="*", bg=c(A=0.25, C=0.25, G=0.25, T=0.25), tags=list(), matrix=matrix(),
                    pseudocounts=numeric(), schneider=logical()){
  new("ICMatrix", ID=ID, name=name, matrixClass=matrixClass, strand=strand, bg=bg,
      tags=tags,
      matrix=matrix, pseudocounts=pseudocounts, schneider=schneider)
}
PWMatrix = function(ID="Unknown", name="Unknown", matrixClass="Unknown",
                    strand="*", bg=c(A=0.25, C=0.25, G=0.25, T=0.25), tags=list(), matrix=matrix(),
                    pseudocounts=numeric()){
  new("PWMatrix", ID=ID, name=name, matrixClass=matrixClass, strand=strand, bg=bg,
      tags=tags,
      matrix=matrix, pseudocounts=pseudocounts)
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
###

setMethod("XMatrixList", "list",
          function(x, use.names=TRUE, type, ...){
            ok = sapply(x, is, "XMatrix")
            if(!all(ok))
              stop("XMatrixList() only accepts XMatrix objects!")
            if(!use.names)
              names(x) = NULL
            IRanges:::newList(type, x)
          }
          )

PFMatrixList = function(..., use.names=TRUE){
  listData = list(...)
  XMatrixList(listData, use.names=use.names, type="PFMatrixList")
}

PWMatrixList = function(..., use.names=TRUE){
  listData = list(...)
  XMatrixList(listData, use.names=use.names, type="PWMatrixList")
}

ICMatrixList = function(..., use.names=TRUE){
  listData = list(...)
  XMatrixList(listData, use.names=use.names, type="ICMatrixList")
}

### ---------------------------------------------------------------------
### Site object: a nucleotide sequence feature object representing (possibly putative) transcription factor binding site. Different from TFBS perl module, here one Site object contains multiple sites.
###

setClass("Site",
         slots=c(views="XStringViews",
                 score="numeric",    # vector
                 strand="character",  ## make it Rle() later.
                 seqname="character", # length 1
                 sitesource="character", # length 1
                 primary="character",  # length 1
                 pattern="PWMatrix"   # length 1
                 )
         )

### -------------------------------------------------------------------
### The Site constructor
###
Site = function(views, score, strand="*",
                seqname="Unknown",
                sitesource="TFBS", primary="TF binding site",
                pattern){
  new("Site", views=views, seqname=seqname, score=score, strand=strand,
      sitesource=sitesource, primary=primary, pattern=pattern)
}


### --------------------------------------------------------------
### SiteList objects
###

setClass("SiteList",
         contains="SimpleList",
         representation(
                        ),
         prototype(
                   elementType="Site"
                   )
         )
### -------------------------------------------------------------
### SiteList() constructor
###
SiteList = function(..., use.names=TRUE){
  listData = list(...)
  if(is(listData[[1]], "list"))
    listData = listData[[1]]
  ok = sapply(listData, is, "Site")
  if(!all(ok))
    stop("SiteList() only accepts Site objects!")
  if(!use.names)
    names(listData) = NULL
  IRanges:::newList("SiteList", listData)
}


### -------------------------------------------------------------------
### SitePair object: a nucleotide sequence feature object representing (possibly putative) transcription factor binding site from A alignment

setClass("SitePair",
         slots=c(site1="Site",
                 site2="Site"
                 )
         )

### -----------------------------------------------------------------
### The SitePair constructor
###
SitePair = function(site1, site2){
  new("SitePair", site1=site1, site2=site2)
}


### ----------------------------------------------------------------
### SitePairList obejct: holds the list of SitePair.
###

setClass("SitePairList",
         contains="SimpleList",
         representation(
                        ),
         prototype(
                   elementType="SitePair"
                   )
         )

### ------------------------------------------------------------
### SitePairList constructor
###
SitePairList = function(..., use.names=TRUE){
  listData = list(...)
  if(is(listData[[1]], "list")) # This is pretty ugly. better solution?
    listData = listData[[1]]
  ok = sapply(listData, is, "SitePair")
  if(!all(ok))
    stop("SitePairList() only accepts SitePair objects!")
  if(!use.names)
    names(listData) = NULL
  IRanges:::newList("SitePairList", listData)
}

### ----------------------------------------------------------------
### The PairwiseAlignmentTFBS object
### Let's define our PairwiseAlignmentTFBS object, based on the PairwiseAlignments object from Biostrings. For simplicity, make it as a slot.
setClass("PairwiseAlignmentTFBS",
         slots=c(alignments="PairwiseAlignments",
                 seqname1="character",
                 seqname2="character",
                 conservation1="numeric",
                 #conservation2="numeric", # because conservation2 is never used.
                 windowSize="integer",
                 cutoff="numeric",
                 seq1length="integer",
                 seq2length="integer"
                 )
         )

### ----------------------------------------------------------------
### The constructor
###
PairwiseAlignmentTFBS = function(pattern, subject, type="global",
                                 substitutionMatrix=NULL, gapOpening=0,
                                 gapExtension=-1,
                                 seqname1="Unknown", seqname2="Unknown",
                                 windowSize=51L, cutoff=0.7){
  alignments = PairwiseAlignments(pattern, subject, type=type,
                                  substitutionMatrix=substitutionMatrix,
                                  gapOpening=gapOpening, gapExtension=gapExtension)
  conservation1 = calConservation(as.character(pattern(alignments)), as.character(subject(alignments)), windowSize=windowSize)
  seq1length = nchar(gsub("(-|_|\\.)", "", as.character(pattern(alignments))))
  seq2length = nchar(gsub("(-|_|\\.)", "", as.character(subject(alignments))))
  new("PairwiseAlignmentTFBS", alignments=alignments, seqname1=seqname1,
      seqname2=seqname2, conservation1=conservation1, windowSize=windowSize,
      cutoff=cutoff, seq1length=seq1length, seq2length=seq2length)
}

##  Is there any way to update the conservation automatically when windowSize is changed???


### ---------------------------------------------------------------
### The "show" method
### Add later... what is the pretty way?



### ---------------------------------------------------------
### Motif and MotifSet class
###
setClass("Motif",
         slot=c(
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
###
Motif = function(motif=GRanges(), motifEvalue=numeric(), subjectSeqs=DNAStringSet()){
  new("Motif", motif=motif, motifEvalue=motifEvalue, subjectSeqs=subjectSeqs)
}


MotifSet = function(motifList=GRangesList(), motifEvalues=numeric(), subjectSeqs=DNAStringSet()){
  new("MotifSet", motifList=motifList, motifEvalues=motifEvalues,
      subjectSeqs=subjectSeqs)
}




### ----------------------------------------------------------------
### The MEME object which holds the result of a MEME run
###
setClass("MEME",
         slots=c(
                 version="character",
                 command="character",
                 motifs="MotifSet"
                 )
         )

### -----------------------------------------------------------------
### The MEME object constructor
###
MEME = function(version=character(), alphabet=c("A", "C", "G", "T"),
                command=character(), motifs){
  new("MEME", version=version, alphabet=alphabet, command=command, motifs=motifs)
}



