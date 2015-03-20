### -----------------------------------------------------------------
### XMatrix Coercion
###
setAs("XMatrix", "matrix", function(from) Matrix(from))
setAs("matrix", "PFMatrix", function(from) PFMatrix(profileMatrix=from))
setAs("matrix", "PWMatrix", function(from) PWMatrix(profileMatrix=from))
setAs("matrix", "ICMatrix", function(from) ICMatrix(profileMatrix=from))
setMethod("as.matrix", "XMatrix",
          function(x){Matrix(x)}
          )

### -----------------------------------------------------------------
### XMatrixList Coercion
setAs("XMatrixList", "data.frame", function(from){})

### -----------------------------------------------------------------
### SiteSet coercion
setAs("SiteSet", "data.frame", function(from){
      if(length(from) == 0L)
        return(data.frame())

      seqs <- DNAStringSet(views(from))
      seqs[strand(from) == "-"] <- reverseComplement(seqs[strand(from) == "-"])
      absScore <- score(from)
      relScore <- relScore(from)
      ans <- data.frame(seqnames=from@seqname,
                        source=from@sitesource,
                        feature=from@sitesource,
                        start=start(views(from)),
                        end=end(views(from)),
                        absScore=absScore,
                        relScore=relScore,
                        strand=strand(from),
                        TF=from@pattern@name,
                        class=from@pattern@matrixClass,
                        sequence=as.character(seqs)
                        )
      return(ans)
          }
)



