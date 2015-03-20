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
      as.data.frame(as(from, "DataFrame"))
          }
)

setMethod("as.data.frame", "SiteSet",
          function(x){
            as(x, "data.frame")
          })

setAs("SiteSet", "DataFrame", function(from){
      if(length(from) == 0L)
        return(DataFrame())

      seqs <- DNAStringSet(views(from))
      seqs[strand(from) == "-"] <- reverseComplement(seqs[strand(from) == "-"])
      absScore <- score(from) 
      relScore <- relScore(from)
      ans <- DataFrame(seqnames=Rle(from@seqname),
                       source=Rle(from@sitesource),
                       feature=Rle(from@sitesource),
                       start=start(views(from)),
                       end=end(views(from)),
                       absScore=absScore,
                       relScore=relScore,
                       strand=Rle(strand(from)),
                       TF=Rle(from@pattern@name),
                       class=Rle(from@pattern@matrixClass),
                       sequence=seqs
                       )
          })


