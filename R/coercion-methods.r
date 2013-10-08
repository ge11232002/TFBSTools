### -----------------------------------------------------------------
### XMatrix Coercion
###
setAs("XMatrix", "matrix", function(from) Matrix(from))
setAs("matrix", "PFMatrix", function(from) PFMatrix(matrix=from))
setMethod("as.matrix", "PFMatrix",
          function(x){Matrix(x)}
          )


### -----------------------------------------------------------------
### XMatrixList Coercion
setAs("XMatrixList", "data.frame", function(from){})

