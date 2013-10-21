### -----------------------------------------------------------------
### XMatrix Coercion
###
setAs("XMatrix", "matrix", function(from) Matrix(from))
setAs("matrix", "PFMatrix", function(from) PFMatrix(matrix=from))
setAs("matrix", "PWMatrix", function(from) PWMatrix(matrix=from))
setAs("matrix", "ICMatrix", function(from) ICMatrix(matrix=from))
setMethod("as.matrix", "XMatrix",
          function(x){Matrix(x)}
          )

### -----------------------------------------------------------------
### XMatrixList Coercion
setAs("XMatrixList", "data.frame", function(from){})

