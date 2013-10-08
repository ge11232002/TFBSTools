### -----------------------------------------------------------------
### Getters
### 
setMethod("name", "XMatrixList", function(x) sapply(x, name))
setMethod("ID", "XMatrixList", function(x) sapply(x, ID))
setMethod("matrixClass", "XMatrixList", function(x) sapply(x, matrixClass))
setMethod("Matrix", "XMatrixList", function(x) lapply(x, Matrix))
setMethod("strand", "XMatrixList", function(x) sapply(x, strand))
setMethod("bg", "XMatrixList", function(x) lapply(x, bg))
setMethod("matrixType", "XMatrixList", function(x) sapply(x, matrixType))
setMethod("pseudocounts", "PWMatrixList", function(x) sapply(x, pseudocounts))
setMethod("schneider", "ICMatrixList", function(x) sapply(x, schneider))
setMethod("tags", "XMatrixList", function(x) lapply(x, tags))

