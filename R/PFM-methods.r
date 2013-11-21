
### -----------------------------------------------------------------
### searchMatrix method. compare two position frequency matrix.
###
compareMatrix = function(pfmSubject, pfmQuery, openPenalty, extPenalty){
  # The true aligning engine. Taking two ordinary matrixs.
  pfmSubject = normargPfm(pfmSubject)
  pfmQuery = normargPfm(pfmQuery)
  ans = .Call("matrixAligner", pfmSubject, pfmQuery, openPenalty, extPenalty)
  return(ans)
}

setMethod("searchMatrix", signature(pfmSubject="matrix", pfmQuery="matrix"),
         function(pfmSubject, pfmQuery, openPenalty=3, extPenalty=0.01){
           score = compareMatrix(pfmSubject, pfmQuery, openPenalty=openPenalty,
                                 extPenalty=extPenalty)
           relScore = 100 * score / max(ncol(pfmSubject), ncol(pfmQuery)) / 2
           return(c(score=score, relScore=relScore))
         }
         )

setMethod("searchMatrix", signature(pfmSubject="PFMatrix", pfmQuery="PFMatrix"),
          function(pfmSubject, pfmQuery, openPenalty=3, extPenalty=0.01){
            ans = searchMatrix(Matrix(pfmSubject), Matrix(pfmQuery), 
                               openPenalty=openPenalty, extPenalty=extPenalty)
            return(ans)
          }
          )

setMethod("searchMatrix", signature(pfmSubject="PFMatrix", pfmQuery="matrix"),
          function(pfmSubject, pfmQuery, openPenalty=3, extPenalty=0.01){
            ans = searchMatrix(Matrix(pfmSubject), pfmQuery,
                               openPenalty=openPenalty, extPenalty=extPenalty)
            return(ans)
          }
          )

setMethod("searchMatrix", signature(pfmSubject="matrix", pfmQuery="PFMatrix"),
          function(pfmSubject, pfmQuery, openPenalty=3, extPenalty=0.01){
            ans = searchMatrix(pfmSubject, Matrix(pfmQuery),
                               openPenalty=openPenalty, extPenalty=extPenalty)
            return(ans)
          }
          )

setMethod("searchMatrix", signature(pfmSubject="PFMatrixList", pfmQuery="matrix"),
          function(pfmSubject, pfmQuery, openPenalty=3, extPenalty=0.01){
            ans = lapply(pfmSubject, searchMatrix, pfmQuery,
                            openPenalty=openPenalty, extPenalty=extPenalty)
            return(ans)
          }
          )

setMethod("searchMatrix", signature(pfmSubject="PFMatrixList", pfmQuery="PFMatrix"),
          function(pfmSubject, pfmQuery, openPenalty=3, extPenalty=0.01){
            ans = lapply(pfmSubject, searchMatrix, pfmQuery,
                         openPenalty=openPenalty, extPenalty=extPenalty)
            return(ans)
          }
          )

setMethod("searchMatrix", signature(pfmSubject="matrix", pfmQuery="character"),
          function(pfmSubject, pfmQuery, openPenalty=3, extPenalty=0.01){
            pfmQuery = strsplit(pfmQuery, "")[[1]]
            if(!all(pfmQuery %in% names(IUPAC_CODE_MAP))){
              stop("All characters must be in IUPAC_CODE_MAP!")
            }
            pfmQueryMatrix = matrix(0L, nrow=4, ncol=length(pfmQuery),
                                    dimnames=list(c("A", "C", "G", "T")))
            for(i in 1:length(pfmQuery)){
              dnaCharacters = strsplit(IUPAC_CODE_MAP[pfmQuery[i]], "")[[1]]
              pfmQueryMatrix[dnaCharacters, i] = 1L
            }
            searchMatrix(pfmSubject, pfmQueryMatrix,
                         openPenalty=openPenalty, extPenalty=extPenalty)
          }
          )

setMethod("searchMatrix", signature(pfmSubject="PFMatrix", pfmQuery="character"),
          function(pfmSubject, pfmQuery, openPenalty=3, extPenalty=0.01){
            searchMatrix(Matrix(pfmSubject), pfmQuery,
                         openPenalty=openPenalty, extPenalty=extPenalty)
          }
          )
                                    
setMethod("searchMatrix", signature(pfmSubject="PFMatrixList", pfmQuery="character"),
          function(pfmSubject, pfmQuery, openPenalty=3, extPenalty=0.01){
            ans = lapply(pfmSubject, searchMatrix, pfmQuery,
                         openPenalty=openPenalty, extPenalty=extPenalty)
            return(ans)
          }
          )


### -----------------------------------------------------------------
### permute the PFM
###
setMethod("permuteMatrix", "matrix",
          function(x, type="intra"){
            if(type == "inter")
              stop("Only permutation within matrix is available for single matrix!")
            type = match.arg(type, c("intra", "inter"))
            x = normargPfm(x)
            index = sample(seq_len(ncol(x)), ncol(x), replace=FALSE)
            x = x[ , index]
            return(x)
          }
          )

setMethod("permuteMatrix", "PFMatrix",
          function(x, type="intra"){
            Matrix(x) = permuteMatrix(Matrix(x), type=type)
            return(x)
          }
          )

setMethod("permuteMatrix", "PFMatrixList",
          function(x, type="intra"){
            type = match.arg(type, c("intra", "inter"))
            if(type == "intra"){
              for(i in seq_len(length(x))){
                x[[i]] = permuteMatrix(x[[i]])
              }
            }else if(type =="inter"){
              allMatrix = do.call(cbind, Matrix(x))
              lengths = sapply(Matrix(x), ncol)
              lengths = c(0, cumsum(lengths))
              index = sample(seq_len(ncol(allMatrix)), ncol(allMatrix), replace=FALSE)
              for(i in seq_len(length(x))){
                Matrix(x[[i]]) = allMatrix[ , index[(lengths[i]+1):lengths[i+1]]]
              }
            }
            return(x)
          }
          )


