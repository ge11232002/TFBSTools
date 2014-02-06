
### ------------------------------------------------------------------------
### The "ICM" generic and methods
setMethod("toICM", "character",
          function(x, pseudocounts=0.8, schneider=FALSE,
                   bg=c(A=0.25, C=0.25, G=0.25, T=0.25)){
            dnaset = DNAStringSet(x)
            toICM(dnaset, schneider=schneider,
                  pseudocounts=pseudocounts,
                  bg=bg)
          }
          )
setMethod("toICM", "DNAStringSet",
          function(x, pseudocounts=0.8, schneider=FALSE,
                   bg=c(A=0.25, C=0.25, G=0.25, T=0.25)){
            if(!isConstant(width(x)))
              stop("'x' must be rectangular (i.e. have a constant width)")
            pfm = consensusMatrix(x)
            toICM(pfm, schneider=schneider, pseudocounts=pseudocounts,
                  bg=bg)
          }
          )
setMethod("toICM", "PFMatrix",
          function(x, pseudocounts=0.8, schneider=FALSE, bg=NULL){
            if(is.null(bg))
              bg = bg(x)
            icmMatrix = toICM(Matrix(x), pseudocounts=pseudocounts,
                              schneider=schneider, bg=bg)
            icm = ICMatrix(ID=ID(x), name=name(x), matrixClass=matrixClass(x),
                           strand=strand(x), bg=bg, 
                           tags=tags(x), matrix=icmMatrix,
                           pseudocounts=pseudocounts, schneider=schneider)
            return(icm)
          }
          )
                   

### Assumes 'x' is a Position *Frequency* Matrix (PFM) and computes the
### corresponding Information *Content* Matrix (ICM).
schneider_Hnb_precomputed = function(colsum){
  if(colsum %% 1 != 0)
    stop("the colsums must be integers")
  if(colsum < 1 || colsum > 30)
    stop("Precomputed params only available for colsums 1 to 30)")
  precomputed = c(0, 0.75, 1.11090234442608, 1.32398964833609, 
                  1.46290503577084, 
                  1.55922640783176, 1.62900374746751, 1.68128673969433,
                  1.7215504663901, 1.75328193031842, 1.77879136615189,
                  1.79965855531179, 1.81699248819687, 1.8315892710679,
                  1.84403166371213, 1.85475371994775, 1.86408383599326,
                  1.87227404728809, 1.87952034817826, 1.88597702438913,
                  1.89176691659196, 1.89698887214968, 1.90172322434865,
                  1.90603586889234, 1.90998133028897, 1.91360509239859,
                  1.91694538711761, 1.92003457997914, 1.92290025302018,
                  1.92556605820924
                  )
  return(precomputed[colsum])
}

schneider_Hnb_exact = function(colsum, bg_probabilities){
## This function is validated with the precomputed above.
  isFlat = length(unique(bg_probabilities)) == 1
  if(colsum == 1)
    return(0)
  ns = c(na=colsum, nc=0, ng=0, nt=0)
  E_Hnb = 0
  while(1){
    Pnb = factorial(colsum) / 
    (factorial(ns["na"]) * factorial(ns["nc"]) * 
     factorial(ns["ng"]) * factorial(ns["nt"])) * 
    prod(bg_probabilities^ns)
    Hnb = -1 * sum((ns/colsum * log2(ns/colsum))[is.finite(log2(ns/colsum))])
    E_Hnb = E_Hnb + Pnb * Hnb

    if(ns["nt"] != 0){
      if(ns["ng"] != 0){
        ns["ng"] = ns["ng"] - 1
        ns["nt"] = ns["nt"] + 1
      }else if(ns["nc"] != 0){
        ns["nc"] = ns["nc"] - 1
        ns["ng"] = ns["nt"] + 1
        ns["nt"] = 0
      }else if(ns["na"] != 0){
        ns["na"] = ns["na"] - 1
        ns["nc"] = ns["nt"] + 1
        ns["nt"] = 0
      }else
        break
    }else{
      if(ns["ng"] != 0){
        ns["ng"] = ns["ng"] - 1
        ns["nt"] = ns["nt"] + 1
      }else if(ns["nc"] != 0){
        ns["nc"] = ns["nc"] - 1
        ns["ng"] = ns["ng"] + 1
      }else{
        ns["na"] = ns["na"] - 1
        ns["nc"] = ns["nc"] + 1
        ns["nt"] = 0
      }
    }
  }
  return(E_Hnb)
}

schneider_Hnb_approx = function(colsum, Hg){
  ans = Hg - 3 / (2 * log(2) * colsum)
  return(ans)
}

schneider_correction = function(x, bg_probabilities){
  EXACT_SCHNEIDER_MAX = 30
  Hg = -sum(bg_probabilities * log2(bg_probabilities))
  isFlat = length(unique(bg_probabilities)) == 1
  saved_Hnb = c()
  for(colsum in unique(colSums(x))){
    if(colsum <= EXACT_SCHNEIDER_MAX){
      if(isFlat)
        Hnb = schneider_Hnb_precomputed(colsum)
      else
        Hnb = schneider_Hnb_exact(colsum, bg_probabilities)
    }else{
      Hnb = schneider_Hnb_approx(colsum, Hg)
    }
    saved_Hnb[as.character(colsum)] = Hnb
  }
  Hnbs = saved_Hnb[as.character(colSums(x))]
  return(-Hg + Hnbs)
}

setMethod("toICM", "matrix",
    ## This is validated by the TFBS perl module implemenation.
          function(x, pseudocounts=0.8, 
                   ## This is the recommended value from 
                   ## http://nar.oxfordjournals.org/content/37/3/939.long.
                   schneider=FALSE,
                   bg=c(A=0.25, C=0.25, G=0.25, T=0.25)){
            x = normargPfm(x)
            ## From here 'x' is guaranteed to have at least 1 column and to have
            ## all its columns sum to the same value.
            ## In fact, these columns sum could be different... 
            ## Modify the .normargPfm a little bit.
            bg= Biostrings:::.normargPriorParams(bg)
            #nseq = sum(x[ ,1L])
            nseq = colSums(x)
            priorN = sum(bg)
            pseudocounts = rep(0, ncol(x)) + pseudocounts
            #if(length(pseudocounts) == 1)
              #p = (x + bg_probabilities*pseudocounts) / (nseq + pseudocounts)
            #  p = sweep(x + bg_probabilities*pseudocounts, MARGIN=2, nseq + pseudocounts, "/")
            #else
              #p = (x + bg_probabilities %*% t(pseudocounts)) / (nseq + pseudocounts)
              p = sweep(x + bg %*% t(pseudocounts), MARGIN=2, 
                        nseq + priorN * pseudocounts, "/")
            D = log2(nrow(x)) + colSums(p * log2(p), na.rm=TRUE)
            #ICMMatrix = t(t(p) * D)
            ICMMatrix = sweep(p, MARGIN=2, D, "*") 
            ## This core function might be better than the operation above
            
            if(schneider){
              correntedColSums = colSums(ICMMatrix) + schneider_correction(x, bg)
              ICMMatrix = sweep(ICMMatrix, MARGIN=2, 
                                correntedColSums/colSums(ICMMatrix), "*")
            }
            return(ICMMatrix)
          }
          )

### --------------------------------------------------------------------
### Plot the seqlogo
###
#setMethod("plotLogo", "ICMatrix",
#          function(x, ic.scale = TRUE, xaxis = TRUE, yaxis = TRUE,
#                   xfontsize = 15, yfontsize = 15){
#            m = Matrix(x)
#            m = sweep(m, MARGIN=2, colSums(m), "/")
#            m = makePWM(m)
#            seqLogo(m, ic.scale = ic.scale, xaxis = xaxis, yaxis = yaxis,
#                             xfontsize = xfontsize, yfontsize = yfontsize)
#          }
#          )
setMethod("seqLogo", "ICMatrix",
          function(x, ic.scale = TRUE, xaxis = TRUE, yaxis = TRUE,
                   xfontsize = 15, yfontsize = 15){
            m = Matrix(x)
            m = sweep(m, MARGIN=2, colSums(m), "/")
            m = makePWM(m)
            seqLogo::seqLogo(m, ic.scale = ic.scale, 
                             xaxis = xaxis, yaxis = yaxis,
                             xfontsize = xfontsize, yfontsize = yfontsize)
          }
          )

### ----------------------------------------------------------------------
### Utilities methods
###
setMethod("totalIC", "ICMatrix",
          function(x) colSums(Matrix(x))
          )

