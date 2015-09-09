### -----------------------------------------------------------------
### Read XML file generated from TFFM python module
### 
readXMLTFFM <- function(fn){
  ghmm <- xmlParse(fn)
  ghmm <- xmlToList(ghmm)
  tffmType <- unname(ghmm$HMM$.attrs["type"])
  emission <- ghmm$HMM[names(ghmm$HMM) == "state"]
  emission <- lapply(lapply(emission, "[[", "discrete"), "[[", "text")
  emission <- lapply(lapply(emission, strsplit, ","), "[[", 1)
  emission <- lapply(emission, as.numeric)

  transition <- ghmm$HMM[names(ghmm$HMM) == "transition"]
  transition <- as.numeric(sapply(transition, "[[", "probability"))

  ans <- TFFMFirst(type=tffmType, emission=emission, transition=transition)
  return(ans)
}

### -----------------------------------------------------------------
### Get the background emission probability: bgEmissionProb
###
setMethod("bgEmissionProb", signature(tffm="TFFMFirst"),
          function(tffm){
            # Retrieve emission proba for the first state which is not
            # 1st-order
            last_emissions <- tffm@emission[[1]]
            names(last_emissions) <- DNA_BASES
            # Compute emission proba for the background
            emissions <- matrix(tffm@emission[[2]], ncol=4) %*% last_emissions
            emissions <- as.numeric(emissions)
            names(emissions) <- DNA_BASES
            return(emissions)
          }
          )

### -----------------------------------------------------------------
### Get the position start of emission: getPosStart
###
setMethod("getPosStart", "TFFMFirst",
          function(tffm){
            # Give the position of the first matching state.
            # Returns: The position of the first matching state of the TFFM.
            return(3L)
          })

### -----------------------------------------------------------------
### Get the position probablity: getPosProb
###
setMethod("getPosProb", "TFFMFirst",
          function(tffm){
            previous_position_proba <- bgEmissionProb(tffm)
            start <- getPosStart(tffm)
            ans <- list()
            i <- 1L
            for(position in start:length(tffm@emission)){
              ans[[i]] <- matrix(tffm@emission[[position]], ncol=4) %*% 
                previous_position_proba
              i <- i + 1L
            }
            ans <- do.call(cbind, ans)
            rownames(ans) <- DNA_BASES
            return(ans)
          })

### -----------------------------------------------------------------
### Get the emission probability at each position: getEmissionProb
###
setMethod("getEmissionProb", "TFFMFirst",
          function(tffm){
            start <- getPosStart(tffm)
            ans <- tffm@emission[start:length(tffm@emission)]
            ans <- do.call(cbind, ans)
            rownames(ans) <- rep(DNA_BASES, 4)
            colnames(ans) <- 1:ncol(ans)
            return(ans)
          })




