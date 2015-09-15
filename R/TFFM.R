### -----------------------------------------------------------------
### Get the background emission probability: bgEmissionProb
###
setMethod("bgEmissionProb", "TFFMFirst",
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

setMethod("bgEmissionProb", "TFFMDetail",
          function(tffm){
            # Compute emission proba for the background

          })

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
              previous_position_proba <- as.numeric(ans[[i]])
              i <- i + 1L
            }
            ans <- do.call(cbind, ans)
            rownames(ans) <- DNA_BASES
            colnames(ans) <- 1:ncol(ans)
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




