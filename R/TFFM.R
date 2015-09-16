### -----------------------------------------------------------------
### getTransition: get the transition probability from i to j. 
### i,j are characters.
### 
setMethod("getTransition", "TFFMFirst",
          function(tffm, i, j){
            tffm@transition[as.character(i),as.character(j)]
          })
setMethod("getTransition", "TFFMDetail",
          function(tffm, i, j){
            tffm@transition[as.character(i),as.character(j)]
          })

### -----------------------------------------------------------------
### ncol: get the length of TFFM, the number of nucleotides in the
### model excluding the background.
### Exported!
setMethod("ncol", "TFFMDetail",
          function(x){
            length(x@emission)/4L - 1L
          })
setMethod("ncol", "TFFMFirst",
          function(x){
            length(x@emission) - 2L
          })

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
            emissions <- t(tffm@transition[as.character(0:3), 
                                           as.character(0:3)]) %*% 
              c(0.25, 0.25, 0.25, 0.25)
            emissions <- drop(emissions)
            names(emissions) <- DNA_BASES
            emissions <- emissions / sum(emissions)
            return(emissions)
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
setMethod("getPosStart", "TFFMDetail",
          function(tffm){
            return(1L)
          })

### -----------------------------------------------------------------
### Get the position probablity: getPosProb
### Exported!
setMethod("getPosProb", "TFFMFirst",
          ## verified!
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

setMethod("getPosProb", "TFFMDetail",
          ## verified!
          function(tffm){
            previous_position_proba <- bgEmissionProb(tffm)
            start <- getPosStart(tffm)
            ans <- list()
            for(position in start:(start+ncol(tffm)-1L)){
              start_state <- 1:4 + (position - 1L) * 4L
              end_state <- 1:4 + position * 4L
              ans[[position]] <- t(tffm@transition[start_state, end_state]) %*%
                previous_position_proba
              previous_position_proba <- as.numeric(ans[[position]])
            }
            ans <- do.call(cbind, ans)
            rownames(ans) <- DNA_BASES
            colnames(ans) <- 1:ncol(ans)
            ans <- sweep(ans, MARGIN=2, colSums(ans), FUN="/")
            return(ans)
          })

### -----------------------------------------------------------------
### Get the emission distribution parameters at each position: getEmissionProb
### Exported!
setMethod("getEmissionProb", "TFFMFirst",
          function(tffm){
            start <- getPosStart(tffm)
            ans <- tffm@emission[start:length(tffm@emission)]
            ans <- do.call(cbind, ans)
            rownames(ans) <- rep(DNA_BASES, 4)
            colnames(ans) <- 1:ncol(ans)
            return(ans)
          })

setMethod("getEmissionProb", "TFFMDetail",
          function(tffm){
            start <- getPosStart(tffm)
            ans <- list()
            for(position in start:(start+ncol(tffm)-1L)){
              start_state <- 1:4 + (position - 1L) * 4L
              end_state <- 1:4 + position * 4L
              emissions <- tffm@transition[start_state, end_state]
              emissions <- sweep(emissions, MARGIN=1, rowSums(emissions), 
                                 FUN="/")
              ans[[position]] <- as.numeric(t(emissions))
            }
            ans <- do.call(cbind, ans)
            rownames(ans) <- rep(DNA_BASES, 4)
            colnames(ans) <- 1:ncol(ans)
            return(ans)
          })

### ----------------------------------------------------------------------
### Information content calculation at each position: totalIC 
### Exported!
setMethod("totalIC", "TFFM",
          function(x){
            pwm <- getPosProb(x)
            ans <- seqLogo:::pwm2ic(pwm)
            return(ans)
          }
          )


