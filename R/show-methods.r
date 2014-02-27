
### -----------------------------------------------------------------------
### The "show" method.
###

setMethod("show", "XMatrix",
          function(object){
            cat("An object of class ", class(object), "\n", sep="")
            cat("ID: ", object@ID, "\n", sep="")
            cat("Name: ", object@name, "\n", sep="")
            cat("Matrix Class: ", object@matrixClass, "\n", sep="")
            cat("strand: ", object@strand, "\n", sep="")
            if(is(object, "PWMatrix")){
              cat("Pseudocounts: ", object@pseudocounts, "\n", sep="")
            }
            if(is(object, "ICMatrix")){
              cat("Schneider correction: ", object@schneider, "\n", sep="")
            }
            #for(i in 1:length(object@tags)){
            #  cat(names(object@tags)[i], ": ", object@tags[[i]], "\n", sep="")
            #}
            cat("Tags: \n")
            print(object@tags)
            cat("Background: ", "\n", sep="")
            print(object@bg)
            cat("Matrix: ", "\n", sep="")
  # add the tags later and print pretty
            print(object@profileMatrix)
          }
          )

### -----------------------------------------------------------------
### The "show" method
### Perhaps it is not a bad idea to show them in gff format.
setMethod("show", "SiteSet",
          function(object){
            cat("An object of class", class(object), "with",
                length(object), "site",
                ifelse(length(object)==1, "sequence", "sequences"))
            cat("\n")
            if(length(object) > 10000)
              object = object[1:10000]
            gff = writeGFF3(object)
            print(gff)
          }
          )


### -----------------------------------------------------------------
### The "show" method
### put them in a extended gff. any good idea?
setMethod("show", "SitePairSet",
          function(object){
            gff1 = writeGFF3(siteset1(object))
            gff2 = writeGFF3(siteset2(object))
            ans = cbind(gff1, gff2)
            cat("An object of class", class(object), "with",
                length(object), "site pair set",
                ifelse(length(object)==1, "sequence", "sequences"))
            cat("\n")
            print(ans)
          }
          )
