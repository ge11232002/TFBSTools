
### -----------------------------------------------------------------------
### The "show" method.
###

setMethod("show", "XMatrix",
          function(object){
            cat("An object of class ", class(object), "\n", sep="")
            printList = list(
                             ID=ID(object),
                             Name=name(object),
                             "Matrix class"=matrixClass(object),
                             Strand=strand(object)
                             )
            if(is(object, "PWMatrix"))
              printList = c(printList, list(Pseudocounts=pseudocounts(object)))
            if(is(object, "ICMatrix"))
              printList = c(printList, list("Schneider correction"=schneider(object)))
            printList = as.data.frame(c(printList, tags(object)))
            print(printList)
            cat("Background:", bg(object), "\n")
            cat("Matrix:", "\n")
  # add the tags later and print pretty
            print(Matrix(object))
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
