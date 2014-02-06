
### -------------------------------------------------------------------
### The SiteSet accessor-like method
###
setMethod("views", "SiteSet", function(x) x@views)

setMethod("score", "SiteSet", function(x) x@score)
setMethod("strand", "SiteSet", function(x) x@strand)

setMethod("seqname", "SiteSet", function(x) x@seqname)

setMethod("sitesource", "SiteSet", function(x) x@sitesource)

setMethod("primary", "SiteSet", function(x) x@primary)

setMethod("pattern", "SiteSet", function(x) x@pattern)

setMethod("length", "SiteSet", function(x) length(views(x)))


### -------------------------------------------------------------------
### The SitePairSet accessor-like method
###
setMethod("siteset1", "SitePairSet", function(x) x@siteset1)

setMethod("siteset2", "SitePairSet", function(x) x@siteset2)

setMethod("length", "SitePairSet", function(x) length(siteset1(x)))

### ------------------------------------------------------------------
### SitePairSet Method
###
setMethod("writeGFF3", "SitePairSet",
          function(x, scoreType=c("absolute", "relative")){
            if(length(x) == 0)
              return(data.frame())
            gff1 = writeGFF3(siteset1(x), scoreType=scoreType)
            gff2 = writeGFF3(siteset2(x), scoreType=scoreType)
            ans = rbind(gff1, gff2)
            return(ans)
          }
          )
setMethod("writeGFF2", "SitePairSet",
          function(x, scoreType=c("absolute", "relative")){
            if(length(x) == 0)
              return(data.frame())
            gff1 = writeGFF2(siteset1(x), scoreType=scoreType)
            gff2 = writeGFF2(siteset2(x), scoreType=scoreType)
            ans = rbind(gff1, gff2)
            return(ans)
          }
          )

### -----------------------------------------------------------
### The SitePairSetList accessor-like methods
###
setMethod("siteset1", "SitePairSetList",
          function(x){
            #ans = SiteSetList(lapply(x, siteset1))
            ans = do.call(SiteSetList, lapply(x, siteset1))
            return(ans)
          }
          )
setMethod("siteset2", "SitePairSetList",
          function(x){
            #ans = SiteSetList(lapply(x, siteset2))
            ans = do.call(SiteSetList, lapply(x, siteset2))
            return(ans)
          }
          )

### ------------------------------------------------------------------
### The getters
###
setMethod("[", "SiteSet",
          function(x, i){
            if(missing(i))
              return(x)
            ans_views = views(x)[i]
            ans_score = score(x)[i]
            ans_strand = strand(x)[i]
            clone(x, views=ans_views, score=ans_score, strand=ans_strand)
          }
          )
### -----------------------------------------------------------------
### Combining
###
setMethod("c", "SiteSet",
          function(x, ...){
            if(missing(x)){
              args = unname(list(...))
              x = args[[1L]]
            }else{
              args = unname(list(x, ...))
            }
            if (length(args) == 1L)
              return(x)
            arg_is_null = sapply(args, is.null)
            if (any(arg_is_null))
              args[arg_is_null] = NULL  
            # remove NULL elements by setting them to NULL!
            if (!all(sapply(args, is, class(x))))
              stop("all arguments in '...' must be ", class(x), 
                   " objects (or NULLs)")
            if(length(unique(sapply(args, slot, "seqname"))) != 1)
              stop("all arguments in '...' must have same seqname ", 
                   x@seqname, "!")
            if(length(unique(sapply(args, slot, "sitesource"))) != 1)
              stop("all arguments in '...' must have same sitesource ", 
                   x@sitesource, "!")
            if(length(unique(sapply(args, slot, "primary"))) != 1)
              stop("all arguments in '...' must have same primary ", 
                   x@primary, "!")
            if(!all(sapply(args, 
                           function(arg, x){
                             identical(arg@pattern, x@pattern)}, x)))
              stop("all arguments in '...' must have same pattern matrix!")
            new_start = unlist(lapply(lapply(args, slot, "views"), start))
            new_end = unlist(lapply(lapply(args, slot, "views"), end))
            new_score = unlist(lapply(args, slot, "score"))
            new_strand = unlist(lapply(args, slot, "strand"))
            ans = update(x, views=Views(subject=subject(x@views), 
                                        start=new_start,
                                        end=new_end),
                         score=new_score, strand=new_strand, seqname=x@seqname,
                         sitesource=x@sitesource, primary=x@primary,
                         pattern=x@pattern
                         )
            ## validObject(ans)
            ## do it when we have setValidity
            return(ans)
          }
          )

### -----------------------------------------------------------------
### Methods
###
setMethod("writeGFF3", "SiteSet",
          function(x, scoreType=c("absolute", "relative")){
            if(length(x) == 0)
              return(data.frame())
            scoreType = match.arg(scoreType, c("absolute", "relative"))
            seqs = DNAStringSet(views(x))
            seqs[strand(x) == "-"] = reverseComplement(seqs[strand(x) == "-"])
            if(scoreType =="absolute"){
              score = score(x)
            }else{
              score = relScore(x)
            }
            gff = list(seqname=seqname(x),
                       source=sitesource(x),
                       feature=sitesource(x),
                       start=start(views(x)),
                       end=end(views(x)),
                       score=score,
                       strand=strand(x),
                       frame=".",
                       attributes=
                        paste(
                              paste("TF", name(pattern(x)), sep="="),
                              paste("class", matrixClass(pattern(x)), sep="="),
                              paste("sequence", seqs, sep="="),
                              sep=";")
                       )
            gff = as.data.frame(gff)
            return(gff)
          }
          )

setMethod("writeGFF2", "SiteSet",
          function(x, scoreType=c("absolute", "relative")){
            if(length(x) == 0)
              return(data.frame())
            scoreType = match.arg(scoreType, c("absolute", "relative"))
            seqs = DNAStringSet(views(x))
            seqs[strand(x) == "-"] = reverseComplement(seqs[strand(x) == "-"])
            if(scoreType == "absolute"){
              score = score(x)
            }else{
              score = relScore(x)
            }
            gff = list(seqname=seqname(x),
                       source=sitesource(x),
                       feature=sitesource(x),
                       start=start(views(x)),
                       end=end(views(x)),
                       score=score,
                       strand=strand(x),
                       frame=".",
                       attributes=
                        paste(
                              paste("TF", 
                                    paste0("\"", name(pattern(x)), "\""), 
                                    sep=" "),
                              paste("class", 
                                    paste0("\"", matrixClass(pattern(x)), "\""),
                                    sep=" "),
                              paste("sequence", paste0("\"", seqs, "\""), 
                                    sep=" "),
                              sep="; ")
                       )
            gff = as.data.frame(gff)
            return(gff)
          }
          )

setMethod("relScore", "SiteSet",
          function(x){
            ans = (score(x) - minScore(Matrix(pattern(x)))) / 
              (maxScore(Matrix(pattern(x))) - minScore(Matrix(pattern(x))))
            return(ans)
          }
          )

setMethod("relScore", "SiteSetList",
          function(x){
            lapply(x, relScore)
          }
          )

### ----------------------------------------------------------------
### SiteSetList Methods
###
setMethod("writeGFF3", "SiteSetList",
          function(x){
            ans = do.call(rbind, lapply(x, writeGFF3))
            return(ans)
          }
          )
setMethod("writeGFF2", "SiteSetList",
           function(x){
             ans = do.call(rbind, lapply(x, writeGFF2))
             return(ans)
           }
           )

