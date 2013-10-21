
### -------------------------------------------------------------------
### The Site accessor-like method
###
setMethod("views", "Site", function(x) x@views)

setMethod("score", "Site", function(x) x@score)
setMethod("strand", "Site", function(x) x@strand)

setMethod("seqname", "Site", function(x) x@seqname)

setMethod("sitesource", "Site", function(x) x@sitesource)

setMethod("primary", "Site", function(x) x@primary)

setMethod("pattern", "Site", function(x) x@pattern)

setMethod("length", "Site", function(x) length(views(x)))


### -------------------------------------------------------------------
### The SitePair accessor-like method
###
setMethod("site1", "SitePair", function(x) x@site1)

setMethod("site2", "SitePair", function(x) x@site2)

setMethod("length", "SitePair", function(x) length(site1(x)))

### ------------------------------------------------------------------
### SitePair Method
###
setMethod("writeGFF3", "SitePair",
          function(x){
            if(length(x) == 0)
              return(data.frame())
            gff1 = writeGFF3(site1(x))
            gff2 = writeGFF3(site2(x))
            ans = rbind(gff1, gff2)
            return(ans)
          }
          )
setMethod("writeGFF2", "SitePair",
          function(x){
            if(length(x) == 0)
              return(data.frame())
            gff1 = writeGFF2(site1(x))
            gff2 = writeGFF2(site2(x))
            ans = rbind(gff1, gff2)
            return(ans)
          }
          )

### -----------------------------------------------------------
### The SitePairList accessor-like methods
###
setMethod("site1", "SitePairList",
          function(x){
            ans = SiteList(lapply(x, site1))
            return(ans)
          }
          )
setMethod("site2", "SitePairList",
          function(x){
            ans = SiteList(lapply(x, site2))
            return(ans)
          }
          )

### ------------------------------------------------------------------
### The getters
###
setMethod("[", "Site",
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
setMethod("c", "Site",
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
              args[arg_is_null] = NULL  # remove NULL elements by setting them to NULL!
            if (!all(sapply(args, is, class(x))))
              stop("all arguments in '...' must be ", class(x), " objects (or NULLs)")
            if(length(unique(sapply(args, slot, "seqname"))) != 1)
              stop("all arguments in '...' must have same seqname ", x@seqname, "!")
            if(length(unique(sapply(args, slot, "sitesource"))) != 1)
              stop("all arguments in '...' must have same sitesource ", x@sitesource, "!")
            if(length(unique(sapply(args, slot, "primary"))) != 1)
              stop("all arguments in '...' must have same primary ", x@primary, "!")
            if(!all(sapply(args, function(arg, x){identical(arg@pattern, x@pattern)}, x)))
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
setMethod("writeGFF3", "Site",
          function(x){
            if(length(x) == 0)
              return(data.frame())
            gff = list(seqname=seqname(x),
                       source=sitesource(x),
                       feature=sitesource(x),
                       start=start(views(x)),
                       end=end(views(x)),
                       score=score(x),
                       strand=strand(x),
                       frame=".",
                       attributes=paste(
                                        paste("TF", name(pattern(x)), sep="="),
                                        paste("class", matrixClass(pattern(x)), sep="="),
                                        paste("sequence", as.character(views(x)), sep="="),
                                        sep=";")
                       )
            gff = as.data.frame(gff)
            return(gff)
          }
          )

setMethod("writeGFF2", "Site",
          function(x){
            if(length(x) == 0)
              return(data.frame())
            gff = list(seqname=seqname(x),
                       source=sitesource(x),
                       feature=sitesource(x),
                       start=start(views(x)),
                       end=end(views(x)),
                       score=score(x),
                       strand=strand(x),
                       frame=".",
                       attributes=paste(
                                        paste("TF", paste0("\"", name(pattern(x)), "\""), sep=" "),
                                        paste("class", paste0("\"", matrixClass(pattern(x)), "\""), sep=" "),
                                        paste("sequence", paste0("\"", as.character(views(x)), "\""), sep=" "),
                                        sep="; ")
                       )
            gff = as.data.frame(gff)
            return(gff)
          }
          )

setMethod("relScore", "Site",
          function(x){
          # Luckliy, the maxScore, minScore implementation is same with TFBS perl module. Validated!
            ans = (score(x) - minScore(Matrix(pattern(x)))) / (maxScore(Matrix(pattern(x))) - minScore(Matrix(pattern(x))))
            return(ans)
          }
          )

### ----------------------------------------------------------------
### SiteList Methods
###
setMethod("writeGFF3", "SiteList",
          function(x){
            ans = do.call(rbind, lapply(x, writeGFF3))
            return(ans)
          }
          )
setMethod("writeGFF2", "SiteList",
           function(x){
             ans = do.call(rbind, lapply(x, writeGFF2))
             return(ans)
           }
           )
