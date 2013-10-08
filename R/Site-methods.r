
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

