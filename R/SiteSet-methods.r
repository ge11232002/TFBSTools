
### -------------------------------------------------------------------
### The SiteSet accessor-like method
###
setMethod("views", "SiteSet", function(x) x@views)
setMethod("start", "SiteSet", function(x) start(x@views))
setMethod("end",   "SiteSet", function(x) end(x@views))

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
            if(length(x) == 0L)
              return(data.frame())
            gff1 = writeGFF3(siteset1(x), scoreType=scoreType)
            gff2 = writeGFF3(siteset2(x), scoreType=scoreType)
            ans = rbind(gff1, gff2)
            return(ans)
          }
          )

setMethod("writeGFF2", "SitePairSet",
          function(x, scoreType=c("absolute", "relative")){
            if(length(x) == 0L)
              return(data.frame())
            gff1 = writeGFF2(siteset1(x), scoreType=scoreType)
            gff2 = writeGFF2(siteset2(x), scoreType=scoreType)
            ans = rbind(gff1, gff2)
            return(ans)
          }
          )

setMethod("writeGFF3", "SitePairSetList",
          function(x, scoreType=c("absolute", "relative")){
            ans = lapply(x, writeGFF3, scoreType)
            ans = do.call(rbind, ans)
            return(ans)
          }
          )

setMethod("writeGFF2", "SitePairSetList",
          function(x, scoreType=c("absolute", "relative")){
            ans = lapply(x, writeGFF2, scoreType)
            ans = do.call(rbind, ans)
            return(ans)
          }
          )

### -----------------------------------------------------------
### The SitePairSetList accessor-like methods
### Exported!
setMethod("siteset1", "SitePairSetList",
          function(x){
            ans <- do.call(SiteSetList, lapply(x, siteset1))
            return(ans)
          }
          )

setMethod("siteset2", "SitePairSetList",
          function(x){
            ans <- do.call(SiteSetList, lapply(x, siteset2))
            return(ans)
          }
          )

setMethod("seqname1", "SitePairSetList",
          function(x){
            ans <- sapply(x, function(x){x@siteset1@seqname})
            return(ans)
          }
          )

setMethod("seqname2", "SitePairSetList",
          function(x){
            ans <- sapply(x, function(x){x@siteset2@seqname})
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
### Can only apply to SiteSet based on same seq, with same pattern matrix.
setMethod("c", "SiteSet",
          function(x, ...){
            if(missing(x) || length(x) == 0L){
              args = unname(list(...))
              x = args[[1L]]
            }else{
              args = unname(list(x, ...))
            }
            if (length(args) == 1L)
              return(x)
            arg_is_null = sapply(args, is.null)
            if (any(arg_is_null)){
              args[arg_is_null] = NULL  
            }
            # remove NULL elements by setting them to NULL!
            if (!all(sapply(args, is, class(x)))){
              stop("all arguments in '...' must be ", class(x), 
                   " objects (or NULLs)")
            }
            if(length(unique(sapply(args, slot, "seqname"))) != 1){
              stop("all arguments in '...' must have same seqname ", 
                   x@seqname, "!")
            }
            if(length(unique(sapply(args, slot, "sitesource"))) != 1){
              stop("all arguments in '...' must have same sitesource ", 
                   x@sitesource, "!")
            }
            if(length(unique(sapply(args, slot, "primary"))) != 1){
              stop("all arguments in '...' must have same primary ", 
                   x@primary, "!")
            }

            if(!all(sapply(args, 
                           function(arg, x){
                             identical(arg@pattern, x@pattern)}, x))){
              stop("all arguments in '...' must have same pattern matrix!")
                           }
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
            if(length(x) == 0L)
              return(data.frame())
            scoreType = match.arg(scoreType, c("absolute", "relative"))
            seqs = DNAStringSet(views(x))
            seqs[strand(x) == "-"] <- 
              reverseComplement(seqs[strand(x) == "-"])
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
            if(length(x) == 0L)
              return(data.frame())
            scoreType = match.arg(scoreType, c("absolute", "relative"))
            seqs = DNAStringSet(views(x))
            seqs[strand(x) == "-"] <- 
              reverseComplement(seqs[strand(x) == "-"])
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
### Exported!
setMethod("writeGFF3", "SiteSetList",
          function(x, scoreType=c("absolute", "relative")){
            ans = do.call(rbind, lapply(x, writeGFF3, scoreType=scoreType))
            return(ans)
          }
          )
setMethod("writeGFF2", "SiteSetList",
           function(x, scoreType=c("absolute", "relative")){
             ans = do.call(rbind, lapply(x, writeGFF2, scoreType=scoreType))
             return(ans)
           }
           )

### -----------------------------------------------------------------
### Get the empirical p-values of the scores.
### Exported!
setMethod("pvalues", "SiteSet",
          function(x, type=c("TFMPvalue", "sampling")){
            pwm = x@pattern@profileMatrix
            bg = x@pattern@bg
            type <- match.arg(type)
            if(type == "TFMPvalue"){
              pvalues <- sapply(x@score, function(x, mat, bg){
                                TFMsc2pv(mat, score=x, bg=bg, type="PWM")
                       }, pwm, bg)
              return(pvalues)
            }else if(type == "sampling"){
              allScores = replicate(1e4, 
                                    sum(apply(pwm, 2, sample, 
                                              size=1, prob=bg)))
              pvalues = sapply(x@score, function(x){
                             sum(x < allScores)/length(allScores)
                                    }
            )
            return(pvalues)
            }
          }
          )

setMethod("pvalues", "SiteSetList",
          function(x, type=c("TFMPvalue", "sampling")){
            ans = lapply(x, pvalues, type)
            return(ans)
          }
          )

#setMethod("pvalues", "SitePairSet",
#          function(x){
#            list(pvalues(siteset1(x)), pvalues(siteset2(x)))
#          }
#          )


### -----------------------------------------------------------------
### get the genomic coordinates in the view of SitePairSet 
### from searchAln for Axt.
### Exported!
setMethod("toGRangesList",
          signature(x="SitePairSetList", axt="Axt"),
          function(x, axt){
            if(length(axt) != length(x)){
              stop("The length of SitePairSetList must be equal to
                   the length of axt object")
            }
            #if(!all(seqname1(x) %in% seqnames(targetBSgenome))){
            #  stop("The target seqnames of SitePairSetList are not subset of 
            #       the seqnames of target BSgenome")
            #}
            #if(!all(seqname2(x) %in% seqnames(queryBSgenome))){
            #  stop("The query seqnames of SitePairSetList are not subset of
            #       the seqnames of query BSgenome")
            #}
            indexNoneZero <- which(sapply(x, length) != 0L)
            x <- x[indexNoneZero]
            axt <- axt[indexNoneZero]
            eachLengths <- lengths(x)
            targetTFBS <- 
              GRanges(seqnames=rep(sapply(x, 
                                     function(x){x@siteset1@seqname}),
                                                eachLengths),
                      ranges=IRanges(start=
                        unlist(lapply(x, 
                          function(x){start(x@siteset1@views)})) + 
                                     rep(start(targetRanges(axt)), 
                                         eachLengths) - 1,
                                     end=
                        unlist(lapply(x,
                          function(x){end(x@siteset1@views)})) +
                                     rep(start(targetRanges(axt)),
                                         eachLengths) - 1
                                     ),
                      strand="+",
                      matrix.ID=rep(sapply(x, 
                                      function(x){x@siteset1@pattern@ID}),
                                    eachLengths),
                      matrix.strand=unlist(lapply(x, 
                                      function(x){x@siteset1@strand})),
                      abs.score=unlist(lapply(x,
                                  function(x){x@siteset1@score})),
                      rel.score=unlist(lapply(x, 
                                  function(x){relScore(x@siteset1)})),
                      sitesSeq=DNAStringSet(unlist(lapply(x, 
                                 function(x){as.character(x@siteset1@views)})))
                      )
            queryTFBS <- 
              GRanges(seqnames=rep(sapply(x,
                                     function(x){x@siteset2@seqname}),
                                                eachLengths),
                      ranges=IRanges(start=
                        unlist(lapply(x,
                          function(x){start(x@siteset2@views)})) +
                                     rep(start(queryRanges(axt)),
                                         eachLengths) - 1,
                                     end=
                        unlist(lapply(x,
                          function(x){end(x@siteset2@views)})) +
                                     rep(start(queryRanges(axt)),
                                         eachLengths) - 1
                                     ),
                      strand=rep(strand(queryRanges(axt)), eachLengths),
                      matrix.ID=rep(sapply(x,
                                      function(x){x@siteset2@pattern@ID}),
                                    eachLengths),
                      matrix.strand=unlist(lapply(x,
                                      function(x){x@siteset2@strand})),
                      abs.score=unlist(lapply(x,
                                      function(x){x@siteset2@score})),
                      rel.score=unlist(lapply(x,
                                      function(x){relScore(x@siteset2)})),
                      sitesSeq=DNAStringSet(unlist(lapply(x, 
                                 function(x){as.character(x@siteset2@views)})))
                      )
            return(GRangesList(targetTFBS=targetTFBS, queryTFBS=queryTFBS))
          }
          )