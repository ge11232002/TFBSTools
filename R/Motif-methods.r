
### -------------------------------------------------------
### Methods
###
setMethod("length", "MotifSet",
          function(x) length(x@motifList)
          )

setMethod("[", "MotifSet",
          function(x, i){
            if(missing(i))
              return(x)
            ans_motifList = x@motifList[i]
            ans_motifEvalues = x@motifEvalues[i]
            ans_subjectSeqs = x@subjectSeqs
            clone(x, motifList=ans_motifList, motifEvalues=ans_motifEvalues, subjectSeqs=ans_subjectSeqs)
          }
          )

setMethod("sitesSeq", "MotifSet",
          function(x, n=10L, type="none"){
            type = match.arg(type, c("all", "left", "right", "none"))
            if(!is(n, "integer"))
              stop("n must be an integer!")
            subjectSeqsAll = x@subjectSeqs
            names(subjectSeqsAll) = sapply(strsplit(names(subjectSeqsAll), "[[:blank:]]+"), "[", 1)
            ans = list()
            for(i in seq_len(length(x))){
              oneRange = x@motifList[[i]]
              motifSeqs = mapply(subseq, subjectSeqsAll[seqnames(oneRange)],
                                 start=start(oneRange), end=end(oneRange))
              motifSeqs = sapply(motifSeqs, as.character)
              leftSeqs = ""
              if(type %in% c("all", "left")){
                leftSeqs = mapply(subseq, subjectSeqsAll[seqnames(oneRange)],
                                  start=pmax(1L, start(oneRange)-n),
                                  end=pmax(1L, start(oneRange)-1L))
                leftSeqs = tolower(sapply(leftSeqs, as.character))
                #leftSeqs = format(leftSeqs, width=n, justify="right")
              }
              rightSeqs = ""
              if(type %in% c("all", "right")){
                rightSeqs = mapply(subseq, subjectSeqsAll[seqnames(oneRange)],
                                   start=pmin(end(oneRange)+1L, width(subjectSeqsAll[seqnames(oneRange)])),
                                   end=pmin(end(oneRange)+n, width(subjectSeqsAll[seqnames(oneRange)])))
                rightSeqs = tolower(sapply(rightSeqs, as.character))
                #rightSeqs = format(rightSeqs, width=n, justify="left")
              }
              #motifSeqs = DNAStringSet(paste0(leftSeqs, motifSeqs, rightSeqs))
              #names(motifSeqs) = seqnames(oneRange)
              motifSeqs = data.frame(leftSeqs=leftSeqs, motifSeqs=motifSeqs, 
                                     rightSeqs=rightSeqs, score=oneRange$score, 
                                     strand=as.character(strand(oneRange)),
                                     stringsAsFactors=FALSE)
              ans[[names(x@motifList)[i]]] = motifSeqs
            }
            #names(ans) = names(x@motifList)
            return(ans)
          }
          )

setMethod("consensusMatrix", "MotifSet",
          function(x, as.prob=FALSE, shift=0L, width=NULL, ...){
            motifSeqs = sitesSeq(x, type="none")
            ans = lapply(motifSeqs, function(x){
                         consensusMatrix(x$motifSeqs, as.prob=as.prob, 
                                         shift=shift, width=NULL, ...)
                                     }
            )
            return(ans)
          }
          )

