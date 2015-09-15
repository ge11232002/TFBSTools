### -----------------------------------------------------------------
### Read meme output file
###
readMEME <- function(fn){
  # get the command
  memeOutput <- readLines(fn)
  oneLine <- memeOutput[grep("^command:", memeOutput)]
  command = gsub("^command: ", "", oneLine)

  # get the motifs information
  revcomp = grepl("revcomp", command)
  oneLines = memeOutput[grep("^MOTIF  \\d+", memeOutput)]
  splittedLines = strsplit(oneLines, "[[:blank:]]+")
  motifWidths = as.integer(sapply(splittedLines, "[", 5))
  motifOccurrences = as.integer(sapply(splittedLines, "[", 8))
  motifEvalues = as.numeric(sapply(splittedLines, "[", 14))

  # get the motifs names
  indexNames = grep("sorted by position p-value", memeOutput)
  oneLines = memeOutput[indexNames]
  motifNames = lapply(strsplit(oneLines, "[[:blank:]]+"), "[", c(2,3))
  motifNames = sapply(motifNames, paste0, collapse=" ")

  # get the motifs ranges
  motifList = list()
  for(i in seq_len(length(indexNames))){
    oneLines = memeOutput[seq(from=indexNames[i]+4,
                              to=indexNames[i]+4+motifOccurrences[i]-1)]
    splittedLines = strsplit(oneLines, "[[:blank:]]+")
    oneRange = GRanges(seqnames=sapply(splittedLines, "[", 1),
                       ranges=IRanges(start=
                                      as.integer(sapply(splittedLines,
                                                        "[",
                                                        ifelse(revcomp, 3, 2))),
                                      width=motifWidths[i]),
                       strand=ifelse(revcomp,
                                     list(sapply(splittedLines, "[", 2)), 
                                     list("+"))[[1]],
                       score=as.numeric(sapply(splittedLines,
                                               "[", ifelse(revcomp, 4, 3))),
                       leftSeqs=sapply(splittedLines, "[", 
                                       ifelse(revcomp, 5, 4)),
                       siteSeqs=sapply(splittedLines, "[",
                                       ifelse(revcomp, 6, 5)),
                       rightSeqs=sapply(splittedLines, "[",
                                        ifelse(revcomp, 7, 6))
                       )

    motifList = c(motifList, oneRange)
  }
  motifList = GRangesList(motifList)
  names(motifList) = motifNames
  ans = list(motifList=motifList, motifEvalues=motifEvalues)
  return(ans)
}


### -----------------------------------------------------------------
### Read XML file generated from TFFM python module
### Exported!
readXMLTFFM <- function(fn, type=c("First", "Detail")){
  type <- match.arg(type)
  ghmm <- xmlParse(fn)
  ghmm <- xmlToList(ghmm)
  tffmType <- unname(ghmm$HMM$.attrs["type"])
  emission <- ghmm$HMM[names(ghmm$HMM) == "state"]
  emission <- lapply(lapply(emission, "[[", "discrete"), "[[", "text")
  emission <- lapply(lapply(emission, strsplit, ","), "[[", 1)
  emission <- lapply(emission, as.numeric)

  transition <- ghmm$HMM[names(ghmm$HMM) == "transition"]
  ## Make sure the name oder is same.
  stopifnot(identical(names(transition[[1]]$.attrs), c("source", "target")))
  dims <- sapply(transition, "[[", ".attrs")
  mode(dims) <- "integer"
  dimsMatrix <- rowMax(dims) - rowMin(dims) + 1L

  transition1 <- matrix(0, ncol=dimsMatrix[2], nrow=dimsMatrix[1])
  colnames(transition1) <- rowMin(dims)[2]:rowMax(dims)[2]
  rownames(transition1) <- rowMin(dims)[1]:rowMax(dims)[1]
  for(eachTransition in transition){
    transition1[eachTransition$.attrs["source"],
                eachTransition$.attrs["target"]] <- 
                  as.numeric(eachTransition$probability)
  }

  if(type == "First"){
    ans <- TFFMFirst(type=tffmType, emission=emission, transition=transition1)
  }else if(type == "Detail"){
    ans <- TFFMDetail(type=tffmType, emission=emission, transition=transition1)
  }else{
    stop("Unsupported type of xml file!")
  }
  return(ans)
}

