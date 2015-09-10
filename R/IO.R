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



