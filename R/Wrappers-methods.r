

### ---------------------------------------------------------------
### The real wrapper function for MEME
###
# inputFastaFn = "/Users/gtan/src/meme_4.9.1/tests/crp0.s"
run_MEME = function(inputFastaFn, binary="meme", seqtype="DNA", 
                    arguments=list()){
  valueArguments = c("-mod", "-nmotifs", "-evt", "-nsites", 
                        "-minsites", "-maxsites", "-wnsites",
                        "-w", "-minw", "-maxw", "-wg", "-ws",
                        "-bfile", "-maxiter", "-distance",
                        "-psp", "-prior", "-b", "-plib",
                        "-spfuzz", "-spmap", "-cons", "-heapsize",
                        "-maxsize", "-p", "-time", "-sf")
  booleanArguments = c("-nomatrim", "-noendgaps", "-revcomp", "-pal",
                       "-x_branch", "-w_branch")
  arguments1 = arguments[names(arguments) %in% valueArguments]
  arguments1 = paste(names(arguments), arguments, collapse=" ")
  arguments2 = arguments[names(arguments) %in% booleanArguments]
  arguments2 = paste(names(arguments2), collapse=" ")
  arguments = paste(arguments1, arguments2)
  cmd = paste(binary, inputFastaFn, "-text", 
              ifelse(seqtype=="DNA", "-dna", "-protein"), 
              arguments, "2>/dev/null")
  memeOutput = my.system(cmd, intern=TRUE)
  #conMemeOutput = pipe(cmd, open="rt")
  #on.exit(close(conMemeOutput))
  
  # get the version of meme used.
  oneLine = memeOutput[grep("^MEME version", memeOutput)]
  version = strsplit(oneLine, " ")[[1]][3]

  # get the command
  oneLine = memeOutput[grep("^command:", memeOutput)]
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
                                     sapply(splittedLines, "[", 2), "+"),
                       score=as.numeric(sapply(splittedLines, 
                                               "[", ifelse(revcomp, 4, 3)))
                       )
    motifList = c(motifList, oneRange)
  }
  motifList = GRangesList(motifList)
  names(motifList) = motifNames
  ans = list(motifList=motifList, motifEvalues=motifEvalues) 
  return(ans)
}



### -------------------------------------------------------------
### The MEME method
###
setMethod("runMEME", "character",
          function(x, binary="meme", seqtype="DNA", arguments=list(), 
                   tmpdir=tempdir()){
            seqtype = match.arg(seqtype, c("DNA", "AA"))
            ans = run_MEME(x, binary=binary, seqtype=seqtype, 
                           arguments=arguments)
            subjectSeqs = switch(seqtype,
                                 "DNA"=readDNAStringSet(filepath=x, 
                                                        format="fasta",),
                                 "AA"=readAAStringSet(filepath=x, 
                                                      format="fasta")
                                 )
            ans = MotifSet(motifList=ans[["motifList"]], 
                           motifEvalues=ans[["motifEvalues"]], 
                           subjectSeqs=subjectSeqs)
            return(ans)
          }
          )

setMethod("runMEME", "DNAStringSet",
          function(x, binary="meme", seqtype="DNA", arguments=list(), 
                   tmpdir=tempdir()){
            tmpFile = tempfile(pattern="MEME_", tmpdir=tmpdir, 
                               fileext = ".fasta")
            writeXStringSet(x, filepath=tmpFile, format="fasta")
            on.exit(unlink(tmpFile))
            ans = run_MEME(tmpFile, binary=binary, 
                           seqtype=seqtype, arguments=arguments)
            ans = MotifSet(motifList=ans[["motifList"]], 
                           motifEvalues=ans[["motifEvalues"]], subjectSeqs=x)
            return(ans)
          }
          )




