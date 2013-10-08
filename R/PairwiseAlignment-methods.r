### ---------------------------------------------------------------
### The PairwiseAlignmentTFBS accessor-like methods
###
setMethod("alignments", "PairwiseAlignmentTFBS",
          function(x) x@alignments)
setMethod("seqname", "PairwiseAlignmentTFBS",
          function(x) c(x@seqname1, x@seqname2))
setMethod("conservation1", "PairwiseAlignmentTFBS",
          function(x) x@conservation1)
setMethod("seqlength", "PairwiseAlignmentTFBS",
          function(x) c(x@seq1length, x@seq2length))
setMethod("alnlength", "PairwiseAlignmentTFBS",
          function(x) nchar(alignments(x)))


calculate_conservation = function(aln1, aln2, windowSize, which="1"){
  ## This function is used to calculate the conservation profiles for a pairwise alignment.
  # x: a DNAStringSet with length 2 holds the alignment
  # windowSize: the smooth window size
  # which: which seq in the alignment is computed.
  windowSize = as.integer(windowSize)
  if(windowSize %% 2 == 0){
    warning("windows size is not even, turned into odd by -1")
    windowSize = windowSize - 1L
  }
  which = match.arg(which, c("1", "2"))
  if(nchar(aln1) != nchar(aln2))
    stop("'aln1' and 'aln2' must have the same number of characters")

  if(which == "1"){
    alignedSeq1 = aln1
    alignedSeq2 = aln2
  }else{
    alignedSeq1 = aln2
    alignedSeq2 = aln1
  }
  alignedSeq1 = strsplit(as.character(alignedSeq1), "")[[1]]
  alignedSeq2 = strsplit(as.character(alignedSeq2), "")[[1]]
  indexGap = alignedSeq1 == "-" | alignedSeq1 == "." | alignedSeq1 == "_"
  alignedSeq1 = alignedSeq1[!indexGap]
  alignedSeq2 = alignedSeq2[!indexGap]
  matches = alignedSeq1 == alignedSeq2
  conservations = runmean(matches, k=windowSize, alg="C", endrule="mean")
  return(conservations)
}

setMethod("calConservation", signature(aln1="DNAString", aln2="DNAString"),
          function(aln1, aln2, windowSize=51L, which="1"){
            calculate_conservation(as.character(aln1), as.character(aln2),
                                   windowSize=windowSize, which=which)
          }
          )
setMethod("calConservation", signature(aln1="DNAStringSet", aln2="missing"),
          function(aln1, aln2, windowSize=51L, which="1"){
            if(length(aln1) != 2)
              stop("'aln1' must be of length 2 when 'aln2' is missing")
            calculate_conservation(as.character(aln1[1]), 
                                   as.character(aln1[2]), 
                                   windowSize=windowSize, which=which)
          }
          )
setMethod("calConservation", signature(aln1="character", aln2="missing"),
          function(aln1, aln2, windowSize=51L, which="1"){
            if(length(aln1) != 2)
              stop("'aln1' must be of length 2 when 'aln2' is missing")
            calculate_conservation(aln1[1], aln1[2], windowSize=windowSize, which=which)
          }
          )
setMethod("calConservation", signature(aln1="character", aln2="character"),
          function(aln1, aln2, windowSize=51L, which="1"){
            calculate_conservation(aln1, aln2, windowSize=windowSize, which=which)
          }
          )

do_sitesearch = function(pwm, aln1, aln2, min.score, windowSize, cutoff, conservation){
# aln1, aln2: characters.
  if(nchar(aln1) != nchar(aln2))
    stop("'aln1' and 'aln2' must have the same number of characters")
  if(cutoff > 1 || cutoff < 0)
    stop("cutoff must be from 0 to 1.")
  seq1 = gsub("(-|_|\\.)", "", aln1)
  seq2 = gsub("(-|_|\\.)", "", aln2) 
  site1 = searchSeq(pwm, seq1, min.score=min.score)
  site2 = searchSeq(pwm, seq2, min.score=min.score)
  siteset1 = views(site1)
  siteset2 = views(site2)
  stopifnot(all(diff(start(siteset1)) >= 1) && all(diff(start(siteset2)) >= 1))
  # not quite sure the views returned by matchPWM is ordered by start, just check here.
  alignedSeq1 = strsplit(aln1, "")[[1]]
  alignedSeq2 = strsplit(aln2, "")[[1]]
  indexGap = alignedSeq1 == "-" | alignedSeq1 == "." | alignedSeq1 == "_"
  seq12aln = seq_len(length(alignedSeq1))[!indexGap]
  indexGap = alignedSeq2 == "-" | alignedSeq2 == "." | alignedSeq2 == "_"
  seq22aln = seq_len(length(alignedSeq2))[!indexGap]

  if(is.null(conservation))
    conservations1 = calConservation(aln1, aln2, windowSize=windowSize, which="1")
  else
    conservations1 = conservation

  pos1_in_aln = seq12aln[start(siteset1)]
  pos2_in_aln = seq22aln[start(siteset2)]
  matchedPairs = match(pos1_in_aln, pos2_in_aln)
  keep = conservations1[start(siteset1)[!is.na(matchedPairs)]] >= cutoff
  #ans_siteset1 = siteset1[(!is.na(matchedPairs))[keep]]
  ans_siteset1 = site1[(!is.na(matchedPairs))[keep]]
  #ans_siteset2 = siteset2[(na.omit(matchedPairs))[keep]]
  ans_siteset2 = site2[(na.omit(matchedPairs))[keep]]
  #return(list(siteset1=ans_siteset1, siteset2=ans_siteset2))
  return(SitePair(site1=ans_siteset1, site2=ans_siteset2))
}

#setMethod("doSiteSearch", signature(aln1="character", aln2="character"),
#          function(pwm, aln1, aln2, min.score="80%", windowSize=51L, cutoff=0.7,
#                   conservation=NULL){
#            do_sitesearch(pwm, aln1, aln2, min.score=min.score, 
#                          windowSize=windowSize, cutoff=cutoff, 
#                          conservation=conservation)
#          }
#          )
#setMethod("doSiteSearch", signature(aln1="character", aln2="missing"),
#          function(pwm, aln1, aln2, min.score="80%", windowSize=51L, cutoff=0.7,
#                   conservation=NULL){
#            if(length(aln1) != 2)
#              stop("'aln1' must be of length 2 when 'aln2' is missing")
#            do_sitesearch(pwm, aln1[1], aln1[2], min.score=min.score, 
#                          windowSize=windowSize, cutoff=cutoff, 
#                          conservation=conservation)
#          }
#          )
#setMethod("doSiteSearch", signature(aln1="DNAStringSet", aln2="missing"),
#          function(pwm, aln1, aln2, min.score="80%", windowSize=51L, cutoff=0.7,
#                   conservation=NULL){
#            if(length(aln1) != 2)
#              stop("'aln1' must be of length 2 when 'aln2' is missing")
#            do_sitesearch(pwm, as.character(aln1[1]), as.character(aln1[2]), 
#                          min.score=min.score, windowSize=windowSize, 
#                          cutoff=cutoff, conservation=conservation)
#          }
#          )
#setMethod("doSiteSearch", signature(aln1="DNAString", aln2="DNAString"),
#          function(pwm, aln1, aln2, min.score="80%", windowSize=51L, cutoff=0.7,
#                   conservation=NULL){
#            do_sitesearch(pwm, as.character(aln1), as.character(aln2),
#                             min.score=min.score, windowSize=windowSize,
#                             cutoff=cutoff, conservation=conservation)
#          }
#          )
#setMethod("doSiteSearch", signature(aln1="PairwiseAlignmentTFBS", aln2="missing"),
#          function(pwm, aln1, aln2, min.score="80%", windowSize=51L, cutoff=0.7,
#                   conservation=NULL){
#            do_sitesearch(pwm, as.character(pattern(alignments(aln1))),
#                          as.character(subject(alignments(aln1))),
#                          min.score=min.score, windowSize=windowSize(aln1),
#                          cutoff=cutoff, conservation=conservation1(aln1))
#          }
#          )



