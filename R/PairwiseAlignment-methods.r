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

do_sitesearchOneStrand = function(pwm, aln1, aln2, 
                                  seqname1="Unknown1", seqname2="Unknown2",
                                  strand, min.score, 
                                  conservation, cutoff, type="any"){
  seq1 = gsub("(-|_|\\.)", "", aln1)
  seq2 = gsub("(-|_|\\.)", "", aln2)
  site1 = searchSeq(pwm, seq1, seqname=seqname1,
                    strand=strand, min.score=min.score)
  site2 = searchSeq(pwm, seq2, seqname=seqname2,
                    strand=strand, min.score=min.score)
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
  pos1_in_aln = seq12aln[start(siteset1)]
  pos2_in_aln = seq22aln[start(siteset2)]
  matchedPairs = match(pos1_in_aln, pos2_in_aln)
  conservations1 = mapply(window, start=start(siteset1), end=end(siteset1), MoreArgs=list(conservation), SIMPLIFY=FALSE)[!is.na(matchedPairs)]
  if(type == "all"){
    keep = sapply(lapply(conservations1, ">=", cutoff), all)
  }else if(type == "any"){
    keep = sapply(lapply(conservations1, ">=", cutoff), any)
  }else{
    stop(type, " is not supported yet!")
  }
  ans_siteset1 = site1[!is.na(matchedPairs)][keep]
  ans_siteset2 = site2[na.omit(matchedPairs)][keep]
  return(list(ans_siteset1=ans_siteset1, ans_siteset2=ans_siteset2))
}

do_sitesearch = function(pwm, aln1, aln2, 
                         seqname1="Unknown1", seqname2="Unknown2", 
                         min.score, windowSize, cutoff, 
                         strand="*", type="any", conservation){
# aln1, aln2: characters.
  strand = match.arg(strand, c("+", "-", "*"))
  type = match.arg(type, c("all", "any"))
  if(nchar(aln1) != nchar(aln2))
    stop("'aln1' and 'aln2' must have the same number of characters")
  if(cutoff > 1 || cutoff < 0)
    stop("cutoff must be from 0 to 1.")
  if(is.null(conservation)){
    conservations1 = calConservation(aln1, aln2, windowSize=windowSize, which="1")
  }else{
    conservations1 = conservation
  }
  sitesetPos = NULL
  sitesetNeg = NULL
  if(strand %in% c("+", "*")){
    sitesetPos = do_sitesearchOneStrand(pwm, aln1, aln2, 
                                        seqname1=seqname1, seqname2=seqname2,
                                        strand="+", min.score=min.score, 
                                        conservation=conservations1, 
                                        cutoff=cutoff, type=type)
  }
  if(strand %in% c("-", "*")){
    sitesetNeg = do_sitesearchOneStrand(pwm, aln1, aln2, 
                                        seqname1=seqname1, seqname2=seqname2,
                                        strand="-", min.score=min.score, 
                                        conservation=conservations1, 
                                        cutoff=cutoff, type=type)
  }
  ans_siteset1 = do.call(c, list(sitesetPos$ans_siteset1, sitesetNeg$ans_siteset1))
  ans_siteset2 = do.call(c, list(sitesetPos$ans_siteset2, sitesetNeg$ans_siteset2))
  return(SitePairSet(siteset1=ans_siteset1, siteset2=ans_siteset2))
}

setMethod("doSiteSearch", signature(aln1="character", aln2="character"),
          function(pwm, aln1, aln2, min.score="80%", windowSize=51L, cutoff=0.7,
                   conservation=NULL){
            do_sitesearch(pwm, aln1, aln2, min.score=min.score, 
                          windowSize=windowSize, cutoff=cutoff, 
                          conservation=conservation)
          }
          )
setMethod("doSiteSearch", signature(aln1="character", aln2="missing"),
          function(pwm, aln1, aln2, min.score="80%", windowSize=51L, cutoff=0.7,
                   conservation=NULL){
            if(length(aln1) != 2)
              stop("'aln1' must be of length 2 when 'aln2' is missing")
            do_sitesearch(pwm, aln1[1], aln1[2], min.score=min.score, 
                          windowSize=windowSize, cutoff=cutoff, 
                          conservation=conservation)
          }
          )
setMethod("doSiteSearch", signature(aln1="DNAStringSet", aln2="missing"),
          function(pwm, aln1, aln2, min.score="80%", windowSize=51L, cutoff=0.7,
                   conservation=NULL){
            if(length(aln1) != 2)
              stop("'aln1' must be of length 2 when 'aln2' is missing")
            do_sitesearch(pwm, as.character(aln1[1]), as.character(aln1[2]), 
                          min.score=min.score, windowSize=windowSize, 
                          cutoff=cutoff, conservation=conservation)
          }
          )
setMethod("doSiteSearch", signature(aln1="DNAString", aln2="DNAString"),
          function(pwm, aln1, aln2, min.score="80%", windowSize=51L, cutoff=0.7,
                   conservation=NULL){
            do_sitesearch(pwm, as.character(aln1), as.character(aln2),
                             min.score=min.score, windowSize=windowSize,
                             cutoff=cutoff, conservation=conservation)
          }
          )
setMethod("doSiteSearch", signature(aln1="PairwiseAlignmentTFBS", aln2="missing"),
          function(pwm, aln1, aln2, min.score="80%", windowSize=51L, cutoff=0.7,
                   conservation=NULL){
            do_sitesearch(pwm, as.character(pattern(alignments(aln1))),
                          as.character(subject(alignments(aln1))),
                          min.score=min.score, windowSize=windowSize(aln1),
                          cutoff=cutoff, conservation=conservation1(aln1))
          }
          )

do_PairBSgenomeSearchPositive = function(pwm, BSgenome1, BSgenome2, chr1, chr2, 
                                 min.score, chain){
  ## I know this is really stupid, but I am almost confused by the strand. So split ito positive and negative.
  ## search with Positive pwm
  ## BSgenome1, BSgenome2 are BSgenome object
  ## chr1, chr2 character
  seq1 = getSeq(BSgenome1, chr1)
  seq2 = getSeq(BSgenome2, chr2)
  site1 = suppressWarnings(searchSeq(pwm, seq1, seqname=chr1, strand="+", min.score=min.score))
  rm(seq1)
  site1GRanges = GRanges(seqnames=chr1, ranges(site1@views), strand="+")
  # we only care about the coordinate based on positive strand, only this coordinate is return by searchSeq.
  site2GRanges = liftOver(site1GRanges, chain)
  # reduce the ranges. can apply on a GRangesList!! Cool!
  site2GRanges = reduce(site2GRanges)
  lengths = sapply(site2GRanges, length)
  site2GRanges = site2GRanges[lengths == 1L] # so far, we drop the region with more ranges. Discuss with Boris for more details.
  # only keep the ranges on chr2
  site2GRanges = site2GRanges[as.character(seqnames(site2GRanges)) == chr2]
  site1 = site1[as.integer(names(site2GRanges))]
  # extend the ranges a bit. Let's use ncol of matrix
  site2GRanges2 = GRanges(seqnames=as.character(seqnames(site2GRanges)), 
                         ranges=IRanges(as.integer(start(site2GRanges)) - ncol(pwm@matrix), 
                                        as.integer(end(site2GRanges)) + ncol(pwm@matrix)
                                        ),
                         strand=as.character(strand(site2GRanges))
                         )
  site2SeqsSet = getSeq(BSgenome2, site2GRanges2)
  site2 = lapply(site2SeqsSet, function(seq1, pwm, min.score){
               searchSeq(pwm, seq1, strand="+", min.score=min.score)
                         }, pwm, min.score)
  lengths = sapply(site2, length)
  site2 = site2[lengths > 0L]
  site1 = site1[lengths > 0L]
  if(length(site2) == 0L){
    site2 = site1
    return(site1=site1, site2=site2)
  }
  site2 = do.call(c, site2)
  site2GRanges2 = site2GRanges2[lengths > 0L]
  # correct the ranges in site2 for negative strand.
  lengthChr2 = length(seq2)
  indexNegative = as.logical(strand(site2GRanges2) == "-")
  ranges(site2GRanges2)[indexNegative] = IRanges(start=lengthChr2-end(site2GRanges2)[indexNegative]+1,
                                                 end=lengthChr2-start(site2GRanges2)[indexNegative]+1
                                                 )
  # build a new site2 with chr2 as subject and new ranges.
  ans_site2 = SiteSet(
                      views=Views(subject=seq2,
                                  start=start(views(site2)) + start(site2GRanges2) - 1,
                                  end=end(views(site2)) + start(site2GRanges2) - 1
                                  ),
                      seqname=chr2,
                      score=score(site2),
                      strand=as.character(strand(site2GRanges2)),
                      sitesource="TFBS", primary="TF binding site",
                      pattern=pwm
                      )
  # correct the strand in site2, some are "-" after liftover
  return(list(site1=site1, site2=ans_site2))
}

do_PairBSgenomeSearchNegative = function(pwm, BSgenome1, BSgenome2, chr1, chr2,
                                         min.score, chain){
  ## deal with the negative strand
  seq1 = getSeq(BSgenome1, chr1)
  seq2 = getSeq(BSgenome2, chr2)
  site1 = suppressWarnings(searchSeq(pwm, seq1, seqname=chr1, strand="-", min.score=min.score))
  rm(seq1)
  site1GRanges = GRanges(seqnames=chr1, ranges(site1@views), strand="+")
  site2GRanges = liftOver(site1GRanges, chain)
  site2GRanges = reduce(site2GRanges)
  lengths = sapply(site2GRanges, length)
  site2GRanges = site2GRanges[lengths == 1L]
  site2GRanges = site2GRanges[as.character(seqnames(site2GRanges)) == chr2]
  site1 = site1[as.integer(names(site2GRanges))]
  site2GRanges2 = GRanges(seqnames=as.character(seqnames(site2GRanges)),
                         ranges=IRanges(as.integer(start(site2GRanges)) - ncol(pwm@matrix),
                                        as.integer(end(site2GRanges)) + ncol(pwm@matrix)
                                        ),
                         strand=as.character(strand(site2GRanges))
                         )
  site2SeqsSet = getSeq(BSgenome2, site2GRanges2)
  site2 = lapply(site2SeqsSet, function(seq1, pwm, min.score){
               searchSeq(pwm, seq1, strand="-", min.score=min.score)
                         }, pwm, min.score)
  lengths = sapply(site2, length)
  site2 = site2[lengths > 0L]
  site1 = site1[lengths > 0L]
  if(length(site2) == 0L){
    site2 = site1
    return(site1=site1, site2=site2)
  }
  site2 = do.call(c, site2)
  site2GRanges2 = site2GRanges2[lengths > 0L]
  # correct the ranges in site2 for negative strand.
  lengthChr2 = length(seq2)
  indexNegative = as.logical(strand(site2GRanges2) == "-")
  ranges(site2GRanges2)[indexNegative] = IRanges(start=lengthChr2-end(site2GRanges2)[indexNegative]+1,
                                                 end=lengthChr2-start(site2GRanges2)[indexNegative]+1
                                                 )
  # build a new site2 with chr2 as subject and new ranges.
  ans_site2 = SiteSet(
                      views=Views(subject=seq2,
                                  start=start(views(site2)) + start(site2GRanges2) - 1,
                                  end=end(views(site2)) + start(site2GRanges2) - 1
                                  ),
                      seqname=chr2,
                      score=score(site2),
                      strand=chartr("+-", "-+", as.character(strand(site2GRanges2))),
                      sitesource="TFBS", primary="TF binding site",
                      pattern=pwm
                      )
  return(list(site1=site1, site2=ans_site2))
}

do_PairBSgenomeSearch = function(pwm, BSgenome1, BSgenome2, chr1, chr2,
                                          strand, min.score, chain){
  strand = match.arg(strand, c("+", "-", "*"))
  sitesetPos = NULL
  sitesetNeg = NULL
  if(strand %in% c("+", "*")){
    sitesetPos = do_PairBSgenomeSearchPositive(pwm, BSgenome1, BSgenome2, chr1, chr2, min.score, chain)
  }
  if(strand %in% c("-", "*")){
    sitesetNeg = do_PairBSgenomeSearchNegative(pwm, BSgenome1, BSgenome2, chr1, chr2, min.score, chain)
  }
  ans_siteset1 = do.call(c, list(sitesetPos$site1, sitesetNeg$site1))
  ans_siteset2 = do.call(c, list(sitesetPos$site2, sitesetNeg$site2))
  return(SitePairSet(siteset1=ans_siteset1, siteset2=ans_siteset2))

}

