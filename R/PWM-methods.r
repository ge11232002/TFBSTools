
### ------------------------------------------------------------------------
### The "PWM" generic and methods. This is a bit different from the implementation of Biostrings.
setMethod("toPWM", "character",
          function(x, type="log2probratio", pseudocounts=0.8, 
                   bg=c(A=0.25, C=0.25, G=0.25, T=0.25)){
            dnaset = DNAStringSet(x)
            toPWM(dnaset, pseudocounts=pseudocounts,
                  bg=bg)
          }
          )
setMethod("toPWM", "DNAStringSet",
          function(x, type="log2probratio", pseudocounts=0.8,
                   bg=c(A=0.25, C=0.25, G=0.25, T=0.25)){
            if(!isConstant(width(x)))
              stop("'x' must be rectangular (i.e. have a constant width)")
            pfm = consensusMatrix(x)
            toPWM(pfm, pseudocounts=pseudocounts,
                  bg=bg)
          }
          )
setMethod("toPWM", "PFMatrix",
          function(x, type="log2probratio", pseudocounts=0.8, bg=NULL){
            if(is.null(bg))
              bg = bg(x)
            pwmMatrix = toPWM(Matrix(x), pseudocounts=pseudocounts,
                              bg=bg)
            pwm = PWMatrix(ID=ID(x), name=name(x), matrixClass=matrixClass(x),
                           strand=strand(x), bg=bg, 
                           tags=tags(x), matrix=pwmMatrix,
                           pseudocounts=pseudocounts)
          }
          )

### Assumes 'x' is a Position *Frequency* Matrix (PFM) and computes the
### corresponding Position *Weight* Matrix (PWM).
setMethod("toPWM", "matrix",
    ## This is validated by the TFBS perl module version.
          function(x, type="log2probratio", pseudocounts=0.8,
                   bg=c(A=0.25, C=0.25, G=0.25, T=0.25)){
            #x = Biostrings:::.normargPfm(x)
            bg = Biostrings:::.normargPriorParams(bg)
            type = match.arg(type, c("log2probratio", "prob"))
            nseq = colSums(x)
            priorN = sum(bg)
            pseudocounts = rep(0, ncol(x)) + pseudocounts
            #if(length(pseudocounts) == 1)
            #  p = sweep(x + bg*pseudocounts, MARGIN=2, nseq + pseudocounts, "/")
              #p = (x + bg_probabilities*pseudocounts) / (nseq + pseudocounts)
            #else
              #p = (x + bg_probabilities %*% t(pseudocounts)) / (nseq + pseudocounts)
              p = sweep(x + bg %*% t(pseudocounts), MARGIN=2, nseq + priorN * pseudocounts, "/")
            if(type == "prob")
              return(p)
            prior.probs = bg / priorN
            #ans = log2(p / prior.probs)
            #Here ans's colSums is 1s. Need to be adapted for seq logo maybe later.
            ans = log2(sweep(p, MARGIN=1, prior.probs, "/"))
            return(ans)
          }
          )
### ---------------------------------------------------------------------
### searchSeq: scans a nucleotide sequence with the pattern represented by the PWM
### Currently we make it as a normal function. Is it necessary to make it a setMethod? Yes. It's necessary to make it a setMethod.
setMethod("searchSeq", "PWMatrix",
# scans a nucleotide sequence with the pattern represented by the PWM.
          function(x, subject, seqname="Unknown", strand="*", min.score="80%"){
            ans_views = matchPWM(unitScale(Matrix(x)), subject, min.score=min.score)
            #score = rep(0, length(ans_views)) # fix the score issue....
            score = PWMscoreStartingAt(unitScale(Matrix(x)), subject(ans_views),
                                       start(ans_views))
            # The score here from PWMscoreStartingAt is the unitscaled score. Let's make it into original one, synced with TFBS module. This is validated!
            score = score * (maxScore(Matrix(x)) - minScore(Matrix(x))) + minScore(Matrix(x))
            stopifnot(strand %in% c("+", "-", "*")) # need to ask Boris strand.
            if(length(strand) == 1)
              strand = rep(strand, length(ans_views))
            stopifnot(length(strand) == length(ans_views))
            ans_site = Site(views=ans_views, seqname=seqname,
                            score=score, strand=strand, 
                            sitesource="TFBS", primary="TF binding site",
                            pattern=x
                            )
          }
          )

setMethod("searchSeq", "PWMatrixList",
# scans a nucleotide sequence with all patterns represented stored in $matrixset;
          function(x, subject, seqname="Unknown", strand="*", min.score="80%"){
            #pwms = lapply(Matrix(x), unitScale)
            #ans = lapply(pwms, matchPWM, subject, min.score)
            ans_list = lapply(x, searchSeq, subject=subject, seqname=seqname, 
                              strand=strand, min.score=min.score)
            ans = SiteList(ans_list)
            return(ans)
          }
          )

### ----------------------------------------------------------------------
### searchAln: Scans a pairwise alignment of nucleotide sequences with the pattern represented by the PWM: it reports only those hits that are present in equivalent positions of both sequences and exceed a specified threshold score in both, AND are found in regions of the alignment above the specified
## Should have a better way for this duplicated code..

setMethod("searchAln", signature(pwm="PWMatrixList", aln1="character", aln2="character"),
          function(pwm, aln1, aln2, min.score="80%", windowSize=51L, cutoff=0.7,
                   conservation=NULL){
            #ans = lapply(x, doSiteSearch, subject, min.score=min.score, windowSize=windowSize, cutoff=cutoff, conservation=conservation)
            ans_list = lapply(pwm, searchAln, aln1, aln2, min.score=min.score, windowSize=windowSize, cutoff=cutoff, conservation=conservation)
            ans = SitePairList(ans_list)
            return(ans)
          }
          )
setMethod("searchAln", signature(pwm="PWMatrixList", aln1="character", aln2="missing"),
          function(pwm, aln1, aln2, min.score="80%", windowSize=51L, cutoff=0.7,
                   conservation=NULL){
            #ans = lapply(x, doSiteSearch, subject, min.score=min.score, windowSize=windowSize, cutoff=cutoff, conservation=conservation)
            ans_list = lapply(pwm, searchAln, aln1, min.score=min.score, windowSize=windowSize, cutoff=cutoff, conservation=conservation)
            ans = SitePairList(ans_list)
            return(ans)
          }
          )
setMethod("searchAln", signature(pwm="PWMatrixList", aln1="DNAStringSet", aln2="missing"),
          function(pwm, aln1, aln2, min.score="80%", windowSize=51L, cutoff=0.7,
                   conservation=NULL){
            #ans = lapply(x, doSiteSearch, subject, min.score=min.score, windowSize=windowSize, cutoff=cutoff, conservation=conservation)
            ans_list = lapply(pwm, searchAln, aln1, min.score=min.score, windowSize=windowSize, cutoff=cutoff, conservation=conservation)
            ans = SitePairList(ans_list)
            return(ans)
          }
          )
setMethod("searchAln", signature(pwm="PWMatrixList", aln1="DNAString", aln2="DNAString"),
          function(pwm, aln1, aln2, min.score="80%", windowSize=51L, cutoff=0.7,
                   conservation=NULL){
            #ans = lapply(x, doSiteSearch, subject, min.score=min.score, windowSize=windowSize, cutoff=cutoff, conservation=conservation)
            ans_list = lapply(pwm, searchAln, aln1, aln2, min.score=min.score, windowSize=windowSize, cutoff=cutoff, conservation=conservation)
            ans = SitePairList(ans_list)
            return(ans)
          }
          )
setMethod("searchAln", signature(pwm="PWMatrixList", aln1="PairwiseAlignmentTFBS", aln2="missing"),
          function(pwm, aln1, aln2, min.score="80%", windowSize=51L, cutoff=0.7,
                   conservation=NULL){
            #ans = lapply(x, doSiteSearch, subject, min.score=min.score, windowSize=windowSize, cutoff=cutoff, conservation=conservation)
            ans_list = lapply(pwm, searchAln, aln1, min.score=min.score, windowSize=windowSize, cutoff=cutoff, conservation=conservation)
            ans = SitePairList(ans_list)
            return(ans)
          }
          )

setMethod("searchAln", signature(pwm="PWMatrix", aln1="character", aln2="character"),
          function(pwm, aln1, aln2, min.score="80%", windowSize=51L, cutoff=0.7,
                   conservation=NULL){
            do_sitesearch(pwm, aln1, aln2, min.score=min.score,
                          windowSize=windowSize, cutoff=cutoff,
                          conservation=conservation)
          }
          )
setMethod("searchAln", signature(pwm="PWMatrix", aln1="character", aln2="missing"),
          function(pwm, aln1, aln2, min.score="80%", windowSize=51L, cutoff=0.7,
                   conservation=NULL){
            if(length(aln1) != 2)
              stop("'aln1' must be of length 2 when 'aln2' is missing")
            do_sitesearch(pwm, aln1[1], aln1[2], min.score=min.score,
                          windowSize=windowSize, cutoff=cutoff,
                          conservation=conservation)
          }
          )
setMethod("searchAln", signature(pwm="PWMatrix", aln1="DNAStringSet", aln2="missing"),
          function(pwm, aln1, aln2, min.score="80%", windowSize=51L, cutoff=0.7,
                   conservation=NULL){
            if(length(aln1) != 2)
              stop("'aln1' must be of length 2 when 'aln2' is missing")
            do_sitesearch(pwm, as.character(aln1[1]), as.character(aln1[2]),
                          min.score=min.score, windowSize=windowSize,
                          cutoff=cutoff, conservation=conservation)
          }
          )
setMethod("searchAln", signature(pwm="PWMatrix", aln1="DNAString", aln2="DNAString"),
          function(pwm, aln1, aln2, min.score="80%", windowSize=51L, cutoff=0.7,
                   conservation=NULL){
            do_sitesearch(pwm, as.character(aln1), as.character(aln2),
                             min.score=min.score, windowSize=windowSize,
                             cutoff=cutoff, conservation=conservation)
          }
          )
setMethod("searchAln", signature(pwm="PWMatrix", aln1="PairwiseAlignmentTFBS", aln2="missing"),
          function(pwm, aln1, aln2, min.score="80%", windowSize=51L, cutoff=0.7,
                   conservation=NULL){
            do_sitesearch(pwm, as.character(pattern(alignments(aln1))),
                          as.character(subject(alignments(aln1))),
                          min.score=min.score, windowSize=windowSize(aln1),
                          cutoff=cutoff, conservation=conservation1(aln1))
          }
          )

