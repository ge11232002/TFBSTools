## Our own normargPfm. The only difference from Biostrings version is we do not require the column sums are identical.
normargPfm = function(x){
    if (!is.matrix(x) || !is.integer(x))
        stop("invalid PFM 'x': not an integer matrix")
    if (is.null(rownames(x)))
        stop("invalid PFM 'x': no row names")
    if (!all(rownames(x) %in% DNA_ALPHABET))
        stop("invalid PFM 'x': row names must be in 'DNA_ALPHABET'")
    if (!all(DNA_BASES %in% rownames(x)))
        stop("invalid PFM 'x': row names must contain A, C, G and T")
    if (any(duplicated(rownames(x))))
        stop("invalid PFM 'x': duplicated row names")
    if (ncol(x) == 0L)
        stop("invalid PFM 'x': no columns")
    if (any(is.na(x)) || any(x < 0L))
        stop("invalid PFM 'x': values cannot be NA or negative")
    if (any(x[!(rownames(x) %in% DNA_BASES), ] != 0L))
        stop("invalid PFM 'x': IUPAC ambiguity letters are represented")
    x <- x[DNA_BASES, , drop = FALSE]  
    x
}

### ------------------------------------------------------------------------
### The "PWM" generic and methods. This is a bit different from the implementation of Biostrings.
###
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
            x = normargPfm(x)
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
            strand = match.arg(strand, c("+", "-", "*"))
            ans_ranges = IRanges()
            ans_score = c()
            ans_strand = c()
            ans_viewsPos = NULL
            ans_viewsNeg = NULL
            if(strand(x)=="+"){
              xPos = x
              xNeg = reverseComplement(x)
            }else{
              xNeg = x
              xPos = reverseComplement(x)
            }
            if(strand %in% c("+", "*")){
              ans_viewsPos = matchPWM(unitScale(Matrix(xPos)), subject, min.score=min.score)
              scorePos = PWMscoreStartingAt(unitScale(Matrix(xPos)), subject(ans_viewsPos),
                                            start(ans_viewsPos))
              # The score here from PWMscoreStartingAt is the unitscaled score. Let's make it into original one, synced with TFBS module. This is validated!
              scorePos = scorePos * (maxScore(Matrix(xPos)) - minScore(Matrix(xPos))) + minScore(Matrix(xPos))
              ans_ranges = c(ans_ranges, ranges(ans_viewsPos))
              ans_score = c(ans_score, scorePos)
              ans_strand = c(ans_strand, rep("+", length(ans_viewsPos)))
            }
            if(strand %in% c("-", "*")){
              ans_viewsNeg = matchPWM(unitScale(Matrix(xNeg)), subject, min.score=min.score)
              scoreNeg = PWMscoreStartingAt(unitScale(Matrix(xNeg)), subject(ans_viewsNeg),
                                            start(ans_viewsNeg))
              scoreNeg = scoreNeg * (maxScore(Matrix(xNeg)) - minScore(Matrix(xNeg))) + minScore(Matrix(xNeg))
              ans_ranges = c(ans_ranges, ranges(ans_viewsNeg))
              ans_score = c(ans_score, scoreNeg)
              ans_strand = c(ans_strand, rep("-", length(ans_viewsNeg)))
            }
            if(!is.null(ans_viewsPos)){
              ans_views = Views(subject=subject(ans_viewsPos), 
                                start=start(ans_ranges),
                                end=end(ans_ranges)
                                )
            }else{
              ans_views = Views(subject=subject(ans_viewsNeg),
                                start=start(ans_ranges),
                                end=end(ans_ranges)
                                )
            }
            stopifnot(isConstant(c(length(ans_strand), length(ans_score),
                                 length(ans_views))))
            ans_site = Site(views=ans_views, seqname=seqname,
                            score=ans_score, strand=ans_strand, 
                            sitesource="TFBS", primary="TF binding site",
                            pattern=xPos
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
                   strand="*", type="any", conservation=NULL){
            #ans = lapply(x, doSiteSearch, subject, min.score=min.score, windowSize=windowSize, cutoff=cutoff, conservation=conservation)
            ans_list = lapply(pwm, searchAln, aln1, aln2, min.score=min.score, 
                              windowSize=windowSize, cutoff=cutoff, 
                              strand=strand, type=type,
                              conservation=conservation)
            ans = SitePairList(ans_list)
            return(ans)
          }
          )
setMethod("searchAln", signature(pwm="PWMatrixList", aln1="character", aln2="missing"),
          function(pwm, aln1, aln2, min.score="80%", windowSize=51L, cutoff=0.7,
                   strand="*", type="any", conservation=NULL){
            #ans = lapply(x, doSiteSearch, subject, min.score=min.score, windowSize=windowSize, cutoff=cutoff, conservation=conservation)
            ans_list = lapply(pwm, searchAln, aln1, min.score=min.score, 
                              windowSize=windowSize, cutoff=cutoff, 
                              strand=strand, type=type, 
                              conservation=conservation)
            ans = SitePairList(ans_list)
            return(ans)
          }
          )
setMethod("searchAln", signature(pwm="PWMatrixList", aln1="DNAStringSet", aln2="missing"),
          function(pwm, aln1, aln2, min.score="80%", windowSize=51L, cutoff=0.7,
                   strand="*", type="any", conservation=NULL){
            #ans = lapply(x, doSiteSearch, subject, min.score=min.score, windowSize=windowSize, cutoff=cutoff, conservation=conservation)
            ans_list = lapply(pwm, searchAln, aln1, min.score=min.score, 
                              windowSize=windowSize, cutoff=cutoff, 
                              strand=strand, type=type,
                              conservation=conservation)
            ans = SitePairList(ans_list)
            return(ans)
          }
          )
setMethod("searchAln", signature(pwm="PWMatrixList", aln1="DNAString", aln2="DNAString"),
          function(pwm, aln1, aln2, min.score="80%", windowSize=51L, cutoff=0.7,
                   strand="*", type="any", conservation=NULL){
            #ans = lapply(x, doSiteSearch, subject, min.score=min.score, windowSize=windowSize, cutoff=cutoff, conservation=conservation)
            ans_list = lapply(pwm, searchAln, aln1, aln2, min.score=min.score, 
                              windowSize=windowSize, cutoff=cutoff, 
                              strand=strand, type=type, 
                              conservation=conservation)
            ans = SitePairList(ans_list)
            return(ans)
          }
          )
#setMethod("searchAln", signature(pwm="PWMatrixList", aln1="PairwiseAlignmentTFBS", aln2="missing"),
#          function(pwm, aln1, aln2, min.score="80%", windowSize=51L, cutoff=0.7,
#                   strand="*", type="any", conservation=NULL){
#            #ans = lapply(x, doSiteSearch, subject, min.score=min.score, windowSize=windowSize, cutoff=cutoff, conservation=conservation)
#            ans_list = lapply(pwm, searchAln, aln1, min.score=min.score, 
#                              windowSize=windowSize, cutoff=cutoff, 
#                              strand=strand, type=type,
#                              conservation=conservation)
#            ans = SitePairList(ans_list)
#            return(ans)
#          }
#          )

setMethod("searchAln", signature(pwm="PWMatrix", aln1="character", aln2="character"),
          function(pwm, aln1, aln2, min.score="80%", windowSize=51L, cutoff=0.7,
                   strand="*", type="any", conservation=NULL){
            do_sitesearch(pwm, aln1, aln2, min.score=min.score,
                          windowSize=windowSize, cutoff=cutoff,
                          strand=strand, type=type,
                          conservation=conservation)
          }
          )
setMethod("searchAln", signature(pwm="PWMatrix", aln1="character", aln2="missing"),
          function(pwm, aln1, aln2, min.score="80%", windowSize=51L, cutoff=0.7,
                   strand="*", type="any", conservation=NULL){
            if(length(aln1) != 2)
              stop("'aln1' must be of length 2 when 'aln2' is missing")
            do_sitesearch(pwm, aln1[1], aln1[2], min.score=min.score,
                          windowSize=windowSize, cutoff=cutoff,
                          strand=strand, type=type,
                          conservation=conservation)
          }
          )
setMethod("searchAln", signature(pwm="PWMatrix", aln1="DNAStringSet", aln2="missing"),
          function(pwm, aln1, aln2, min.score="80%", windowSize=51L, cutoff=0.7,
                   strand="*", type="any", conservation=NULL){
            if(length(aln1) != 2)
              stop("'aln1' must be of length 2 when 'aln2' is missing")
            do_sitesearch(pwm, as.character(aln1[1]), as.character(aln1[2]),
                          min.score=min.score, windowSize=windowSize,
                          cutoff=cutoff, strand=strand,
                          type=type, conservation=conservation)
          }
          )
setMethod("searchAln", signature(pwm="PWMatrix", aln1="DNAString", aln2="DNAString"),
          function(pwm, aln1, aln2, min.score="80%", windowSize=51L, cutoff=0.7,
                   strand="*", type="any", conservation=NULL){
            do_sitesearch(pwm, as.character(aln1), as.character(aln2),
                             min.score=min.score, windowSize=windowSize,
                             cutoff=cutoff, strand=strand, 
                             type=type, conservation=conservation)
          }
          )
#setMethod("searchAln", signature(pwm="PWMatrix", aln1="PairwiseAlignmentTFBS", aln2="missing"),
#          function(pwm, aln1, aln2, min.score="80%", windowSize=51L, cutoff=0.7,
#                   strand="*", type="any", conservation=NULL){
#            do_sitesearch(pwm, as.character(pattern(alignments(aln1))),
#                          as.character(subject(alignments(aln1))),
#                          min.score=min.score, windowSize=windowSize(aln1),
#                          cutoff=cutoff, strand=strand, 
#                          type=type, conservation=conservation1(aln1))
#          }
#          )

