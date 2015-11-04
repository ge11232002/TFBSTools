
### -----------------------------------------------------------------
### a better system call
### Not Exported!
my.system = function(cmd, echo=TRUE, intern=FALSE, ...){
  if (echo){
    message(cmd)
  }
  res = system(cmd, intern=intern, ...)
  if (!intern){
    stopifnot(res == 0)
  }
  return(res)
}

### -----------------------------------------------------------------
### findOverlapsBases
### Not exported!
findLargestOverlaps = function(query, subject){
  hits = findOverlaps(query, subject, select="all")
  if(length(hits) == 0L){
    ## if no overlaps, return hits with length 0.
    return(hits)
  }
  hitsQuery = query[queryHits(hits)]
  hitsSubject = subject[subjectHits(hits)]
  overlapLength = sapply(mapply(intersect, hitsQuery, hitsSubject, 
                                SIMPLIFY=FALSE), length)
  splittedLength = split(overlapLength, queryHits(hits))
  groupMax = sapply(splittedLength, max)
  stopifnot(length(groupMax) == length(splittedLength))
  for(i in 1:length(splittedLength)){
    splittedLength[[i]] = splittedLength[[i]] == groupMax[i]
  }
  return(hits[unname(unlist(splittedLength))])
}

### -----------------------------------------------------------------
### shannon.entropy
### Exported!
shannon.entropy <- function(p)
{
  if (min(p) < 0 || sum(p) <= 0)
    return(NA)
  p.norm <- p[p>0]/sum(p)
  -sum(log2(p.norm)*p.norm)
}

### -----------------------------------------------------------------
### IUPAC_CODE_MAP to the matrix conversion
### Exported!
IUPAC2Matrix <- function(x){
  x <- as.character(x)
  x <- strsplit(x, "")[[1]]
  if(!all(x %in% names(IUPAC_CODE_MAP))){
    stop("All characters must be in IUPAC_CODE_MAP!")
  }
  ans <- matrix(0L, nrow=4, ncol=length(x),
                dimnames=list(c("A", "C", "G", "T")))
  for(i in 1:length(x)){
    dnaCharacters <- strsplit(IUPAC_CODE_MAP[x[i]], "")[[1]]
    ans[dnaCharacters, i] <- 1L
  }
  return(ans)
}


### -----------------------------------------------------------------
### Sampling ranges from areas of subject sequence based on input ranges
### sampleRanges exported!
sampleRangesOneStrand <- function(inputGRanges, subjectGRanges){
  if(length(inputGRanges) == 0L){
    return(GRanges())
  }
  if(length(subjectGRanges) == 0L){
    return(GRanges())
  }
  widthsInput <- width(inputGRanges)
  widthsSubject <- width(subjectGRanges)
  indexAll <- lapply(widthsInput, 
                     function(x, widthsSubject){which(x<=widthsSubject)},
                     widthsSubject)
  indexSampling <- sapply(indexAll, sample, size=1L)
  selectedSubjectGRanges <- subjectGRanges[indexSampling]
  sampledStart <- sapply(end(selectedSubjectGRanges) - 
                         width(selectedSubjectGRanges) + 1L,
                         function(x){sample(1L:x, size=1L)})
  sampledGRanges <- GRanges(seqnames=seqnames(selectedSubjectGRanges),
                      ranges=IRanges(start=start(selectedSubjectGRanges)+
                                     sampledStart-1L,
                                     width=width(inputGRanges)),
                            strand="*",
                            seqinfo=seqinfo(selectedSubjectGRanges))
  stopifnot(length(sampledGRanges) == length(inputGRanges))
  stopifnot(all(width(sampledGRanges) <= width(selectedSubjectGRanges)))
  return(sampledGRanges)
}

sampleRanges <- function(inputGRanges, subjectGRanges, ignore.strand=TRUE){
  if(ignore.strand){
    ans <- sampleRangesOneStrand(inputGRanges, subjectGRanges)
  }else{
    orderPostive <- which(strand(inputGRanges)=="+")
    sampledGRangesPostive <- sampleRangesOneStrand(
                               inputGRanges[orderPostive],
                               subjectGRanges[strand(subjectGRanges)=="+"])
    strand(sampledGRangesPostive) <- "+"
    orderNegative <- which(strand(inputGRanges)=="-")
    sampledGRangesNegative <- sampleRangesOneStrand(
                                inputGRanges[orderNegative],
                                subjectGRanges[strand(subjectGRanges)=="-"])
    strand(sampledGRangesNegative) <- "-"
    orderUnknow <- which(strand(inputGRanges)=="*")
    sampledGRangesUnknown <- sampleRangesOneStrand(
                               inputGRanges[orderUnknow],
                               subjectGRanges)
    ans <- c(sampledGRangesPostive, sampledGRangesNegative, 
             sampledGRangesUnknown)[order(c(orderPostive, orderNegative, orderUnknow))]
  }
  return(ans)
}


