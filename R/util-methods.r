
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



