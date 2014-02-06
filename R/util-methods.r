
### -----------------------------------------------------------------
### a bettern system call
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

