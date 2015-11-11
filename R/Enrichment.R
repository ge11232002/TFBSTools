
tfbsEnrichment <- function(pwm, query, background, cutoff=0.8){
  queryHits <- searchSeq(pwm, query, min.score=cutoff)
  foo1 <- as(queryHits, "GRanges")
  backgroundHits <- searchSeq(pwm, background, min.score=cutoff)
  foo2 <- as(backgroundHits, "GRanges")
}
