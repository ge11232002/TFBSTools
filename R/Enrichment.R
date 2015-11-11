
tfbsEnrichment <- function(pwm, query, background, cutoff=0.8){
  queryHits <- searchSeq(pwm, query, strand="+", min.score=cutoff)
  foo1 <- as(queryHits, "GRanges")
  backgroundHits <- searchSeq(pwm, background, strand="+", min.score=cutoff)
  foo2 <- as(backgroundHits, "GRanges")
}
