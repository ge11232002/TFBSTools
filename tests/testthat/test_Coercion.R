test_that("SiteSetAsDataFrame", {
  data(MA0003.2)
  pwm1 <- toPWM(MA0003.2)
  seq1 <- "GAATTCTCTCTTGTTGTAGCATTGCCTCAGGGCACACGTGCAAAATG"
  siteset <- searchSeq(pwm1, seq1, seqname="seq1", strand="+", min.score="80%")
  expect_equal(ncol(as(siteset, "DataFrame")), 12L)
})
