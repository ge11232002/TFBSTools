test_that("test_searchSeq", {
    library(IRanges)
    library(Biostrings)
    data(MA0003.2)
    data(MA0004.1)
    pwm1 <- toPWM(MA0003.2)
    pwm2 <- toPWM(MA0004.1)
    pwmList <- PWMatrixList(pwm1 = pwm1, pwm2 = pwm2)
    seq1 <- "GAATTCTCTCTTGTTGTAGCATTGCCTCAGGGCACACGTGCAAAATG"
    seq2 <- "GTTTCACCATTGCCTCAGGGCATAAATATATAAAAAAATATAATTTTCATC"
    siteset <- searchSeq(pwm1, seq1, seqname = "seq1", strand = "*", 
        min.score = "80%")
    expect_identical(3L, length(siteset))
    siteset <- searchSeq(pwm1, seq1, seqname = "seq1", strand = "+", 
        min.score = "80%")
    expect_identical(1L, length(siteset))
    siteset <- searchSeq(pwm1, seq1, seqname = "seq1", strand = "-", 
        min.score = "80%")
    expect_identical(2L, length(siteset))
    seqs <- DNAStringSet(c(seq1 = seq1, seq2 = seq2))
    sitesetList <- searchSeq(pwmList, seqs, min.score = "80%")
    expect_identical(IRanges(start = c(20L, 22L, 23L, 8L, 10L, 
        11L, 35L, 35L), end = c(34L, 36L, 37L, 22L, 24L, 25L, 
        40L, 40L)), ranges(unname(as(sitesetList, "GRanges"))))
})

