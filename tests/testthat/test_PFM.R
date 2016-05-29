test_that("test_PFMSimilarity", {
  library(Biostrings)
  library(JASPAR2016)
  profileMatrix <- matrix(as.integer(
    c(13, 13,  3,  1, 54,  1,  1,  1,  0,  3,  2,  5,
      13, 39,  5, 53,  0,  1, 50,  1,  0, 37,  0, 17,
      17,  2, 37,  0,  0, 52,  3,  0, 53,  8, 37, 12,
      11,  0,  9,  0,  0,  0,  0, 52,  1,  6, 15, 20)),
    nrow=4, byrow=TRUE, dimnames=list(DNA_BASES))
  pfmQuery <- PFMatrix(profileMatrix=profileMatrix)
  pfmSubjects <- getMatrixSet(JASPAR2016,
                              opts=list(ID=c("MA0500", "MA0499", "MA0521",
                                             "MA0697")))
  scores <- PFMSimilarity(pfmSubjects, pfmQuery)
  expect_equal(unname(scores[[1]]), c(20.65002, 93.86375), tolerance=0.01)
  expect_equal(unname(scores[[2]]), c(20.52325, 85.51356), tolerance=0.01)
  expect_equal(unname(scores[[3]]), c(20.44110, 92.91409), tolerance=0.01)
  expect_equal(unname(scores[[4]]), c(20.04603, 83.52514), tolerance=0.01)
}
)