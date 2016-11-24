test_that("parseMEMEOutput", {
  memeOutput <- file.path(system.file("extdata", package="TFBSTools"),
                          "meme.output")
  ans <- parseMEMEOutput(memeOutput)
  expect_equal(length(ans$motifList), 5L)
  expect_equivalent(lengths(ans$motifList), c(17L, 5L, 2L, 5L, 2L))
  expect_equal(ans$motifEvalues, c(4.1e-09, 2.8e+02, 1.2e+03, 1.4e+03, 1.8e+03))
})
