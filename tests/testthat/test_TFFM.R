test_that("test_readXMLTFFM", {
    xmlFirst <- file.path(system.file("extdata", package = "TFBSTools"), 
        "tffm_first_order.xml")
    tffmFirst <- readXMLTFFM(xmlFirst, type = "First")
    expect_identical(17L, length(tffmFirst@emission))
    expect_identical(c(17L, 16L), dim(tffmFirst@transition))
    expect_equal(8.744361, sum(totalIC((tffmFirst))), tolerance=1e-5)
    xmlDetail <- file.path(system.file("extdata", package = "TFBSTools"), 
        "tffm_detailed.xml")
    tffmDetail <- readXMLTFFM(xmlDetail, type = "Detail")
    expect_identical(64L, length(tffmDetail@emission))
    expect_identical(c(64L, 64L), dim(tffmDetail@transition))
    expect_equal(7.745002, sum(totalIC((tffmDetail))), tolerance=1e-5)
})

