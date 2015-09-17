
test_readXMLTFFM <- function(){
  # TFFMFirst
  xmlFirst <- file.path(system.file("extdata", package="TFBSTools"),
                        "tffm_first_order.xml")
  tffmFirst <- readXMLTFFM(xmlFirst, type="First")
  ## check the length of emission distribution parameters
  checkIdentical(length(tffmFirst@emission), 17L)
  ## check the dimensions of transition probabilities.
  checkIdentical(dim(tffmFirst@transition), c(17L, 16L))
  ## check the information content of TFFMFirst
  checkEquals(sum(totalIC((tffmFirst))), 8.744361, tolerance=1e-5)
  
  # TFFMDetail
  xmlDetail <- file.path(system.file("extdata", package="TFBSTools"),
                         "tffm_detailed.xml")
  tffmDetail <- readXMLTFFM(xmlDetail, type="Detail")
  ## check the length of emission distribution parameters
  checkIdentical(length(tffmDetail@emission), 64L)
  ## check the dimensions of transition probabilities.
  checkIdentical(dim(tffmDetail@transition), c(64L, 64L))
  ## check the information content of TFFMDetail
  checkEquals(sum(totalIC((tffmDetail))), 7.745002, tolerance=1e-5)
}

