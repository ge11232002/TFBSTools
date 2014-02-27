
test_XMatrixConstructor <- function(){
  # just check this does throw an error
  checkException(PFMatrix(pseudocounts=0.8), silent=TRUE)
  # just check this does throw an error
  checkException(PWMatrix(schneider=TRUE), silent=TRUE)

  # check the wrong bg
  ## wrong base
  checkException(PFMatrix(bg=c(A=0.25, C=0.25, G=0.25, U=0.25)), silent=TRUE)
  ## negative frequency
  checkException(PFMatrix(bg=c(A=-0.5, C=0.25, G=0.25, T=0.25)), silent=TRUE)
  ## wrong oder
  checkException(PFMatrix(bg=c(A=-0.5, C=0.25, T=0.25, G=0.25)), silent=TRUE)

  # check the length of XMatrix's arguments
  checkException(PFMatrix(ID=c("12", "34")), silent=TRUE)
  checkException(PFMatrix(name=c("12", "34")), silent=TRUE)
  checkException(PFMatrix(matrixClass=c("12", "34")), silent=TRUE)
  checkException(PFMatrix(strand=c("+", "-")), silent=TRUE)

  # check the strand's value
  checkException(PFMatrix(strand=c("a")), silent=TRUE)

  # check the profileMatrix
  checkException(PFMatrix(profileMatrix=matrix(1, ncol=2, nrow=2)), silent=TRUE)
}

test_XMatrixListConstructor <- function(){
  ## PFMatrixList constructor
  checkException(PFMatrixList(PWMatrix(), PWMatrix()), silent=TRUE)
  checkException(PWMatrixList(PFMatrix(), PFMatrix()), silent=TRUE)
  checkException(ICMatrixList(PWMatrix(), PWMatrix()), silent=TRUE)
}

test_Accessor <- function(){
  pfm <- PFMatrix(ID="MA0004.1", name="Arnt", matrixClass="Zipper-Type", 
                  strand="+",
                  bg=c(A=0.25, C=0.25, G=0.25, T=0.25),
                  tags=list(family="Helix-Loop-Helix", 
                            species="10090", tax_group="vertebrates",
                            medline="7592839", type="SELEX", ACC="P53762", 
                            pazar_tf_id="TF0000003",
                            TFBSshape_ID="11", TFencyclopedia_ID="580"),
                  profileMatrix=matrix(c(4L,  19L, 0L,  0L,  0L,  0L,
                                         16L, 0L,  20L, 0L,  0L,  0L,
                                         0L,  1L,  0L,  20L, 0L,  20L,
                                         0L,  0L,  0L,  0L,  20L, 0L),
                                       byrow=TRUE, nrow=4, 
                                       dimnames=list(c("A", "C", "G", "T")))
                  )
  checkIdentical(ID(pfm), "MA0004.1")
  checkIdentical(name(pfm), "Arnt")
  checkIdentical(matrixClass(pfm), "Zipper-Type")

}

