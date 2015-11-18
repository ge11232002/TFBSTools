test_that("test_Accessor", {
    pfm <- PFMatrix(ID = "MA0004.1", name = "Arnt", matrixClass = "Zipper-Type", 
        strand = "+", bg = c(A = 0.25, C = 0.25, G = 0.25, T = 0.25), 
        tags = list(family = "Helix-Loop-Helix", species = "10090", 
            tax_group = "vertebrates", medline = "7592839", type = "SELEX", 
            ACC = "P53762", pazar_tf_id = "TF0000003", TFBSshape_ID = "11", 
            TFencyclopedia_ID = "580"), profileMatrix = matrix(c(4L, 
            19L, 0L, 0L, 0L, 0L, 16L, 0L, 20L, 0L, 0L, 0L, 0L, 
            1L, 0L, 20L, 0L, 20L, 0L, 0L, 0L, 0L, 20L, 0L), byrow = TRUE, 
            nrow = 4, dimnames = list(c("A", "C", "G", "T"))))
    expect_identical("MA0004.1", ID(pfm))
    expect_identical("Arnt", name(pfm))
    expect_identical("Zipper-Type", matrixClass(pfm))
})

test_that("test_XMatrixConstructor", {
    expect_error(PFMatrix(pseudocounts = 0.8))
    expect_error(PWMatrix(schneider = TRUE))
    expect_error(PFMatrix(bg = c(A = 0.25, C = 0.25, G = 0.25, 
        U = 0.25)))
    expect_error(PFMatrix(bg = c(A = -0.5, C = 0.25, G = 0.25, 
        T = 0.25)))
    expect_error(PFMatrix(bg = c(A = -0.5, C = 0.25, T = 0.25, 
        G = 0.25)))
    expect_error(PFMatrix(ID = c("12", "34")))
    expect_error(PFMatrix(name = c("12", "34")))
    expect_error(PFMatrix(strand = c("+", "-")))
    expect_error(PFMatrix(strand = c("a")))
    expect_error(PFMatrix(profileMatrix = matrix(1, ncol = 2, 
        nrow = 2)))
})

test_that("test_XMatrixListConstructor", {
    expect_error(PFMatrixList(PWMatrix(), PWMatrix()))
    expect_error(PWMatrixList(PFMatrix(), PFMatrix()))
    expect_error(ICMatrixList(PWMatrix(), PWMatrix()))
})

