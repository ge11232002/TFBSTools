### -----------------------------------------------------------------
### make FlatFileDir from JASPAR database: including *.pfm and matrix_list.txt
### The matrix_list.txt has columns: ID, InformationContent, Name, Class, Tags
### Exported!
makeFlatFileDir <- function(JASPAR){
  outputDir <- "FlatFileDir"
  dir.create(outputDir)
  
  # Make the matrix_list.txt
  pfms <- getMatrixSet(JASPAR, opts=list(all=TRUE))
  pfmIDs <- ID(pfms)
  pfmICs <- sapply(pfms, function(x){sum(totalIC(toICM(x)))})
  pfmNames <- name(pfms)
  ## There are cases that have more than two classes
  pfmClasses <- sapply(matrixClass(pfms), paste, collapse="/")
  pfmTags <- lapply(pfms, function(x){sapply(x@tags, paste, 
                                              collapse="/")})
  pfmTags <- sapply(pfmTags, function(x){x[order(names(x))]})
  pfmTags <- sapply(pfmTags, function(x){paste(";", paste0(names(x), " \"", x, 
                                                "\" ", collapse="; "))})
  ans <- cbind(data.frame(IDs=pfmIDs,
                          ICs=pfmICs,
                          Names=pfmNames,
                          Classes=pfmClasses,
                          Tags=pfmTags
                          ))
  write.table(ans, file=file.path(outputDir, "matrix_list.txt"),
              quote=FALSE, sep="\t",
              col.names = FALSE, row.names = FALSE)
  
  # Make the *.pfm
  lapply(pfms, function(x){write.table(x@profileMatrix,
                                       file=file.path(outputDir,
                                                      paste0(x@ID, ".pfm")),
                                       quote=FALSE, sep="\t",
                                       col.names = FALSE, row.names = FALSE)})
  return("success")
}

### -----------------------------------------------------------------
### readJASPARMatrix: read the jaspar format PFM in txt file
### "individual" format:
### >MA0001.1 AGL3
### A  [ 0  3 79 40 66 48 65 11 65  0 ]
### C  [94 75  4  3  1  2  5  2  3  3 ]
### G  [ 1  0  3  4  1  0  5  3 28 88 ]
### T  [ 2 19 11 50 29 47 22 81  1  6 ]
### "all" format: multiple "individual" matrices and seperated with a blank line
### Exported
.processJASPARText <- function(text, matrixClass){
  ID <- sub("^>", "", strsplit(text[1], "\t")[[1]][1])
  name <- strsplit(text[1], "\t")[[1]][2]
  if(is.na(name))
    name <- ID
  if(!identical(substr(text[2:5], 1, 1), DNA_BASES)){
    stop("The second to fifth lines of the file must start with",
         "`A`, `C`, `G`, `T`, respectively.")
  }
  profileMatrix <- do.call(rbind, strsplit(sub(" *]$", "", 
                                               sub("^(A|C|G|T)( |\t)+\\[ *", "",
                                                   text[2:5])), "( |\t)+"))
  rownames(profileMatrix) <- DNA_BASES
  if(matrixClass == "PFM"){
    mode(profileMatrix) <- "integer"
    ans <- PFMatrix(ID=ID, name=name, profileMatrix=profileMatrix)
  }else if(matrixClass == "PWM"){
    mode(profileMatrix) <- "numeric"
    ans <- PWMatrix(ID=ID, name=name, profileMatrix=profileMatrix)
  }else if(matrixClass == "PWMProb"){
    mode(profileMatrix) <- "numeric"
    bg <- c(A=0.25, C=0.25, G=0.25, T=0.25)
    if(any(profileMatrix==0)){
      profileMatrix <- profileMatrix + 1e-4
    }
    profileMatrix <- log2(sweep(profileMatrix, MARGIN=1, bg, "/"))
    ans <- PWMatrix(ID=ID, name=name, profileMatrix=profileMatrix)
  }
  return(ans)
}

readJASPARMatrix <- function(fn, matrixClass=c("PFM", "PWM", "PWMProb")){
  matrixClass <- match.arg(matrixClass)
  
  text <- readLines(fn)
  text <- text[text != ""]
  
  
  if(length(text) >= 5L){
    if(length(text) %% 5 != 0L && length(text) %% 5 != 0L){
      stop("The `all` format is supposed to have a number of lines",
           "mutipled by 5!")
    }
    text2 <- split(text, rep(1:(length(text)/5), rep(5, length(text)/5)))
    ans <- lapply(text2, .processJASPARText, matrixClass)
    names(ans) <- lapply(ans, name)
    ans <- switch(matrixClass,
                  "PFM"=do.call(PFMatrixList, ans),
                  "PWM"=do.call(PWMatrixList, ans),
                  "PWMProb"=do.call(PWMatrixList, ans))
  }else{
    stop("Fewer than 5 lines in input file!")
  }
  return(ans)
}
