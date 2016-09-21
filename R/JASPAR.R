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
.processJASPARText <- function(text){
  ID <- sub("^>", "", strsplit(text[1], "\t")[[1]][1])
  name <- strsplit(text[1], "\t")[[1]][2]
  if(!identical(substr(text[2:5], 1, 1), DNA_BASES)){
    stop("The second to fifth lines of the file must start with",
         "`A`, `C`, `G`, `T`, respectively.")
  }
  profileMatrix <- do.call(rbind, strsplit(sub(" *]$", "", 
                                               sub("^(A|C|G|T)  \\[ *", "",
                                                   text[2:5])), " +"))
  mode(profileMatrix) <- "integer"
  rownames(profileMatrix) <- DNA_BASES
  ans <- PFMatrix(ID=ID, name=name, profileMatrix=profileMatrix)
}

readJASPARMatrix <- function(fn, type=c("individual", "all")){
  type <- match.arg(type)
  text <- readLines(fn)
  if(type == "individual"){
    if(length(text) != 5L){
      stop("The `individual` format is supposed to have 5 lines!")
    }
    ans <- .processJASPARText(text)
  }else{
    if(length(text) %% 6 != 0L){
      stop("The `all` format is supposed to have a number of lines",
           "mutipled by 6!")
    }
    text2 <- split(text, rep(1:(length(text)/6), rep(6, length(text)/6)))
    ans <- lapply(text2, .processJASPARText)
    ans <- do.call(PFMatrixList, ans)
  }
  return(ans)
}