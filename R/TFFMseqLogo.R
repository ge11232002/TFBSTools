######################
##
## plot sequence logo
## Based on the code from seqLogo package. Add the transparency parameter: alpha
######################

letterA <- function(x.pos,y.pos,ht,wt,alpha=1,id=NULL){

  x <- c(0,4,6,10,8,6.8,3.2,2,0,3.6,5,6.4,3.6)
  y <- c(0,10,10,0,0,3,3,0,0,4,7.5,4,4)
  x <- 0.1*x
  y <- 0.1*y

  x <- x.pos + wt*x
  y <- y.pos + ht*y

  if (is.null(id)){
    id <- c(rep(1,9),rep(2,4))
  }else{
    id <- c(rep(id,9),rep(id+1,4))
  }

  fill <- c("green","white")
  alpha <- c(alpha, 1)
  list(x=x,y=y,id=id,fill=fill,alpha=alpha)
}

## T
letterT <- function(x.pos,y.pos,ht,wt,alpha=1,id=NULL){

  x <- c(0,10,10,6,6,4,4,0)
  y <- c(10,10,9,9,0,0,9,9)
  x <- 0.1*x
  y <- 0.1*y

  x <- x.pos + wt*x
  y <- y.pos + ht*y

  if (is.null(id)){
    id <- rep(1,8)
  }else{
    id <- rep(id,8)
  }

  fill <- "red"

  list(x=x,y=y,id=id,fill=fill,alpha=alpha)
}

## C
letterC <- function(x.pos,y.pos,ht,wt,alpha=1,id=NULL){
  angle1 <- seq(0.3+pi/2,pi,length=100)
  angle2 <- seq(pi,1.5*pi,length=100)
  x.l1 <- 0.5 + 0.5*sin(angle1)
  y.l1 <- 0.5 + 0.5*cos(angle1)
  x.l2 <- 0.5 + 0.5*sin(angle2)
  y.l2 <- 0.5 + 0.5*cos(angle2)

  x.l <- c(x.l1,x.l2)
  y.l <- c(y.l1,y.l2)

  x <- c(x.l,rev(x.l))
  y <- c(y.l,1-rev(y.l))

  x.i1 <- 0.5 +0.35*sin(angle1)
  y.i1 <- 0.5 +0.35*cos(angle1)
  x.i1 <- x.i1[y.i1<=max(y.l1)]
  y.i1 <- y.i1[y.i1<=max(y.l1)]
  y.i1[1] <- max(y.l1)

  x.i2 <- 0.5 +0.35*sin(angle2)
  y.i2 <- 0.5 +0.35*cos(angle2)

  x.i <- c(x.i1,x.i2)
  y.i <- c(y.i1,y.i2)

  x1 <- c(x.i,rev(x.i))
  y1 <- c(y.i,1-rev(y.i))

  x <- c(x,rev(x1))
  y <- c(y,rev(y1))

  x <- x.pos + wt*x
  y <- y.pos + ht*y

  if (is.null(id)){
    id <- rep(1,length(x))
  }else{
    id <- rep(id,length(x))
  }

  fill <- "blue"

  list(x=x,y=y,id=id,fill=fill,alpha=alpha)
}

## G
letterG <- function(x.pos,y.pos,ht,wt,alpha=1,id=NULL){
  angle1 <- seq(0.3+pi/2,pi,length=100)
  angle2 <- seq(pi,1.5*pi,length=100)
  x.l1 <- 0.5 + 0.5*sin(angle1)
  y.l1 <- 0.5 + 0.5*cos(angle1)
  x.l2 <- 0.5 + 0.5*sin(angle2)
  y.l2 <- 0.5 + 0.5*cos(angle2)

  x.l <- c(x.l1,x.l2)
  y.l <- c(y.l1,y.l2)

  x <- c(x.l,rev(x.l))
  y <- c(y.l,1-rev(y.l))

  x.i1 <- 0.5 +0.35*sin(angle1)
  y.i1 <- 0.5 +0.35*cos(angle1)
  x.i1 <- x.i1[y.i1<=max(y.l1)]
  y.i1 <- y.i1[y.i1<=max(y.l1)]
  y.i1[1] <- max(y.l1)

  x.i2 <- 0.5 +0.35*sin(angle2)
  y.i2 <- 0.5 +0.35*cos(angle2)

  x.i <- c(x.i1,x.i2)
  y.i <- c(y.i1,y.i2)

  x1 <- c(x.i,rev(x.i))
  y1 <- c(y.i,1-rev(y.i))

  x <- c(x,rev(x1))
  y <- c(y,rev(y1))

  h1 <- max(y.l1)
  r1 <- max(x.l1)

  h1 <- 0.4
  x.add <- c(r1,0.5,0.5,r1-0.2,r1-0.2,r1,r1)
  y.add <- c(h1,h1,h1-0.1,h1-0.1,0,0,h1)

  if (is.null(id)){
    id <- c(rep(1,length(x)),rep(2,length(x.add)))
  }else{
    id <- c(rep(id,length(x)),rep(id+1,length(x.add)))
  }

  x <- c(rev(x),x.add)
  y <- c(rev(y),y.add)

  x <- x.pos + wt*x
  y <- y.pos + ht*y


  fill <- c("orange","orange")
  alpha <- c(alpha, alpha)
  list(x=x,y=y,id=id,fill=fill,alpha=alpha)

}

addLetter <- function(letters,which,x.pos,y.pos,ht,wt,alpha=1){

  if (which == "A"){
    letter <- letterA(x.pos,y.pos,ht,wt,alpha)
  }else if (which == "C"){
    letter <- letterC(x.pos,y.pos,ht,wt,alpha)
  }else if (which == "G"){
    letter <- letterG(x.pos,y.pos,ht,wt,alpha)
  }else if (which == "T"){
    letter <- letterT(x.pos,y.pos,ht,wt,alpha)
  }else{
    stop("which must be one of A,C,G,T")
  }

  letters$x <- c(letters$x,letter$x)
  letters$y <- c(letters$y,letter$y)

  lastID <- ifelse(is.null(letters$id),0,max(letters$id))
  letters$id <- c(letters$id,lastID+letter$id)
  letters$fill <- c(letters$fill,letter$fill)
  letters$alpha <- c(letters$alpha, letter$alpha)
  letters
}

seqLogoTFFM <- function(x, xfontsize=10, yfontsize=10, xaxis=TRUE, yaxis=TRUE){
  emissions <- getEmissionProb(x)
  posProb <- getPosProb(x)
  bgProb <- bgEmissionProb(x)

  posProb <- cbind(bgProb, posProb[ ,-ncol(posProb)])
  letters <- list(x=NULL,y=NULL,id=NULL,fill=NULL,alpha=NULL)

  npos <- ncol(posProb)
  wt <- 1
  facs <- 2
  x.pos <- 2
  chars <- DNA_BASES 
  for(j in 1:npos){
    y.pos <- 8
    for(k in 1:4){
      column <- emissions[((k-1)*4+1):((k-1)*4+4),j]
      hts <- 0.95*column*facs
      letterOrder <- order(hts)
      y.posNew <- y.pos - 2

      for (i in 1:4){
        letter <- chars[letterOrder[i]]
        ht <- hts[letterOrder[i]]
        if (ht>0) letters <- addLetter(letters,letter,x.pos,y.pos,ht,wt,
                                       posProb[k,j])
        y.pos <- y.pos + ht + 0.01
      }
      y.pos <- y.posNew
    }
    x.pos <- x.pos + wt
  }

  ## Add the left letters
  letters <- addLetter(letters, "A", 0, 8, ht=1.9, wt=1.7)
  letters <- addLetter(letters, "C", 0, 6, ht=1.9, wt=1.7)
  letters <- addLetter(letters, "G", 0, 4, ht=1.9, wt=1.7)
  letters <- addLetter(letters, "T", 0, 2, ht=1.9, wt=1.7)

  ## Add the information content logos
  x.pos <- 2
  ylim <- 2
  pwm <- getPosProb(x)
  ylab <- "Information content"
  facs <- seqLogo:::pwm2ic(pwm)
  for(j in 1:npos){
    column <- pwm[,j]
    hts <- 0.95*column*facs[j]
    letterOrder <- order(hts)
    y.pos <- 0 
    for(i in 1:4){
      letter <- chars[letterOrder[i]]
      ht <- hts[letterOrder[i]]
      if (ht>0) letters <- addLetter(letters,letter,x.pos,y.pos,ht,wt)
      y.pos <- y.pos + ht + 0.01
    }
    x.pos <- x.pos + wt
  }

  grid.newpage()
  bottomMargin = ifelse(xaxis, 2 + xfontsize/3.5, 2)
  leftMargin = ifelse(yaxis, 2 + yfontsize/3.5, 2)
  pushViewport(plotViewport(c(bottomMargin,2,5,3)))
  pushViewport(dataViewport(0:ncol(pwm),0:ylim*4,name="vp1"))
  grid.polygon(x=unit(letters$x,"native"), y=unit(letters$y,"native"),
               id=letters$id,
               gp=gpar(fill=letters$fill,col="transparent", 
                       alpha=letters$alpha))
  grid.xaxis(at=seq(2.5,ncol(pwm)+2-0.5),label=1:ncol(pwm), 
             gp=gpar(fontsize=xfontsize))
  grid.text("Position",y=unit(-3,"lines"), gp=gpar(fontsize=xfontsize))
  #grid.yaxis(gp=gpar(fontsize=yfontsize))
  #grid.text(ylab,x=unit(-3,"lines"),rot=90, gp=gpar(fontsize=yfontsize))

  popViewport()
  popViewport()
  par(ask=FALSE)
}

setMethod("seqLogo", "TFFM",
          function(x, ic.scale = TRUE, xaxis = TRUE, yaxis = TRUE,
                   xfontsize = 10, yfontsize = 10){
            seqLogoTFFM(x, xfontsize=xfontsize, yfontsize=yfontsize,
                        xaxis=xaxis, yaxis=yaxis)
          })


