if ( !require(pracma) )
{
  install.packages("pracma")
  library(pracma)
}


if ( !require(funHDDC) )
{
  install.packages("funHDDC")
  library(funHDDC)
}


if ( !require(plotrix) )
{
  install.packages("plotrix")
  library(plotrix)
}


if ( !require(coda) )
{
  install.packages("coda")
  library(coda)
}





majority <- function(v) {
  if(length(as.numeric(names(which(table(v)==max(table(v)))))) > 1) {
    a <- 1
  }
  else if (length(as.numeric(names(which(table(v)==max(table(v)))))) == 1) 
    a <- as.numeric(names(which(table(v)==max(table(v)))))
  return(a)
}


euc.dist <- function(x1, x2) sqrt(sum((x1 - x2) ^ 2))

hilbert_dist <- function(x1,x2,t1,x3,x4,t2) sqrt(trapz(t1,(x1-x2)^2)+trapz(t2,(x3-x4)^2))



diagnostic_performance <- function(y) {
  CC <- (y[1,1] +y[2,2]) / sum(y)
  sen <- y[2,2] / (y[1,2]+y[2,2])
  spe <- y[1,1] / (y[1,1]+y[2,1])
  ppv <- y[2,2] / (y[2,1]+y[2,2])
  npv <- y[1,1] / (y[1,1]+y[1,2])
  
  return(list("CC" = CC, "Sensitivity"= sen, 
              "Specificity" = spe, "PPV" = ppv, "NPV" = npv))
}

diagnostic_performance2 <- function(y) {
  CC <- (y[1,1] +y[3,2]) / sum(y)
  sen <- y[3,2] / (y[1,2]+y[2,2]+y[3,2])
  spe <- y[1,1] / (y[1,1]+y[2,1]+y[3,1])
  ppv <- y[3,2] / (y[3,1]+y[3,2])
  npv <- y[1,1] / (y[1,1]+y[1,2])
  
  return(list("CC" = CC, "Sensitivity"= sen, 
              "Specificity" = spe, "PPV" = ppv, "NPV" = npv))
}










make_basis <- function(p,nu,psi, t1,t2) {
  m1 <- length(t1)
  m2 <- length(t2)
  bs1 <- matrix(rep(0,(p-1)*m1),nrow=m1,ncol=p-1)
  for (i in 1:m1) {
    for (j in 1:(p-1)) {
      bs1[i,j] <- exp(-nu*(t1[i]-j*psi)^2) 
    }
  }
  b1 <- cbind(rep(1,m1),bs1)
  
  
  bs2 <- matrix(rep(0,(p-1)*m2),nrow=m2,ncol=p-1)
  for (i in 1:m2) {
    for (j in 1:(p-1)) {
      bs2[i,j] <- exp(-nu*(t2[i]-j*psi)^2) 
    }
  }
  b2 <- cbind(rep(1,m2),bs2)
  return(list('b1'=b1, 'b2'=b2))
  
}