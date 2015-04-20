## Tue Mar 18 2015
## Original file Copyright © 2015 S. Chen
## This file is part of the R package kdfe.
## Su Chen <su.chen@memphis.edu> 
## Modified by: Danny Arends <Danny.Arends@gmail.com>
## Department of Mathmematical Sciences
## Faculty of Statistics
## University of Memphis
## Memphis, TN
## USA
##############################################################################

rot.kdfe <- function(x, type = "s", mu = mean(x)) {
  if (length(x) < 2L) stop("need at least 2 data points")
  if (is.na(type)) stop("Please specify type: 's' or 'l'")        # Danny: here we need to check if the type is correct, this just checks for NA

  r   <- quantile(x,c(0.25,0.75))
  h   <- (r[2L] - r[1L]) / 1.34
  sig <- min(sd(x),h)
  if (type == "s") {
    return(1.515717 * sig * length(x)^(-2/5))
  }else if (type == "l") {
    return((4*(2+(sig/mu)^2))^{-1/5}*sig*length(x)^(-2/5))
  }else stop("Note that type must be 's' or 'l'")
}

# Preparation for Direct Plug-in 
d2k <- function(u){
  (1/2)*exp(-(1/2)*u^2)*sqrt(2)*(u^2-1)/sqrt(pi)  # second derivative of Guassian N(0,1)
}

#Direct Plug-in bandwidth for Kernel Functional Estimation (KFE): location and scale
# dpi.kdfe(x=data, type=c("s","l"),pilot=c("rot.d","ucv.d","bcv.d","dpi.d","ste.d","rot.kdfe","manual")) 
# By default, pilot is chosen as "rot.d".
dpi.kdfe <- function(x, type = "s", pilot = "rot.d", g = NA){
  if (length(x) < 2L) { stop("need at least 2 data points") }
  if (is.na(type)) { stop("Please specify type: 's' or 'l'") }

  if (pilot == "rot.d"){
    g = bw.nrd(x)
  }else if (pilot == "ucv.d"){ 
    g = bw.ucv(x, nb = 1000, lower = 0.1, upper = 5,tol = 0.01)
  }else if (pilot == "bcv.d"){ 
    g = bw.bcv(x, nb = 1000, lower = 0.1, upper = 5,tol = 0.01)
  }else if (pilot == "dpi.d"){
    g = bw.SJ(x, nb = 1000, lower = 0.1, upper =5, method = "dpi", tol = 0.01)
  }else if (pilot == "ste.d"){
    g = bw.SJ(x, nb = 1000, lower = 0.1, upper =5, method = "ste", tol = 0.01)
  }else if (pilot == "rot.kdfe"){
    g = rot.kdfe(x, type)
  }else if (pilot=="manual"){ 
    g = g 
  }

  if(is.na(g)) { stop("Please specify pilot type or input the bandwidth parameter g") }

  matrix1 <- dnorm(outer(x,x,"-")/g)/g
  matrix2 <- d2k(outer(x,x,"-")/g)/(g^3)

  if (type=="s"){  
    Rf      <- mean(matrix1[upper.tri(matrix1, diag = FALSE)])
    phi     <- mean(matrix2[row(matrix2)!=col(matrix2)]) 
  } else if (type=="l") {
    matrix3 <- outer(x^2,x^2,"+")*matrix1/2
    matrix4 <- diag(x) %*% matrix2
    Rf      <- mean(matrix3[upper.tri(matrix3, diag = FALSE)])
    phi     <- mean(matrix4[row(matrix4)!=col(matrix4)])
  } else { stop("Note that type must be 's' or 'l' and pilot must be 'rot.d','ucv.d','bcv.d','dpi.d','ste.d' or 'rot.kdfe'") }
  return((Rf/(sqrt(pi)*(phi^2)))^(1/5)*(length(x)^(-2/5)))
}

