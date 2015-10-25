data(uu.picks, envir=environment())
data(cu.picks, envir=environment())
data(adja, envir=environment())
data(h95, envir=environment())

#Unbiased Estimate of the third central moment (Kleijnen, Kloppenburg, and Meeuwsen 1986)
u3 <- function(x) {
  n <- length(x)
  m <- sum(x)/n
  sum((x - m)^3) / ((n - 1) * (n - 2))
}

skew <- function(x) length(x) * u3(x) / sd(x)^3

# Takes an mns object and uses it to interpolate values using natural splines
mns.interpolate <- function(mns, new.x) {
  
  # This ridiculous looking thing takes each new.x value, matches it up with its
  # corresponding ns object in the mns object, performs a prediction, then
  # outputs a vector containing all the values with names appropriate for
  # comparing to the coefficients present in the mns object... or something
  p <- do.call(c, lapply(seq(length(new.x)), 
                         function(i, x, mns) {
                           a <- as.vector(predict(mns[["ns"]][[names(x[i])]], x[i]))
                           names(a) <- paste0(names(x[i]), seq(length(a)))
                           return(a)
                         }, 
                         x = new.x, mns = mns))
  
  #Calculate the new design matrix
  dm <- sapply(mns$dmo, function(x, p) prod(p[x], na.rm = TRUE), p = p)
  dm[1] = 1 # set the intercept to 1
  pred <- as.vector(dm %*% mns$coef) #Matrix multiple the result with the beta coeff
  #   if(pred > a) pred = a
  #   if(pred < 0.00) pred = 0
  return(pred)
  
}
