# Creates generic function "ucl"
# ucl <- function(x, type, confidence = 0.95, N = 2000, d, ...) UseMethod("ucl")

# This is the default function for the generic function 'ucl'. It requires a 
# numeric vector 'x' and optionally a logical vector 'd' for censored datasets 
# indicating whether the values in 'x' are detections. 'type' allows the user to
# pick a specific type of ucl to be calculated. If missing, it will use the ucl
# picking algorithm to select an appropriate ucl. Alternative confidence and N 
# values can be provided, but if missing will default to 0.95 and 2000 
# respectively. If 'fast' = FALSE then a rucl object is returned, which will 
# allow for additional graphing and summarizing of the data. This can be slow, 
# especially with large datasets. Setting 'fast' = TRUE will only calculate the 
# ucl(s) that the packages deems to be appropriate and return a numeric vector 
# with those values. The names of that vector will contain the ucl type 
# contained in that cell.
ucl <- function(x, type = "fast", confidence = 0.95, N = 2000, d, ...) {

  if(!is.numeric(x)) stop("x must be a numeric vector")
  
  if(length(unique(x)) < 2) {
    return(NA)
  }
  
  i <- !is.na(x)
  
  if(!missing(d)) {
    d <- as.logical(d)
    if(length(x) != length(d)) {
      stop("vectors x and d must be of same length")
    }
    d <- d[i]
  }
  
  x <- x[i]
    
  x[x <= 0] <- 9e-6 # Prevents errors in log and gamma calculations
  
  suppressWarnings(if(type == "detailed") {
    
    if(missing(d)) {
      u <- u.ucl(x, confidence, N, ...)
    } else {
      u <- c.ucl(x, d, confidence, N, ...)
    }
    
  } else if(type == "fast") {
    if(missing(d)) {
      u <- u.ucl.fast(x, confidence, N, ...)
    } else {
      u <- cu.ucl.fast(x, d, confidence, N, ...)
    }
    
  } else {
    
    u <- singleUCL(x, d, type, confidence, N, ...)
    
  })

  u
  
}

#This function is called when the type arguement in UCL indicates that the user 
#wants a specific UCL calculated. It will take the list of values in type and 
#calculate the UCL that corresponds to each one. As long as one value matches a
#UCL it will complete sucessfully. If it cannot find any UCLs to calculate it
#will stop with an error.
singleUCL <- function(x, d, type, confidence, N, ...) {
  
  fm <- function(x, n) sum(x) / n
  fsd <- function(x, n) sqrt(sum((x - fm(x, n)) ^ 2) / (n - 1))
  fse <- function(x, n) (sqrt(sum((x - fm(x, n)) ^ 2) / (n - 1))) / sqrt(n)
  hdl <- function(x, d) {
    x[!d] = x/2
    x
  }
  
  
  if(missing(d)) d <- rep(TRUE, length(x))
  
  n <- length(x)
  y <- log(x)
  
  # Unparsed expressions containing code to calculate parameters necessary to
  # calculate UCLs
  parameters <- list(
    n.tucl = list("m = fm(x, n)", "se = fse(x, n)"),
    n.modz = list("m = fm(x, n)", "se = fse(x, n)", "sk = skew(x)"),
    n.modt = list("m = fm(x, n)", "u3 = u3(x)", "se = fse(x, n)"),
    l.hucl = list("my = fm(y, n)", "sdy = fsd(y, n)"),
    l.mvue95 = list("my = fm(y, n)", "vy = fsd(y, n)^2"),
    l.mvue975 = list("my = fm(y, n)", "vy = fsd(y, n)^2"),
    l.mvue99 = list("my = fm(y, n)", "vy = fsd(y, n)^2"),
    g.appgam = list("m = fm(x, n)", "k = gshape(x)"),
    g.adjgam = list("m = fm(x, n)", "k = gshape(x)"),
    o.zucl = list("m = fm(x, n)", "se = fse(x, n)"),
    o.zboot = list("B = u.bootstrap(x, N)"),
    o.tboot = list("m = fm(x, n)", "se = fse(x, n)"),
    o.hallboot = list("m = fm(x, n)", "sk = skew(x)", "se = fse(x, n)"),
    o.pboot = list("B = u.bootstrap(x, N)"),
    o.bcaboot = list("B = u.bootstrap(x, N)", "m = fm(x, n)"),
    o.cheb95 = list("m = fm(x, n)", "se = fse(x, n)"),
    o.cheb975 = list("m = fm(x, n)", "se = fse(x, n)"),
    o.cheb99 = list("m = fm(x, n)", "se = fse(x, n)"),
    n.halfdl.t = list("dnm = fm(hdl(x, d), n)", "dnse = fse(hdl(x,d), n)"),
    n.ros.t = list("r = ros(x, d)", "nrm = fm(r$norm, n)", "nrs = fse(r$norm, n)"),
    n.ros.z = list("r = ros(x, d)", "nrm = fm(r$norm, n)", "nrs = fse(r$norm, n)"),
    ln.halfdl.h = list("dlm = fm(log(hdl(x, d)), n)", "dlsd = fsd(log(hdl(x,d)), n)"),
    ln.ros.t = list("r = ros(x, d)", "lrm = fm(r$ln, n)", "lrs = fse(r$ln, n)"),
    ln.ros.pboot = list("r = ros(x, d)", "lrB = u.bootstrap(r$ln, N)"),
    ln.ros.bcaboot = list("r = ros(x, d)", "lrB = u.bootstrap(r$ln, N)", "lrm = fm(r$ln, n)"),
    ln.ros.h = list("r = ros(x, d)", "llrm = fm(log(r$ln), n)", "llrs = fse(log(r$ln), n)"),
    g.ros.appgamma = list("r = ros(x, d)", "grm = fm(r$gamma, n)", "rk = gshape(r$gamma)"),
    g.ros.adjgamma = list("r = ros(x, d)", "grm = fm(r$gamma, n)", "rk = gshape(r$gamma)"),
    o.km.t = list("km = ple(x, d)"),
    o.km.z = list("km = ple(x, d)"),
    o.km.tboot = list("km = ple(x, d)"),
    o.km.bcaboot = list("kB = c.bootstrap(x, d, N)", "km = ple(x, d)"),
    o.km.pboot = list("kB = c.bootstrap(x, d, N)"),
    o.km.cheb95 = list("km = ple(x, d)"),
    o.km.cheb975 = list("km = ple(x, d)"),
    o.km.cheb99 = list("km = ple(x, d)")
  )
  
  #Unparsed expressions containing code to calculate the UCLs needed
  ucls <- list(
    n.tucl = "b.tucl(m, se, n, confidence)",
    n.modz = "b.modz(m, se, sk, n, confidence)",
    n.modt = "b.modt(m, u3, se, n, confidence)",
    l.hucl = "b.hucl(my, sdy, n, confidence)",
    l.mvue95 = "b.mvue(my, vy, n, 0.95)",
    l.mvue975 = "b.mvue(my, vy, n, 0.975)",
    l.mvue99 = "b.mvue(my, vy, n, 0.99)",
    g.appgam = "b.appgamma(m, k, n, confidence)",
    g.adjgam = "b.adjgamma(m, k, n, confidence)",
    o.zucl = "b.zucl(m, se, confidence)",
    o.zboot = "b.zboot(B, confidence)",
    o.tboot = "u.tboot(x, m, se, confidence, N)",
    o.hallboot = "b.hallboot(x, m, n, sk, se, confidence, N)",
    o.pboot = "b.pboot(B, confidence)",
    o.bcaboot = "u.bcaboot(B, x, m, n, confidence, N)",
    o.cheb95 = "b.cheb(m, se, 0.95)",
    o.cheb975 = "b.cheb(m, se, 0.975)",
    o.cheb99 = "b.cheb(m, se, 0.99)",
    n.halfdl.t = "b.tucl(dnm, dnse, n, confidence)",
    n.ros.t = "b.tucl(nrm, nrs, n, confidence)",
    n.ros.z = "b.zucl(nrm, nrs, confidence)",
    ln.halfdl.h = "b.hucl(dlm, dlsd, n, confidence)",
    ln.ros.t = "b.tucl(lrm, lrs, n, confidence)",
    ln.ros.pboot = "b.pboot(lrB, confidence)",
    ln.ros.bcaboot = "u.bcaboot(lrB, r$ln, lrm, n, confidence, N)",
    ln.ros.h = "b.hucl(llrm, llrs, n, confidence)",
    g.ros.appgamma = "b.appgamma(grm, rk, n, confidence)",
    g.ros.adjgamma = "b.adjgamma(grm, rk, n, confidence)",
    o.km.t = "b.tucl(km$mean, km$se, n, confidence)",
    o.km.z = "b.zucl(km$mean, km$se, confidence)",
    o.km.tboot = "c.tboot(x, d, km$mean, km$se, confidence, N)",
    o.km.bcaboot = "c.bcaboot(kB, x, d, km$mean, n, confidence, N)",
    o.km.pboot = "b.pboot(kB, confidence)",
    o.km.cheb95 = "b.cheb(km$mean, km$se, 0.95)",
    o.km.cheb975 = "b.cheb(km$mean, km$se, 0.975)",
    o.km.cheb99 = "b.cheb(km$mean, km$se, 0.99)"
  )
  
  p <- unique(unlist(parameters[unlist(type)]))
  if(length(p) == 0) stop("Cannot find any matching UCLs")
  for(i in p) eval(parse(text = i))
  
  v <- sapply(type, function(i) eval(parse(text = ucls[i])))
  names(v) <- type
  
  v
    
}