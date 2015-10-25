plot.rucl <- function(x, ...) {
  
  barplot(x, ...)
  readline("Press <Enter> for the next plot")
  qqnorm(x, ...)
  readline("Press <Enter> for the next plot")
  hist(x, ...)
  readline("Press <Enter> for the next plot")
  boxplot(x, ...)
  readline("Press <Enter> for the next plot")
  timeplot(x, ...)
  
}

timeplot <- function(x, ...) {
  
  if(x$type == "censored") {
    
    col = rep("black", x$summary.statistics$n)
    col[!x$detects] = "red"
     
  } else {
    
    col = "black"    
    
  }
  
  plot(x$data, pch=19, col=col, main = paste0("Time Series of ", x$data.name))
  lines(x$data, lty=3)
  tl <- lm(x$data~seq(x$summary.statistics$n))
  lines(tl$fitted.values, lwd=2)
  
}


boxplot.rucl <- function(x, ...) {
  
  if(x$type == "uncensored") {
    boxplot(x$data, main = paste0("Boxplot of ", x$data.name), ylab = "Concentration", ...)
  } else {
    
    col <- rep("wheat", 5)
    d <- x$distribution.tests
    if(d$Normal == TRUE) col[3] = "lightgreen"
    if(d$Lognormal == TRUE) col[4] = "lightgreen"
    if(d$Gamma == TRUE) col[5] = "lightgreen"
    
    boxplot(list(x$data, x$data[x$detects], x$ros$norm, x$ros$ln, x$ros$gamma), 
            names = c("With Non-Detects", "Detects Only", "Normal ROS", "Lognormal ROS", "Gamma ROS"),
            main = paste0("Boxplot of ", x$data.name), ylab = "Concentration", col=col, ...)
    
  }
  
}

hist.rucl <- function(x, ...) {
  
  if(x$type == "censored") {
  
    t <- hist(x$data, plot = FALSE, ...)
    nd <- t
    for(i in seq(length(nd$breaks)-1)) {
      nd$counts[i] = sum(!x$detects[x$data > nd$breaks[i] & x$data <= nd$breaks[i+1]])
    }
    plot(t, col = "wheat", main = paste0("Histogram of ", x$data.name))
    plot(nd, col = "pink", add=TRUE)
    
  } else {
    
    hist(x$data, main = paste0("Histogram of ", x$data.name), ...)
    
  }
  
}

barplot.rucl <- function(height, ...) {
  
  x <- height
  
  par(mfrow=c(1,1))
  if(x$type == "censored") {
    m <- x$summary.statistics$km.mean
    mx <- x$summary.statistics$detect.max
  } else {
    m <- x$summary.statistics$mean
    mx <- x$summary.statistics$max
  }
  
  xmax = length(x$ucl.desc)+3
  col = rep("wheat", length(x$ucls))
  col[names(x$ucls) %in% names(x$rec)] = "lightgreen"
  
  barplot(unlist(x$ucls), 
          main = "UCL Comparison", space=0, col = col, border="grey50",
          ylim = c(0,mx*1.1), ylab = "Concentration", 
          xlim = c(0.5, xmax), xlab = "UCLs", xaxt = "n", ...)
  text(seq(length(x$ucl.desc))-0.75, mx/50, x$ucl.desc, srt=90, pos = 4, cex=0.7)
  abline(m, 0, lty=2, lwd = 2, col="grey50")
  text(xmax-1, m, "Mean", pos=3, col="grey30")
  abline(mx, 0, lty=1, lwd = 2, col="grey50")
  text(xmax-1, mx, "Max Value", pos=3, col="grey30")
  
}

qqnorm.rucl <- function(y, ...) {
  
  if(y$type == "censored") {
    z <- sort(y$data[y$detects])
    n <- length(z)
    j <- sort(ros.pp(y$data, y$detects)[y$detects])
    m <- y$summary.statistics$km.mean
    s <- y$summary.statistics$km.sd
    lm <- log(m)
    ls <- s
    k <- y$summary.statistics$g.ros.shape
    r <- y$summary.statistics$g.ros.rate
  } else {
    z <- sort(y$data)
    n <- length(z)
    j <- (seq(n) - 3/8)/(n+1/4)
    m <- y$summary.statistics$mean
    s <- y$summary.statistics$sd
    lm <- y$summary.statistics$ln.mean
    ls <- y$summary.statistics$ln.sd
    k <- y$summary.statistics$g.shape
    r <- y$summary.statistics$g.rate
  }
  
  d <- y$distribution.tests
  par(mfrow=c(3,1))
    
  nx <- pnorm(z, m, s)
  ny <- j
  par(mar=c(1,5,4,2))
  plot(ny,nx, main = "Distribution Testing",
       xlab = "", 
       ylab = "Emperical Percentiles",
       xlim = c(0, 1),
       ylim = c(0, 1), ...)
  if(d$Normal) lcol = "lightgreen" else lcol = "pink"
  abline(0, 1, lwd = 2, col = lcol)
  text(-0.02, 0.9, "Normal Distribution", pos = 4, cex = 1.35)
  
  lx <- plnorm(z, lm, ls)
  ly <- j
  par(mar=c(2.5,5,2.5,2))
  plot(ly, lx,
       xlab = "", 
       ylab = "Emperical Percentiles",
       xlim = c(0, 1),
       ylim = c(0, 1), ...)
  if(d$Lognormal) lcol = "lightgreen" else lcol = "pink"
  abline(0, 1, lwd = 2, col = lcol)
  text(-0.02, 0.9, "Lognormal Distribution", pos = 4, cex = 1.35)
  
  gx <- pgamma(z, k, r)
  gy <- (seq(n)-1/2)/n
  par(mar=c(4,5,1,2))
  plot(gy, gx, 
       xlab = "Theoretical Percentiles", 
       ylab = "Emperical Percentiles",
       xlim = c(0, 1),
       ylim = c(0, 1), ...)
  if(d$Gamma) lcol = "lightgreen" else lcol = "pink"
  abline(0, 1, lwd = 2, col = lcol)
  text(-0.02, 0.9, "Gamma Distribution", pos = 4, cex = 1.35)
  
  par(mfrow=c(1,1))
  
}