summary.rucl <- function(object, ...) {
  

}

print.rucl.summary <- function(x, ...) {
  
  cat("Dataset:", x$data.name, "\n")
  cat("Dataset Type:", x$type, "\n")
  cat("Number of Samples:", x$summary.statistics$n, "\n")
  cat("Recommended UCL(s)\n")
  cat("------------------\n")
  l <- max(nchar(x$rec.desc))
  desc <- sapply(x$rec.desc, 
                 function(v) paste0(paste0(rep(" ", length.out = l - nchar(v)), collapse = ""), v, collapse = ""))
  values <- sapply(x$rec, function(v) {
    t <- as.character(signif(v, 3))
    if(nchar(t) > 6) t <- substr(t, 1, 6)
    paste0(t, paste0(rep(" ", length.out =  6 - nchar(t)), collapse = ""), collapse = "")
  })
  for(i in seq(length(desc))) cat(paste0(desc[i], ":"), values[i], "\n")  
    
}