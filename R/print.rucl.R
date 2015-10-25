print.rucl <- function(x, ...) {
  
  sf <- 3
  
  if(x$type == "uncensored") {
  
    ss.desc <- c("Sample Size", "Minimum Value", "Maximum Value", "Arithmetic Mean", 
                 "Geometric Mean", "Median", "Standard Deviation", "Standard Error", 
                 "Coefficient of Variance", "Skewness", "Minimum Value of Log Data", 
                 "Maximum Value of Log Data", "Mean of Log Data", 
                 "Standard Deviation of Log Data", "MLE of Gamma Shape",
                 "MLE of Gamma Scale", "MLE of Gamma Rate",
                 "MLE of Gamma Mean", "MLE of Gamma SSD")
    
    ucl.desc <- c("Student's t UCL", "Adjusted Central Limit Theorem UCL", 
                  "Modified Student's t UCL", "Land's H UCL", 
                  "95% Chebyshev UCL (MVUE)", "97.5% Chebyshev UCL (MVUE)", 
                  "99% Chebyshev UCL (MVUE)", "Approximate Gamma UCL",
                  "Adjusted Gamma UCL", "Central Limit Theorem UCL", 
                  "Standard Bootstrap UCL", "Bootstrap t UCL",
                  "Hall's Bootstrap UCL", "Simple Percentile Bootstrap UCL", 
                  "BCA Percentile Bootstrap", "95% Chebyshev UCL (mean, sd)", 
                  "97.5% Chebyshev UCL (mean, sd)", "99% Chebyshev UCL (mean, sd)"
    )
    
  } else {
    
    ss.desc <- list("Sample Size", "# of Distinct Detections", "# of Detections", 
                    "# of Non-Detects", "Detection Rate", "Minimum Detection", 
                    "Maximum Detection", "Arith. Mean of Detections", 
                    "SD of Detections", "Arith. Mean using 1/2 DL", 
                    "SD using 1/2 DL", "SE using 1/2 DL", "Normal ROS Mean",
                    "Normal ROS SD", "Normal ROS SE", "Kaplan-Meier Mean", 
                    "Kaplan-Meier SD", "Kaplan-Meier SE", 
                    "Minimum of Log Detections", "Maximum of Log Detections", 
                    "Mean of Log Detections", "SD of Log Detections", 
                    "Mean of Log Data using 1/2 DL", "SD of Log Data using 1/2 DL",
                    "SE of Log Data using 1/2 DL", "Log ROS-based Mean - Log",
                    "Log ROS-based SD - Log", "Log ROS-based Mean",
                    "Log ROS-based SD", "Log ROS-based SE", "MLE of Gamma Shape - Detects",
                    "MLE of Gamma Scale - Detects", "MLE of Gamma Rate - Detects",
                    "MLE of Gamma Mean - ROS", "MLE of Gamma Shape - ROS",
                    "MLE of Gamma Scale - ROS", "MLE of Gamma Rate - ROS")
    
    
    ucl.desc <- list("Student's t UCL (1/2 DL)", "Student's t UCL (Normal ROS)",
                     "CLT UCL (Normal ROS)", "Land's H UCL (1/2 DL)", 
                     "Student's t UCL (LN ROS)", "Percentile Bootstrap UCL (LN ROS)", 
                     "BCA Bootstrap (LN ROS)", "Land's H UCL (LN ROS)", 
                     "Approximate Gamma UCL (Gamma ROS)", "Adjusted Gamma UCL (Gamma ROS)",
                     "Student's t UCL (KM)", "CLT UCL (KM)", "Bootstrap-t UCL (KM)",
                     "BCA Percentile Bootstrap UCL (KM)", "Percentile Bootstrap UCL (KM)",
                     "95% Chebyshev UCL (KM)", "97.5% Chebyshev UCL (KM)",
                     "99% Chebyshev UCL (KM)")
  
  }
  
  cat("\n")
  cat("--------------------------------------------------------------------------------\n")
  cat("  Summary for", x$type, "dataset", paste0("'", x$data.name, "'") , "\n")
  cat("--------------------------------------------------------------------------------\n")
  cat("    Confidence Coefficient:", paste0(x$confidence*100, "%"), "\t\t\t\t\t\t")
  cat("Number of Bootstrap Operations(N):", x$N, "\n")
  cat("\n")
  
  cat("--------------------------------------------------------------------------------\n")
  cat("  Summary Statistics\n")
  cat("--------------------------------------------------------------------------------\n")
  l <- max(nchar(ss.desc))
  desc <- sapply(ss.desc, 
                 function(v) paste0(paste0(rep(" ", length.out = 30 - nchar(v)), collapse = ""), v, collapse = ""))
  values <- sapply(x$summary.statistics, function(v) {
    t <- as.character(signif(v, sf))
    if(nchar(t) > 6) t <- substr(t, 1, 6)
    paste0(t, paste0(rep(" ", length.out =  6 - nchar(t)), collapse = ""), collapse = "")
  })
  hc <- ceiling(length(desc)/2)
  for(i in seq(hc)) {if(i+hc <= length(desc)) {
    cat(paste0(desc[i], ":"), values[i], "|", paste0(desc[i+hc], ":"), values[i+hc], "\n")
  } else {
    cat(paste0(desc[i], ":"), values[i], "|\n")
  }
  }
  cat("\n")
  
  cat("--------------------------------------------------------------------------------\n")
  cat("  Distribution Testing\n")
  cat("--------------------------------------------------------------------------------\n")
  cat("    Normal:", x$dist.test$Normal, "          ")
  cat("Lognormal:", x$dist.test$Lognormal, "          ")
  cat("Gamma:", x$dist.test$Gamma, "\n")
  cat("\n")
  cat("    Recommendation:", x$dist.test$Recommend, "\n")
  cat("\n")
  
  cat("--------------------------------------------------------------------------------\n")
  cat("  Upper Confidence Limits of the Mean (UCLs)\n")
  cat("--------------------------------------------------------------------------------\n")
  l <- max(nchar(ucl.desc))
  m <- rep("    ", length(x$ucls))
  m[names(x$ucls) %in% names(x$rec)] = "****"
  desc <- sapply(ucl.desc, 
                 function(v) paste0(paste0(rep(" ", length.out = l - nchar(v)), collapse = ""), v, collapse = ""))
  values <- sapply(x$ucls, function(v) {
    t <- as.character(signif(v, sf))
    if(nchar(t) > 6) t <- substr(t, 1, 6)
    paste0(t, paste0(rep(" ", length.out =  6 - nchar(t)), collapse = ""), collapse = "")
  })
  for(i in seq(length(desc))) cat(paste0(desc[i], ":"), values[i], m[i], "\n")  
  cat("\n****: Recommended UCL(s)")
  cat("\n\n")
  
}