#!/usr/bin/env Rscript

# Test minimal case
cat("Testing minimal case...\n")

# Source model_classes first
source("hpolcc/R/model_classes.R")

# Test if classes are defined
cat("iid class defined:", exists("iid"), "\n")
cat("iwp class defined:", exists("iwp"), "\n")

# Try to create a simple setMethod call
tryCatch({
  setGeneric("design", function(term, data) standardGeneric("design"))
  setMethod("design", "iid", function(term, data) {
    cat("iid method called\n")
    NULL
  })
  cat("Simple setMethod works\n")
}, error = function(e) {
  cat("Error with simple setMethod:", e$message, "\n")
})

# Now try to source the problematic part
cat("\nTrying to source effects.R...\n")
tryCatch({
  source("hpolcc/R/effects.R")
  cat("Success!\n")
}, error = function(e) {
  cat("Error:", e$message, "\n")
  # Print the exact line that's causing the issue
  lines <- readLines("hpolcc/R/effects.R")
  cat("Line 332:", lines[332], "\n")
  cat("Line 333:", lines[333], "\n")
})
