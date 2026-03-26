#!/usr/bin/env Rscript

# Load the package
library(hpolcc)

# Test 1: Create iid model
cat("Test 1: Creating iid model...\n")
iid_term <- iid("x")
cat("iid class:", class(iid_term), "\n")
cat("iid slots:", slotNames(iid_term), "\n")

# Test 2: Create fpoly model
cat("\nTest 2: Creating fpoly model...\n")
fpoly_term <- fpoly("x", p = 2)
cat("fpoly class:", class(fpoly_term), "\n")
cat("fpoly slots:", slotNames(fpoly_term), "\n")

# Test 3: Create rpoly model
cat("\nTest 3: Creating rpoly model...\n")
rpoly_term <- rpoly("x", p = 2, ref_value = 0)
cat("rpoly class:", class(rpoly_term), "\n")
cat("rpoly slots:", slotNames(rpoly_term), "\n")

# Test 4: Create iwp model
cat("\nTest 4: Creating iwp model...\n")
iwp_term <- iwp("x", p = 2, ref_value = 0, knots = c(0, 1, 2))
cat("iwp class:", class(iwp_term$iwp), "\n")
cat("iwp slots:", slotNames(iwp_term$iwp), "\n")

# Test 5: Create hiwp model
cat("\nTest 5: Creating hiwp model...\n")
hiwp_term <- hiwp("x", p = 2, ref_value = 0, knots = c(0, 1, 2), group_var = "group")
cat("hiwp class:", class(hiwp_term$hiwp), "\n")
cat("hiwp slots:", slotNames(hiwp_term$hiwp), "\n")

# Test 6: Test coercion from hiwp to iwp
cat("\nTest 6: Testing hiwp to iwp coercion...\n")
hiwp_obj <- hiwp_term$hiwp
iwp_from_hiwp <- as(hiwp_obj, "iwp")
cat("Original hiwp class:", class(hiwp_obj), "\n")
cat("Coerced to iwp class:", class(iwp_from_hiwp), "\n")

# Test 7: Test design matrix creation
cat("\nTest 7: Testing design matrix creation...\n")
test_data <- data.frame(x = 1:10)
tryCatch({
  iwp_design <- design(iwp_term$iwp, test_data)
  cat("iwp design matrix created successfully, dimensions:", dim(iwp_design), "\n")
}, error = function(e) {
  cat("Error creating iwp design matrix:", e$message, "\n")
})

cat("\nAll tests completed!\n")