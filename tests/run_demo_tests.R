#!/usr/bin/env Rscript

# =============================================================================
# Test Runner: Runs all scripts in demo mode to verify functionality
# Usage: Rscript tests/run_demo_tests.R [--quick] [--script SCRIPT_NUM]
# =============================================================================

suppressPackageStartupMessages({
  library(optparse)
})

# Parse arguments
option_list <- list(
  make_option("--quick", action = "store_true", default = FALSE,
              help = "Quick test (scripts 1-2 only)"),
  make_option("--script", type = "integer", default = NULL,
              help = "Run specific script number (1-6)"),
  make_option("--output_base", type = "character", default = "test_output",
              help = "Base output directory [default: %default]"),
  make_option("--cleanup", action = "store_true", default = FALSE,
              help = "Remove output after successful tests")
)

parser <- OptionParser(
  usage = "Usage: %prog [options]",
  option_list = option_list,
  description = "Test runner for Seurat vignette scripts"
)

args <- parse_args(parser)

# Define test cases
tests <- list(
  list(
    name = "01_basic_analysis",
    script = "scripts/01_basic_analysis.R",
    args = c("--demo"),
    description = "Standard workflow with PBMC 3K data"
  ),
  list(
    name = "02_sctransform", 
    script = "scripts/02_sctransform.R",
    args = c("--demo"),
    description = "SCTransform normalization"
  ),
  list(
    name = "03_integration",
    script = "scripts/03_integration.R", 
    args = c("--demo"),
    description = "Data integration with IFNB data"
  ),
  list(
    name = "04_differential_expression",
    script = "scripts/04_differential_expression.R",
    args = c("--demo"),
    description = "Differential expression analysis"
  ),
  list(
    name = "05_cell_cycle",
    script = "scripts/05_cell_cycle.R",
    args = c("--demo"),
    description = "Cell cycle scoring"
  ),
  list(
    name = "06_visualization",
    script = "scripts/06_visualization.R",
    args = c("--demo"),
    description = "Visualization suite"
  )
)

# Filter tests
if (!is.null(args$script)) {
  if (args$script < 1 || args$script > 6) {
    stop("Script number must be 1-6")
  }
  tests <- tests[args$script]
} else if (args$quick) {
  tests <- tests[1:2]
}

# Get script directory
get_script_dir <- function() {
  # Method 1: Rscript --file argument
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg) > 0) {
    return(dirname(normalizePath(sub("^--file=", "", file_arg))))
  }
  # Method 2: source() - search through frames for ofile
  for (i in seq_len(sys.nframe())) {
    ofile <- sys.frame(i)$ofile
    if (!is.null(ofile)) {
      return(dirname(normalizePath(ofile)))
    }
  }
  # Fallback to current directory
  return(".")
}

script_dir <- get_script_dir()
project_dir <- normalizePath(file.path(script_dir, ".."))

# Run tests
cat("\n")
cat("=" , rep("=", 60), "\n", sep = "")
cat("  SEURAT CLI - DEMO TEST RUNNER\n")
cat("=", rep("=", 60), "\n", sep = "")
cat("\nProject directory:", project_dir, "\n")
cat("Running", length(tests), "test(s)\n\n")

results <- list()
start_time <- Sys.time()

for (i in seq_along(tests)) {
  test <- tests[[i]]
  output_dir <- file.path(args$output_base, test$name)
  
  cat("-", rep("-", 60), "\n", sep = "")
  cat("Test", i, "/", length(tests), ":", test$name, "\n")
  cat("Description:", test$description, "\n")
  cat("-", rep("-", 60), "\n", sep = "")
  
  # Build command
  script_path <- file.path(project_dir, test$script)
  cmd_args <- c(test$args, "--output", output_dir, "--quiet")
  cmd <- paste("Rscript", shQuote(script_path), paste(cmd_args, collapse = " "))
  
  cat("Command:", cmd, "\n\n")
  
  test_start <- Sys.time()
  
  # Run the test
  exit_code <- system(cmd)
  
  test_end <- Sys.time()
  test_duration <- round(difftime(test_end, test_start, units = "secs"), 1)
  
  # Check results
  if (exit_code == 0) {
    # Check for expected output files
    expected_files <- c("parameters.csv", "session_info.txt")
    found_files <- file.exists(file.path(output_dir, expected_files))
    
    if (all(found_files)) {
      cat("\n✓ PASSED (", test_duration, "s)\n\n", sep = "")
      results[[test$name]] <- list(status = "PASSED", duration = test_duration)
    } else {
      cat("\n✗ FAILED - Missing output files\n\n")
      results[[test$name]] <- list(status = "FAILED", duration = test_duration, 
                                    reason = "Missing output files")
    }
  } else {
    cat("\n✗ FAILED - Exit code:", exit_code, "\n\n")
    results[[test$name]] <- list(status = "FAILED", duration = test_duration,
                                  reason = paste("Exit code:", exit_code))
  }
}

# Summary
end_time <- Sys.time()
total_duration <- round(difftime(end_time, start_time, units = "mins"), 1)

cat("\n")
cat("=", rep("=", 60), "\n", sep = "")
cat("  TEST SUMMARY\n")
cat("=", rep("=", 60), "\n\n", sep = "")

passed <- sum(sapply(results, function(x) x$status == "PASSED"))
failed <- sum(sapply(results, function(x) x$status == "FAILED"))

for (name in names(results)) {
  r <- results[[name]]
  status_icon <- if (r$status == "PASSED") "✓" else "✗"
  cat(sprintf("  %s %-35s %s (%ss)\n", status_icon, name, r$status, r$duration))
}

cat("\n")
cat("Total:", length(results), "| Passed:", passed, "| Failed:", failed, "\n")
cat("Duration:", total_duration, "minutes\n")

# Cleanup
if (args$cleanup && failed == 0) {
  cat("\nCleaning up test output...\n")
  unlink(args$output_base, recursive = TRUE)
  cat("Done.\n")
}

# Exit with appropriate code
if (failed > 0) {
  quit(status = 1)
}

cat("\nAll tests passed! ✓\n\n")
