# Contributing to Seurat CLI

Thank you for your interest in contributing to this project! This document provides guidelines for contributing.

## How to Contribute

### Reporting Bugs

If you find a bug, please open an issue with:

1. A clear, descriptive title
2. Steps to reproduce the problem
3. Expected behavior vs. actual behavior
4. Your environment details:
   - R version (`R --version`)
   - Seurat version (`packageVersion("Seurat")`)
   - Operating system
5. Any relevant error messages or logs

### Suggesting Features

Feature requests are welcome! Please open an issue with:

1. A clear description of the feature
2. The use case / motivation
3. Example of how it might work

### Submitting Changes

1. **Fork the repository** and create your branch from `main`
2. **Make your changes** following the code style guidelines below
3. **Test your changes** with demo mode and your own data if possible
4. **Update documentation** if you're adding or changing functionality
5. **Submit a pull request** with a clear description of changes

## Code Style Guidelines

### R Code Style

We follow the [tidyverse style guide](https://style.tidyverse.org/) with these specifics:

```r
# Use snake_case for variables and functions
variable_name <- "value"
my_function <- function(arg1, arg2) { }

# Use 2-space indentation
if (condition) {
  do_something()
}

# Comment your code meaningfully
# Calculate QC metrics for filtering
qc_metrics <- calculate_qc_metrics(seurat_obj)

# Use explicit returns at end of functions
my_function <- function(x) {
  result <- x * 2
  return(result)
}
```

### Script Structure

When adding or modifying scripts, follow this structure:

```r
#!/usr/bin/env Rscript

# =============================================================================
# Script Name: XX_script_name.R
# Description: Brief description of what the script does
# =============================================================================

# Load common utilities
source(file.path(dirname(sys.frame(1)$ofile), "../utils/common.R"))

# --- Argument Parsing ---
# Define and parse arguments

# --- Main Function ---
main <- function(args) {
  # Implementation
}

# --- Execute ---
if (!interactive()) {
  args <- parse_args(parser)
  main(args)
}
```

### Documentation

- Add help text for all new command-line arguments
- Update README.md if adding new features or parameters
- Include examples in help text where appropriate

## Testing

Before submitting:

1. **Test with demo mode:**
   ```bash
   Rscript scripts/your_script.R --demo --output ./test_output
   ```

2. **Test help output:**
   ```bash
   Rscript scripts/your_script.R --help
   ```

3. **Test with different input formats** if applicable

4. **Verify output files** are created correctly

## Pull Request Process

1. Ensure your code follows the style guidelines
2. Update the README.md with details of changes if applicable
3. Add an entry to CHANGELOG.md under "Unreleased"
4. The PR will be reviewed and merged once approved

## Questions?

Feel free to open an issue for any questions about contributing!
