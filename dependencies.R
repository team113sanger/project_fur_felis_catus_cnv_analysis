#!/usr/bin/env Rscript

# Simple R package dependencies installer
# Uses a fixed repository snapshot for compatibility with R 4.2.1

#########################################################################
# CONFIGURATION SECTION - Modify these values as needed                 #
#########################################################################

# Options for pak and repository configuration
OPTIONS <- list(
  # Repository settings
  repos_snapshot_date = "2022-10-28",
  distribution = "jammy",

  # Use Bioconductor
  use_bioconductor = TRUE,
  bioc_version = "3.16",

  # Pak configuration
  pkg.r_versions = as.character(getRversion()),
  pkg.sysreqs = TRUE,
  pkg.sysreqs_platform = paste0("ubuntu-", "22.04"),
  pkg.sysreqs_verbose = TRUE
)

# Packages to install
PACKAGES <- c(
    # Bioconductor packages
    "bioc::DNAcopy@1.72",
    "bioc::ComplexHeatmap@2.14",
    # Target packages
    "ggplot2",
    "optparse",
    "dplyr",
    "tidyr",
    "circlize",
    "randomcoloR"
)

#########################################################################
# SCRIPT FUNCTIONS                                                      #
#########################################################################

# Function to initialize repositories
setup_repositories <- function(options) {
  # Construct repository URL with specific snapshot date
  base_url <- sprintf("https://p3m.dev/cran/__linux__/%s/%s",
                     options$distribution,
                     options$repos_snapshot_date)

  message(sprintf("Using CRAN snapshot from %s for %s",
                 options$repos_snapshot_date,
                 options$distribution))

  # Set repositories for both pak and install.packages
  repos <- c(CRAN = base_url)

  # Configure R options for repositories
  options(repos = repos)
  options(pkg.repos = repos)

  # Set pak-specific options
  for (option_name in names(options)) {
    if (startsWith(option_name, "pkg.")) {
      options(setNames(list(options[[option_name]]), option_name))
    }
  }

  # Configure BioConductor if needed
  if (options$use_bioconductor) {
    message("Setting up Bioconductor repositories...")

    # Install BiocManager if needed
    if (!requireNamespace("BiocManager", quietly = TRUE)) {
      install.packages("BiocManager", repos = repos)
    }

    # Initialize Bioconductor repositories with specific version
    BiocManager::repositories(version = options$bioc_version)
  }

  return(repos)
}

# Function to parse package identifier
parse_package_name <- function(pkg_identifier) {
  # Check if it's a Bioconductor package
  is_bioc <- grepl("^bioc::", pkg_identifier)

  # Remove the "bioc::" prefix if present
  pkg_clean <- ifelse(is_bioc, sub("^bioc::", "", pkg_identifier), pkg_identifier)

  # Extract package name and version if specified (e.g., "package@version")
  pkg_parts <- strsplit(pkg_clean, "@")[[1]]
  pkg_name <- pkg_parts[1]
  pkg_version <- if (length(pkg_parts) > 1) pkg_parts[2] else NULL

  return(list(
    name = pkg_name,
    version = pkg_version,
    is_bioc = is_bioc,
    original = pkg_identifier
  ))
}

# Global variable to track critical warnings
critical_warnings <- character(0)

# Function to install packages using pak
install_packages <- function(packages, repos) {
  # Setup warning handler to capture non-zero exit status warnings
  withCallingHandlers({
    # Install pak if not available
    if (!requireNamespace("pak", quietly = TRUE)) {
      message("Installing pak...")
      install.packages("pak", repos = repos)
    }

    message("\n=== Installing packages with pak ===")
    message("Packages: ", paste(packages, collapse = ", "))

    # Install all packages using pak
    pak_success <- tryCatch({
      pak::pkg_install(packages, ask = FALSE)
      message("\n✓ All packages installed successfully with pak")
      TRUE
    }, error = function(e) {
      message("\n✗ Error installing packages with pak: ", e$message)
      message("\nFalling back to direct installation methods")
      FALSE
    })

    if (!pak_success) {
      # Try installing packages one by one with appropriate method
      results <- sapply(packages, function(pkg_identifier) {
        # Parse the package identifier
        pkg_info <- parse_package_name(pkg_identifier)

        if (pkg_info$is_bioc) {
          # Bioconductor package installation
          tryCatch({
            message("Installing Bioconductor package: ", pkg_info$name)
            if (!requireNamespace("BiocManager", quietly = TRUE)) {
              install.packages("BiocManager", repos = repos)
            }

            # Determine install parameters
            install_args <- list(pkgs = pkg_info$name, update = FALSE, ask = FALSE)
            if (!is.null(pkg_info$version)) {
              message("  (requesting version ", pkg_info$version, ")")
              # Note: BiocManager doesn't directly support version in the install call
              # Version constraints would need to be handled differently if needed
            }

            do.call(BiocManager::install, install_args)
            message("✓ Successfully installed ", pkg_info$name)
            TRUE
          }, error = function(e) {
            message("✗ Failed to install ", pkg_info$name, ": ", e$message)
            FALSE
          })
        } else {
          # Regular CRAN package installation
          tryCatch({
            message("Installing CRAN package: ", pkg_identifier)
            install.packages(pkg_identifier, repos = repos)
            message("✓ Successfully installed ", pkg_identifier)
            TRUE
          }, error = function(e) {
            message("✗ Failed to install ", pkg_identifier, ": ", e$message)
            FALSE
          })
        }
      })

      # Return overall success status
      all_succeeded <- all(results)
      if (all_succeeded) {
        message("\n✓ All fallback installations completed successfully")
      } else {
        message("\n✗ Some fallback installations failed")
      }
      return(all_succeeded)
    }
    return(TRUE)
  }, warning = function(w) {
    # Capture critical warnings
    if (grepl("non-zero exit status", w$message)) {
      pkg_name <- sub(".*package '([^']+)'.*", "\\1", w$message)
      critical_warnings <<- c(critical_warnings, sprintf("Package '%s': %s", pkg_name, w$message))
    }
    # Continue execution
    invokeRestart("muffleWarning")
  })
}

#########################################################################
# MAIN EXECUTION                                                        #
#########################################################################

main <- function() {
  # Set up repositories and options
  repos <- setup_repositories(OPTIONS)

  # Install packages
  success <- install_packages(PACKAGES, repos)

  # Report result
  if (success) {
    message("\n=== Installation completed successfully ===")
  } else {
    message("\n=== Installation completed with errors ===")
  }

  # Check for critical warnings
  if (length(critical_warnings) > 0) {
    message("\n=== CRITICAL WARNINGS DETECTED ===")
    for (warning in critical_warnings) {
      message("- ", warning)
    }
    message("\nExiting with error status due to critical warnings.")
    # Exit with an error code
    quit(status = 1)
  }
}

# Run the script
main()
