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

# Function to install packages using pak
install_packages <- function(packages, repos) {
  # Install pak if not available
  if (!requireNamespace("pak", quietly = TRUE)) {
    message("Installing pak...")
    install.packages("pak", repos = repos)
  }

  message("\n=== Installing packages with pak ===")
  message("Packages: ", paste(packages, collapse = ", "))

  # Install all packages using pak
  tryCatch({
    pak::pkg_install(packages, ask = FALSE)
    message("\n✓ All packages installed successfully")
    return(TRUE)
  }, error = function(e) {
    message("\n✗ Error installing packages: ", e$message)
    message("\nFalling back to direct installation using install.packages()")

    # Try installing packages one by one with install.packages()
    results <- sapply(packages, function(pkg) {
      tryCatch({
        install.packages(pkg, repos = repos)
        message("✓ Successfully installed ", pkg)
        return(TRUE)
      }, error = function(e) {
        message("✗ Failed to install ", pkg, ": ", e$message)
        return(FALSE)
      })
    })

    return(all(results))
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
}

# Run the script
main()
