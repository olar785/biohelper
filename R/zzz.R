.onLoad <- function(libname, pkgname) {
  # Disable automatic Miniconda configuration
  Sys.setenv(RETICULATE_MINICONDA_ENABLED = "FALSE")
  Sys.setenv(RETICULATE_AUTOCONFIGURE = "0")

  # Define the Python environment
  python_env <- "~/.virtualenvs/my_env/bin/python"

  if (file.exists(python_env)) {
    Sys.setenv(RETICULATE_PYTHON = python_env)  # Force reticulate to use this Python version
  } else {
    Sys.setenv(RETICULATE_PYTHON = "/usr/bin/python")  # Fallback system Python
  }

  # Load reticulate and apply settings BEFORE it initializes
  library(reticulate, quietly = TRUE, warn.conflicts = FALSE)

  # Ensure reticulate uses the defined Python
  reticulate::py_discover_config()  # Forces reticulate to apply the settings early
}
