.onLoad <- function(libname, pkgname) {
  # Disable Miniconda to prevent prompts
  Sys.setenv(RETICULATE_MINICONDA_ENABLED = "FALSE")
  Sys.setenv(RETICULATE_AUTOCONFIGURE = "0")

  # Define Python environment path (Modify based on your system)
  python_env <- "~/.virtualenvs/my_env/bin/python"

  # Check if the specified Python environment exists
  if (file.exists(python_env)) {
    Sys.setenv(RETICULATE_PYTHON = normalizePath(python_env))
  } else {
    system_python <- Sys.which("python3")  # Find system Python

    if (system_python != "") {
      Sys.setenv(RETICULATE_PYTHON = system_python)
    } else {
      stop("No valid Python environment found. Please install Python and set RETICULATE_PYTHON manually.")
    }
  }

  # Load reticulate AFTER setting the environment
  library(reticulate, quietly = TRUE, warn.conflicts = FALSE)

  # Ensure Python is set up correctly
  reticulate::py_config()
}
