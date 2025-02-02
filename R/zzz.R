.onLoad <- function(libname, pkgname) {
  # Prevent reticulate from prompting for Miniconda installation
  Sys.setenv(RETICULATE_MINICONDA_ENABLED = "FALSE")
  Sys.setenv(RETICULATE_AUTOCONFIGURE = "0")  # Disable auto-configuration prompt

  # Define the Python environment
  python_env <- "~/.virtualenvs/my_env/bin/python"

  if (file.exists(python_env)) {
    reticulate::use_python(python_env, required = TRUE)
  } else {
    reticulate::use_python("/usr/bin/python", required = FALSE)
  }

  # Ensure reticulate initializes immediately
  reticulate::py_config()
}

