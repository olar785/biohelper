.onLoad <- function(libname, pkgname) {
  Sys.setenv(RETICULATE_MINICONDA_ENABLED = "FALSE")

  # Set Python environment
  python_env <- "~/.virtualenvs/my_env/bin/python"

  if (file.exists(python_env)) {
    reticulate::use_python(python_env, required = TRUE)
  } else {
    reticulate::use_python("/usr/bin/python", required = FALSE)
  }
}
