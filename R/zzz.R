.onLoad <- function(libname, pkgname) {
  op <- options()
  op.devtools <- list(
    devtools.path = "~/dev/R-sandbox",
    devtools.install.args = "",
    devtools.name = "M. Drew LaMar",
    devtools.desc.author = '"M. Drew LaMar <drew.lamar@gmail.com> [aut, cre]"',
    devtools.desc.license = "MIT",
    devtools.desc.suggests = NULL,
    devtools.desc = list()
  )
  toset <- !(names(op.devtools) %in% names(op))
  if(any(toset)) options(op.devtools[toset])

  dirs <- rappdirs::app_dir("Milkweed", "LaMar")
  op.milkweed <- list(
    milkweed.cache = file.path(dirs$data(), "calculated")
  )
  toset <- !(names(op.milkweed) %in% names(op))
  if(any(toset)) options(op.milkweed[toset])

  # Create directory for cached calculations
  dir.create(file.path(dirs$data(), "calculated"),
             recursive = TRUE,
             showWarnings = FALSE)

  invisible()
}
