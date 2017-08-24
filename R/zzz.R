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
  
  op.milkweed <- list(
    milkweed.dirs = app_dir("Milkweed", "LaMar")
  )
  toset <- !(names(op.milkweed) %in% names(op))
  if(any(toset)) options(op.milkweed[toset])
  
  invisible()
}