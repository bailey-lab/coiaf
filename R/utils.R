# Silence text output
quiet <- function(x) {
  sink(nullfile())
  on.exit(sink())
  invisible(force(x))
}
