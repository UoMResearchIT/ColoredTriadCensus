# Fix uniTri
# This is a temporary solution, better would be to use memoization
# https://github.com/UoMResearchIT/ColoredTriadCensus/issues/4

devtools::load_all()

#' Update cached `uniTri[["n type"]]`
#' for `n` in `ncolors` and `type` in `types`
update_uniTri <- function(types, ncolors = 1:10, dry.run = FALSE) {

  any.changes <- FALSE

  for (type in types) {
    stopifnot(type %in% names(triad.types))

    for (n in ncolors) {

      # from uniTri[[n type]]
      old <- unique.triad.index(type, n, use.cache = TRUE)
      stopifnot(is.null(old) | length(old) == n^3)

      # re-calculataed from triad.types[[type]] triangle
      new <- unique.triad.index(type, n, use.cache = FALSE)
      stopifnot(length(new) == n^3)

      if (!isTRUE(all.equal(old, new))) {
        ndiffs <- ifelse(is.null(old), length(new), sum(old != new))
        cat("Replacing", ndiffs, "counts from", n, type,"\n")
        uniTri[[paste(n, type, sep = " ")]] <- new
        any.changes <- TRUE
      } else {
        cat("Matching counts on", n, type,"\n")
      }
    }
  }

  if (!any.changes | dry.run) return()

  # move R/sysdata.rda to R/sysdata.rda.bak
  if (file.exists("R/sysdata.rda")) {
    file.rename("R/sysdata.rda", "R/sysdata.rda.bak")
  }

  usethis::use_data(uniTri, internal = TRUE, overwrite = TRUE)
}

restore_uniTri <- function() {
  if (file.exists("R/sysdata.rda.bak")) {
    file.rename("R/sysdata.rda.bak", "R/sysdata.rda")
  } else {
    stop("No backup of sysdata.rda found.")
  }
}
