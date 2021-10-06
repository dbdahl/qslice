#' @docType package
#' @usage NULL
#' @useDynLib cucumber, .registration = TRUE
NULL

.Kall <- function(...) {
  x <- .Call(...)
  if ( inherits(x,"error") ) stop(x) else x
}
