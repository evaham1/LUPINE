
#' Progress bar for long running tasks
#'
#' @param fn function to wrap
#' @param total total number of iterations
#'
#' @return wrapped function
#' @export
#'
completed_counter <- function(fn, total) {
  completed <- 0
  function(...) {
    completed <<- completed + 1
    # Use \r to overwrite the current line with updated progress
    cat(sprintf("Completed: %d out of %d\r", completed, total))
    flush.console()  # Ensure the console output is flushed immediately
    fn(...)
  }
}
