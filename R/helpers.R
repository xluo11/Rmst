#' Helper functions for MST Assembly and Simulation
#' @description helper functions for MST assembly and simulation
#' @name helpers
NULL


#' @rdname helpers
#' @description \code{mst_get_form_index} maps the input form indices to the actual indices
#' @param x the MST object
#' @param indices the input form indices
#' @param method the assembly method: topdown or bottomup
#' @keywords internal
mst_get_form_index <- function(x, indices, method){
  if(is.null(method))
    method <- x$method
  if(method == 'topdown'){
    if(is.null(indices))
      indices <- 1:x$n_routes
    indices <- x$route[x$route$index %in% indices, 1:x$n_stages, drop=FALSE]
  } else if(method == 'bottomup') {
    if(is.null(indices))
      indices <-1:x$n_modules
    indices <- x$module[x$module$index %in% indices, 'index', drop=FALSE]
  } else {
    stop('Invalid method')
  }
  indices
}
