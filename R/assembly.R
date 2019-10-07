#' Assemble Computerized Adaptive Multistage Testing
#' @name assembly
#' @examples
#' ## generate item pool
#' items <- Rirt::model_mixed_gendata(1, n_3pl=200)$items
#'
#' ## Ex. 1: 1-2-2 MST, 2 panels, 20 items, topdown
#' ## maximize info. at -1 and 1 for easy and hard routes
#' x <- mst(items, "1-2-2", n_panels=2, method='topdown', test_len=10, max_use=1)
#' x <- mst_objective(x, -1, indices=1:2)
#' x <- mst_objective(x, 1, indices=3:4)
#' x <- mst_assemble(x, 'lpsolve', time_limit=2)
#' plot(x, byroute=TRUE, label=TRUE)
#'
#' \donttest{
#' ## Ex. 2: 1-2-3 MST, 2 panels, bottomup,
#' ## remove two routes with large theta change: 1-2-6, 1-3-4
#' ## 10 items in each module, content= and 3 items in content area 1 in each module
#' ## maximize info. at -1, 0 and 1 for easy, medium, and hard modules
#' x <- mst(items, "1-2-3", 1, 'bottomup', len=10, max_use=1)
#' x <- mst_route(x, c(1, 2, 6), "-")
#' x <- mst_route(x, c(1, 3, 4), "-")
#' x <- mst_objective(x,  0, indices=c(1, 5))
#' x <- mst_objective(x, -1, indices=c(2, 4))
#' x <- mst_objective(x,  1, indices=c(3, 6))
#' x <- mst_assemble(x, timeout=30)
#' plot(x, byroute=FALSE)
#' plot(x, byroute=TRUE)
#' }
NULL


#' @rdname assembly
#' @description \code{mst} creates a multistage (MST) assembly model
#' @param pool the item pool (a list of '3pl', 'gpcm', and 'grm' items)
#' @param design the MST design (string): e.g., "1-3", "1-2-2", "1-2-3"
#' @param n_panels the number of panels (integer)
#' @param method the design method (string): 'topdown' or 'bottomup'
#' @param test_len the module/route length (integer)
#' @param max_use the maximum selection of items (integer)
#' @param ... additional arguments
#' @return \code{mst} returns a \code{mst} object.
#' @details
#' A \code{mst} object stores the definitions of the MST. When \code{mst_assemble}
#' is called, definitions are converted to a real mixed integer programming model for
#' assembly. If the model is solved, assembled items are appended to the original object. \cr
#' The bottom-up approach adds objectives and constraints on individual modules, and the top-down
#' approach adds objectives and constraints on routes.
#' @importFrom Rata ata
#' @export
mst <- function(pool, design, n_panels=1, method=c('topdown', 'bottomup'), test_len=NULL, max_use=NULL, ...){
  method <- match.arg(method)
  design <- as.integer(unlist(strsplit(design, split="[-, /]")))
  n_stages <- length(design)
  n_modules <- sum(design)

  # model labels
  labels <- list('1'=c('M'),
                 '2'=c('E', 'H'),
                 '3'=c('E', 'M', 'H'),
                 '4'=c('XE', 'E', 'H', 'XH'),
                 '5'=c('XE', 'E', 'M', 'H', 'XH'))

  # module-index map
  module <- data.frame()
  if(max(design) <= length(labels)) {
    for(s in 1:n_stages)
      module <- rbind(module, data.frame(stage=s, module=1:design[s], label=labels[[design[s]]], stringsAsFactors=FALSE))
  } else {
    for(s in 1:n_stages)
      module <- rbind(module, cbind(stage=s, module=1:design[s]))
  }
  module$index <- 1:nrow(module)

  # route-index map
  route <- expand.grid(split(module$index, module$stage))
  colnames(route) <- paste('stage', 1:n_stages, sep='')
  priority <- apply(route[, 1:n_stages], 1, function(r) sum(r * 10^(n_stages:1)))
  route <- route[order(priority), ]
  route$index <- 1:nrow(route)
  n_routes <- nrow(route)

  # ata
  x <- list(pool=pool, design=design, method=method, module=module, route=route,
            n_panels=n_panels, n_stages=n_stages, n_modules=n_modules, n_routes=n_routes,
            ata=ata(pool, n_forms=n_panels*n_modules, max_use=max_use, ...))
  class(x) <- "mst"

  # constraint: test length
  if(!is.null(test_len)) {
    if(length(test_len) > 2)
      stop("Invalid test_len: expect 1 or 2 elements")
    if(length(test_len) == 1)
      test_len <- rep(test_len, 2)
    x <- mst_constraint(x, 1, min=test_len[1], max=test_len[2])
  }

  # constraint: at least 1 item per stage
  x <- mst_stage_length(x, 1:n_stages, min=1)

  # constraint: no item reuse in the same route
  mat <- matrix(0, nrow=x$n_panels*x$n_routes*x$ata$n_items, ncol=x$ata$n_lpvar)
  for(p in 1:x$n_panels)
    for(r in 1:x$n_routes) {
      idx <- outer(1:x$ata$n_items, (unlist(x$route[r, 1:x$n_stages]) - 1) * x$ata$n_items, "+")
      idx <- idx + (p-1) * x$n_modules * x$ata$n_items
      for(j in 1:ncol(idx))
        mat[cbind(1:x$ata$n_items+(r-1)*x$ata$n_items+(p-1)*x$n_routes*x$ata$n_items, idx[, j])] <- 1
    }
  x$ata$mat <- rbind(x$ata$mat, mat)
  x$ata$dir <- c(x$ata$dir, rep("<=", nrow(mat)))
  x$ata$rhs <- c(x$ata$rhs, rep(1, nrow(mat)))

  x
}


#' @rdname assembly
#' @description \code{mst_route} adds/removes a route to/from the assembly model
#' @param x the MST object
#' @param route a MST route (a vector of module index)
#' @param op \code{"+"} to add a route and \code{"-"} to remove a route
#' @export
mst_route <- function(x, route, op=c("+", "-")){
  if(class(x) != "mst")
    stop("not a 'mst' object: ", class(x))

  op <- match.arg(op)
  idx <- apply(x$route[, 1:x$n_stages], 1, function(r) all(r == route))

  if(op == "+") {
    if(any(idx))
      stop("The route already exists")
    if(!all(route %in% 1:x$n_modules))
      stop("Invalid route: module index is out of bound.")
    x$route <- rbind(x$route, c(route, NA))
  } else if(op == "-") {
    if(!any(idx))
      stop("The route hasn't been added yet")
    x$route <- x$route[!idx, ]
  }

  # re-index
  priority <- apply(x$route[, 1:x$n_stages], 1, function(r) sum(r * 10^(x$n_stages:1)))
  x$route <- x$route[order(priority), ]
  x$route$index <- 1:nrow(x$route)
  x$n_routes <- nrow(x$route)

  x
}


#' @rdname assembly
#' @description \code{mst_objective} adds an objective to the assembly model
#' @param coef the coefficients (numeric vector or string)
#' @param mode the optimization direction: \code{"max"} or \code{"min"}
#' @param indices the indices of the route (topdown) or the module (bottomup) where the objective is added
#' @param target the target values of the absolute objectives, \code{NULL} for the relative objective
#' @details
#' \code{coef} in \code{mst_objective} can be a vector of theta points where TIFs are optimized, or
#' a continuous variable in the pool where the item attribute is optimized, or a numeric value with
#' the same length of the pool at either item or group level.
#' @importFrom Rata ata_relative_objective ata_absolute_objective
#' @export
mst_objective <- function(x, coef, mode=c("max", "min"), indices=NULL, target=NULL, method=NULL, ...) {
  if(class(x) != "mst")
    stop("not a 'mst' object: ", class(x))

  indices <- mst_get_form_index(x, indices, method)
  for(i in 1:x$n_panels) {
    for(j in 1:nrow(indices)) {
      f <- unlist(indices[j, ]) + (i - 1) * x$n_modules
      if(is.null(target) || is.na(target)) {
        x$ata <- ata_relative_objective(x$ata, coef, mode=mode, forms=f, collapse=TRUE, ...)
      } else {
        x$ata <- ata_absolute_objective(x$ata, coef, target=target, forms=f, collapse=TRUE, ...)
      }
    }
  }

  x
}

#' @rdname assembly
#' @description \code{mst_constraint} adds constraints to the assembly model
#' @param level the constrained level of categorical item attribute,
#' \code{NULL} for continuous item attributes
#' @param min the lower bound of the constraint
#' @param max the upper bound of the constraint
#' @importFrom Rata ata_constraint
#' @export
mst_constraint <- function(x, coef, min=NA, max=NA, level=NULL, indices=NULL, method=NULL){
  if(class(x) != "mst")
    stop("not a 'mst' object: ", class(x))

  indices <- mst_get_form_index(x, indices, method)
  for(i in 1:x$n_panels){
    for(j in 1:nrow(indices)){
      f <- unlist(indices[j,] + (i - 1) * x$n_modules)
      x$ata <- ata_constraint(x$ata, coef, min, max, level, forms=f, collapse=TRUE)
    }
  }

  x
}


#' @rdname assembly
#' @description \code{mst_stage_length} sets length limits on stages
#' @param stages the stage indices
#' @importFrom Rata ata_constraint
#' @export
mst_stage_length <- function(x, stages, min=NA, max=NA){
  if(class(x) != "mst")
    stop("not a 'mst' object: ", class(x))

  if(length(min) == 1)
    min <- rep(min, length(stages))
  if(length(max) == 1)
    max <- rep(max, length(stages))
  if(length(stages) != length(min) || length(stages) != length(max))
    stop("Different lengths of stage, min and max")

  for(i in 1:length(stages)){
    if(!stages[i] %in% 1:x$n_stages)
      stop("Invalid stage: out of bound")
    f <- x$module[x$module$stage == stages[i], "index"]
    f <- as.vector(outer(f, (1:x$n_panels-1)*x$n_modules, "+"))
    x$ata <- ata_constraint(x$ata, 1, min[i], max[i], forms=f, collapse=FALSE)
  }

  x
}


#' @rdname assembly
#' @description \code{mst_rdp} anchors the routing decision point (rdp) between adjacent modules
#' @param tol tolerance parameter (numeric)
#' @importFrom stats aggregate
#' @export
mst_rdp <- function(x, theta, indices, tol=.5) {
  if(class(x) != "mst")
    stop("Not a 'mst' object: ", class(x))

  if(length(theta) != 1)
    stop("Set rdp on one theta point each time")

  if(length(indices) != 2 || abs(indices[1] - indices[2]) != 1)
    stop("Modules are not adjacent")

  items <- Map(function(x) if(nrow(x) == 0) NULL else x, x$pool)
  coef <- model_mixed_info(theta, items, D=x$ata$opts$D, combine=TRUE)
  mat <- matrix(0, nrow=x$n_panels, ncol=x$ata$n_lpvar)
  for(i in 1:x$n_panels) {
    idx <- 1:x$ata$n_items + (indices[1] - 1) * x$ata$n_items + (i - 1) * x$n_modules * x$ata$n_items
    mat[i, idx] <- coef
    idx <- 1:x$ata$n_items + (indices[2] - 1) * x$ata$n_items + (i - 1) * x$n_modules * x$ata$n_items
    mat[i, idx] <- -1 * coef
  }
  x$ata$mat <- rbind(x$ata$mat, mat, mat)
  x$ata$dir <- c(x$ata$dir, rep(">=", x$n_panels), rep("<=", x$n_panels))
  x$ata$rhs <- c(x$ata$rhs, rep(-tol, x$n_panels), rep(tol, x$n_panels))

  x
}


#' @rdname assembly
#' @description \code{mst_module_info} sets the information requirements for modules
#' @param theta the theta point where TIF is controlled
#' @importFrom stats aggregate
#' @importFrom Rirt model_mixed_info
#' @export
mst_module_info <- function(x, theta, min=NA, max=NA, indices) {
  if(class(x) != "mst")
    stop("Not a 'mst' object: ", class(x))

  if(is.na(min) && is.na(max))
    return(x)

  if(any(indices < 1 | indices > x$n_modules))
    stop("Invalid module index")

  if(length(theta) != 1)
    stop("Control module information function at one theta point each time")

  items <- Map(function(x) if(nrow(x) == 0) NULL else x, x$pool)
  coef <- model_mixed_info(theta, items, D=x$ata$opts$D, combine=TRUE)
  mat <- matrix(0, nrow=x$n_panels*length(indices), ncol=x$ata$n_lpvar)
  for(p in 1:x$n_panels) {
    idx <- outer(1:x$ata$n_items, (indices - 1) * x$ata$n_items, "+")
    idx <- idx + (p - 1) * x$n_modules * x$ata$n_items
    for(i in 1:ncol(idx))
      mat[i + (p - 1) * length(indices), idx[,i]] <- coef
  }

  if(!is.na(min) && is.na(max)) {
    x$ata$mat <- rbind(x$ata$mat, mat)
    x$ata$dir <- c(x$ata$dir, rep(">=", x$n_panels * length(indices)))
    x$ata$rhs <- c(x$ata$rhs, rep(min, x$n_panels * length(indices)))
  } else if(is.na(min) && !is.na(max)) {
    x$ata$mat <- rbind(x$ata$mat, mat)
    x$ata$dir <- c(x$ata$dir, rep("<=", x$n_panels * length(indices)))
    x$ata$rhs <- c(x$ata$rhs, rep(max, x$n_panels * length(indices)))
  } else {
    x$ata$mat <- rbind(x$ata$mat, mat, mat)
    x$ata$dir <- c(x$ata$dir, rep(">=", x$n_panels * length(indices)))
    x$ata$rhs <- c(x$ata$rhs, rep(min, x$n_panels * length(indices)))
    x$ata$dir <- c(x$ata$dir, rep("<=", x$n_panels * length(indices)))
    x$ata$rhs <- c(x$ata$rhs, rep(max, x$n_panels * length(indices)))
  }

  x
}


#' @rdname assembly
#' @description \code{mst_assemble} tries to solve the assembly model
#' @param solver the MIP solver: \code{"lpsolve"} or \code{"glpk"}
#' @param silent \code{TRUE} to mute solving status
#' @param time_limit the time limit for solving the model in seconds
#' @param message \code{TRUE} to print messages from the solver
#' @importFrom Rata ata_solve
#' @export
mst_assemble <- function(x, solver=c('lpsolve', 'glpk'), silent=FALSE, time_limit=30, message=FALSE, ...){
  if(class(x) != "mst")
    stop("not a 'mst' object: ", class(x))

  x$ata <- ata_solve(x$ata, solver=solver, return_format='model', silent=silent, time_limit=time_limit, message=message, ...)
  if(!is.null(x$ata$items))
    x$items <- Map(function(f) {
      f$panel <- ceiling(f$form / x$n_modules)
      f$index <- ((f$form - 1) %% x$n_modules) + 1
      f <- merge(f, x$module, by='index', all.x=TRUE)
      f$form <- NULL
      f
    }, x$ata$items)

  x
}


#' @rdname assembly
#' @description \code{mst_get_items} retrieves items from the assembly results
#' @param panel_ix the panel index (int vector)
#' @param stage_ix the stage index (int vector)
#' @param module_ix the module index (int vector)
#' @param route_ix the route index (int vector)
#' @return \code{mst_get_items} returns the assembled forms in a list of 3pl, gpcm, and grm items
#' @export
mst_get_items <- function(x, panel_ix=NULL, stage_ix=NULL, module_ix=NULL, route_ix=NULL){
  if(class(x) != "mst")
    stop("Not a 'mst' object: ", class(x))

  if(is.null(x$items))
    stop('The mst has not been assembled yet.')

  Map(function(items) {
    if(!is.null(panel_ix))
      items <- subset(items, items$panel %in% panel_ix)
    if(!is.null(stage_ix))
      items <- subset(items, items$stage %in% stage_ix)
    if(!is.null(module_ix))
      items <- subset(items, items$index %in% module_ix)
    if(!is.null(route_ix))
      items <- subset(items, items$index %in% unlist(x$route[x$route$index == route_ix, 1:x$n_stages]))
    if(nrow(items) == 0)
      items <- NULL
    items
  }, x$items)
}

