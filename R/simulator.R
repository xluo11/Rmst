#' Simulate the Administration of Multistage Tests
#' @name sim
#' @examples
#' ## assemble a MST
#' items <- Rirt::model_mixed_gendata(1, n_3pl=150)$items
#' x <- mst(items, "1-3", 2, 'topdown', len=20, max_use=1)
#' x <- mst_objective(x, -1, indices=1)
#' x <- mst_objective(x,  0, indices=2)
#' x <- mst_objective(x,  1, indices=3)
#' x <- mst_stage_length(x, 1:2, min=5)
#' x <- mst_assemble(x, 'lpsolve', time_limit=2)
#'
#' sim1 <- mst_sim(x, true=.5)
#' print(sim1)
#' plot(sim1)
#'
#' sim2 <- mst_sim(x, true=-.5, rdp=list('stage2'=c(-.44, .44)))
#' print(sim2)
#' plot(sim2)
NULL

#' @rdname sim
#' @description \code{mst_sim} simulates the administration of the assembled MST panel(s)
#' @param x the assembled MST object
#' @param true the true theta parameter (numeric)
#' @param rdp routing decision points (list)
#' @param estimator the estimator of the ability parameter (function)
#' @param ... additional option/control parameters
#' @details
#' Use \code{theta} to set the initial theta, \code{panel} to select the MST panel,
#' \code{prior} to set the prior for theta estimation, \code{bounds} to set the
#' bounds of theta estimation, and \code{D} to set the scaling constant.
#' @importFrom stats rnorm rbinom rmultinom
#' @importFrom Rirt model_mixed_eap model_mixed_prob model_mixed_info freq
#' @export
mst_sim <- function(x, true, rdp=NULL, estimator=model_mixed_eap, ...){
  if(is.null(x$items))
    stop("the mst has not been assembled yet")

  opts <- list(...)
  if(is.null(opts$D)) D <- 1.702 else D <- opts$D
  if(is.null(opts$prior)) prior <- c(0, 1) else prior <- opts$prior
  if(is.null(opts$bounds)) bounds <- c(-4, 4) else bounds <- opts$bounds
  if(is.null(opts$panel)) panel <- sample(x$n_panels, 1) else panel <- opts$panel
  if(is.null(opts$theta)) theta <- rnorm(1, 0, .1) else theta <- opts$theta

  admin <- list('3pl'=NULL, 'gpcm'=NULL, 'grm'=NULL)
  stats <- matrix(nrow=x$n_stages, ncol=4, dimnames=list(NULL, c("route", "t", "info", "se")))

  if(!is.null(rdp)) {
    if(length(rdp) != x$n_stages - 1)
      stop("Invalid routing decision points.")
    rdp <- lapply(rdp, function(x) data.frame(lower=c(-Inf, x), upper=c(x, Inf)))
    rdp <- Reduce(rbind, rdp)
    rdp$index <- 2:x$n_modules
  }

  # MST administration
  for(i in 1:x$n_stages){
    # select module
    if(i == 1) {
      next_module <- unique(x$route[, i])
      next_module <- sample(next_module, 1)
    } else {
      next_module <- x$route[, i-1] == stats[i-1, 'route']
      next_module <- x$route[next_module, i]
      next_module <- sort(unique(next_module))
      if(is.null(rdp)) {
        info <- model_mixed_info(theta, mst_get_items(x, panel_ix=panel), D=D, combine=TRUE)[1,]
        index <- unlist(Map(function(x) x$index, mst_get_items(x, panel_ix=panel)))
        info <- aggregate(info, by=list(module=index), sum)
        info <- info[info$module %in% next_module, ]
        next_module <- info$module[which.max(info$x)]
      } else {
        next_module <- rdp[rdp$index %in% next_module, ]
        next_module$lower[1] <- -Inf
        next_module$upper[nrow(next_module)] <- Inf
        next_module <- min(next_module[theta < next_module$upper, 'index'])
      }
    }

    # generate responses
    items <- mst_get_items(x, panel_ix=panel, stage_ix=i, module_ix=next_module)
    rsp <- Map(function(x) {
      if(is.null(x))
        return(NULL)
      if(length(dim(x)) == 2)
        return(rbinom(ncol(x), 1, x[1,]))
      if(length(dim(x)) == 3)
        return(apply(x, 1:2, function(x) which.max(rmultinom(1, 1, x)[,1]) - 1)[1,])
    }, model_mixed_prob(true, items, D=D))
    admin <- Map(function(x, y, z) rbind(z, cbind(x, rsp=y)), items, rsp, admin)

    # estimate ability
    u <- matrix(Reduce(c, Map(function(x) x$rsp, admin)), nrow=1)
    u <- rbind(u, u)
    theta <- estimator(u, admin, D=D, priors=prior, bounds_t=bounds)$t[1]
    info <- sum(model_mixed_info(theta, admin, D=D, combine=TRUE))
    se <- 1/sqrt(info)
    stats[i, c('route', 't', 'info', 'se')] <- c(next_module, theta, info, se)
  }

  stats <- data.frame(stats, stringsAsFactors=FALSE)
  stats$n_items <- freq(unlist(sapply(admin, function(x) x$index)), stats$route)$freq
  rs <- list(admin=admin, stats=stats, true=true, theta=theta)
  class(rs) <- "mst_sim"
  rs
}




