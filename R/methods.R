#' @rdname assembly
#' @export
print.mst <- function(x, ...){
  cat("The MST design has", x$n_stages, "stages,", x$n_modules, "modules, and", x$n_routes, "routes.\n")

  cat("\nRoutes:\n")
  cat('----------------------------------------------\n')
  print(x$route)

  cat("\nModules:\n")
  cat('----------------------------------------------\n')
  print(x$module)

  cat("\nATA Results:\n")
  cat('----------------------------------------------\n')
  print(x$ata)

  invisible(x)
}


#' @rdname assembly
#' @details
#' \code{plot.mst} draws module information functions when \code{byroute=FALSE}
#' and route information functions when \code{byroute=TRUE}. Use \code{label=TRUE}
#' to put labels on routes and modules.
#' @import ggplot2
#' @export
plot.mst <- function(x, ...){
  if(class(x) != "mst")
    stop("Not a 'mst' object: ", class(x))

  if(is.null(x$items))
    stop('The MST has not been assembled yet.')

  opts <- list(...)
  if(is.null(opts$byroute))
    opts$byroute <- FALSE
  if(is.null(opts$theta))
    opts$theta <- round(seq(-3, 3, .1), 1)
  if(is.null(opts$label))
    opts$label <- FALSE

  data <- NULL
  if(opts$byroute) {
    # compute the route information functions
    for(i in 1:x$n_routes)
      for(j in 1:x$n_panels){
        items <- mst_get_items(x, panel_ix=j, route_ix=i)
        info <- model_mixed_info(opts$theta, items, D=x$ata$opts$D, combine=TRUE)
        info <- data.frame(panel=j, index=i, theta=opts$theta, info=rowSums(info, na.rm=TRUE))
        data <- rbind(data, info)
      }
    # format
    data$panel <- paste("Panel", data$panel)
    data <- merge(data, x$route, by='index', all.x=TRUE)
    if(opts$label)
      for(k in 1:x$n_stages) {
        label <- x$module$label[match(data[,paste('stage', k, sep='')], x$module$index)]
        data[,paste('stage', k, sep='')] <- paste(k, label, sep='')
      }
    data$route <- apply(data[, paste('stage', 1:x$n_stages, sep='')], 1, paste, collapse='-')
    # produce ggplot
    g <- ggplot(data, aes_string(x="theta", y="info", color="route")) +
      geom_line() + xlab(expression(theta)) + ylab("Route Information Functions") +
      theme_bw() + theme(legend.key=element_blank()) +
      guides(color=guide_legend("Routes")) +
      facet_grid(. ~ panel)
  } else {
    # compute the module information functions
    for(i in 1:x$n_panels)
      for(j in 1:x$n_modules){
        items <- mst_get_items(x, panel_ix=i, module_ix=j)
        info <- model_mixed_info(opts$theta, items, D=x$ata$opts$D, combine=TRUE)
        info <- data.frame(panel=i, index=j, theta=opts$theta, info=rowSums(info, na.rm=TRUE))
        data <- rbind(data, info)
      }
    data <- merge(data, x$module, by='index', all.x=TRUE)
    if(opts$label) {
      data$module <- paste(data$stage, data$label, sep="")
    } else {
      data$module <- factor(as.integer(data$index))
    }
    data$panel <- paste("Panel", data$panel)
    data$stage <- paste("Stage", data$stage)
    g <- ggplot(data, aes_string(x="theta", y="info", color="module")) +
      geom_line() + xlab(expression(theta)) + ylab("Module Information Functions") +
      theme_bw() + theme(legend.key=element_blank()) +
      guides(color=guide_legend("Modules")) +
      facet_grid(panel ~ stage)
  }

  g
}


#' @rdname sim
#' @export
print.mst_sim <- function(x, ...){
  cat("true=", round(x$true, 2), ", est.=", round(x$theta, 2), ":\n", sep="")
  print(round(x$stats, 2))
  cat("Call x$admin to see administered items ('x' is the mst_sim object).\n")
}


#' @rdname sim
#' @importFrom stats qnorm
#' @importFrom Rirt model_mixed_info
#' @import ggplot2
#' @export
plot.mst_sim <- function(x, ...) {
  opts <- list(...)
  if(is.null(opts$ci_width)) opts$ci_width <- qnorm(.975)
  if(is.null(opts$ylim)) opts$ylim <- c(-3, 3)

  if(!is.null(x$admin$'grm')){
    cols <- grep('^b[0-9]+$', colnames(x$admin$'grm'), value=TRUE)
    x$admin$'grm'$b <- rowMeans(x$admin$'grm'[,cols], na.rm=TRUE)
    x$admin$'grm'[,cols] <- x$admin$'grm'[,cols] - x$admin$'grm'$b
  }
  data <- Reduce(rbind, Map(function(x) x[, c('stage', 'rsp', 'a', 'b')], x$admin))
  data <- data[order(data$stage), ]
  data$Position <- 1:nrow(data)
  data$Scores <- factor(data$rsp, levels=sort(unique(data$rsp)))
  x$stats$lb <- x$stats$t - opts$ci_width * x$stats$se
  x$stats$ub <- x$stats$t + opts$ci_width * x$stats$se
  x$stats$position <- cumsum(x$stats$n_items)

  ggplot(data, aes_string(x="Position", y="b")) +
    geom_point(aes_string(size="a", color="Scores")) +
    geom_pointrange(data=x$stats, aes_string(x="position", y="t", ymin="lb", ymax="ub"), lty=2, pch=4, col="coral") +
    xlab("Position") + ylab("Item Difficulty") + guides(size=F, fill=F) +
    coord_cartesian(ylim=opts$ylim) + scale_size_continuous(range=c(1, 3)) + theme_bw()
}
