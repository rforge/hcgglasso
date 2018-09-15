#' @title Print Values
#'
#' Print a \code{\link{MLGL}} object
#'
#'  
#' @param x \code{\link{MLGL}} object
#' @param ... Not used.
#' 
#' 
#' @method print MLGL
#' 
#' @seealso \link{MLGL} \link{summary.MLGL}
#' 
#' @export
print.MLGL <- function(x, ...)
{
  cat("$lambda\n")
  print(x$lambda)
  cat("$nVar\n")
  print(x$nVar)
  cat("$nGroup\n")
  print(x$nGroup)
}

#' @title Object Summaries
#'
#' Summary of a \code{\link{MLGL}} object
#'  
#'
#' @param object \code{\link{MLGL}} object
#' @param ... Not used.
#' 
#' @return A matrix with fitted values or estimated coefficients for given values of s.
#' 
#' @method summary MLGL
#' 
#' @seealso \link{MLGL} \link{print.MLGL}
#' 
#' @export
summary.MLGL <- function(object, ...)
{
  cat("#### MLGL\n")
  cat("## Data \n")
  cat("Number of individuals:", object$dim[1], "\n")
  cat("Number of variables:", object$dim[2], "\n")
  cat("\n")
  cat("## Hierarchical clustering \n")
  cat("HC proveded by user:", "hc"%in%names(object$call), "\n")
  cat("Time:", object$time[1],"s\n")
  cat("\n")
  cat("## Group-lasso\n")
  cat("Loss:", object$loss, "\n")
  cat("Intercept:", object$intercept,"\n")
  cat("Number of lambda:", length(object$lambda), "\n")
  cat("Number of selected variables:", head(object$nVar), "...\n")
  cat("Number of selected groups:", head(object$nGroup), "...\n")
  cat("Time:", object$time[2],"s\n")
  cat("\n")
  cat("Total elapsed time:", sum(object$time, na.rm = TRUE),"s\n")
}

#' @title Print Values
#'
#' Print a \code{\link{fullProcess}} object
#'
#'  
#' @param x \code{\link{fullProcess}} object
#' @param ... Not used.
#' 
#' 
#' @method print fullProcess
#' 
#' @seealso \link{fullProcess} \link{summary.fullProcess}
#' 
#' @export
print.fullProcess <- function(x, ...)
{
  cat("Group-lasso\n")
  cat("$res$lambda\n")
  print(x$res$lambda)
  cat("$res$nVar\n")
  print(x$resnVar)
  cat("$res$nGroup\n")
  print(x$res$nGroup)
  cat("Test output")
  cat("$lambdaOpt")
  cat(x$lambdaOpt)  
  cat("$selectedGroups")
  cat(x$selectedGroups)
}

#' @title Object Summaries
#'
#' Summary of a \code{\link{fullProcess}} object
#'  
#'
#' @param object \code{\link{fullProcess}} object
#' @param ... Not used.
#' 
#' @return A matrix with fitted values or estimated coefficients for given values of s.
#' 
#' @method summary fullProcess
#' 
#' @seealso \link{fullProcess} \link{print.fullProcess}
#' 
#' @export
summary.fullProcess <- function(object, ...)
{
  summary(object$reds)
  cat("#### Multiple Hierarchical testing\n")
  cat("## Data \n")
  cat("alpha:", object$alpha, "\n")
  cat("control:", object$control, "\n")
  cat("optimal lambda:", object$lambdaOpt, "\n")
  cat("Selected groups:", object$selectedGroups, "\n")
  cat("Selected variables:", object$var, "\n")
  cat("Time:", object$time[3],"s\n")
  cat("\n")
  cat("Total elapsed time:", sum(object$time, na.rm = TRUE),"s\n")
}