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
  cat("Total elapsed time:", sum(object$time),"s\n")
}