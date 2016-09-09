#'
#' Obtain a sparse matrix of the coefficients of the path
#' 
#' @param x \code{\link{MLGL}} object
#' @param row "lambda" or "covariates". If row="covariates", each row of the output matrix represents a covariates else ff row="lambda", it represents a value of lambda.
#'
#' @return a sparse matrix containing the estimated coefficients for different values of lambda
#'
#' @details This functions can be used with a \code{\link{MLGL}} object to obtain a matrix with all estimated coefficients for the p original variables.
#' In case of overlapping groups, coefficients from repeated variables are summed. 
#'
#' @examples 
#' # Simulate gaussian data with block-diagonal variance matrix containing 12 blocks of size 5
#' X <- simuBlockGaussian(50, 12, 5, 0.7)
#' # Generate a response variable
#' y <- drop(X[,c(2,7,12)]%*%c(2,2,-2)+rnorm(50,0,0.5))
#' # Apply MLGL method
#' res <- MLGL(X,y)
#' # Convert output in sparse matrix format
#' beta <- listToMatrix(res)
#'
#' @seealso \link{MLGL}, \link{overlapgglasso}
#'
#' @export
listToMatrix <- function(x, row = c("covariates","lambda"))
{
  row = match.arg(row)
  if(row == "covariates")
  {
    bet <- Matrix(0, ncol = length(x$lambda), nrow = x$dim[2])
    for(i in 1:length(x$lambda))
    {
      dup <- duplicated(x$var[[i]])
      bet[x$var[[i]][!dup], i] = bet[x$var[[i]][!dup],i] + x$beta[[i]][!dup]
      bet[x$var[[i]][dup] , i] = bet[x$var[[i]][dup] ,i] + x$beta[[i]][dup]
    }
    
    return(bet)
  }
  else
  {
    bet <- Matrix(0, nrow = length(x$lambda), ncol = x$dim[2])
    for(i in 1:length(x$lambda))
    {
      dup <- duplicated(x$var[[i]])
      bet[i,x$var[[i]][!dup]] = bet[i,x$var[[i]][!dup]] + x$beta[[i]][!dup]
      bet[i,x$var[[i]][dup]]  = bet[i,x$var[[i]][dup]]  + x$beta[[i]][dup]
    }
    
    return(bet)
  }

}

#'
#' Plot the path obtained from \code{\link{MLGL}} function
#'
#' @param x \code{\link{MLGL}} object
#' @param log.lambda If TRUE, use log(lambda) instead of lambda in abscissa
#' @param lambda.lines if TRUE, add vertical lines at lambda values
#' @param ... Other parameters for plot function
#' 
#' @examples
#' \dontrun{
#' # Simulate gaussian data with block-diagonal variance matrix containing 12 blocks of size 5
#' X <- simuBlockGaussian(50, 12, 5, 0.7)
#' # Generate a response variable
#' y <- drop(X[,c(2,7,12)]%*%c(2,2,-2)+rnorm(50,0,0.5))
#' # Apply MLGL method
#' res <- MLGL(X,y)
#' # Plot the solution path
#' plot(res) 
#' }
#' 
#' @method plot MLGL
#' 
#' @seealso \link{MLGL}
#' 
#' @export
plot.MLGL <- function(x, log.lambda = FALSE, lambda.lines = FALSE,...)
{
  #check log
  if(length(log.lambda)!=1)
    stop("log must be a boolean.")
  if(!is.logical(log.lambda))
    stop("log must be a boolean.")
  #check log
  if(length(lambda.lines)!=1)
    stop("lambda.lines must be a boolean.")
  if(!is.logical(lambda.lines))
    stop("lambda.lines must be a boolean.") 
  
  bet <- listToMatrix(x,"lambda")
  
  #abscissa : log or not ?
  absc <- x$lambda
  if(log.lambda)
    absc = log(absc)
  
  #plot
  matplot(absc, bet, type = "l", lty = 1, xlab = ifelse(log.lambda, expression(paste("log(",lambda,")")), expression(lambda)),
          ylab = "Coefficients",...)
  
  #add vertical lines for lambda values
  if(lambda.lines)
    abline(v = absc, col = "blue", lty = "dotted", lwd = 0.5)
  
}


#'
#' Plot the cross-validation obtained from \code{\link{cv.MLGL}} function
#'
#' @param x \code{\link{cv.MLGL}} object
#' @param log.lambda If TRUE, use log(lambda) instead of lambda in abscissa
#' @param ... Other parameters for plot function
#' 
#' @examples
#' \dontrun{
#' set.seed(42)
#' # Simulate gaussian data with block-diagonal variance matrix containing 12 blocks of size 5
#' X <- simuBlockGaussian(50, 12, 5, 0.7)
#' # Generate a response variable
#' y <- drop(X[,c(2,7,12)]%*%c(2,2,-2)+rnorm(50,0,0.5))
#' # Apply cv.MLGL method
#' res <- cv.MLGL(X,y)
#' # Plot the cv error curve
#' plot(res) 
#' }
#' 
#' @method plot cv.MLGL
#' 
#' @seealso \link{cv.MLGL}
#' 
#' @export
plot.cv.MLGL <- function(x, log.lambda = FALSE,...)
{
  #check log
  if(length(log.lambda)!=1)
    stop("log must be a boolean.")
  if(!is.logical(log.lambda))
    stop("log must be a boolean.")
  
  #abscissa : log or not ?
  absc <- x$lambda
  lam <- c(x$lambda.min,x$lambda.1se)
  if(log.lambda)
  {
    absc = log(absc)
    lam = log(lam)
  }
  
  #plot
  matplot(absc, cbind(x$cvm,x$cvupper,x$cvlower), type = "l", lty = c(1,2,2), col = c(1,2,2),
          xlab = ifelse(log.lambda, expression(paste("log(",lambda,")")), expression(lambda)), ylab = "Error",...)
  abline(v = lam, col = "blue", lty = "dashed")
}


#'
#' Plot the stability path obtained from \code{\link{stability.MLGL}} function
#'
#' @param x \code{\link{stability.MLGL}} object
#' @param log.lambda If TRUE, use log(lambda) instead of lambda in abscissa
#' @param threshold Threshold for selection frequency
#' @param ... Other parameters for plot function
#' 
#' @return A list containing :
#' \describe{
#' \item{var}{Index of selected variables for the given threshold.}
#' \item{group}{Index of the associated group.}
#' \item{threshold}{Value of threshold}
#' } 
#' 
#' @examples
#' \dontrun{
#' set.seed(42)
#' # Simulate gaussian data with block-diagonal variance matrix containing 12 blocks of size 5
#' X <- simuBlockGaussian(50, 12, 5, 0.7)
#' # Generate a response variable
#' y <- drop(X[,c(2,7,12)]%*%c(2,2,-2)+rnorm(50,0,0.5))
#' # Apply stability.MLGL method
#' res <- stability.MLGL(X,y)
#' selected <- plot(res)
#' print(selected)
#' }
#' 
#' @method plot stability.MLGL
#' 
#' @seealso \link{stability.MLGL}
#' 
#' @export
plot.stability.MLGL <- function(x, log.lambda = FALSE, threshold = 0.75,...)
{
  #check log
  if(length(log.lambda)!=1)
    stop("log must be a boolean.")
  if(!is.logical(log.lambda))
    stop("log must be a boolean.")
  #threshold
  if( !is.numeric(threshold) | (length(threshold)!=1) )
    stop("threshold must be a positive real lesser than 1.")
  if( (threshold<0) | (threshold>1) ) 
    stop("threshold must be a positive real lesser than 1.")
  
  #abscissa : log or not ?
  absc <- x$lambda
  lam <- c(x$lambda.min,x$lambda.1se)
  if(log.lambda)
  {
    absc = log(absc)
    lam = log(lam)
  }
  
  #determine color according to threshold
  col <- apply(x$stability, 2, FUN=function(x){ifelse(any(x>threshold), 2, 1)})
    
  #plot
  matplot(absc, x$stability, type = "l", lty = 1,col = col,
          xlab = ifelse(log.lambda, expression(paste("log(",lambda,")")), expression(lambda)), ylab = "Probability selection",...)
  abline(v = lam,col = "blue",lty = "dashed")
  
  #determine selected groups and variables
  selectedGroup <- which(col==2)
  indsel <- x$group%in%selectedGroup
  return(list(var = x$var[indsel], group = x$group[indsel], threshold = threshold))
}

