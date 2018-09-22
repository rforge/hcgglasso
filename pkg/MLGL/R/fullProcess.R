#' Run hierarchical clustering following by a group-lasso on all the different partition and a hierarchical testing procedure.
#'
#' @title Full process of MLGL
#' 
#' @author Quentin Grimonprez
#' @param X matrix of size n*p
#' @param y vector of size n. If loss = "logit", elements of y must be in {-1,1} 
#' @param hc output of \code{\link{hclust}} function. If not provided, \code{\link{hclust}} is run with ward.D2 method
#' @param control either "FDR" or "FWER"
#' @param alpha control level for testing procedure
#' @param test test used in the testing procedure. Default is partialFtest for loss = "ls" and partialChisqTest for loss = "logit"
#' @param loss a character string specifying the loss function to use, valid options are: "ls" least squares loss (regression) and "logit" logistic loss (classification)
#' @param fractionSampleMLGL a real between 0 and 1 : the fraction of individuals to use in the sample for MLGL (see Details).
#' @param ... Others parameters for MLGL
#'
#' @return a list containing :
#' \describe{
#'   \item{res}{output of \link{MLGL} function}
#'   \item{lambdaOpt}{lambda values maximizing the number of rejects}
#'   \item{var}{A vector containing the index of selected variables for the first \code{lambdaOpt} value}
#'   \item{group}{A vector containing the values index of selected groups for the first \code{lambdaOpt} value}
#'   \item{selectedGroups}{Selected groups for the first \code{lambdaOpt} value}
#'   \item{reject}{Selected groups for all lambda values}
#'   \item{alpha}{Control level}
#'   \item{test}{Test used in the testing procedure}
#'   \item{control}{"FDR" or "FWER"}
#'   \item{time}{Elapsed time}
#' } 
#'
#' @details
#' Divide the n individuals in two samples. Then the three following steps are done :
#' 1) Hierarchical CLustering of the variables of X based on the first sample of individuals
#' 2) MLGL on the second sample of individuals
#' 3) Hierarchical testing procedure on the first sample of individuals.
#'
#' @examples
#' # least square loss
#' set.seed(42)
#' X <- simuBlockGaussian(50, 12, 5, 0.7)
#' y <- drop(X[,c(2,7,12)] %*% c(2,2,-2) + rnorm(50, 0, 0.5))
#' res <- fullProcess(X, y)
#' 
#' # Logistic loss
#' y <- 2*(rowSums(X[,1:4]) > 0) - 1
#' res <- fullProcess(X, y, loss = "logit", test = partialChisqtest)
#'
#'
#' @seealso \link{MLGL}, \link{hierarchicalFDR}, \link{hierarchicalFWER}, \link{selFDR}, \link{selFWER}
#'
#' @export
fullProcess <- function(X, y, control = c("FWER", "FDR"), alpha = 0.05, test = partialFtest, hc = NULL, loss = c("ls", "logit"), fractionSampleMLGL = 1/2, ...)
{
  loss = match.arg(loss)
  if(loss == "logit" & identical(test, partialFtest))
    test = partialChisqtest
  .checkFullProcess(X, y, hc, alpha, test, fractionSampleMLGL, loss)
    
  n <- nrow(X)
  
  # Split the data in 2
  ind2 <- sample(n, floor(n * fractionSampleMLGL))
  ind1 <- (1:n)[-ind2]

  ##### part 1 : hierarchical clustering with half of the data
  if(is.null(hc))
  {
    
    # center variables and sd = 1
    Xb <- scale(X, center = TRUE, scale = FALSE)
    Xb = scale(Xb, center = FALSE, scale = sqrt(colSums(Xb^2)/n))
    
    # enclidian distance
    d <- dist(t(Xb[ind1,]))
    
    # hierarchical clustering
    hc = fastcluster::hclust(d, method = "ward.D2")
      
  }
  
  
  ##### part 2 : group-lasso
  res <- MLGL(X[ind2,], y[ind2], hc = hc, loss = loss, ...)
  
  
  ##### part 3 : testing procedure
  outTest <- HMT(res, X[ind1,], y[ind1], control, alpha, test)

  outObj <- c(list(res = res), outTest)
  class(outObj) = "fullProcess"
  
  return(outObj)
}


#' @param formula an object of class "formula" (or one that can be coerced to that class): a symbolic description of the model to be fitted. 
#' @param data an optional data frame, list or environment (or object coercible by as.data.frame to a data frame) containing the variables in the model. If not found in data, the variables are taken from environment(formula)
#'
#'
#' @rdname fullProcess
#' 
#' @export
fullProcess.formula <- function(formula, data, control = c("FWER", "FDR"), alpha = 0.05, test = partialFtest, hc = NULL, loss = c("ls", "logit"), fractionSampleMLGL = 1/2, ...)
{
  cl <- match.call()
  mf <- match.call(expand.dots = FALSE)
  m <- match(c("formula", "data"), names(mf), 0L)
  mf <- mf[c(1L, m)]
  mf$drop.unused.levels <- TRUE
  mf[[1L]] <- quote(stats::model.frame)
  mf <- eval(mf, parent.frame())
  mt <- attr(mf, "terms")
  
  y <- model.response(mf, "numeric")
  X <- model.matrix(mt, mf)
  X = as.matrix(X)
  
  res <- fullProcess(X, y, control, alpha, test, hc, loss, fractionSampleMLGL, ...)
  
  return(res)
}

# check parameters of MLGL function
.checkFullProcess <- function(X, y, hc, alpha, test, fractionSampleMLGL, loss)
{
  #check X
  if(!is.matrix(X)) 
    stop("X has to be a matrix.")
  if(any(is.na(X))) 
    stop("Missing values in X not allowed.")
  if(!is.numeric(X))
    stop("X has to be a matrix of real.")
  
  #check y
  if(!is.numeric(y))
    stop("y has to be a vector of real.")
  if(any(is.na(y))) 
    stop("Missing values in y not allowed.")
  if (loss == "logit" && any(y %in% c(-1, 1) == FALSE)) 
    stop("Classification method requires the response y to be in {-1,1}")
  
  #check if X and y are compatible
  if(nrow(X)!=length(drop(y)))
    stop("The length of y and the number of rows of X don't match.")
  
  #check hc
  if(!is.null(hc))
  {
    #check if hc is a hclust object
    if(class(hc)!="hclust")
      stop("hc must be an hclust object.")
    #check if hc and X are compatible
    if(length(hc$order)!=ncol(X))
      stop("hc is not a clustering of the p covariates of X.")
    
  }
  
  #alpha
  if(length(alpha)!=1)
    stop("alpha must be a real between 0 and 1.")
  if((alpha <=0) || (alpha>1))
    stop("alpha must be a real between 0 and 1.")
  
  # check if test is a function
  if(!is.function(test))
    stop("test must be a function.")
  
  # fractionSampleMLGL
  if(length(fractionSampleMLGL)!=1)
    stop("fractionSampleMLGL must be a real between 0 and 1.")
  if((fractionSampleMLGL <=0) || (fractionSampleMLGL>=1))
    stop("fractionSampleMLGL must be a real between 0 and 1.")
  
  invisible(return(NULL))
}
