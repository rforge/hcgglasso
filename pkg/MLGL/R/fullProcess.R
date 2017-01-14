#' Run hierarchical clustering following by a group-lasso on all the different partition and a hierarchical testing procedure.
#'
#' @title Full process of MLGL
#' 
#' @author Quentin Grimonprez
#' @param X matrix of size n*p
#' @param y vector of size n
#' @param hc output of \code{\link{hclust}} function. If not provided, \code{\link{hclust}} is run with ward.D2 method
#' @param control either "FDR" or "FWER"
#' @param alpha control elvel for testing procedure
#' @param test test used in the testing procedure. Default is partialFtest
#' @param plot If TRUE plot the number of groups selected before and after the testing procedure
#' @param fractionSampleMLGL a real between 0 and 1 : the fraction of individuals to use in the sample for MLGL (see Details).
#' @param ... Others parameters 
#'
#' @return a list containing :
#' \describe{
#'   \item{res}{output of \link{MLGL} function}
#'   \item{lambdaOpt}{lambda values maximizing the number of rejects}
#'   \item{var}{A vector containing the index of selected variables for \code{lambdaOpt}}
#'   \item{group}{A vector containing the values index of selected groups for \code{lambdaOpt}}
#' } 
#'
#' @details
#' Divide the n individuals in two samples. Then the three following steps are done :
#' 1) Hierarchical CLustering of the variables of X based on the first sample of individuals
#' 2) MLGL on the second sample of individuals
#' 3) Hierarchical testing procedure on the first sample of individuals.
#'
#' @examples
#' \dontrun{
#' set.seed(42)
#' X <- simuBlockGaussian(50,12,5,0.7)
#' y <- drop(X[,c(2,7,12)]%*%c(2,2,-2)+rnorm(50,0,0.5))
#' res <- fullProcess(X, y)
#'}
#'
#' @seealso \link{MLGL}, \link{hierarchicalFDR}, \link{hierarchicalFWER}, \link{selFDR}, \link{selFWER}
#'
#' @export
fullProcess <- function(X, y, control = c("FWER", "FDR"), alpha = 0.05, test = partialFtest, hc = NULL, plot = TRUE, fractionSampleMLGL = 1/2,...)
{
  control = match.arg(control)
  .checkFullProcess(X, y, hc, plot, alpha, test, fractionSampleMLGL)
    
  n <- nrow(X)
  
  # Split the data in 2
  ind2 <- sample(n, floor(n * fractionSampleMLGL))
  ind1 <- (1:n)[-ind2]

  ##### part 1 : hierarchical clustering with half of the data
  if(is.null(hc))
  {
    
    # center variables and sd = 1
    Xb <- scale(X, center = TRUE, scale = FALSE)
    Xb = scale(Xb, center = FALSE, scale = sqrt(colSums(X^2)/n))
    
    # enclidian distance
    d <- dist(t(Xb[ind1,]))
    
    # hierarchical clustering
    hc = fastcluster::hclust(d, method = "ward.D2")
      
  }
  
  
  ##### part 2 : group-lasso
  res <- MLGL(X[ind2,], y[ind2], hc = hc)
  
  
  ##### part 3 : testing procedure
  
  # choose the right function
  hierTestFunction <- hierarchicalFWER
  selFunction <- selFWER
  if(control == "FDR")
  {
    hierTestFunction = hierarchicalFDR
    selFunction = selFDR  
  }
    
  # testing procedure for each lambda
  REJECT <- list()
  nbReject <- rep(0, length(res$lambda))
  prevSelGroup = selGroup <- c()
  for(i in 1:length(res$lambda))
  {
    # if no groups are selected we do nothing
    if(length(res$group[[i]])>0)
    {
      
      selGroup = unique(res$group[[i]])
      
      # if the selected groups have not changed compared with the last iteration, we copy the result
      if(setequal(prevSelGroup, selGroup))
      {
        REJECT[[i]] = REJECT[[i-1]]
        nbReject[i] = nbReject[i-1]
      }
      else
      {
        # hierarchical testing and selection
        resTest <- hierTestFunction(X[ind1,], y[ind1], res$group[[i]], res$var[[i]], test)
        resSel <- selFunction(resTest, alpha, ...)
        
        # keep outerNode (need for FDR outer = FALSE, do not change in other cases)
        groupSel <- outerNode(resSel$toSel, resTest$hierMatrix)
        # Id of rejected groups
        REJECT[[i]] = (resSel$groupId[resSel$toSel])[groupSel] 
        
        # number of rejects for the lambda value
        nbReject[i] = length(REJECT[[i]])
      }

    }# end if no selection
    
    prevSelGroup = selGroup
    
  }# end for lambda
  
  # indice of optimal lambda : the one with the gretest number of reject
  indLambdaOpt <- which.max(nbReject)
  
  
  # plot the number of groups selevted by MLGL before and after the testing procedure
  if(plot)
  {
    matplot(res$lambda, cbind(res$nGroup, nbReject), type ="l", lwd = 1.5, xlab = expression(lambda), ylab = "Number of groups")
    # vertical line for optimal lambda
    abline(v = res$lambda[indLambdaOpt], lwd = 1, lty = "dotted")
    legend("topright", legend = c("groups",  "rejected groups"), lty = 1:2, lwd = 1.5, col = 1:2)
  }
  
  # group selected for the lambda optimal
  indGroupSel <- res$group[[indLambdaOpt]]%in% REJECT[[indLambdaOpt]]
  group <- res$group[[indLambdaOpt]][indGroupSel]
  var   <- res$var[[indLambdaOpt]][indGroupSel]
  
  
  return(list(res = res, lambdaOpt = res$lambda[indLambdaOpt], selectedGroups = REJECT[[indLambdaOpt]], group = group, var = var))
  
}


# check parameters of MLGL function
.checkFullProcess <- function(X, y, hc, plot, alpha, test, fractionSampleMLGL)
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
  
  #check if plot is a boolean
  if(length(plot)!=1)
    stop("plot must be a boolean.")
  if(!is.logical(plot))
    stop("plot must be a boolean.")

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
