#
# compute the principal component of each group and perform a OLS with the principal component as variables
#
# X : design matrix
# y : response
# group & var : vectors describing the group structure. group[i] contains the number of the group associated to the variable var[i] 
#
# return a list containing :
# - lm the summary of OLS
# - newdata matrix containing the principal components
# - y response (same as parameter)
# - group name of the columns of newdata
#
acpOLStest <- function(X, y, group, var)
{
  
  #acp + ols
  newdata = do.call(cbind, tapply(var, group, FUN = function(indvar)
  {
    respca=PCA(X[,indvar], scale.unit=TRUE, ncp=1)
    return(respca$ind$coord)
  }))
  #colnames(newdata) = unique(group)
  colnames(newdata) = unique(group)[order(unique(tapply(var, group, NULL)))]
  
  
  #ols on new data
  reslm = lm(y ~ newdata)
  sumlm = summary(reslm)
  
  rownames(sumlm$coefficients)[-1] = colnames(newdata)
  
  #return pval
  return(list(lm = sumlm, newdata = newdata, y = y, group = as.numeric(colnames(newdata))))
}


#'
#' Perform a partial F-test
#'
#' @title Partial F-test
#'
#' @param X design matrix of size n*p
#' @param y response vector of length n 
#' @param varToTest vector containing the index of the column of X to test
#'
#' @return a vector of the same length as varToTest containing the p-values of the test.
#' 
#' @details 
#' y = X * beta + epsilon
#' 
#' null hypothesis : beta[varToTest] = 0
#' alternative hypothesis : it exists an index k in varToTest such that beta[k] != 0
#' 
#' The test statistic is based on a full and a reduced model.
#' full : y = X * beta + epsilon
#' reduced : y = X * beta[-varToTest] + epsilon
#' 
#' @seealso \link{partialFtest}
#' 
#' @export
partialFtest <- function(X, y, varToTest)
{
  reduced=c()
  if(length(varToTest)==ncol(X))
    reduced = lm(y ~ 1)
  else
    reduced = lm(y ~ X[,-varToTest])
  
  full = lm(y ~ X)
  
  outPartialFtest = anova(reduced, full)
  
  return(outPartialFtest[['Pr(>F)']][2])
}


#'
#' Perform a F-test
#'
#' @title F-test
#'
#' @param X design matrix of size n*p
#' @param y response vector of length n 
#' @param varToTest vector containing the index of the column of X to test
#'
#' @return a vector of the same length as varToTest containing the p-values of the test.
#' 
#' @details 
#' y = X * beta + epsilon
#' 
#' null hypothesis : beta[varToTest] = 0
#' alternative hypothesis : it exists an index k in varToTest such that beta[k] != 0
#' 
#' The test statistic is based on a full and a reduced model.
#' full : y = X * beta[varToTest]  + epsilon
#' reduced : the null model
#' 
#' @seealso \link{partialFtest}
#' 
#' @export
Ftest <- function(X, y, varToTest)
{
  anova(lm(y~X[,varToTest]))[["Pr(>F)"]][1]
}  

# #'
# #' Perform a score test
# #'
# #' @title Score test
# #'
# #' @param X design matrix of size n*p
# #' @param y response vector of length n 
# #' @param varToTest vector containing the index of the column of X to test
# #'
# #' @return a vector of the same length as varToTest containing the p-values of the test.
# #' 
# #' @details 
# #' y = X * beta + epsilon
# #' 
# #' null hypothesis : beta[varToTest] = 0
# #' alternative hypothesis : it exists an index k in varToTest such that beta[k] != 0
# #' 
# #' see the reference for more details.
# #' 
# #' @seealso \link{partialFtest}, \link{Ftest}
# #' 
# #' @references Goeman, J. J., Van De Geer, S. A. and Van Houwelingen, H. C. (2006), Testing against a high dimensional alternative. Journal of the Royal Statistical Society: Series B (Statistical Methodology), 68: 477-493. doi:10.1111/j.1467-9868.2006.00551.x
# #' 
# #' @export
# scoreTest <- function(X, y, varToTest)
# {
#   globaltest::p.value(gt(y, X[,varToTest], model = "linear"))
# }
