#' @import gglasso MASS Matrix fastcluster FactoMineR
#' 
#' @title MLGL
#' @docType package
#' @aliases MLGL-package
#' @name MLGL-package
#' @description  
#' Group-Lasso with Hierarchical Clustering
#'
#' 
#' @details
#' 
#'   \tabular{ll}{
#' Package: \tab HCgglasso\cr
#' Type: \tab Package\cr
#' Version: \tab 0.3\cr
#' Date: \tab 2016-09-09\cr
#' License: \tab GPL (>=2) \cr
#' }
#' 
#' 
#' This package presents a method combining Hierarchical Clustering and Group-lasso. Usually, a single partition of the covariates is used in the group-lasso.
#' Here, we provides several partition from the hierarchical tree.
#' 
#' A post-treatment method based on statistical test (with FWER and FDR control) for selecting the regularization parameter and the optimal group for this value is provided.
#' This method can be applied for the classical group-lasso and our method.  
#' 
#' 
#' @author Quentin Grimonprez 
#' 
#' Maintainer: Quentin Grimonprez  <quentin.grimonprez@@inria.fr>
#'  
#' 
#' @examples 
#' # Simulate gaussian data with block-diagonal variance matrix containing 12 blocks of size 5
#' X <- simuBlockGaussian(50, 12, 5, 0.7)
#' # Generate a response variable
#' y <- drop(X[,c(2,7,12)]%*%c(2,2,-2)+rnorm(50,0,0.5))
#' # Apply MLGL method
#' res <- MLGL(X,y)
#' 
#' @seealso \link{MLGL}, \link{cv.MLGL}
#' 
#' @keywords package
NULL