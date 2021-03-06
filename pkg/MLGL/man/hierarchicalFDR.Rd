% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/hierarchicalFDR.R
\name{hierarchicalFDR}
\alias{hierarchicalFDR}
\title{Hierachical testing with FDR control}
\usage{
hierarchicalFDR(X, y, group, var, test = partialFtest)
}
\arguments{
\item{X}{original data}

\item{y}{associated response}

\item{group}{vector with index of groups. group[i] contains the index of the group of the variable var[i].}

\item{var}{vector whith the variables contained in each group. group[i] contains the index of the group of the variable var[i].}

\item{test}{function for testing the nullity of a group of coefficients in linear regression. 3 parameters : X : design matrix, y response and varToTest : vector of variables to test; return a pvalue}
}
\value{
a list containing :
\describe{
  \item{pvalues}{pvalues of the different test (without correction)}
  \item{adjPvalues}{adjusted pvalues}
  \item{groupId}{Index of the group}
  \item{hierMatrix}{Matrix describing the hierarchical tree.}
  }
}
\description{
Apply hierarchical test for each hierarchy, and test external variables for FDR control at level alpha
}
\details{
Version of the hierarchical testing procedure of Yekutieli for MLGL output. You can use th \link{selFDR} function to select groups
at a desired level alpha.
}
\examples{
set.seed(42)
X = simuBlockGaussian(50,12,5,0.7)
y = X[,c(2,7,12)]\%*\%c(2,2,-2) + rnorm(50,0,0.5)
res = MLGL(X,y)
test = hierarchicalFDR(X, y, res$group[[20]], res$var[[20]])


}
\references{
Yekutieli, Daniel. "Hierarchical False Discovery Rate-Controlling Methodology." Journal of the American Statistical Association 103.481 (2008): 309-16.
}
\seealso{
\link{selFDR}, \link{hierarchicalFWER}
}
