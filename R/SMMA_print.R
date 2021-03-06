#
#     Description of this R script:
#     R interface for SMMA routines.
#
#     Intended for use with R.
#     Copyright (C) 2017 Adam Lund
# 
#     This program is free software: you can redistribute it and/or modify
#     it under the terms of the GNU General Public License as published by
#     the Free Software Foundation, either version 3 of the License, or
#     (at your option) any later version.
# 
#     This program is distributed in the hope that it will be useful,
#     but WITHOUT ANY WARRANTY; without even the implied warranty of
#     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#     GNU General Public License for more details.
# 
#     You should have received a copy of the GNU General Public License
#     along with this program.  If not, see <http://www.gnu.org/licenses/>
#

#' @title Print Function for objects of Class SMMA
#'
#' @description This function will print some information about the SMMA object.
#'  
#' @param x a SMMA object
#' @param ... ignored
#' 
#' @examples  
#' 
#' ##size of example 
#' n1 <- 65; n2 <- 26; n3 <- 13; p1 <- 13; p2 <- 5; p3 <- 4
#' 
#' ##marginal design matrices (Kronecker components)
#' X1 <- matrix(rnorm(n1 * p1, 0, 0.5), n1, p1) 
#' X2 <- matrix(rnorm(n2 * p2, 0, 0.5), n2, p2) 
#' X3 <- matrix(rnorm(n3 * p3, 0, 0.5), n3, p3) 
#' X <- list(X1, X2, X3)
#' 
#' component <- rbinom(p1 * p2 * p3, 1, 0.1) 
#' Beta1 <- array(rnorm(p1 * p2 * p3, 0, .1) + component, c(p1 , p2, p3))
#' Beta2 <- array(rnorm(p1 * p2 * p3, 0, .1) + component, c(p1 , p2, p3))
#' mu1 <- RH(X3, RH(X2, RH(X1, Beta1)))
#' mu2 <- RH(X3, RH(X2, RH(X1, Beta2)))
#' Y1 <- array(rnorm(n1 * n2 * n3, mu1), dim = c(n1, n2, n3))
#' Y2 <- array(rnorm(n1 * n2 * n3, mu2), dim = c(n1, n2, n3))
#' 
#' Y <- array(NA, c(dim(Y1), 2))
#' Y[,,, 1] <- Y1; Y[,,, 2] <- Y2;
#' 
#' fit <- softmaximin(X, Y, zeta = 10, penalty = "lasso", alg = "npg")
#' fit
#' 
#' 
#' @method print SMMA
# @S3method print SMMA
#' @author Adam Lund
#' @export
 
print.SMMA <- function(x, ...) {
out <- data.frame(Df = x$df, lambda = x$lambda)
print(out)	
}
