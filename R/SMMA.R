#
#     Description of this R script:
#     R interface/wrapper for the Rcpp function pga in the SMMA package.
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

#' @name SMMA
#' @aliases pga
#' @title Maximin Estimation for Large Scale Array Data with Known Groups
#' 
#' @description  Efficient design matrix free procedure for solving a soft maximin problem for  
#' large scale array-tensor structured models, see  \cite{Lund et al., 2017}. 
#' Currently Lasso and SCAD penalized estimation is implemented.
#'  
#' @usage  softmaximin(X, 
#'             Y, 
#'             penalty = c("lasso", "scad"),
#'             nlambda = 30,
#'             lambda.min.ratio = 1e-04,
#'             lambda = NULL,
#'             penalty.factor = NULL,
#'             reltol = 1e-05,
#'             maxiter = 15000,
#'             steps = 1,
#'             btmax = 100,
#'             zeta = 2, 
#'             c = 0.001,
#'             Delta0 = 1,
#'             nu = 1,
#'             alg = c("npg", "mfista"),
#'             log = TRUE)
#'            
#' @param X A list containing the Kronecker components (2 or 3) of the Kronecker design matrix.
#'  These are  matrices of sizes \eqn{n_i \times p_i}.
#' @param Y The response values, an array of size \eqn{n_1 \times\cdots\times n_d \times G}. 
#' @param penalty A string specifying the penalty. Possible values 
#' are \code{"lasso", "scad"}.
#' @param nlambda The number of \code{lambda} values.
#' @param lambda.min.ratio The smallest value for \code{lambda}, given as a fraction of 
#' \eqn{\lambda_{max}}; the (data dependent) smallest value for which all coefficients are zero.
#' @param lambda The sequence of penalty parameters for the regularization path.
#' @param penalty.factor An array of size \eqn{p_1 \times \cdots \times p_d}. Is multiplied 
#' with each element in \code{lambda} to allow differential shrinkage on the coefficients. 
#' @param reltol The convergence tolerance for the inner loop.
#' @param maxiter The maximum number of  iterations allowed for each \code{lambda}
#' value, when  summing over all outer iterations for said \code{lambda}.
#' @param steps The number of steps used in the multi-step adaptive lasso algorithm for non-convex penalties. 
#' Automatically set to 1 when \code{penalty = "lasso"}.
#' @param btmax Maximum number of backtracking steps allowed in each iteration. Default is \code{btmax = 100}.
#' @param zeta Constant controlling  the softmax apprximation accuracy. Must be strictly positive. Default is \code{zeta = 2}.
#' @param c constant used in the NPG algorithm. Must be strictly positive. Default is \code{c = 0.001}.
#' @param Delta0 constant used to bound the stepsize. Must be strictly positive. Default is \code{Delta0 = 1}.
#' @param nu constant used to control the stepsize. Must be positive. A small value gives a big stepsize. Default is \code{nu = 1}.
#' @param alg string indicating which algortihm to use. Possible values are \code{"npg", "mfista"}.
#' @param log logical variable indicating wheter to use log-loss to or not.  TRUE is default and yields the problem described below.
#'
#' @details We consider the mixed model setup from \cite{Meinshausen and B{u}hlmann, 2015} for array data with a known fixed group structure 
#' and tensor structured design matrix. 
#' 
#' For  \eqn{g \in \{1,\ldots,G\}} let \eqn{n} be the number of observations in each group. 
#' With  \eqn{Y_g:=(y_{i},\ldots,y_{i_{n}})^\top}  the   group-specific \eqn{n_1\times \cdots \times n_d} response array  
#' and \eqn{X :=(x_{i}\mid\ldots\mid x_{i_{n}})^\top} a \eqn{n\times p} design matrix, with tensor structure
#'  \deqn{X = \bigotimes_{i=1}^d X_i,} 
#'  where for \eqn{d =2,3}, \eqn{X_1,\ldots, X_d} are the marginal \eqn{n_i\times p_i} design matrices (Kronecker components). 
#'  
#' We use the array model framework, see  \cite{Currie et al., 2006}, to write the model equation as
#'  \deqn{Y_g = \rho(X_d,\rho(X_{d-1},\ldots,\rho(X_1,B_g))) + E,}
#' where \eqn{\rho} is the so called rotated \eqn{H}-transfrom, \eqn{B_g} for each \eqn{g} is a random \eqn{p_1\times\cdots\times p_d} parameter array
#'  and \eqn{E}  is  \eqn{n_1\times \cdots \times n_d} error array uncorrelated with \eqn{X}.
#'  
#' In \cite{Meinshausen and B{u}hlmann, 2015} it is suggested to maximize the minimal  empirical explained variance 
#' \deqn{\hat V_g(\beta):=\frac{1}{n}(2\beta^\top X^\top y_g-\beta^\top X^\top X\beta),}
#' where \eqn{y_g:=vec(Y_g)}. In \cite{Lund et al., 2017} a soft version of this problem, 
#' the soft maximin problem, given as
#' \deqn{\min_{\beta}\log\bigg(\sum_{g=1}^G \exp(-\zeta \hat V_g(\beta))\bigg) + \lambda J (\beta),} 
#'  is suggested, for \eqn{J} a proper and convex function and \eqn{\zeta > 0}. 
#'  
#'  For \eqn{d=2,3} and using only the marginal matrices \eqn{X_1,X_2,\ldots}, the function \code{softmaximin} 
#' solves the soft maximin problem for a sequence of penalty parameters \eqn{\lambda_{max}>\ldots >\lambda_{min}>0}. The underlying algorithm is based on a non-monotone  
#' proximal gradient method. We note that if \eqn{J} is not  convex, as with the SCAD penalty,
#' we use the multiple step adaptive lasso procedure to loop over the proximal algorithm, see \cite{Lund et al., 2017} for more details.
#'   
#' @return An object with S3 Class "SMMA". 
#' \item{spec}{A string indicating the array dimension (2 or 3) and the penalty.}  
#' \item{coef}{A \eqn{p_1\cdots p_d \times} \code{nlambda} matrix containing the estimates of 
#' the model coefficients (\code{beta}) for each \code{lambda}-value.}
#' \item{lambda}{A vector containing the sequence of penalty values used in the estimation procedure.}
#' \item{Obj}{A matrix containing the objective values for each iteration and each model.}
#' \item{df}{The number of nonzero coefficients for each value of \code{lambda}.}	
#' \item{dimcoef}{A vector giving the dimension of the model coefficient array \eqn{\beta}.}
#' \item{dimobs}{A vector giving the dimension of the observation (response) array \code{Y}.}
#' \item{Iter}{A list with 4 items:  
#' \code{bt_iter}  is total number of backtracking steps performed,
#' \code{bt_enter} is the number of times the backtracking is initiated,
#' and \code{iter_mat} is a vector containing the  number of  iterations for each \code{lambda} value 
#' and  \code{iter} is total number of iterations i.e. \code{sum(Iter)}.}  
#'  
#' @author  Adam Lund
#' 
#' Maintainer: Adam Lund, \email{adam.lund@@math.ku.dk}
#' 
#' @references 
#' Lund, A., S. W. Mogensen and N. R. Hansen  (2017). Estimating Soft Maximin Effects in Heterogeneous Large-scale Array Data.
#' \emph{Preprint}. 
#' 
#' Meinshausen, N and P. B{u}hlmann (2015). Maximin effects in inhomogeneous large-scale data.
#' \emph{The Annals of Statistics}. 43, 4, 1801-1830.
#' 
#' Currie, I. D., M. Durban, and P. H. C. Eilers (2006). Generalized linear
#' array models with applications to multidimensional
#' smoothing. \emph{Journal of the Royal Statistical Society. Series B}. 68, 259-280.
#' 
#' @keywords package 
#'
#' @examples 
#' 
#' ##size of example 
#' n1 <- 65; n2 <- 26; n3 <- 13; p1 <- 13; p2 <- 5; p3 <- 4
#' 
#' ##marginal design matrices (Kronecker components)
#' X1 <- matrix(rnorm(n1 * p1), n1, p1) 
#' X2 <- matrix(rnorm(n2 * p2), n2, p2) 
#' X3 <- matrix(rnorm(n3 * p3), n3, p3) 
#' X <- list(X1, X2, X3)
#' 
#' component <- rbinom(p1 * p2 * p3, 1, 0.1) 
#' Beta1 <- array(rnorm(p1 * p2 * p3, 0, 0.1) + component, c(p1 , p2, p3))
#' mu1 <- RH(X3, RH(X2, RH(X1, Beta1)))
#' Y1 <- array(rnorm(n1 * n2 * n3), dim = c(n1, n2, n3)) + mu1
#' Beta2 <- array(rnorm(p1 * p2 * p3, 0, 0.1) + component, c(p1 , p2, p3))
#' mu2 <- RH(X3, RH(X2, RH(X1, Beta2)))
#' Y2 <- array(rnorm(n1 * n2 * n3), dim = c(n1, n2, n3)) + mu2
#' Beta3 <- array(rnorm(p1 * p2 * p3, 0, 0.1) + component, c(p1 , p2, p3))
#' mu3 <- RH(X3, RH(X2, RH(X1, Beta3)))
#' Y3 <- array(rnorm(n1 * n2 * n3), dim = c(n1, n2, n3)) + mu3
#' Beta4 <- array(rnorm(p1 * p2 * p3, 0, 0.1) + component, c(p1 , p2, p3))
#' mu4 <- RH(X3, RH(X2, RH(X1, Beta4)))
#' Y4 <- array(rnorm(n1 * n2 * n3), dim = c(n1, n2, n3)) + mu4
#' Beta5 <- array(rnorm(p1 * p2 * p3, 0, 0.1) + component, c(p1 , p2, p3))
#' mu5 <- RH(X3, RH(X2, RH(X1, Beta5)))
#' Y5 <- array(rnorm(n1 * n2 * n3), dim = c(n1, n2, n3)) + mu5
#' 
#' Y <- array(NA, c(dim(Y1), 5))
#' Y[,,, 1] <- Y1; Y[,,, 2] <- Y2; Y[,,, 3] <- Y3; Y[,,, 4] <- Y4; Y[,,, 5] <- Y5;
#'
#' fit <- softmaximin(X, Y, penalty = "lasso", alg = "npg")
#' Betafit <- fit$coef
#' 
#' modelno <- 15
#' m <- min(Betafit[ , modelno], c(component))
#' M <- max(Betafit[ , modelno], c(component))
#' plot(c(component), type="l", ylim = c(m, M))
#' lines(Betafit[ , modelno], col = "red")
#' 

#' @export
#' @useDynLib SMMA, .registration = TRUE
#' @importFrom Rcpp evalCpp
softmaximin <-function(X,
               Y, 
               penalty = c("lasso", "scad"),
               nlambda = 30,
               lambda.min.ratio = 1e-04,
               lambda = NULL,
               penalty.factor = NULL,
               reltol = 1e-05,
               maxiter = 15000,
               steps = 1,
               btmax = 100,
               zeta = 2, 
               c = 0.001,
               Delta0 = 1,
               nu = 1,
               alg = c("npg", "mfista"),
               log = TRUE) {
  
##get dimensions of problem
dimglam <- length(X)

if (dimglam < 2 || dimglam > 3){
  
stop(paste("the dimension of the model must be 2 or 3!"))
  
}else if (dimglam == 2){X[[3]] <- matrix(1, 1, 1)} 

X1 <- X[[1]]
X2 <- X[[2]]
X3 <- X[[3]]

dimX <- rbind(dim(X1), dim(X2), dim(X3))

n1 <- dimX[1, 1]
n2 <- dimX[2, 1]
n3 <- dimX[3, 1]
p1 <- dimX[1, 2]
p2 <- dimX[2, 2]
p3 <- dimX[3, 2]
n <- prod(dimX[,1])
p <- prod(dimX[,2])
G <- dim(Y)[length(dim(Y))]

####reshape Y into matrix
Z <- array(NA, c(n1, n2 * n3, G))
for(i in 1:G){Z[, , i] <- matrix(Y[, , , i], n1, n2 * n3)}

####check on algorithm
if(sum(alg == c("npg", "mfista")) != 1){stop(paste("algorithm must be correctly specified"))}

if(alg == "npg"){npg <- 1}else{npg <- 0}

####check on loss type
if(log == TRUE){ll <- 1}else{ll <- 0}


####check on constants
if(c <= 0){stop(paste("c must be strictly positive"))}

if(Delta0 <= 0){stop(paste("Delta0 must be strictly positive"))}

if(zeta <= 0){stop(paste("zeta must be strictly positive"))}

####check on penalty
if(sum(penalty == c("lasso", "scad")) != 1){stop(paste("penalty must be correctly specified"))}

if(penalty == "lasso"){steps <- 1}

# ####check on weights 
# if(is.null(Weights)){
#   
# Weights <- Z * 0 + 1
# 
# }else{
# 
# if(min(Weights) < 0){stop(paste("only positive weights allowed"))}    
# 
# W <- array(NA, c(n1, n2 * n3, G))
# for(i in 1:G){W[, , i] <- matrix(Weights[, , , i], n1, n2 * n3)}
# 
# }

####check on lambda 
if(is.null(lambda)){
  
makelamb <- 1
lambda <- rep(NA, nlambda)
  
}else{
  
if(length(lambda) != nlambda){

stop(
paste("number of elements in lambda (", length(lambda),") is not equal to nlambda (", nlambda,")", sep = "")
)

}
  
makelamb <- 0

}  

####check on penalty.factor
if(is.null(penalty.factor)){
  
penalty.factor <- matrix(1, p1, p2 * p3)
  
}else if(prod(dim(penalty.factor)) != p){
  
stop(
paste("number of elements in penalty.factor (", length(penalty.factor),") is not equal to the number of coefficients (", p,")", sep = "")
)
  
}else {
  
if(min(penalty.factor) < 0){stop(paste("penalty.factor must be positive"))}    
  
penalty.factor <- matrix(penalty.factor, p1, p2 * p3)
  
}  

##run  algorithm
res <- pga(X1, X2, X3,
           Z, 
           penalty,
           zeta, #approx accu
           c, #bactrack  par
           lambda, nlambda, makelamb, lambda.min.ratio,
           penalty.factor,
           reltol, 
           maxiter,
           steps,
           btmax,
           Delta0,
           nu,
           npg,
           ll)

####checks
if(res$Stops[2] == 1){

warning(paste("program exit due to maximum number of iterations (",maxiter,") reached for model no ",res$endmodelno,""))

}

if(res$Stop[3] == 1){

warning(paste("program exit due to maximum number of backtraking steps reached for model no ",res$endmodelno,""))

}

endmodelno <- res$endmodelno
Iter <- res$Iter

maxiterpossible <- sum(Iter > 0)
maxiterreached <- sum(Iter >= (maxiter - 1)) 

if(maxiterreached > 0){

warning(
paste("maximum number of inner iterations (",maxiter,") reached ",maxiterreached," time(s) out of ",maxiterpossible," possible")
)

}

out <- list()

class(out) <- "SMMA"

out$spec <- paste("", dimglam,"-dimensional ", penalty," penalized model") 
out$coef <- res$Beta[ , 1:endmodelno]
out$lambda <- res$lambda[1:endmodelno] 
out$df <- res$df[1:endmodelno]
out$dimcoef <- c(p1, p2, p3)[1:dimglam]
out$dimobs <- c(n1, n2, n3)[1:dimglam]
out$Obj <- res$Obj

Iter <- list()
Iter$bt_enter <- res$btenter
Iter$bt_iter <- res$btiter
Iter$iter_mat <- res$Iter[1:endmodelno, ]
Iter$iter <- sum(Iter$iter_mat, na.rm = TRUE)

out$Iter <- Iter

return(out)

}