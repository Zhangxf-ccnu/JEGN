#' Joint gene network inference across multiple subpopulations and data types
#'
#' The complete procedure for reconstructin gene networks from gene expression data using JEGN (Algorithm 2
#' in the Supplementary Information). For details, refer to Supplementary Section S3.4.
#' @param X A matrix (\eqn{K \times G}) of list of data matrices (\eqn{n_g \times p}), where K is
#' the number of data types and G is the number of subpopulations. The (k,g)-th element is a
#' \eqn{n_g \times p} data matrix corresponding to the k-th data type and the g-th subpopulation.
#' @param lambda The tuning parameter for controlling the level of sparsity of the estimated networks.
#' @param alpha The tuning parameter for controlling the extend of similarity among the estimated
#' subpopulation-specific networks.
#' @param model A character string indicating which model is used to fit the data. "Gaussian" (default)
#' and "nonparanormal" can be used. If model == "Gaussian", the Gaussian graphical model will be used.
#' If model == "nonparanormal", the nonparanormal graphical model will be used.
#' @param weights Determines the putative sample size of each subpopulation's data.  Allowed values:
#' a vector with length equal to the number of subpopulations; "equal", giving each subpopulation
#'  weight 1; "sample.size", giving each subpopulation weight corresponding to its sample size.
#' @param penalize.diagonal Determines whether the sparsity penalty is applied to the diagonal.
#'
#' @return a list with the following components
#' \item{Omega.hat}{A matrix (\eqn{K \times G}) of list of estimated precision matrices.}
#' \item{R.hat}{A matrix (\eqn{K \times 1}) of list of estimated common components.}
#' \item{M.hat}{A matrix (\eqn{K \times G}) of list of estimated subpopulation-unique components.}
#' \item{Omega.bar}{A matrix (\eqn{1 \times G}) of list of estimated gene networks.}
#' \item{R.bar}{A matrix of estimated common subnetworks.}
#' \item{M.bar}{A matrix (\eqn{1 \times G}) of list of estimated subpopulation-unique subnetworks.}
#' @export
#' @importFrom stats cor cov
#' @importFrom Matrix nearPD
#' @author Xiao-Fei Zhang  <zhangxf@mail.ccnu.edu.cn>
#'
#' @references Xiao-Fei Zhang, Le Ou-Yang, Ting Yan, Xiaohua Hu and Hong Yan (2019),
#' A joint graphical model for inferring gene networks across multiple subpopulations and data types,
#' @examples
# TCGA breast cancer data
#'data("TCGA.BRCA")
#'result = JEGN(TCGA.BRCA$X, 0.95, 0.4, model = "nonparanormal", weights = "equal")

JEGN = function(X, lambda, alpha, model = "Gaussian", weights = "equal", penalize.diagonal = FALSE){

  K = dim(X)[1]
  G = dim(X)[2]
  p = dim(X[[1,1]])[2]
  n = rep(0, G)
  for (g in 1:G){
    n[g] =  dim(X[[1,g]])[1]
  }

  S = matrix(list(), K, G)

  if (model == "Gaussian"){
    for (k in 1:K){
      for (g in 1:G){
        S[[k,g]] = Corr(X[[k,g]], method = "pearson")
      }
    }
  }

  if (model == "nonparanormal"){
    for (k in 1:K){
      for (g in 1:G){
        S[[k,g]] = Corr(X[[k,g]], method = "kendall")
      }
    }
  }

  if (weights=="equal"){
    n = rep(1,G)
  }else if (weights=="sample.size"){
    n = n
  }else{
    n = weights
  }


  result = JEGN.admm(S, lambda, alpha, n = n, penalize.diagonal = penalize.diagonal)

  result$Omega.bar = matrix(list(), 1, G)
  result$M.bar = matrix(list(), 1, G)
  result$R.bar = result$R.hat[[1]]!=0
  diag(result$R.bar) <- 0
  for (g in 1:G){
    result$Omega.bar[[g]] =  result$Omega.hat[[1,g]]!=0
    diag(result$Omega.bar[[g]]) <- 0
    result$M.bar[[g]] =  (result$Omega.hat[[1,g]]!=0)&(result$R.bar==0)
    diag(result$M.bar[[g]]) <- 0

  }

  result
}




# Algorithm 1



#' ADMM algorithm for JEGN
#'
#' ADMM algorithm for JEGN (Algorithm 1 in the Supplementary Information).
#' For details, refer to Supplementary Section S3.3.
#' @param S A matrix (\eqn{K \times G}) of list of sample covariance matrices (\eqn{p \times p}),
#'  where K is the number of data types and G is the number of subpopulations. The (k,g)-th
#'  element is a \eqn{p \times p} sample covariance matrix corresponding to the k-th data type
#'  and the g-th subpopulation.
#' @param lambda The tuning parameter for controlling the level of sparsity of networks.
#' @param alpha The tuning parameter for controlling the extend of similarity among subpopulation-specific networks.
#' @param n The sample size. A vector with length equal to the number of subpopulations.
#' @param penalize.diagonal Determines whether the sparsity penalty is applied to the diagonal.
#' @param epsilon The tolerance parameter for convergence criteria.
#' @param maxiter The maximum number of iterations for the ADMM algorithm.
#' @param rho The penalty parameter in the ADMM algorithm.
#' @param rho.incr The increase step parameter for varying penalty parameter rho.
#' @param rho.max The maximum value of rho.
#'
#' @return a list with the following components
#' \item{Omega.hat}{A matrix (\eqn{K \times G}) of list of estimated precision matrices.}
#' \item{R.hat}{A matrix (\eqn{K \times 1}) of list of estimated common components.}
#' \item{M.hat}{A matrix (\eqn{K \times G}) of list of estimated subpopulation-unique components.}
#'
#' @export
#'
#' @author Xiao-Fei Zhang  <zhangxf@mail.ccnu.edu.cn>
#'
#' @references Xiao-Fei Zhang, Le Ou-Yang, Ting Yan, Xiaohua Hu and Hong Yan (2019),
#' A joint graphical model for inferring gene networks across multiple subpopulations and data types,
#' @examples
# TCGA breast cancer data
#'data("TCGA.BRCA")
#'result = JEGN(TCGA.BRCA$Sigma, 0.95, 0.4, model = "nonparanormal", weights = "equal")


JEGN.admm = function(S, lambda, alpha, n = NULL, penalize.diagonal = FALSE, epsilon = 1e-5,
                maxiter = 500, rho = 0.1,  rho.incr = 1.2, rho.max = 1e10){
    #
    K = dim(S)[1]
    G = dim(S)[2]
    p = dim(S[[1,1]])[1]


    if(is.null(n)){
      n = rep(1, G)
    }

    # initialize:
    Omega = matrix(list(), K, G)
    R = matrix(list(), K, 1)
    M = matrix(list(), K, G)
    Q = matrix(list(), K, G)

    for (k in 1:K){
      R[[k]] = diag(p)
      for (g in 1:G){
        M[[k,g]] = diag(p)
        Omega[[k,g]] = diag(p)
        Q[[k,g]] = matrix(0, p, p)
      }
    }


    for (i in 1:maxiter){

      # updata Omega according to EQ (12)
      Omega_prev = Omega
      for (k in 1:K){
        for (g in 1:G){
          Omega[[k,g]] = expand(R[[k]] + M[[k,g]] - (1/rho)*(n[g]*S[[k,g]] + Q[[k,g]]), rho, n[g])
        }
      }

      # updata R according to EQ (16)
      A = matrix(list(), K, 1)
      for (k in 1:K){
        A[[k]] = matrix(0, p, p)
        for (g in 1:G){
          A[[k]] = A[[k]] + (Omega[[k,g]] - M[[k,g]] + Q[[k,g]]/rho)
        }
        A[[k]] = A[[k]]/G
      }
      R = group_lasso(A, lambda*alpha/rho, penalize.diagonal)


      # updata M according to EQ (19)
      for (g in 1:G){
        B = matrix(list(), K, 1)
        for (k in 1:K){
          B[[k]] = Omega[[k,g]] - R[[k]] + Q[[k,g]]/rho
        }
        M[,g] = group_lasso(B, lambda*(1-alpha)/rho, penalize.diagonal)
      }

      # updata dual variables Q
      for (k in 1:K){
        for (g in 1:G){
          Q[[k,g]] =  Q[[k,g]] + rho*(Omega[[k,g]] - (R[[k]] + M[[k,g]]))
        }
      }

      # Check the convergence condition;
      diff_value_1 = 0
      diff_value_2 = 0
      norm_value = 0
      for(k in 1:K){
        for (g in 1:G){
          diff_value_1 = diff_value_1 + sum(abs(Omega[[k,g]] - Omega_prev[[k,g]]))
          diff_value_2 = diff_value_2 + sum(abs(Omega[[k,g]] - (R[[k]] + M[[k,g]])))
          norm_value  = norm_value + sum(abs(Omega[[k,g]]))
        }
      }
      if (max(diff_value_1, diff_value_2) <= norm_value*epsilon){
        break
      }

      rho = min(rho*rho.incr,rho.max)


    }

    for(k in 1:K){
      for(g in 1:G){
        Omega[[k,g]] = R[[k]] + M[[k,g]]
      }
    }

    out = list(Omega.hat = Omega, R.hat = R, M.hat = M)
  }






# Compute sample covariance matrices
Corr <- function(X, method = "pearson"){
  n = dim(X)[1]
  if(method=="pearson"){
    S = (n-1)*stats::cov(X, method = "pearson")/n
  }

  if (method == "kendall"){
    S = stats::cor(X, method = "kendall")
    S = sin(S*pi/2)
    S[is.na(S)] = 0
    diag(S) <- 1
    S = as.matrix(Matrix::nearPD(S, corr = TRUE)$mat)
  }
  S
}


# Expand operator, EQ (9)
expand = function(A, rho, n){
  edecomp <- eigen(A)
  D <- edecomp$values
  U <- edecomp$vectors
  D2 <- 0.5*(D + sqrt(D^2 + 4*n/rho ))
  Omega <- U %*% diag(D2) %*% t(U)
  Omega
}




# group lasso, EQ (13) and (14)
group_lasso =  function(A, lambda, penalize.diagonal=FALSE){

    K = length(A)
    normA = A[[1]]*0
    for(k in 1:K){
      normA = normA + (A[[k]])^2
    }
    normA = sqrt(normA)

    notshrunk = (normA>lambda)*1
    # reset 0 elements of normA to 1 so we don't get NAs later.
    normA = normA + (1-notshrunk)

    out = A
    for(k in 1:K)
    {
      out[[k]] = A[[k]]*(1-lambda/normA)
      out[[k]] = out[[k]]*notshrunk
      if(!penalize.diagonal){
        diag(out[[k]]) <- diag(A[[k]])
      }
    }
    out
  }




