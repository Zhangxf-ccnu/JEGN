require("igraph")
require("MASS")
require("Matrix")

generate_data <-
  function(p = 100, n = 100, tau = 0.5, K = 3, G = 4, model = "ER", umin = 0.5, umax = 1){

    if (model == "ER"){
      Net.Comm = as_adjacency_matrix(sample_gnp(p, 0.02), type = "both", sparse=FALSE)
    }

    if (model == "SF"){
      Net.Comm = as_adjacency_matrix(sample_pa(p,directed = FALSE), type = "both", sparse=FALSE)
    }

    NO.Comm = sum(Net.Comm)/2
    Net.Spec = list()
    for (g in 1:G){
      Net.Spec[[g]] = as_adjacency_matrix(sample_gnm(p, floor(NO.Comm*tau)), type = "both", sparse=FALSE)
    }

    M = matrix(list(), K, G)
    Omega = matrix(list(), K, G)
    R = matrix(list(), K, 1)


    for (k in 1:K){

      temp_comm = matrix(runif(p*p, min = umin, max = umax)*(2*rbinom(p*p, 1, 0.5) - 1), p, p)
      R[[k]] = Net.Comm*temp_comm
      R[[k]] = R[[k]]*upper.tri(R[[k]])
      R[[k]] = R[[k]] + t(R[[k]])

      eigen.min = c()
      for (g in 1:G){
        temp_spec = matrix(runif(p*p, min = umin, max = umax)*(2*rbinom(p*p, 1, 0.5) - 1), p, p)
        M[[k,g]] = Net.Spec[[g]]*temp_spec
        M[[k,g]] = M[[k,g]]*upper.tri(M[[k,g]])
        M[[k,g]] = M[[k,g]] + t(M[[k,g]])
        Omega[[k,g]] = R[[k]] + M[[k,g]]
        eigen.min[g] = eigen(Omega[[k,g]])$values[p]
      }

      eigen_min = min(eigen.min)

      for (g in 1:G){
        Omega[[k,g]] = R[[k]] + M[[k,g]] + (abs(eigen_min) + 0.1)*diag(p)
      }
    }


    X = matrix(list(), K, G)
    Sigma = matrix(list(), K, G)
    for (k in 1:K){
      for (g in 1:G){
        X[[k,g]] = mvrnorm(n, rep(0,p), solve(Omega[[k,g]]))

         Sigma[[k,g]] = cov(X[[k,g]])*(n-1)/n
    #    Sigma[[k,g]] =   (1/n)*t(X[[k,g]])%*%X[[k,g]]
     #   Sigma[[k,g]] = as.matrix(nearPD((1/n)*t(X[[k,g]])%*%X[[k,g]])$mat)

      }
    }

    result = list(Sigma = Sigma, Omega = Omega, R = R, M = M, X = X)
  }
