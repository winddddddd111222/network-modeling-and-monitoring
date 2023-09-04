library(doParallel)
library(parallel)
library(LaplacesDemon)

## This script include all the functions for the VEM algorithm

series <- function(k, delta){
  return((delta[2] / (delta[2] + k)) ** delta[1] / k)
}
# This function is for the E-step corresponding to network topology of the algorithm
E.step.topology <- function(i, r, g, adjacent.tensor, adjacent.mat, beta1, tau1, rs.mat){
  community.num <- dim(beta1$nu)[2]
  node.num <- dim(adjacent.tensor)[1]
  value <- 0
  for (s in seq(community.num)){
    for (j in seq(node.num)){
      if (i != j){
        if(adjacent.mat[i, j] == 0){
          summationij <- 0
        }else{
          summationij <- sum(log(seq(adjacent.mat[i, j])))
        }
        if(adjacent.mat[j, i] == 0){
          summationji <- 0
        }else{
          summationji <- sum(log(seq(adjacent.mat[j, i])))
        }
        value <- value + tau1[j, s] * ((adjacent.mat[i, j] == 0)*(digamma(beta1$delta.pi[r, s, 1, g]) - digamma(sum(beta1$delta.pi[r, s, , g]))) + 
              (adjacent.mat[i, j] > 0) * ((digamma(beta1$delta.pi[r, s, 2, g]) - digamma(sum(beta1$delta.pi[r, s, , g]))) + adjacent.mat[i, j] *
              (digamma(beta1$delta.theta[r, s, 1, g]) - log(beta1$delta.theta[r, s, 2, g])) - beta1$delta.theta[r, s, 1, g] / beta1$delta.theta[r, s, 2, g] - 
              summationij + rs.mat[r, s] + sum(adjacent.tensor[i, j,] * (digamma(beta1$delta.variphi[r, s, , g]) - digamma(sum(beta1$delta.variphi[r, s, , g]))))) + 
              (adjacent.mat[j,i] == 0) * (digamma(beta1$delta.pi[s, r, 1, g]) - digamma(sum(beta1$delta.pi[s, r, , g]))) + 
              (adjacent.mat[j,i] > 0) * ((digamma(beta1$delta.pi[s, r, 2, g]) - digamma(sum(beta1$delta.pi[s, r, , g]))) + adjacent.mat[j,i] *
              (digamma(beta1$delta.theta[s, r, 1, g]) - log(beta1$delta.theta[s, r, 2, g])) - beta1$delta.theta[s, r, 1, g] / beta1$delta.theta[s, r, 2, g] - 
               summationji + rs.mat[s, r] + sum(adjacent.tensor[j, i, ] * 
               (digamma(beta1$delta.variphi[s, r, , g]) - digamma(sum(beta1$delta.variphi[s, r, , g]))))))
      } 
    }
  }
  return(value)
}

# This function is for the E-step corresponding to the node attributes of the algorithm
E.step.attribute <- function(i, r, g, attribute.tensor, d.mu, weight){
  attribute.num <- dim(attribute.tensor)[2]
  value <- 0
  for (t in seq(attribute.num)) {
    value <- value + weight * attribute.tensor[i, t, g] * (digamma(d.mu[r, t]) - digamma(sum(d.mu[r , ])))
  }
  return(value)
}

E.step.remain <- function(r, g, nu){
  return(digamma(nu[g, r]) - digamma(sum(nu[g, ])))
}
# This function is for one iteration of E-step of the algorithm 
E.step.iter <- function(g, adjacent.tensor, attribute.tensor, beta1, tau1, rs.mat, weight){
  node.num <- dim(adjacent.tensor)[1]
  community.num <- dim(rs.mat)[1]
  adjacent.mat <- apply(adjacent.tensor[, , , g], c(1, 2), sum)
  log.tau <- matrix(NA, nrow = node.num, ncol = community.num)
  for (i in seq(node.num)) {
    for (r in seq(community.num)) {
      log.tau[i, r] <- E.step.topology(i, r, g, adjacent.tensor[, , , g], adjacent.mat, beta1, tau1, rs.mat) +
        E.step.attribute(i, r, g, attribute.tensor, beta1$delta.mu, weight) +
        E.step.remain(r, g, beta1$nu)
    }
  }
  return(log.tau)
}

rsmat.cal <- function(d.theta, K){
  community.num <- dim(d.theta)[1]
  rs.mat <- matrix(NA, nrow = community.num, ncol = community.num)
  for (r in seq(community.num)) {
    for (s in seq(community.num)) {
      rs.mat[r, s] <- sum(sapply(1:K, series, delta = d.theta[r, s, ]))
    }
  }
  return(rs.mat)
}
# This function is for the E-step of the algorithm 
E.step.beta1 <- function(adjacent.tensor, attribute.tensor, beta1, max.iter, tol = 1e-8, no.cluster, weight){
  network.num <- dim(adjacent.tensor)[4]
  par.fun <- function(g){
    node.num <- dim(adjacent.tensor)[1]
    difference <- 10
    iter.index <- 0
    tau1 <- beta1$tau[, , g]
    tau2 <- tau1
    rs.mat <- matrix(0, nrow = 3, ncol = 3)
    while((iter.index <= max.iter) & (difference > tol)){
      log.tau <- E.step.iter(g, adjacent.tensor, attribute.tensor, beta1, tau1, rs.mat, weight)
      for (i in seq(node.num)) {
        tau2[i, ] <- exp(log.tau[i, ] - max(log.tau[i, ])) / sum(exp(log.tau[i, ] - max(log.tau[i, ])))
      }
      difference <- max(abs(tau1 - tau2))
      tau1 <- tau2
      iter.index <- iter.index + 1
    }
    return(tau1)
  }
  cl <- makeCluster(no.cluster)
  clusterExport(cl, c('max.iter', 'tol', 'beta1', 'weight', 'adjacent.tensor', 'attribute.tensor'), envir = environment())
  clusterExport(cl, c('E.step.topology', 'E.step.attribute', 'E.step.remain', 'E.step.iter', 'series', 'rsmat.cal'))
  res <- parLapply(cl, 1:network.num, par.fun)
  stopCluster(cl)
  for (g in seq(network.num)) {
    beta1$tau[, , g] <- res[[g]]
  }
  beta1$tau[which(beta1$tau == 0)] <- 1e-10
  new.tau <- beta1$tau
  for (g in seq(network.num)){
    for (i in seq(node.num)){
      beta1$tau[i, , g] <- beta1$tau[i, , g] / sum(new.tau[i, , g])
    }
  }
  return(beta1$tau)
}
# This function is used for update nu in the M-step
M.step.nu <- function(nu, tau, alpha) {
  network.num <- dim(tau)[3]
  community.num <- dim(tau)[2]
  for (g in seq(network.num)) {
    for (r in seq(community.num)) {
      nu[g, r] <- sum(tau[, r, g]) + alpha[r]
    }
  }
  return(nu)
}

# This function is used for update delta^mu in the M-step
M.step.dmu <- function(d.mu, tau, eta, attribute.tensor, weight){
  community.num <- dim(tau)[2]
  attribute.num <- dim(d.mu)[2]
  attribute.dim <- c()
  for (r in seq(community.num)) {
    for (t in seq(attribute.num)) {
      d.mu[r, t] <- weight * sum(tau[ , r, ] * attribute.tensor[, t, ]) + eta[r, t]
    }
  }
  return(d.mu)
}

# This function is used for update delta^pi in the M-step
M.step.dpi <- function(d.pi, tau, M, adjacent.tensor){
  community.num <- dim(tau)[2]
  adjacent.mat <- apply(adjacent.tensor, c(1, 2), sum)
  node.num <- dim(adjacent.tensor)[1]
  d.pi[, , 1] <- M[, , 1]
  d.pi[, , 2] <- M[, , 2]
  d.pi[, , 1] <- M[, , 1] + t(tau) %*% ((adjacent.mat == 0) - diag(1, nrow = node.num, ncol = node.num)) %*% tau
  d.pi[, , 2] <- M[, , 2] + t(tau) %*% (adjacent.mat > 0) %*% tau
  return(d.pi)
}

# This function is used for update delta^varphi in the M-step
M.step.dvariphi <- function(d.variphi, tau, omega, adjacent.tensor) {
  community.num <- dim(tau)[2]
  node.num <- dim(adjacent.tensor)[1]
  layer.num <- dim(adjacent.tensor)[3]
  for (l in seq(layer.num)) {
    d.variphi[, , l] <- omega[, , l]
    d.variphi[, , l] <- omega[, , l] + t(tau) %*% adjacent.tensor[, , l] %*% tau
  }
  return(d.variphi)
}

# This function defined the part of ELBO related to delta^theta
obj.theta <- function(x, tau, N, adjacent.mat, logf.mat){
  K <- 500
  node.num <- dim(adjacent.mat)[1]
  temp <- 0
  value <- 0
  mid.mat <- (adjacent.mat > 0) *
    (adjacent.mat * (digamma(x[1]) -  log(x[2])) - x[1] / x[2] -
       logf.mat + temp)
  value <- tau[, 1] %*% mid.mat %*% tau[, 2]
  value <- value + (N[1]-1) * (digamma(x[1]) - log(x[2])) - N[2] * 
    x[1] / x[2] + N[1] * log(N[2]) - lgamma(N[1]) -
    ((x[1] - 1) * (digamma(x[1]) - log(x[2])) - x[1] + x[1] * log(x[2]) - lgamma(x[1])) 
  return(- value)
}

# This function defined the gradient of the part of ELBO related to delta^theta
gra.theta <- function(x, tau, N, adjacent.mat, logf.mat){
  K <- 500
  temp1 <- temp2 <- 0
  mid.mat <- (adjacent.mat > 0) *
    (adjacent.mat * trigamma(x[1]) - 1 / x[2]  + temp1)
  value1 <- tau[, 1] %*% mid.mat %*% tau[, 2]
  value1 <- value1 + (N[1]-1) * trigamma(x[1]) - N[2] / x[2] -
    ((digamma(x[1]) - log(x[2])) + (x[1] - 1) * trigamma(x[1]) - 1 + log(x[2]) - digamma(x[1]))
  mid.mat <- (adjacent.mat > 0) *
    (adjacent.mat * (- 1 / x[2]) + x[1] / x[2]**2 +  + temp2)
  value2 <- tau[, 1] %*% mid.mat %*% tau[, 2]
  value2 <- value2 + (N[1]-1) * (- 1 / x[2]) + N[2] * x[1] / x[2]**2 -
    ((x[1] - 1) * (- 1 / x[2]) + x[1] / x[2])
  return(c(-value1, -value2))
}

# This function is used for update delta^theta in the M-step
M.step.dtheta <- function(d.theta, tau, N, adjacent.tensor){
  community.num <- dim(d.theta)[1]
  adjacent.mat <- apply(adjacent.tensor, c(1, 2), sum)
  node.num <- dim(adjacent.mat)[1]
  logf <- matrix(0, node.num, node.num)
  for (i in seq(node.num)) {
    for (j in seq(node.num)) {
      if(adjacent.mat[i,j] > 0){
        logf[i,j] <- sum(log(seq(adjacent.mat[i,j])))
      }
    }
  }
  for (r in seq(community.num)) {
    for (s in seq(community.num)) {
      results <- constrOptim(d.theta[r, s, ], f = obj.theta, grad = gra.theta, tau = tau[, c(r, s)], N = N[r, s, ], adjacent.mat = adjacent.mat, logf.mat = logf,
                             ui = rbind(c(1, 0), c(0, 1)), ci = c(0, 0), method = 'BFGS', control = list(maxit = 800))
      d.theta[r, s, ] <- results$par
      if (results$convergence == 1){
        warning('reached max iteration number')
      }
    }
  }
  return(d.theta)
}

# This function is for the M-step of the algorithm 
M.step.beta1 <- function(adjacent.tensor, attribute.tensor, beta1, beta2, epsilon = 1e-10, no.cluster, weight){
  network.num <- dim(attribute.tensor)[3]
  beta1$nu <- M.step.nu(beta1$nu, beta1$tau, beta2$alpha)
  beta1$nu[which(beta1$nu == 0)] <- epsilon
  beta1$delta.mu <- M.step.dmu(beta1$delta.mu, beta1$tau, beta2$eta, attribute.tensor, weight)
  beta1$delta.mu[which(beta1$delta.mu == 0)] <- epsilon
  par.fun <- function(g){
    d.pi <- M.step.dpi(beta1$delta.pi[, , , g], beta1$tau[, , g], beta2$M, adjacent.tensor[, , , g])
    d.variphi <- M.step.dvariphi(beta1$delta.variphi[, , , g], beta1$tau[, , g], beta2$omega, adjacent.tensor[, , , g])
    d.theta <- M.step.dtheta(beta1$delta.theta[, , , g], beta1$tau[, , g], beta2$N, adjacent.tensor[, , , g])
    return(list(d.pi = d.pi, d.variphi = d.variphi, d.theta = d.theta))
  }
  cl <- makeCluster(no.cluster)
  clusterExport(cl, c('beta1', 'beta2', 'adjacent.tensor'), envir = environment())
  clusterExport(cl, c('M.step.dpi', 'M.step.dvariphi', 'obj.theta', 'gra.theta', 'M.step.dtheta'))
  result <- parLapply(cl, 1:network.num, par.fun)
  stopCluster(cl)
  for (g in seq(network.num)) {
    beta1$delta.pi[, , , g] <- result[[g]]$d.pi
    beta1$delta.variphi[, , , g] <- result[[g]]$d.variphi
    beta1$delta.theta[, , , g] <- result[[g]]$d.theta
  }
  beta1$delta.pi[which(beta1$delta.pi == 0)] <- epsilon
  beta1$delta.variphi[which(beta1$delta.variphi == 0)] <- epsilon
  beta1$delta.theta[which(beta1$delta.theta == 0)] <- epsilon
  return(beta1)
}
# This is the main function for the # This function is for the E-step of the algorithm VEM algorithm
# This function get the adjacent tensor and attribute matrix and return all the variational parameters and tau
opt.beta1 <- function(adjacent.tensor, attribute.tensor, beta1, beta2, max.iter, tol = 1e-2, no.cluster, weight){
  iter.index <- 0
  old.beta1 <- beta1
  difference.indicator <- 1
  tau.diff <- 10
  while ((iter.index <= max.iter) & (max(tau.diff) > tol)) {
    print('beta1')
    print(iter.index)
    beta1$tau <- E.step.beta1(adjacent.tensor = adjacent.tensor, attribute.tensor = attribute.tensor, beta1 = old.beta1, max.iter = 100, no.cluster = no.cluster, weight = weight)
    beta1 <- M.step.beta1(adjacent.tensor = adjacent.tensor, attribute.tensor = attribute.tensor, beta1 = beta1, beta2 = beta2, no.cluster = no.cluster, weight = weight)
    iter.index <- iter.index + 1
    tau.diff <- c()
    nu.diff <- c()
    d.pi.diff <- c()
    d.theta.diff <- c()
    d.variphi.diff <- c()
    for (g in seq(dim(adjacent.tensor)[4])) {
      tau.diff[g] <- sum(abs(old.beta1$tau[, , g]-beta1$tau[, , g]))
      nu.diff[g] <- sum(abs(old.beta1$nu[g, ] - beta1$nu[g, ]))
      d.pi.diff[g] <- sum(abs(old.beta1$delta.pi[, , , g] - beta1$delta.pi[, , , g]))
      d.theta.diff[g] <- sum(abs(old.beta1$delta.theta[, , , g] - beta1$delta.theta[, , , g]))
      d.variphi.diff[g] <- sum(abs(old.beta1$delta.variphi[, , , g] - beta1$delta.variphi[, , , g]))
    }
    old.beta1 <- beta1
  }
  return(beta1 = beta1)
}
