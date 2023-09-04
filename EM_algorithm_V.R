library(doParallel)
library(parallel)
library(LaplacesDemon)

## This script is for the VEM algorithm of sequential-NGM, all the functions here are similar to 
## those in EM_algorithm.R

series <- function(k, delta){
  return((delta[2] / (delta[2] + k)) ** delta[1] / k)
}

E.step.topology <- function(i, r, adjacent.tensor, adjacent.mat, beta1, tau1, rs.mat){
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
        value <- value + tau1[j, s] * ((adjacent.mat[i, j] == 0)*(digamma(beta1$delta.pi[r, s, 1]) - digamma(sum(beta1$delta.pi[r, s, ]))) + 
                                         (adjacent.mat[i, j] > 0) * ((digamma(beta1$delta.pi[r, s, 2]) - digamma(sum(beta1$delta.pi[r, s, ]))) + adjacent.mat[i, j] *
                                         (digamma(beta1$delta.theta[r, s, 1]) - log(beta1$delta.theta[r, s, 2])) - beta1$delta.theta[r, s, 1] / beta1$delta.theta[r, s, 2] - 
                                         summationij + rs.mat[r, s] + sum(adjacent.tensor[i, j,] * 
                                         (digamma(beta1$delta.variphi[r, s, ]) - digamma(sum(beta1$delta.variphi[r, s, ]))))) + 
                                         (adjacent.mat[j,i] == 0) * (digamma(beta1$delta.pi[s, r, 1]) - digamma(sum(beta1$delta.pi[s, r, ]))) + 
                                         (adjacent.mat[j,i] > 0) * ((digamma(beta1$delta.pi[s, r, 2]) - digamma(sum(beta1$delta.pi[s, r, ]))) + adjacent.mat[j,i] *
                                         (digamma(beta1$delta.theta[s, r, 1]) - log(beta1$delta.theta[s, r, 2])) - beta1$delta.theta[s, r, 1] / beta1$delta.theta[s, r, 2] - 
                                         summationji + rs.mat[s, r] + sum(adjacent.tensor[j, i, ] * 
                                         (digamma(beta1$delta.variphi[s, r, ]) - digamma(sum(beta1$delta.variphi[s, r, ]))))))
      } 
    }
  }
  return(value)
}

E.step.attribute.online <- function(i, r, attribute.tensor, mu, weight){
  attribute.num <- dim(attribute.tensor)[2]
  value <- 0
  for (t in seq(attribute.num)) {
    value <- value + weight * attribute.tensor[i, t] * log(mu[r, t])
  }
  return(value)
}

E.step.remain <- function(r, nu){
  return(digamma(nu[r]) - digamma(sum(nu)))
}

E.step.iter.online <- function(adjacent.tensor, attribute.tensor, beta1, tau1, rs.mat, mu, weight){
  node.num <- dim(adjacent.tensor)[1]
  community.num <- dim(rs.mat)[1]
  log.tau <- matrix(NA, nrow = node.num, ncol = community.num)
  adjacent.mat <- apply(adjacent.tensor, c(1, 2), sum)
  for (i in seq(node.num)) {
    for (r in seq(community.num)) {
      log.tau[i, r] <- E.step.topology(i, r, adjacent.tensor,  adjacent.mat, beta1, tau1, rs.mat) +
        E.step.attribute.online(i, r, attribute.tensor, mu, weight) +
        E.step.remain(r, beta1$nu)
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

E.step.beta1.online <- function(adjacent.tensor, attribute.tensor, beta1, mu, max.iter, tol = 1e-8, no.cluster, weight){
  node.num <- dim(adjacent.tensor)[1]
  difference <- 10
  iter.index <- 0
  tau1 <- beta1$tau
  tau2 <- tau1
  rs.mat <- matrix(0, nrow = 3, ncol = 3)
  while((iter.index <= max.iter) & (difference > tol)){
    log.tau <- E.step.iter.online(adjacent.tensor, attribute.tensor, beta1, tau1, rs.mat, mu, weight)
    for (i in seq(node.num)) {
      tau2[i, ] <- exp(log.tau[i, ] - max(log.tau[i, ])) / sum(exp(log.tau[i, ] - max(log.tau[i, ])))
    }
    difference <- max(abs(tau1 - tau2))
    tau1 <- tau2
    iter.index <- iter.index + 1
  }
  tau1[which(tau1 == 0)] <- 1e-10
  beta1$tau<- tau1
  new.tau <- beta1$tau
  for (i in seq(node.num)){
    beta1$tau[i, ] <- beta1$tau[i, ] / sum(new.tau[i, ])
  }
  return(beta1$tau)
}


M.step.nu.online <- function(nu, tau, alpha){
  community.num <- dim(tau)[2]
  for (r in seq(community.num)) {
    nu[r] <- sum(tau[, r]) + alpha[r]
  }
  return(nu)
}

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


M.step.beta1.online <- function(adjacent.tensor, attribute.tensor, beta1, beta2, epsilon = 1e-10, no.cluster){
  beta1$nu <- M.step.nu.online(beta1$nu, beta1$tau, beta2$alpha)
  beta1$delta.pi <- M.step.dpi(beta1$delta.pi, beta1$tau, beta2$M, adjacent.tensor)
  beta1$delta.variphi <- M.step.dvariphi(beta1$delta.variphi, beta1$tau, beta2$omega, adjacent.tensor)
  beta1$delta.theta <- M.step.dtheta(beta1$delta.theta, beta1$tau, beta2$N, adjacent.tensor)
  return(beta1)
}


opt.beta1.online <- function(adjacent.tensor, attribute.tensor, beta1, mu, beta2, max.iter, tol = 1e-2, no.cluster, weight){             # optimize variation parameters
  iter.index <- 0
  old.beta1 <- beta1
  tau.diff <- 10
  difference.indicator <- 1
  ELBO.old.list <- c()
  while ((iter.index <= max.iter) & ((tau.diff) > tol)) {
    print('beta1')
    print(iter.index)
    beta1$tau <- E.step.beta1.online(adjacent.tensor = adjacent.tensor, attribute.tensor = attribute.tensor, 
                                     beta1 = old.beta1, mu = mu, max.iter = 100, no.cluster = no.cluster, weight = weight)
    beta1 <- M.step.beta1.online(adjacent.tensor = adjacent.tensor, attribute.tensor = attribute.tensor, 
                                 beta1 = beta1, beta2 = beta2, no.cluster = no.cluster)
    iter.index <- iter.index + 1
    tau.diff <- sum(abs(old.beta1$tau-beta1$tau))
    print('tau.diff')
    print(tau.diff)
    old.beta1 <- beta1
  }
  return(beta1 = beta1)
}

# This function generates multilayer attributed network using the input model parameters
# The output is adjacent tensor, attribute matrix and community ID matrix
net.ger <- function(node.num, layer.num, attribute.num, gamma, mu, pi, theta, varphi){
  adjacent.tensor <- array(data = 0, dim = c(node.num, node.num, layer.num))
  attribute.tensor <- array(data = 0, dim = c(node.num, attribute.num))
  node.community <- matrix(data = NA, nrow = 1, node.num, byrow = TRUE)
  for (i in seq(node.num)){
    node.community[i] <- which(rmultinom(1, 1, prob = gamma) == 1, arr.ind = TRUE)[1,1]
    for (t in seq(attribute.num)) {
      if(mu[node.community[i], t] == 0.1){
        attribute.tensor[i, t] <- rbern(1, 0.9)
      }
      else{
        attribute.tensor[i, t] <- rbern(1, 0.1)
      }
    }
  }
  for (i in seq(node.num)){
    for (j in seq(node.num)){
      if(i != j){
        indicator <- rbern(1, 1 - pi[node.community[i], node.community[j]])
        if(indicator){
          total.edge <- rpois(1, lambda = theta[node.community[i], node.community[j]])
          while(total.edge == 0){
            total.edge <- rpois(1, lambda = theta[node.community[i], node.community[j]])
          }
          adjacent.tensor[i, j, ] <- rmultinom(1, total.edge, prob = varphi[node.community[i], node.community[j],])
        }
      }
    }
  }
  return(list(attribute.tensor = attribute.tensor, adjacent.tensor = adjacent.tensor, node.community = node.community))
}