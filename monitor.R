library(doParallel)
library(parallel)
library(LaplacesDemon)
#source("ELBO_calculation.R")
source("EM_algorithm_V.R")

repli.time <- 10000
no.cluster <- 100

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
        attribute.tensor[i, t] <- rbern(1, 0.9/8)
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

fun <- function(replication){
  max.iter1 <- 1000
  max.iter2 <- 1000
  node.num <- 45
  layer.num <- 3
  attribute.num <- 24
  community.num <- 3
  epsilon <- 1e-10
  rho <- 0.05
  weight <- 3 * node.num / attribute.num
  # U.gamma <- 7.6404
  # U.pi <- 15.8270
  # U.theta <- 15.8160
  # U.varphi <- 24.1968
  U.gamma <- 12.9370
  U.pi <- 25.9151
  U.theta <- 25.8597
  U.varphi <- 39.4552
  # adjacent.tensor <- array(data = 0, dim = c(node.num, node.num, layer.num, network.num))
  # attribute.tensor <- array(data = 0, dim = c(node.num, attribute.num, network.num))
  mu <- matrix(0.2/16, nrow = community.num, ncol = attribute.num)
  for (r in seq(community.num)) {
    mu[r, ((8 * (r - 1) + 1) : (8 * r))] <- 0.1
  }
  nu <- matrix(data = NA, nrow = 1, ncol = community.num)
  gamma0 <- c(1/3, 1/3, 1/3)
  gamma <- gamma0
  pi0 <- c(0.2, 0.5, 0.5, 0.5, 0.23, 0.5, 0.5, 0.5, 0.2)
  pi0 <- matrix(data = pi0, nrow = community.num, ncol = community.num, byrow = T)
  pi <- pi0
  theta0 <- c(20, 10, 10, 10, 25, 10, 10, 10, 24)
  theta0 <- matrix(data = theta0, nrow = community.num,  ncol = community.num, byrow = T)
  theta <- theta0
  varphi0 <- c(1/7, 1/3, 1/3, 1/3, 3/5, 1/3, 1/3, 1/3, 3/5, 4/7, 1/3, 1/3, 1/3, 1/10, 1/3, 1/3, 1/3, 1/10, 2/7, 1/3, 1/3, 1/3, 3/10, 1/3, 1/3, 1/3, 3/10)
  varphi0 <- array(data = varphi0, dim = c(community.num, community.num, layer.num))
  varphi <- varphi0
  tau <- array(data = 0, dim = c(node.num, community.num))
  delta.pi <- array(data = NA, dim = c(community.num, community.num, 2))
  delta.theta <- array(data = 3000, dim = c(community.num, community.num, 2))
  delta.theta[, , 2] <- 200
  delta.variphi <- array(data = NA, dim = c(community.num, community.num, layer.num))
  # variation parameter initialization
  alpha <- rho * node.num * gamma0
  # the parameter of the prior distribution of attribute
  M <- array(data = 1, dim = c(community.num, community.num, 2))
  N <- array(data = 1, dim = c(community.num, community.num, 2))
  omega <- array(data = 1, dim = c(community.num, community.num, layer.num))
  for (r in seq(community.num)) {
    for (s in seq(community.num)) {
      M[r, s, 1] <- rho * node.num * (node.num - 1) * gamma0[r] * gamma0[s] * pi0[r, s]
      M[r, s, 2] <- rho * node.num * (node.num - 1) * gamma0[r] * gamma0[s] * (1 - pi0[r, s])
      N[r, s, 1] <- rho * node.num * (node.num - 1) * gamma0[r] * gamma0[s] * (1 - pi0[r, s]) * theta0[r, s]
      N[r, s, 2] <- rho * node.num * (node.num - 1) * gamma0[r] * gamma0[s] * (1 - pi0[r, s])
      omega[r, s, ] <- rho * node.num * (node.num - 1) * gamma0[r] * gamma0[s] * (1 - pi0[r, s]) * theta0[r, s] * varphi0[r, s, ]
    }
  }
  beta1 <- list(nu = nu, tau = tau, delta.pi = delta.pi, delta.theta = delta.theta, delta.variphi = delta.variphi)
  beta2 <- list(alpha = alpha, omega = omega, M = M, N = N)
  kl.pi <- array(NA, dim = c(community.num, community.num))
  kl.theta <- array(NA, dim = c(community.num, community.num))
  kl.varphi <- array(NA, dim = c(community.num, community.num))
  flag <- 1
  RL <- 0
  abnormal.chart <- 0
  change.point <- 5
  t <- 0
  while(flag == 1) {
    t <- t + 1
    if(t > change.point){
      #abnormal.gamma <- gamma0 + c(-0.1, 0, 0.1)
      # abnormal.pi <- pi0
      # abnormal.pi[1,1] <- abnormal.pi[1,1] + 0.01
      network.set <- net.ger(node.num = node.num, layer.num = layer.num, 
                             attribute.num = attribute.num, gamma = gamma0, mu = mu, pi = pi0, 
                             theta = theta0, varphi = abnormal.varphi) 
    }
    else{
      network.set <- net.ger(node.num = node.num, layer.num = layer.num, 
                             attribute.num = attribute.num, gamma = gamma0, mu = mu, pi = pi0, 
                             theta = theta0, varphi = varphi0) 
    }
    adjacent.tensor <- network.set$adjacent.tensor
    attribute.tensor <- network.set$attribute.tensor
    node.community <- network.set$node.community
    beta1$tau[ ,] <- 1 / community.num
    # for (i in seq(node.num)) {
    #   beta1$tau[i , node.community[i]] <- 1
    # }
    beta1$nu <- M.step.nu.online(beta1$nu, beta1$tau, beta2$alpha)
    beta1$delta.pi <- M.step.dpi(beta1$delta.pi, beta1$tau, beta2$M, adjacent.tensor)
    beta1$delta.variphi <- M.step.dvariphi(beta1$delta.variphi, beta1$tau, beta2$omega, adjacent.tensor)
    beta1$delta.theta <- M.step.dtheta(beta1$delta.theta, beta1$tau, beta2$N, adjacent.tensor)
    beta1 <- opt.beta1.online(adjacent.tensor = adjacent.tensor, attribute.tensor = attribute.tensor,
                              beta1 = beta1, mu = mu, beta2 = beta2, max.iter = 1000, tol = 1e-15,
                              no.cluster = network.num, weight = weight)
    gamma <- beta1$nu / sum(beta1$nu)
    pi <- beta1$delta.pi[, , 1] / apply(beta1$delta.pi, c(1, 2), sum)
    theta <- beta1$delta.theta[, , 1] / beta1$delta.theta[, , 2]
    for (l in seq(layer.num)) {
      varphi[, , l] <- beta1$delta.variphi[, , l] / apply(beta1$delta.variphi, c(1, 2), sum)
    }
    kl.gamma <- sum(gamma * log(gamma / gamma0))
    Q.gamma <- 2 * (1 + rho) * node.num * kl.gamma
    Q.pi <- Q.theta <- Q.varphi <- 0
    for (r in seq(community.num)) {
      for (s in seq(community.num)) {
        kl.pi[r, s] <- pi[r, s] * log(pi[r, s]/pi0[r, s]) +(1 - pi[r, s]) * log((1 - pi[r, s])/(1 - pi0[r, s]))
        #Q.pi[g] <- Q.pi[g] + 2 * node.num * (node.num - 1) * gamma[r] * gamma[s] * kl.pi[r, s, g]
        Q.pi <- Q.pi + 2 * (beta1$tau[, r] %*% (matrix(1, nrow = node.num, ncol = node.num) - diag(1, nrow = node.num, ncol = node.num)) 
                            %*% beta1$tau[, s] + sum(beta2$M[r, s, ])) * kl.pi[r, s]
        kl.theta[r, s] <- theta[r, s] * log(theta[r, s]/ theta0[r, s]) + theta0[r, s] - theta[r, s]
        # Q.theta[g] <- Q.theta[g] + 2 * node.num * (node.num - 1) * gamma[r] * gamma[s] * (1 - pi[r, s]) * kl.theta[r, s, g]
        Q.theta <- Q.theta +2 * (beta1$tau[, r] %*% (apply(adjacent.tensor, c(1, 2), sum) > 0) %*% beta1$tau[, s] + 
                                         beta2$N[r, s, 2]) * kl.theta[r, s]
        kl.varphi[r, s] <- sum(varphi[r, s, ] * log(varphi[r, s, ]/varphi0[r, s, ]))
        # Q.varphi[g] <- Q.varphi[g] + 2 * node.num * (node.num - 1) * gamma[r] * gamma[s] * (1 - pi[r, s]) * theta[r, s] * kl.varphi[r, s, g]
        Q.varphi <- Q.varphi + 2 * (beta1$tau[, r] %*% (apply(adjacent.tensor, c(1, 2), sum)) %*% beta1$tau[, s]+ sum(beta2$omega[r, s, ])) * kl.varphi[r, s]
      }
    }
    beta2$alpha <- rho * node.num * gamma
    for (r in seq(community.num)) {
      for (s in seq(community.num)) {
        beta2$M[r, s, 1] <- rho * node.num * (node.num - 1) * gamma[r] * gamma[s] * pi[r, s]
        beta2$M[r, s, 2] <- rho * node.num * (node.num - 1) * gamma[r] * gamma[s] * (1 - pi[r, s])
        beta2$N[r, s, 1] <- rho * node.num * (node.num - 1) * gamma[r] * gamma[s] * (1 - pi[r, s]) * theta[r, s]
        beta2$N[r, s, 2] <- rho * node.num * (node.num - 1) * gamma[r] * gamma[s] * (1 - pi[r, s])
        beta2$omega[r, s, ] <- rho * node.num * (node.num - 1) * gamma[r] * gamma[s] * (1 - pi[r, s]) * theta[r, s] * varphi[r, s, ]
      }
    }
    rm(network.set)
    rm(adjacent.tensor)
    rm(attribute.tensor)
    if(t > change.point){
      if(Q.gamma > U.gamma){
        abnormal.chart <- 1
        flag <- 0
        break
      }else if(Q.pi > U.pi){
        abnormal.chart <- 2
        flag <- 0
        break
      }else if(Q.theta > U.theta){
        abnormal.chart <- 3
        flag <- 0
        break
      }else if(Q.varphi > U.varphi){
        abnormal.chart <- 4
        flag <- 0
        break
      }else{
        RL <- RL + 1
      }
    }
  }
  return(list(RL = RL, abnormal.chart = abnormal.chart))
}
abnormal.list <- c(0, 0.001 ,0.002, 0.003, 0.004, 0.005, 0.006, 0.007, 0.008, 0.009, 0.010)
for (abnormal.value in abnormal.list) {
  print(abnormal.value)
  gamma0 <- c(1/3, 1/3, 1/3)
  pi0 <- c(0.2, 0.5, 0.5, 0.5, 0.23, 0.5, 0.5, 0.5, 0.2)
  pi0 <- matrix(data = pi0, nrow = community.num, ncol = community.num, byrow = T)
  theta0 <- c(20, 10, 10, 10, 25, 10, 10, 10, 24)
  theta0 <- matrix(data = theta0, nrow = community.num,  ncol = community.num, byrow = T)
  varphi0 <- c(1/7, 1/3, 1/3, 1/3, 3/5, 1/3, 1/3, 1/3, 3/5, 4/7, 1/3, 1/3, 1/3, 1/10, 1/3, 1/3, 1/3, 1/10, 2/7, 1/3, 1/3, 1/3, 3/10, 1/3, 1/3, 1/3, 3/10)
  varphi0 <- array(data = varphi0, dim = c(community.num, community.num, layer.num))
  abnormal.varphi <- varphi0
  abnormal.varphi[1,1,] <- abnormal.varphi[1, 1,] + c(-abnormal.value, 0, abnormal.value)
  cl <- makeCluster(no.cluster)
  clusterExport(cl, c('series', 'E.step.topology', 'E.step.attribute.online', 'E.step.remain', 'abnormal.varphi',
                      'E.step.iter.online', 'rsmat.cal', 'E.step.beta1.online', 'M.step.nu.online', 'M.step.dpi',
                      'M.step.dvariphi', 'obj.theta', 'M.step.dtheta', 'M.step.beta1.online', 'opt.beta1.online', 'net.ger', 'gra.theta'))
  clusterEvalQ(cl, {library(LaplacesDemon)})
  result <- parLapply(cl, 1:repli.time, fun)
  stopCluster(cl)
  RL <- c()
  abnormal.chart <- c()
  for (index in seq(repli.time)) {
    RL[index] <- result[[index]]$RL
    abnormal.chart[index] <- result[[index]]$abnormal.chart
  }
  path <- paste('./data&results/data-three community-online-UCL/RL-varphi-', abnormal.value, '.csv', sep = '')
  write.table(RL, path, sep = ',', row.names = FALSE, col.names = FALSE)
  path <- paste('./data&results/data-three community-online-UCL/abnormal_chart-varphi-', abnormal.value,'.csv', sep = '')
  write.table(abnormal.chart, path, sep = ',', row.names = FALSE, col.names = FALSE)
}
