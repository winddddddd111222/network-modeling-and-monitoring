
library(doParallel)
library(parallel)
library(LaplacesDemon)
source("ELBO_calculation.R")
source("EM_algorithm.R")

## This script provide an example of using VEM perform static-NGM parameters estimation

# data preparation. Please note that all three .csv files are outputs of simulation_data_generation script
# Thus, all the basic settings, network number, node number, community number, attribute number, should be the 
# same as those in generation process
# Also, here we allow parallel computation. If one does not need, please set no.cluster to 1
adjacent.matrix <- read.csv(file = './adjacent_tensor.csv', sep = ',', header = F)
adjacent.matrix <- as.matrix(adjacent.matrix)
attribute.matrix <- read.csv(file = './attribute_tensor.csv', sep = ',', header = F)
attribute.matrix <- as.matrix(attribute.matrix)
node.community <- read.csv(file = './node_community.csv', sep = ',', header = F)
node.community <- as.matrix(node.community)
max.iter1 <- 1000
max.iter2 <- 1000
network.num <- 100 # number of networks
node.num <- 80 # node number in each network
layer.num <- 3 # layer number in each network
attribute.num <- 24 # total attribute number
community.num <- 3 # community number in each network
epsilon <- 1e-10
sigma <- 10 # the value of weighted coefficient sigma
no.cluster <- 50 # for parallel computation, define the kernal number used
adjacent.tensor <- array(data = NA, dim = c(node.num, node.num, layer.num, network.num))
attribute.tensor <- array(data = NA, dim = c(node.num, attribute.num, network.num))
for (g in seq(network.num)){
  for (l in seq(layer.num)){
    adjacent.tensor[ , , l, g] <- adjacent.matrix[(node.num * (l - 1) + 1): (node.num * l), (node.num * (g - 1) + 1): (node.num *g)]
  }
}
for (g in seq(network.num)){
  attribute.tensor[ , , g] <- attribute.matrix[, (attribute.num * (g - 1) + 1) :(attribute.num *g)] 
}

# define variation parameters & hyperparameters
nu <- matrix(data = NA, nrow = network.num, ncol = community.num)
tau <- array(data = NA, dim = c(node.num, community.num, network.num))
delta.mu <- matrix(data = NA, nrow = community.num, ncol = attribute.num, byrow = TRUE)
delta.pi <- array(data = NA, dim = c(community.num, community.num, 2, network.num))
delta.theta <- array(data = 1, dim = c(community.num, community.num, 2, network.num))
delta.variphi <- array(data = NA, dim = c(community.num, community.num, layer.num, network.num))
# variation parameter initialization
tau[ , ,] <- 1 / community.num
# hyperparameters setting

alpha <- c(0.05 * node.num / community.num, 0.05 * node.num / community.num, 
           0.05 * node.num / community.num)

eta <- matrix(0.02, nrow = community.num, ncol = attribute.num)

eta[1, 1:8] <- 0.05
eta[2, 9:16] <- 0.05
eta[3, 17:24] <- 0.05

temp <- node.num * (node.num - 1) / community.num ** 2

M <- c( 0.02, 0.05, 0.05, 
        0.05, 0.02, 0.05, 
        0.05, 0.05, 0.02,
        0.08, 0.05, 0.05, 
        0.05, 0.08, 0.05,
        0.05, 0.05, 0.08) * temp 
M <- array(data = M, dim = c(community.num, community.num, 2))
N <- c(0.80, 0.25, 0.25,  
       0.25, 0.80, 0.25,  
       0.25, 0.25, 0.80,  
       0.08, 0.05, 0.05, 
       0.05, 0.08, 0.05,
       0.05, 0.05, 0.08) * temp
N <- array(data = N, dim = c(community.num, community.num, 2))
omega <- array(data = 1, dim = c(community.num, community.num, layer.num))
nu <- M.step.nu(nu, tau, alpha)

delta.mu <- M.step.dmu(delta.mu, tau, eta, attribute.tensor, weight = sigma)

par.fun <- function(g){
  d.pi <- M.step.dpi(delta.pi[, , , g], tau[, , g], M, adjacent.tensor[, , , g])
  d.variphi <- M.step.dvariphi(delta.variphi[, , , g], tau[, , g], omega, adjacent.tensor[, , , g])
  d.theta <- M.step.dtheta(delta.theta[, , , g], tau[, , g], N, adjacent.tensor[, , , g])
  return(list(d.pi = d.pi, d.variphi = d.variphi, d.theta = d.theta))
}
cl <- makeCluster(no.cluster)
clusterExport(cl, c('delta.pi', 'tau', 'delta.variphi', 'M', 'omega', 'N', 'delta.theta', 'adjacent.tensor'))
clusterExport(cl, c('M.step.dpi', 'M.step.dvariphi', 'gra.theta','obj.theta', 'M.step.dtheta'))
result <- parLapply(cl, 1:network.num, par.fun)
stopCluster(cl)
for (g in seq(network.num)) {
  delta.pi[, , , g] <- result[[g]]$d.pi
  delta.variphi[, , , g] <- result[[g]]$d.variphi
  delta.theta[, , , g] <- result[[g]]$d.theta
}
beta1 <- list(nu = nu, tau = tau, delta.mu = delta.mu, delta.pi = delta.pi, delta.theta = delta.theta, delta.variphi = delta.variphi)
beta2 <- list(alpha = alpha, eta = eta, omega = omega, M = M, N = N)

# Run the VEM algorithm to estimate model parameter
start_time <- Sys.time()
results <- opt.beta1(adjacent.tensor = adjacent.tensor, attribute.tensor = attribute.tensor, 
                     beta1 = beta1, beta2 = beta2, max.iter = 1000, tol = 1e-10, 
                     no.cluster = no.cluster, weight = sigma)
end_time <- Sys.time()

# inference results output, note that all the output are variational parameters, and the estimated model parameters
# can be abtained by calculate the mean of variational distributions
write.table(results$nu, './nu.csv', sep = ',', row.names = FALSE, col.names = FALSE)
tau.exp <- results$tau[, , 1]
for (g in seq(network.num-1)) {
  tau.exp <- cbind(tau.exp, results$tau[, , g + 1])
}
write.table(tau.exp, './tau.csv', sep = ',', row.names = FALSE, col.names = FALSE)
d.pi.exp <- rbind(results$delta.pi[, , 1, 1], results$delta.pi[, , 2, 1])
for (g in seq(network.num-1)) {
  d.pi.exp <- cbind(d.pi.exp, rbind(results$delta.pi[, , 1, g + 1], results$delta.pi[, , 2, g + 1]) )
}
write.table(d.pi.exp, './delta_pi.csv', sep = ',', row.names = FALSE, col.names = FALSE)
d.theta.exp <- rbind(results$delta.theta[, , 1, 1], results$delta.theta[, , 2, 1])
for (g in seq(network.num-1)) {
  d.theta.exp <- cbind(d.theta.exp, rbind(results$delta.theta[, , 1, g + 1], results$delta.theta[, , 2, g + 1]) )
}
write.table(d.theta.exp, './delta_theta.csv', sep = ',', row.names = FALSE, col.names = FALSE)
d.variphi.exp <- results$delta.variphi[, , 1, 1]
for (l in seq(layer.num - 1)) {
  d.variphi.exp <- rbind(d.variphi.exp, results$delta.variphi[, , l + 1, 1])
}
for (g in seq(network.num - 1)) {
  d.variphi.layer <- results$delta.variphi[, , 1, g + 1]
  for (l in seq(layer.num - 1)) {
    d.variphi.layer <- rbind(d.variphi.layer, results$delta.variphi[, , l + 1, g + 1])
  }
  d.variphi.exp <- cbind(d.variphi.exp, d.variphi.layer)
}
write.table(d.variphi.exp, './delta_variphi.csv', sep = ',', row.names = FALSE, col.names = FALSE)
d.mu <- results$delta.mu
write.table(d.mu, './delta_mu.csv', sep = ',', row.names = FALSE, col.names = FALSE)

