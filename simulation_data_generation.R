library(LaplacesDemon)

## This R script is used to generate simulated network data according to the setting parameters, and the output file 
## record flattened adjacent tensor, node attribute matrix and node community ID

network.num <- 2 # simulated network number
node.num <- 45 # node number for each network
layer.num <- 3 # layer number for each network
attribute.num <- 24 # total attribute number
community.num <- 3 # community number for each network
adjacent.tensor <- array(data = NA, dim = c(node.num, node.num, layer.num, network.num))
attribute.tensor <- array(data = 0, dim = c(node.num, attribute.num, network.num))
node.community <- matrix(data = NA, nrow = network.num, node.num, byrow = TRUE)

pagamma <- c(1/3, 1/3, 1/3)

# model parameter setting for simulated networks
mu <- matrix(0.2/16, nrow = community.num, ncol = attribute.num)
for (r in seq(community.num)) {
  mu[r, ((8 * (r - 1) + 1) : (8 * r))] <- 0.1
}
pi <- c(0.20, 0.50, 0.50, 
        0.50, 0.23, 0.50, 
        0.50, 0.50, 0.20)
pi <- matrix(data = pi, nrow = community.num, ncol = community.num, byrow = T)
theta <- c(20, 10, 10, 
           10, 25, 10, 
           10, 10, 24) # the parameter of the prior distribution of theta 
theta <- matrix(data = theta, nrow = community.num,  ncol = community.num, byrow = T)
varphi <- c(1/7,  1/3,  2/3, 
            1/3,  3/5,  1/3, 
            1/3,  1/6,  3/5, 
            4/7,  1/3,  1/6, 
            1/3, 1/10,  1/3, 
            1/3,  2/3, 1/10,
            2/7, 1/3, 1/3,
            1/3, 1/10, 1/3,
            1/3, 1/3, 1/10)
varphi <- array(data = varphi, dim = c(community.num, community.num, layer.num))
# network data generation
for (g in seq(network.num)){
  print(g)
  # attribute data generation
  for (i in seq(node.num)){
    node.community[g, i] <- which(rmultinom(1, 1, prob = pagamma) == 1, arr.ind = TRUE)[1,1]
    for (t in seq(attribute.num)) {
      if(mu[node.community[g, i], t] == 0.1){
        attribute.tensor[i, t, g] <- rbern(1, 0.9)
      }
      else{
        attribute.tensor[i, t, g] <- rbern(1, 0.1)
      }
    }
  }
  # topology data generation
  for (i in seq(node.num)){
    for (j in seq(node.num)){
      if(i != j){
        indicator <- rbern(1, 1 - pi[node.community[g, i], node.community[g, j]])
        if(indicator){
          total.edge <- rpois(1, lambda = theta[node.community[g, i], node.community[g,j]])
          while(total.edge == 0){
            total.edge <- rpois(1, lambda = theta[node.community[g, i], node.community[g, j]])
          }
          adjacent.tensor[i , j, , g] <- rmultinom(1, total.edge, prob = varphi[node.community[g, i], node.community[g, j],])
        }
      }
    }
  }
}
adjacent.tensor[is.na(adjacent.tensor)] <- 0

# tensor expansion
attribute.exp <- attribute.tensor[ , , 1]
for (g in seq(network.num - 1)){
  attribute.exp <- cbind(attribute.exp, attribute.tensor[ , ,g+1])
}
adjacent.exp <- adjacent.tensor[ , , 1, 1]
for (l in seq(layer.num -1)){
  adjacent.exp <- rbind(adjacent.exp, adjacent.tensor[ , , l + 1, 1])
}
for (g in seq(network.num - 1)){
  adjacent.layer <- adjacent.tensor[, , 1, g + 1]
  for (l in seq(layer.num - 1)){
    adjacent.layer <- rbind(adjacent.layer, adjacent.tensor[ , , l + 1, g + 1])
  }
  adjacent.exp <- cbind(adjacent.exp, adjacent.layer)
}

#simulation data output
path1 <- './adjacent_tensor.csv'
path2 <- './attribute_tensor.csv'
path3 <- './node_community.csv'
write.table(adjacent.exp, path1, sep = ',', row.names = FALSE, col.names = FALSE)
write.table(attribute.exp, path2, sep = ',', row.names = FALSE, col.names = FALSE)
write.table(node.community, path3, sep = ',', row.names = FALSE, col.names = FALSE)



