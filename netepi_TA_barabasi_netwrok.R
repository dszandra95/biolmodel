####### Barabási-Albert network modell #########

library(igraph)
library(ggplot2)

# sTARTING POINT
N <- 10000 # Number of individualsz
S <- N - 1 # Number of susceptible individuals
I <- N - S # Number of infected individuals
R <- N - S- I # Number of individuals who had recovered from the diseas
Conn <- 0.3
links.per.step <- 3

generate.network.B <- function(N, links.per.step){
  L <- matrix(nrow = 0, ncol = 2)
  deg <- integer(N)
  for (i in 2:N){
    n.new <- min(links.per.step, i-1)
    linkto <- sample(i-1, n.new, prob = deg[1:(i - 1) + 1])
    L <- rbind(L, newlinks)
    deg[i] = deg[i] + n.new
    deg[linkto] = deg[linkto] + 1
  }
  colnames(L) <- NULL
  L
}


simlength <- 15
p.t <- 0.3 
plot.spread <- TRUE
links <- generate.network.B(N, 4)
infected <- logical(N)
patientzero <- sample(N, 1)
infected[patientzero] <- TRUE

if (plot.spread){
  netwotk.i <- graph.edgelist(links, directed = FALSE)
  fixlayout <- layout.kamada.kawai(network.i)
  node.colour <- rep("SkyBlue", N)
    node.colour[patientzero] <- "red"
    plot(netwotk.i, layout=fixlayout, main = "Time = 0", vertex.color = node.colour)
}

for (i in 1:simlength){
  discordant.links <- which(xor(infected[links[, 1]], infected[,2]))
  transmit <- rbinom(length(discordant.links), 1, p.t)
  nodes.of.transmitter.links <- unique(as.vector(links[transmitter.links, 1:2]))
  infected[nodes.of.transmitter.links] <- TRUE
  if (plot.spread){
    node.colour[infected] <- "red"
    plot(network.i, layout = fixlayout, main = paste("Time = 0"), vertex.colour = node.colour)
  }
}

