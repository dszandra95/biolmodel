library(igraph)
library(lsa)
plots <- c()
# generate.network.B <- function(N, links.per.step)
# Use: generates a random network by adding nodes one by one and linking the new nodes to old ones such that old nodes with more existing links
#      are more likely to receive the new links (preferential attachment). This is simple example of a Barabási random graph.
# Input: 
#      N: the number of nodes in the network to be generated
#      links.per.step: the number of new links to be added after each node
# Output:
#      L: a 2-column matrix of all links with each row defining one link (i.e. incidence list). Each link occurs only once.

# Define metadata for individuals
define_metadata <- function(N){
  L <- as.data.frame(matrix(nrow=N,ncol=3)) # initialize matrix with zero rows
  colnames(L) <- c('ID', 'age', 'sex')
  for (i in 1:N){
    # Define sex 
    node_sex <- sample(c('F', 'M'), size = 1, prob = c(0.5, 0.5))
    # Define age group
    node_age <- sample(c('child', 'adult', 'elderly'), size = 1, prob = c(0.24, 0.63, 0.13))
    L[i, 1] <- i
    L[i, 2] <- node_age
    L[i, 3] <- node_sex
    
  }
  L
}
# Get metadata information
metadata <- define_metadata(N)
print(metadata)

generate.network.B <- function(N,links.per.step){
  L <- matrix(nrow=0,ncol=2) # initialize matrix with zero rows
  deg <- integer(N) # initialize vector for degrees of the nodes
  for (i in 2:N) { # note that we start assigning edges from the second node only
    
    # Make connections more probable between similar age groups
    # See if age group is the same - results in true/false values
    is_age_same <- metadata[i,2] == metadata[, 2]
    
    # If you lower the added value, the difference will be more drastic
    # with 0.5 value you get a 1:3 ratio
    prob_times <- is_age_same + 0.5
    
    n.new <- min(links.per.step,i-1) # this is needed because at the beginning there might be too few old nodes to select from
    # Define probability of new connection
    myprob <- (deg[1:(i-1)]+ 1) * prob_times[1:(i-1)]
    
    # Select new nodes
    linkto <- sample(i-1,n.new,prob=myprob) # in this example, the "attractiveness of isolated nodes is 1, and it increases linearly with degree
    newlinks <- cbind(rep(i,n.new),linkto)
    L <- rbind(L,newlinks) # append the new links to the matrix of all links
    deg[i] = deg[i] + n.new
    deg[linkto] = deg[linkto] + 1
  }
  colnames(L) <- NULL # remove the column name automatically added by cbind
  L
}


# Simulate epidemics with barabasi graph

N = 200 # number of individuals
simlength <- 30 # number of time steps

plot.spread <- TRUE # switch to set whether you would like to plot the spreading of the epidemic over the network (slows down the simulation when TRUE)


links <- generate.network.B(N,2) # generate an edge list for a Barabási graph
infected <- logical(N) # initialize infection status
patientzero <- sample(N,1) # select 'patient zero'

infected[patientzero] <- TRUE


if (plot.spread) {
  network.i <- graph.edgelist(links,directed=FALSE)
  fixlayout <- layout.kamada.kawai(network.i)  # store a fixed layout for the graph
  node.colour <- rep("SkyBlue2",N) # initialize node colours (SkyBlue2 is also the default node colour in igraph)
    node.colour[patientzero] <- "red" # infected nodes will be coloured red
  for (person in 1:N){
    if (metadata[person,2] == 'elderly'){
    node.colour[person] <- "green"
    }
    else if (metadata[person,2] == 'child'){
      node.colour[person] <- 'pink'
    }
    else{
      node.colour[person] <- "blue"
    }
  }
  plot(network.i,layout=fixlayout, main="Time = 0", vertex.color=node.colour)
}

# Define transmission probability based on age group for all connections based on transmitter
# Index of value in vector will correspond to row of interaction in links
transm_prob_vector = c()
for (n in 1:nrow(links)){
  if (metadata[links[n,][1], 2] == 'child'){
    transm_prob_vector[n] <- 0.8
  }
  else{
    transm_prob_vector[n] <- 0.2
  }
}

color_vector <- c('orange', 'blue', 'red', 'green', 'purple', 'pink', 'yellow')

for (nn in 1:1){
  time <- c()
  num_infected <- c()
  
  links <- generate.network.B(N,2) # generate an edge list for a Barabási graph
  infected <- logical(N) # initialize infection status
  patientzero <- sample(N,1) # select 'patient zero'
  
  infected[patientzero] <- TRUE
  for (i in 1:simlength) {
    discordant.links <- which(xor(infected[links[,1]],infected[links[,2]])) # find the indeces of links that connect an infected individual to an uninfected
    
    # Determine which links of infected individual (discordant links) will transmit the disease
    # Define counter to follow index of does_transmit vector
    # Does_transmit: same dimension as links, 1 if will get infected, 0 if not
    counter <- 1
    tarnsmit <- c()
    
    for (conn in discordant.links){
      # generate random number between 1 and 0
      r <- runif(1)
      # If random number is smaller than transmitting probability, the link will be infectious
      # Transmitting probability is defined based on the transmitter's age group (see above)
      if (r < transm_prob_vector[conn]){
        transmit[counter] <- 1
        counter <- counter + 1
      }
      else{
        transmit[counter] <- 0
        counter <- counter + 1
      }
    }
    
    # let me update the infection vector in three steps to make it easier to read:
    transmitter.links <- discordant.links[transmit==1]
    nodes.of.transmitter.links <- unique(as.vector(links[transmitter.links,1:2])) # gets both nodes of the transmitter links into a single vector; unique just filters out repetitions
    infected[nodes.of.transmitter.links] <- TRUE # here I simply set both nodes to TRUE (although the transmitter already had 'TRUE'). In more complex models, you might want to do a further check here and overwrite only the newly infected nodes.
    
    num_infected[i] <- sum(infected, na.rm = TRUE)
    time[i] <- i
  
  }
  
  if (nn == 1){
    plot(time, num_infected, type = 'b', pch=16, col = color_vector[nn], main=' Same probs: 33% child') 
  }
  else{
    lines(time, num_infected, type = 'b', pch=16, col = color_vector[nn])
  }
}


