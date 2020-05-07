library(igraph)

generate.network.B <- function(N,links.per.step){ 
  L <- matrix(nrow=0,ncol=2) # initialize matrix with zero rows
  deg <- integer(N) # initialize vector for degrees of the nodes
  for (i in 2:N) { # note that we start assigning edges from the second node only
    n.new <- min(links.per.step,i-1) # this is needed because at the beginning there might be too few old nodes to select from
    linkto <- sample(i-1,n.new,prob=deg[1:(i-1)]+1) # in this example, the "attractiveness of isolated nodes is 1, and it increases linearly with degree
    newlinks <- cbind(rep(i,n.new),linkto)
    L <- rbind(L,newlinks) # append the new links to the matrix of all links
    deg[i] = deg[i] + n.new
    deg[linkto] = deg[linkto]+1
  }
  colnames(L) <- NULL # remove the column name automatically added by cbind
  L
}

time_infected <- function(N, p.t, plot.spread=F, connectivity=2, ){ # N = number of nodes, p.t = probability of transmitting the disease to an uninfected node, along a connector link, connectivity = network parameter
  
  links <- generate.network.B(N,connectivity) # generate an edge list for a Barabási graph
  infected <- logical(N) # initialize infection status
  patientzero <- sample(N,1) # select 'patient zero'
  infected[patientzero] <- TRUE
  
  if (plot.spread) {
    network.i <- graph.edgelist(links,directed=FALSE)
    fixlayout <- layout.kamada.kawai(network.i)  # store a fixed layout for the graph
    node.colour <- rep("SkyBlue2",N) # initialize node colours (SkyBlue2 is also the default node colour in igraph)
    node.colour[patientzero] <- "red" # infected nodes will be coloured red
    plot(network.i,layout=fixlayout, main="Time = 0", vertex.color=node.colour)
  }
  
  c=0 # counter
  while (length(infected[infected == F]) != 0) {
    discordant.links <- which(xor(infected[links[,1]],infected[links[,2]])) # find the indeces of links that connect an infected individual to an uninfected
    transmit <- rbinom(length(discordant.links),1,p.t) # determine randomly which of the discordant links transmit the disease
    # let me update the infection vector in three steps to make it easier to read:
    transmitter.links <- discordant.links[transmit==1]
    nodes.of.transmitter.links <- unique(as.vector(links[transmitter.links,1:2])) # gets both nodes of the transmitter links into a single vector; unique just filters out repetitions
    infected[nodes.of.transmitter.links] <- TRUE # here I simply set both nodes to TRUE (although the transmitter already had 'TRUE'). In more complex models, you might want to do a further check here and overwrite only the newly infected nodes.
    c = c + 1
    if (plot.spread) {
      node.colour[infected] <- "red"
      # readline() # waits for the user to press <ENTER> before proceeding; you need to switch to the console to do this
      plot(network.i,layout=fixlayout, main=paste("Time =", c), vertex.color=node.colour)
    }
    
  }
  c} # c = time steps

repN = 1000 # number of simulations
small = numeric(repN) # inital states
for (rep in 1:repN){small[rep] = time_infected(N=50, p.t=0.2, plot.spread = F, connectivity = 2)}

s2 =numeric(repN)
for (rep in 1:repN){s2[rep] = time_infected(N=100, p.t=0.2, plot.spread = F, connectivity = 2)}

s3 =numeric(repN)
for (rep in 1:repN){s3[rep] = time_infected(N=200, p.t=0.2, plot.spread = F, connectivity = 2)}

s4 =numeric(repN)
for (rep in 1:repN){s4[rep] = time_infected(N=300, p.t=0.2, plot.spread = F, connectivity = 2)}

s5 =numeric(repN)
for (rep in 1:repN){s5[rep] = time_infected(N=400, p.t=0.2, plot.spread = F, connectivity = 2)}

medium = numeric(repN)
for (rep in 1:repN){medium[rep] = time_infected(N=500, p.t=0.2, plot.spread = F, connectivity = 2)}

s6 =numeric(repN)
for (rep in 1:repN){s6[rep] = time_infected(N=1000, p.t=0.2, plot.spread = F, connectivity = 2)}

s7 =numeric(repN)
for (rep in 1:repN){s7[rep] = time_infected(N=2000, p.t=0.2, plot.spread = F, connectivity = 2)}

s8 =numeric(repN)
for (rep in 1:repN){s8[rep] = time_infected(N=3000, p.t=0.2, plot.spread = F, connectivity = 2)}

s9 =numeric(repN)
for (rep in 1:repN){s9[rep] = time_infected(N=4000, p.t=0.2, plot.spread = F, connectivity = 2)}

large = numeric(repN)
for (rep in 1:repN){large[rep] = time_infected(N=5000, p.t=0.2, plot.spread = F, connectivity = 2)}

min=c(min(small), min(s2), min(s3), min(s4), min(s5), min(medium), min(s6), min(s7), min(s8), min(s9), min(large))
median=c(median(small), median(s2), median(s3), median(s4), median(s5), median(medium), median(s6), median(s7), median(s8), median(s9), median(large))
sizepop=c(50,100,200,300,400,500,1000,2000,3000,4000,5000)
max=c(max(small), max(s2), max(s3), max(s4), max(s5), max(medium), max(s6), max(s7), max(s8), max(s9), max(large))
Ndf <- data.frame(sizepop, min, median, max)
Ndf


###incorrect predictions, keep trying

logsizepop=c(50,500,5000)
logmedian=c(median(small), median(medium), median(large))
logdf = data.frame(logsizepop, logmedian)
model <- lm(median ~ sizepop, data=Ndf)
pred.df <- data.frame(sizepop=c(50000,500000, 5000000))
prediction <- predict(model, newdata=pred.df)
prediction
predicted.sizepop = c(50,100,200,300,400,500,1000,2000,3000,4000,5000, 50000, 500000)
predicted.median = c(median(small), median(s2), median(s3), median(s4), median(s5), median(medium), median(s6), median(s7), median(s8), median(s9), median(large), 155.0961,  1353.8504)
pred.median2 = c(median(small), median(medium), median(large), )
predicteddf<-  data.frame (predicted.sizepop, predicted.median)
plot(predicteddf)  #probably incorrect

#plots
plot(Ndf, xlab= "Size of population", ylab="Median of time steps", main= "Time needed until everyone is infected")

plot(sizepop,median)
boxplot(small, s2,s3,s4,s5, medium,s6,s7,s8,s9, large,huge, main= "Time needed to be infected in different populations", xlab="Size of population", ylab="Time", names=c("50","100","200","300","400", "500","1000","2000","3000","4000", "5000", "10000"))

hist(small,main="N=50", xlab="Time steps", ylab="Frequency of infected nodes", col="blue")
hist(medium,main="N=500", xlab="Time steps", ylab="Frequency of infected nodes", col="orange", breaks=50)
hist(large,main="N=5000", xlab="Time steps", ylab="Frequency of infected nodes", col="red")

### changing connectivity in large populations
s9con =numeric(repN)
for (rep in 1:repN){s9con[rep] = time_infected(N=4000, p.t=0.2, plot.spread = F, connectivity = 3)}

largecon = numeric(repN)
for (rep in 1:repN){largecon[rep] = time_infected(N=5000, p.t=0.2, plot.spread = F, connectivity = 4)}

boxplot(s9, s9con, main= "Time needed to be infected in different populations (N=4000)", xlab="Connectivity", ylab="Time", names=c("2","3"))

boxplot(large, largecon, main= "Time needed to be infected in different populations (N=5000)", xlab="Connectivity", ylab="Time", names=c("2","4"))
