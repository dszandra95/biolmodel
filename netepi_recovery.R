# 4/23/20
# Benedek Dankó

library(igraph)
library(ggplot2)
library(ggpubr)


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


################

# MAIN FUNCTION:

# network epidemics model with recovery
# N = number of nodes, p.t. = probability of transmission, recovery.time = minimal time needed for recovery, rec.prob = probability of recovery (for rbinom function)

pandemic.simulation <- function(N, p.t, connectivity, recovery.time, rec.prob){

  
infected <- numeric(N) # initialize infection status, 0: susceptible, 1: infected, 2: recovered
patientzero <- sample(N,1) # select 'patient zero'
infected[patientzero] <- 1
links <- generate.network.B(N,connectivity)

populations.df <- data.frame(time=integer(), susceptible=integer(), infected=integer(), recovered=integer()) # stores the group numbers per time step

infection.df <- data.frame(node=c(1:N), time=rep(Inf,N)) # to store the infection time of all nodes
infection.df[patientzero,2] <- 0
previous.inf.nodes <- c()

time=0 # time step

while (length(infected[infected == 1]) > 0) { # run simulation until all nodes become recovered
  inf <- as.logical(infected)
  discordant.links <- which(xor(inf[links[,1]],inf[links[,2]])) # find the indeces of links that connect an infected individual to an uninfected
  transmit <- rbinom(length(discordant.links),1,p.t) # determine randomly which of the discordant links transmit the disease
  transmitter.links <- discordant.links[transmit==1]
  
  total.infected.nodes <- which(infected %in% 1) # get all infected nodes
  recovered.nodes <- which(infected %in% 2) # all recovered nodes
  
  nodes.of.transmitter.links <- unique(as.vector(links[transmitter.links,1:2])) # gets both nodes of the transmitter links into a single vector; unique just filters out repetitions
  nodes.of.transmitter.links <- nodes.of.transmitter.links[!nodes.of.transmitter.links %in% recovered.nodes] # to avoid infecting recovered nodes
  
  if (length(as.vector(nodes.of.transmitter.links)) >= 0) {new.inf.nodes <- setdiff(total.infected.nodes, previous.inf.nodes) # get new infected nodes
  previous.inf.nodes <- total.infected.nodes
  
  if (length(new.inf.nodes) > 0){infection.df[new.inf.nodes, 2] <- time}} # store the infection time of new infected nodes
  infected[nodes.of.transmitter.links] <- 1 # here I simply set both nodes to TRUE (although the transmitter already had 'TRUE'). In more complex models, you might want to do a further check here and overwrite only the newly infected nodes.
  max.recover <- time-recovery.time # minimal time after infection, that is needed to recover (minimal length of the disease)
  potential.recovers <- as.vector(infection.df[infection.df$time <= max.recover,]$node)
  if (max.recover >= 0){
    recovered.nodes <- rbinom(length(potential.recovers), 1, rec.prob)
    recovered.nodes <- potential.recovers[recovered.nodes == 1]
    infected[recovered.nodes] <- 2} # recovered nodes
  
  time = time + 1
  current.df <- data.frame(time, length(infected[infected == 0]), length(infected[infected == 1]), length(infected[infected == 2]))
  names(current.df) <- c("time", "susceptible", "infected", "recovered")  
  populations.df <- rbind(populations.df, current.df)
}


return (populations.df) # susceptible, infected, recovered states in time
}

#test:
#pandemic.simulation(100, 0.2, 3, 3, 0.4)

#########################################

### Plot simulations in 2x2 grid

simulations <- list()
params <- c(0.2, 0.4, 0.6, 0.8) # p.t. parameters


for (i in 1:4) {
  sim <- pandemic.simulation(1000, params[i], 4, 3, 0.4)
  simulations[[i]] <- sim  # create and add new data frame
  }

plots <- list() # store plots

for (i in (1:length(simulations))){
df <- simulations[[i]]
max.time <- df$time[df$infected == max(df$infected)][1]
p <-ggplot(data=df) +
  geom_line(mapping=aes(y=susceptible,x=time,color="Susceptible"),size=1) +
  geom_line(mapping=aes(y=infected,x=time,color="Infected"),size=1) +
  geom_line(mapping=aes(y=recovered,x=time,color="Recovered"),size=1) +
  scale_color_manual(values = c(
    'Susceptible' = 'blue',
    'Infected' = 'red',
    'Recovered' = "green")) +
  labs(color = 'Groups') +
  ggtitle(paste("Network epidemic simulation", as.character(i))) +
  theme(plot.title = element_text(hjust = 0.5)) + labs(x = "Time step", y="Number of individuals") +
  geom_vline(xintercept = max.time, linetype="dotted", 
                color = "red", size=0.7) +
  scale_x_continuous(n.breaks=max(df$time))
plots[[i]] <- p
}


labels <- sapply('p.t.', paste, params, sep=" ")
grid_plot <- ggarrange(plots[[1]],plots[[2]], plots[[3]], plots[[4]], ncol=2, nrow=2, common.legend = TRUE, legend="bottom", labels=labels)
grid_plot
