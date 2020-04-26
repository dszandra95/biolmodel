# 4/26/20
# Benedek Dankó


library(igraph)
library(ggplot2)
library(ggpubr)
library(reshape2)


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


##### MAIN FUNCTION #####

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
  return (time)
  }
  

#test:
#pandemic.simulation(50, 0.2, 3, 0, 0.4)


#### Simulations #####

params <- c(0, 2, 4, 6, 8, 10) # different recovery time parameters
repetitions <- 50 # repeat each simulation 50 times

df <- data.frame('1'=numeric(repetitions), '2'=numeric(repetitions), '3'=numeric(repetitions), # initialize empty data frame
                 '4'=numeric(repetitions), '5'=numeric(repetitions), '6'=numeric(repetitions))


for (i in (1:length(params))){ # run simulations and fill up data frame
values <- numeric(repetitions)
for (rep in (1:repetitions)){
values[rep] <- pandemic.simulation(50, 0.2, 3, params[i], 0.4)
}
df[paste0("X", i)] <- values
}

# transform dataframe for boxplots:
s_df <- stack(df)
s_df <- transform(s_df, group = substring(ind, 1, 1),
                             obs = substring(ind, 2))

# plot the result:
plt <- ggplot(s_df, aes(x = ind, y = values)) +
  stat_boxplot(geom ='errorbar', width = 0.6) +
  geom_boxplot() +
  stat_summary(fun = mean, geom = "errorbar", aes(ymax = ..y.., ymin = ..y..),
               width = .75, linetype = "dashed") + 
  ggtitle("Time steps needed to become all nodes\n recovered after infection") +
  xlab("Min. recovery time") + ylab("Time step") +
  scale_x_discrete(labels=c("X1" = "0", "X2" = "2",
                            "X3" = "4", "X4" = "6", 
                            "X5" = "8", "X6"="10")) +
  theme(axis.text=element_text(size=15), axis.title=element_text(size=15), 
        plot.title = element_text(size=16, hjust = 0.5))
plt # boxplots with increasing recovery time parameter

# To save plot:
#ggsave("/path/boxplot.png")
#plt
#dev.off()
