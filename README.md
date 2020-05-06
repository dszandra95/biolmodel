<center> <h1> Network models of epidemics </h1> </center>
<center> <h2> Computer modeling in biology course 2019/2020/2 </h2> </center>
https://tb.ethz.ch/education/learningmaterials/modelingcourse/level-2-modules/network.html
<br/><br/>

Luca_epimodel_agegroup:
  - Add metadata dataframe for each individual specifying age group and sex based on demography data
  - Set up network so that connections within age groups is more likely than outside age groups
  - Make children transmit the disease with a higher probability
  - Plots number of infected people at each time point (**would it be better to plot the newly infected at each timepoint?**) while changeing the percentage of age groups in the population
    - As the children transmit the disease with bigger probability I would expect the higher children percentage populations to have a steeper curve (also it would rise earlier on)
    - **However from the current figures I can only see the shift on x axis and not in steepness**
  <br/>
  
netepi_recovery (Benedek):
  - Recovery added to the basic epidemic network model: at each time step recovered individuals are chosen based on               recovery.time and rec.prob parameters. Recovery time indicates the required minimal elapsed time to recover. 
  - Network dynamics: swapping edges in each time step, keep node degree the same. 
  - Group info in each time step is stored in a dataframe (recovered, infected, susceptible groups).
  - Plots: epidemic line plot of change of groups (S, I, R) in time; boxplot of repeated epidemic simulation durations; histogram of the distribution of epidemic simulations
  <br/>
