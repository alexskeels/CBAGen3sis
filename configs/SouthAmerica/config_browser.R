########################
### General settings ###
########################
### browser() #first,browse#
### NOTE THAT THIS SCRIPT IS DETERMINISTIC!



# set the random seed for the simulation
random_seed = NA

# set the starting time step or leave NA to use the earliest/highest timestep
start_time = NA

# set the end time step or leave as NA to use the lates/lowest timestep (0)
end_time = 59

# maximum total number of species in the simulation before it is aborted
max_number_of_species = 10000

# maximum number of species within one cell before the simulation is aborted
max_number_of_coexisting_species = 1000

# a list of traits to include with each species
# a "dispersion" trait is implictly added in any case
#trait_names = c("t_min", "a_min", "competition", "dispersion")
trait_names = c("temp", "niche_width", "dispersal") # "prec",

# ranges to scale the input environemts with:
# not listed variable:         no scaling takes place
# listed, set to NA:           the environmental variable will be scaled from [min, max] to [0, 1]
# lsited with a given range r: the environmental variable will be scaled from [r1, r2] to [0, 1]

# environmental_ranges = list("temp" = c(9, 31))

# a place to inspect the internal state of the simulation and collect additional information if desired
end_of_timestep_observer = function(data, vars, config){
  oldpar <- par(no.readonly = TRUE)
  on.exit(par(oldpar))
  #browser() #first,browse#
  par(mfrow=c(2,2))
  plot_richness(data$all_species, data$landscape)
  plot_species_presence(data$all_species[[1]], data$landscape)
  plot_species_presence(data$all_species[[2]], data$landscape)
  plot_species_presence(data$all_species[[3]], data$landscape)
  
  #save
  # save_species()
  # save_landscape()
}



######################
### Initialization ###
######################
# the intial abundace of a newly colonized cell, both during setup and later when colonizing a cell during the dispersal
initial_abundance = 1
# place species within rectangle, our case entire globe
create_ancestor_species <- function(landscape, config) {
  # browser() #first,browse# 
  range <- c(-180, 180, -90, 90)
  co <- landscape$coordinates
  selection <- co[, "x"] >= range[1] &
    co[, "x"] <= range[2] &
    co[, "y"] >= range[3] &
    co[, "y"] <= range[4]
  # Hard codded initial conditions for 3 new species
  opt_temps <- c(22, 25, 28) # hist(landscape$environment[,"temp"])
  nich_widths <- c(0.5, 0.75, 1) # trait trade-off based on competitive interactions assumptions
  disps <- c(200, 100, 60) # traits trade-off based on disp v.s. stress tolerance v.s. competitive ability
  # create 3 new species
  new_species <- list()
  for(i in 1:3){
    initial_cells <- rownames(co)[selection]
    #initial_cells <- sample(initial_cells, 1)
    new_species[[i]] <- create_species(initial_cells, config)
    #set local adaptation to max optimal temp equals local temp
    new_species[[i]]$traits[ , "temp"] <- opt_temps[i]
    new_species[[i]]$traits[, "niche_width"] <- nich_widths[i]
    new_species[[i]]$traits[ , "dispersal"] <- disps[i] 
    #plot_species_presence(landscape, species=new_species[[i]])
  }
  return(new_species)
}


#################
### Dispersal ###
#################
# returns n dispersal values
get_dispersal_values <- function(n, species, landscape, config) {
  browser() #first,browse# 
  # if(landscape$timestep==stop_time){browser()} #second,browse# 
  # hint, look at #78 dispersal.R num_draws <- length(free_cells) * length(presence_spi_ti). Start looking for: get_dispersal_values With Ctrl+Shift+F
  # at the species level... e.g. plot_species_presence(species, landscape)
  mean_trait_disp <- mean(species$traits[, "dispersal"])
  values <- rep(mean_trait_disp, n)
  # values <- rexp(n, rate=1/11) ### VARY
  return(values)
}


##################
### Speciation ###
##################
# threshold for genetic distance after which a speciation event takes place
divergence_threshold=4

# factor by which the divergence is increased between geographicaly isolated population
# can also be a matrix between the different population clusters
get_divergence_factor <- function(species, cluster_indices, landscape, config) {
  return(1)
}


################
### Mutation ###
################
# mutate the traits of a species and return the new traits matrix
apply_evolution <- function(species, cluster_indices, landscape, config) {
  # browser() #first,browse# 
  # trait_evolutionary_power <-0.005 ### VARY
  traits <- species[["traits"]] # no evolution
  # option: homogenize traits in case they change:
  # cells <- rownames(traits)
  # #homogenize trait based on abundance
  # for(cluster_index in unique(cluster_indices)){
  #   # cluster_index <- 1
  #   cells_cluster <- cells[which(cluster_indices == cluster_index)]
  #   # hist(traits[cells_cluster, "temp"], main="before")
  #   mean_abd <- mean(species$abundance[cells_cluster])
  #   weight_abd <- species$abundance[cells_cluster]/mean_abd
  #   traits[cells_cluster, "temp"] <- mean(traits[cells_cluster, "temp"]*weight_abd)
  #   # hist(traits[cells_cluster, "temp"], main="after")
  # }
  # option: add random mutations:
  # mutation_deltas <-rnorm(length(traits[, "temp"]), mean=0, sd=trait_evolutionary_power)
  # traits[, "temp"] <- traits[, "temp"] + mutation_deltas
  return(traits)
}


###############
### Ecology ###
###############
# called for every cell with all occuring species, this functin calculates the who survives in the current cells
# returns a vector of abundances
# set the abundance to 0 for every species supposed to die
apply_ecology <- function(abundance, traits, landscape, config, abundance_scale = 10, abundance_threshold = 1) {
  # browser() #first,browse# 
  
  # for arbitrary real constants a, b and non zero c. 
  # It is named after the mathematician Carl Friedrich Gauss. 
  # The graph of a Gaussian is a characteristic symmetric "bell curve" shape. 
  # The parameter a is the height of the curve's peak, b is the position of the center of the peak 
  # and c (the standard deviation, sometimes called the Gaussian RMS width) controls the width of the "bell". 
  #
  # fg <- function(x,a,b,c){
  #   v <- a*exp(-((x-b)^2/(2*c^2)))
  #   return(v)
  # }
  # this functions scales the abundance
  # fgs <- function(x,a,b,c){
  #   v <- (a/c)*exp(-((x-b)^2/(2*c^2)))
  #   return(v)
  # 
  # plot(fg(x=20, a=10, b=seq(10,31,1), c=1), type='l', xaxt="n") # c ranges from 0.001 to 0.3 (very wide niche)
  # axis(1, at=1:22, labels =seq(10,31,1) )
  # abline(h=1)
  # abline(v=8.3)
  # gaussian
  abundance <- ((abundance_scale)*exp(-((traits[, "temp"] - landscape[, "temp"])**2/(2*traits[, "niche_width"]**2)))) ### VARRY niche breath
  
  ############# INTERACT ################
  # To understand this... play with the ecology function and see what is happening, for example...
  # within the temperature gradient:
  # abundance_scale <- 1
  # traits <- seq(0,1,0.01)
  # landscape_temp <- 0.5
  # niche_width <- 0.05
  # plot((abundance_scale /niche_width)*exp(-((traits - landscape_temp)**2/(2*niche_width**2))), type='l')
  # abline(h=1, col="red")
  # text(10,1,"all bellow die =(...", adj=0)
  
  #abundance thhreashold
  abundance[abundance<abundance_threshold] <- 0
  return(abundance)
}
