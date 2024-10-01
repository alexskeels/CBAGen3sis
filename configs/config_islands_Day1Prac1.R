########################
### General settings ###
########################

random_seed =76985
start_time = 250
end_time = NA
max_number_of_species =50000
max_number_of_coexisting_species =20000
initial_abundance =  1

# set for first time
# assign("ss_eff", ss_eff_emp, envir = .GlobalEnv)
# a list of traits to include with each species, traits with ss_eff_ are hack traits to extract in cell processes
trait_names = c("dispersal", "mean_temp", "temp_width")

# ranges to scale the input environemts with:
# not listed variable:         no scaling takes place
# listed, set to NA:           the environmental variable will be scaled from [min, max] to [0, 1]
# listed with a given range r: the environmental variable will be scaled from [r1, r2] to [0, 1]

# environmental_ranges = list("mean_temp"=c(9,26), "min_temp"=c(9,26),  "max_temp"=c(9,26))

#########################
### Observer Function ###
#########################

end_of_timestep_observer = function(data, vars, config){
  plot_richness(data$all_species, data$landscape) # plot richness
  save_species() # saves a species and landscape objects for desired timesteps
}


######################
### Initialization ###
######################

create_ancestor_species <- function(landscape, config) {
  #### BROWSER ! ----------
  browser()
  co <- landscape$coordinates
  new_species <- create_species(rownames(co), config)
  new_species$traits[ , "dispersal"] <- 5 # denominator of exponential distribution
  new_species$traits[ , "mean_temp"] <- 20
  new_species$traits[ , "temp_width"] <- 2
  return(list(new_species))
}


#################
### Dispersal ###
#################

# returns n dispersal values (prob distrib function)
get_dispersal_values <- function(n, species, landscape, config){
  values <- rexp(n, rate = 1/species$traits[,"dispersal"])
  return(values)
}


##################
### Speciation ###
##################

# threshold for genetic distance after which a speciation event takes place
divergence_threshold=65 # between 10 and 50 ? as 0.1 to 0.5 Myrs or 100 - 500 kyrs

# adds a value of 1 to each geographic population cluster
get_divergence_factor <- function(species, cluster_indices, landscape, config) {
  return(1)
}

#######################
### Trait Evolution ###
#######################

apply_evolution <- function(species, cluster_indices, landscape, config) {
  # cell names
  trait_evolutionary_power <-0.01
  traits <- species[["traits"]]
  cells <- rownames(traits)
  #homogenize trait based on abundance
  # Homogenize all traits by weighted abundance, attention, exclude here any trait that should be more neutral
  trn <- config$gen3sis$general$trait_names
  for(cluster_index in unique(cluster_indices)){
    # cluster_index <- 1
    cells_cluster <- cells[which(cluster_indices == cluster_index)]
    # hist(traits[cells_cluster, "temp"], main="before")
    mean_abd <- mean(species$abundance[cells_cluster])
    weight_abd <- species$abundance[cells_cluster]/mean_abd
    for (ti in trn){
      traits[cells_cluster, ti] <- mean(traits[cells_cluster, ti]*weight_abd)
    }
    # hist(traits[cells_cluster, "temp"], main="after")
  }
  
  # mutate mean temperature
  mutation_deltas <-rnorm(length(traits[, ti]), mean=0, sd=trait_evolutionary_power)
  traits[, "mean_temp"] <- traits[, "mean_temp"] + mutation_deltas
  return(traits)
}

#################################################
### Environmental and Ecological Interactions ###
#################################################

apply_ecology <- function(abundance, traits, landscape, config) {
  # get the difference between species and site mean temp
  diff <- abs(traits[, "mean_temp"]-landscape[,"mean_temp"])
  # set the abundance of species with a difference in mean temp larger than the species temp width to zero and the ones below to one
  abundance[diff>traits[,"temp_width"]] <- 0
  abundance[diff<=traits[,"temp_width"]] <- 1
  return(abundance)
}
