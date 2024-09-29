########################
### General settings ###
########################

random_seed = 663
start_time = 50
end_time = NA
max_number_of_species =50000
max_number_of_coexisting_species =20000
initial_abundance =  1

# set for first time
# assign("ss_eff", ss_eff_emp, envir = .GlobalEnv)
# a list of traits to include with each species, traits with ss_eff_ are hack traits to extract in cell processes
trait_names = c("dispersal", "temp_niche_centre", "temp_niche_width", "start_island")

# ranges to scale the input environemts with:
# not listed variable:         no scaling takes place
# listed, set to NA:           the environmental variable will be scaled from [min, max] to [0, 1]
# listed with a given range r: the environmental variable will be scaled from [r1, r2] to [0, 1]

# environmental_ranges = list("mean_temp"=c(9,26), "min_temp"=c(9,26),  "max_temp"=c(9,26))

#########################
### Observer Function ###
#########################

end_of_timestep_observer = function(data, vars, config){  
  
  plot_richness(data$all_species, data$landscape)
  
  save_species() # saves a species and landscape objects for desired timesteps
  
}


######################
### Initialization ###
######################

create_ancestor_species <- function(landscape, config) {
  
  co <- landscape$coordinates
  
  # however, we made thinks simpler to you.
  # get patches vector
  pv <- landscape$environment[, "patch"]

  new_species <- list()
  
  manual_traits <- list(
    "dispersal" = c(10, 7.5, 5, 2.5),
    "temp_niche_width" = c(0.5, 1, 2, 4)
  )
  
  for (sp_i in 1:4){
    # create a species empty object
    new_species[[sp_i]] <- create_species(names(pv[pv==sp_i]), config)
    # set manual dispersal and niche width
    new_species[[sp_i]]$traits[ , "dispersal"] <- manual_traits$"dispersal"[sp_i]
    new_species[[sp_i]]$traits[ , "temp_niche_width"] <- manual_traits$"temp_niche_width"[sp_i]
    # set species 1 start to island 1, species 2 to island 2, etc...
    new_species[[sp_i]]$traits[ , "start_island"] <- unique(landscape$environment[pv==sp_i, "patch"]) # this is just a sanity check
    # set species mean temp to the mean temp of the patches
    new_species[[sp_i]]$traits[ , "temp_niche_centre"] <- mean(landscape$environment[pv==sp_i, "mean_temp"])
  }
  
  return(new_species)
}



#################
### Dispersal ###
#################

# returns n dispersal values (prob distrib function)
get_dispersal_values <- function(n, species, landscape, config){
  #browser()
  values <- rexp(n, rate = 1/species$traits[,"dispersal"])
  return(values)
}


##################
### Speciation ###
##################

# threshold for genetic distance after which a speciation event takes place
divergence_threshold=5 # between 10 and 50 ? as 0.1 to 0.5 Myrs or 100 - 500 kyrs

# adds a value of 1 to each geographic population cluster
get_divergence_factor <- function(species, cluster_indices, landscape, config) {
  return(1)
}

#######################
### Trait Evolution ###
#######################

apply_evolution <- function(species, cluster_indices, landscape, config) {
  # cell names
  trait_evolutionary_power <-0.05
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
  traits[, "temp_niche_centre"] <- traits[, "temp_niche_centre"] + mutation_deltas
  return(traits)
}
#################################################
### Environmental and Ecological Interactions ###
#################################################


inflection <- 0.2 # extinction parameter: inflexion point of extinction probability curve
decay <- 100  # extinction parameter: decay rate of extinction probability curve
K_max <- 2

apply_ecology <- function(abundance, traits, landscape, config) {
  
  temp_omega <- traits[, "temp_niche_width"]
  
  # define the traits
  temp_opt <- traits[, "temp_niche_centre"]
  temp_site <- landscape[, "mean_temp"]

  # define K with environment
  K_site <- K_max 
  K_c <- 1
  # determine abundance
  abundance <- K_c * exp(-((temp_opt-temp_site)^2/(2*(temp_omega^2)))) 
  # determine whether species go extinct
  prob_extinction <- (1/(1+exp(-decay*(inflection - abundance))))
  abundance[which(sapply(prob_extinction, FUN=function(x){sample(c(0,1),1, prob= c(x, 1-(x)))})==0)] <- 0
  abundance_sum <- sum(abundance)
  
  # if the total abundance in a grid cell is above K, then sequentially reduce the abudnace of species until it isn't anymore
  while_i <- 0.01
  while_sd  <- 0.2
  
  
  while(abundance_sum > K_site){
    
    #number of species to reduce # can go slow when doing one at time if lots of species in the cell
    # select 10% of species
    n_species_to_reduce <- pmax(1, round(0.1*length(which(abundance>0))))
    
    reduce_prob <- rep(1, length(abundance))
    reduce_prob[which(abundance==0)] <- 0
    
    # choose a species inversely proportional to how well adapated they are
    species_to_reduce <- sample(1:length(abundance), n_species_to_reduce, prob=reduce_prob)
    
    # reduce abundance by a value drawn from a normal distribution with mean 0
    abundance[species_to_reduce] <- abundance[species_to_reduce] - abs(rnorm(1, 0, sd=while_sd+while_i))
    
    # recalculate extinction probabilities
    prob_extinction <- (1/(1+exp(-decay*(inflection - abundance))))
    abundance[which(sapply(prob_extinction, FUN=function(x){sample(c(0,1),1, prob= c(x, 1-(x)))})==0)] <- 0
    
    # condition to stop the while loop and reduce all species
    if(while_i >= 1){
      #print("breaking while loop")
      abundance[order(abundance, decreasing=T)[which(cumsum(abundance[order(abundance, decreasing=T)]) > K_site)]] <- 0
    }
    # sum again abundances
    abundance_sum <- sum(abundance, na.rm=T)
    while_i <- while_i + 0.01
  }
  
  return(abundance)
  
  
  
  
}
