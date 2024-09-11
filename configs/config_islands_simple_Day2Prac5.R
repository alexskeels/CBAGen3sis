########################
### General settings ###
########################

random_seed = 357
start_time = NA
end_time = NA
max_number_of_species =50000
max_number_of_coexisting_species =20000
initial_abundance =  1

# set for first time
# assign("ss_eff", ss_eff_emp, envir = .GlobalEnv)
# a list of traits to include with each species, traits with ss_eff_ are hack traits to extract in cell processes
trait_names = c("dispersal", "temp_niche_centre", "temp_width", "start_island")


### HI, I'm a custom function of a gaussian
# ecological environmental function
fg <- function(x,a,b,c){
  a <- a/c
  v <- a*exp(-((x-b)^2/(2*c^2)))
  return(v)
}
# plot(fg(20:30, 10, 25, 5), type="l")
#########################
### Observer Function ###
#########################

end_of_timestep_observer = function(data, vars, config){
  # browser()
  save_phylogeny()
  par(mfrow=c(1,2))
  plot_richness(data$all_species, data$landscape)
  plot(read.nexus(file.path(config$directories$output, "phylogeny", paste0("phylogeny_t_", vars$ti, 
                                                                         ".nex"))))
  save_species() # saves a species and landscape objects for desired timesteps
}


######################
### Initialization ###
######################

create_ancestor_species <- function(landscape, config) {
  # browser()
  co <- landscape$coordinates
  
  # If we wouldn't have passe the patches though the environment, we could have done this:
  # sp1 <- co[which(co[,1]<20&co[,2]<30),]
  # sp2 <- co[which(co[,1]>20&co[,2]<30),]
  # sp3 <- co[which(co[,1]<20&co[,2]>30),]
  # sp4 <- co[which(co[,1]>20&co[,2]>30),]
  
  # however, we made thinks simpler to you.
  # get patches vector
  pv <- landscape$environment[, "island_id"]
  
  # species_coords <- list(sp1,
  #                        sp2,
  #                        sp3,
  #                        sp4)
  new_species <- list()

  
  n_ids <- length(unique(pv))
  
  manual_traits <- list(
    "dispersal" = rep(5, n_ids),
    "temp_width" = rep(3, n_ids),
    "temp_niche_centre" = rep(27, n_ids)
  )
  
  for (sp_i in 1:n_ids){
    # sp_i <- 1
    # create a species empty object
    new_species[[sp_i]] <- create_species(names(pv[pv==sp_i]), config)
    # set manual dispersal and niche width
    new_species[[sp_i]]$traits[ , "dispersal"] <- manual_traits$"dispersal"[sp_i]
    new_species[[sp_i]]$traits[ , "temp_width"] <- manual_traits$"temp_width"[sp_i]
    # set species 1 start to island 1, species 2 to island 2, etc...
    new_species[[sp_i]]$traits[ , "start_island"] <- unique(landscape$environment[pv==sp_i, "island_id"]) # this is just a sanity check
    # set species mean temp to the mean temp of the patches
    new_species[[sp_i]]$traits[ , "temp_niche_centre"] <- manual_traits$"temp_niche_centre"[sp_i]
  }
  
  return(new_species)
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
divergence_threshold=10 # between 10 and 50 ? as 0.1 to 0.5 Myrs or 100 - 500 kyrs

# adds a value of 1 to each geographic population cluster
get_divergence_factor <- function(species, cluster_indices, landscape, config) {
  return(1)
}

#######################
### Trait Evolution ###
#######################

apply_evolution <- function(species, cluster_indices, landscape, config) {
  # browser()
  # cell names
  trait_evolutionary_power <-0.01
  traits <- species[["traits"]]
  cells <- rownames(traits)
  #homogenize trait based on abundance
  # Homogenize all traits by weighted abundance, attention, exclude here any trait that should be more neutral
  trn <- config$gen3sis$general$trait_names # get trait names
  trn <- trn[!trn %in% c("start_island")] # exclude start_island
  for(cluster_index in unique(cluster_indices)){
    # cluster_index <- 1
    cells_cluster <- cells[which(cluster_indices == cluster_index)]
    # hist(traits[cells_cluster, "temp"], main="before")
    mean_abd <- mean(species$abundance[cells_cluster])
    weight_abd <- species$abundance[cells_cluster]/mean_abd
    for (ti in trn){
      traits[cells_cluster, ti] <- mean(traits[cells_cluster, ti]*weight_abd)
    }
    # hist(traits[cells_cluster, "temp"], main="after") # usefull to see spread
  }
  
  # mutate mean temperature
  mutation_deltas <-rnorm(length(traits[, ti]), mean=0, sd=trait_evolutionary_power)
  traits[, "temp_niche_centre"] <- traits[, "temp_niche_centre"] + mutation_deltas
  return(traits)
}

#################################################
### Environmental and Ecological Interactions ###
#################################################

apply_ecology <- function(abundance, traits, landscape, config) {
  # browser()
  
  #old# get the difference between species and site mean temp
  #old# diff <- abs(traits[, "temp_niche_centre"]-landscape[,"temp"])
  
  diff <- fg(landscape[,"temp"], a=10, b=traits[,"temp_niche_centre"], c=traits[,"temp_width"])
  
  # first set abundances that have diff bellow 1 to zero
  abundance[diff<1] <- 0
  
  diff[diff<1] <- 1
  
  # local extinction probability based on env. fitness
  p_ext <- 1/diff
  
  # generate a vector of 0s and 1s based on the probabilities
  # 0 indicates the value will be set to zero, 1 indicates the value remains unchanged
  set_to_zero <- rbinom(length(abundance), size = 1, prob = 1 - p_ext)
  
  # Set abundances to zero where set_to_zero is 0
  abundance <- abundance * set_to_zero
  return(abundance)
}
