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

fun_1 <- function(x){x+1}
create_ancestor_species <- function(landscape, config) {
  #browser()
  co <- landscape$coordinates
  
  # If we wouldn't have passe the patches though the environment, we could have done this:
  # sp1 <- co[which(co[,1]<20&co[,2]<30),]
  # sp2 <- co[which(co[,1]>20&co[,2]<30),]
  # sp3 <- co[which(co[,1]<20&co[,2]>30),]
  # sp4 <- co[which(co[,1]>20&co[,2]>30),]
  
  # however, we made thinks simpler to you.
  # get patches vector
  pv <- landscape$environment[, "patch"]
  
  # species_coords <- list(sp1,
  #                        sp2,
  #                        sp3,
  #                        sp4)
  new_species <- list()
  

  manual_traits <- list(
    "dispersal" = c(5, 5, 5, 5),
    "temp_niche_width" = c(1, 1, 1, 1)
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
  # get the difference between species and site mean temp
  diff <- abs(traits[, "temp_niche_centre"]-landscape[,"mean_temp"])
  # set the abundance of species with a difference in mean temp larger than the species temp width to zero and the ones below to one
  abundance[diff>traits[,"temp_niche_width"]] <- 0
  abundance[diff<=traits[,"temp_niche_width"]] <- 1
  return(abundance)
}
