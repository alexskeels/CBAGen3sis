########################
### General settings ###
########################

random_seed = 666
start_time = 50
end_time = NA
max_number_of_species =50000
max_number_of_coexisting_species =20000
initial_abundance =  1

# set for first time
# assign("ss_eff", ss_eff_emp, envir = .GlobalEnv)
# a list of traits to include with each species, traits with ss_eff_ are hack traits to extract in cell processes
trait_names = c("dispersal", "temp_niche_centre", "temp_niche_width", "start_island", "body_size", "metabolic_rate", "abundance", "mass_specific_metabolic_rate")

# ranges to scale the input environemts with:
# not listed variable:         no scaling takes place
# listed, set to NA:           the environmental variable will be scaled from [min, max] to [0, 1]
# listed with a given range r: the environmental variable will be scaled from [r1, r2] to [0, 1]

# environmental_ranges = list("mean_temp"=c(9,26), "min_temp"=c(9,26),  "max_temp"=c(9,26))

#########################
### Observer Function ###
#########################

end_of_timestep_observer = function(data, vars, config){  
  
  #plot_richness(data$all_species, data$landscape)
  
  save_species() # saves a species and landscape objects for desired timesteps
  
}


######################
### Initialization ###
######################
E <- 0.65 # activation energy (eV)
k <- 8.62 * (10^-5) # boltzmann constant

create_ancestor_species <- function(landscape, config) {
 #browser()
  co <- landscape$coordinates
  sp1 <- co[which(co[,1]>20&co[,2]>30),]
  sp2 <- co[which(co[,1]<20&co[,2]<30),]
  sp3 <- co[which(co[,1]>20&co[,2]<30),]
  sp4 <- co[which(co[,1]<20&co[,2]>30),]
  
  species_coords <- list(sp1,
                         sp2,
                         sp3,
                         sp4)
  new_species <- list()
  
  # Species 1 - high dispersal, low niche breadth
  new_species[[1]] <- create_species(rownames(species_coords[[1]]), config)
  new_species[[1]]$traits[ , "dispersal"] <- 10
  new_species[[1]]$traits[ , "temp_niche_centre"] <- landscape$environment[rownames(species_coords[[1]]), "mean_temp"]
  new_species[[1]]$traits[ , "temp_niche_width"] <- 0.5
  new_species[[1]]$traits[ , "start_island"] <- unique(landscape$environment[rownames(species_coords[[1]]), "patch"])
  new_species[[1]]$traits[ , "body_size"] <- 1
  new_species[[1]]$traits[ , "abundance"] <- 1
  new_species[[1]]$traits[ , "metabolic_rate"] <- (new_species[[1]]$traits[ , "body_size"]^(3/4)) *(exp(-E/(k*mean((new_species[[1]]$traits[ , "temp_niche_centre"]+273.15)))))
  new_species[[1]]$traits[ , "mass_specific_metabolic_rate"] <- new_species[[1]]$traits[ , "metabolic_rate"]/new_species[[1]]$traits[ , "body_size"]
  
  # Species 2 - medium-high dispersal, medium-low niche breadth
  new_species[[2]] <- create_species(rownames(species_coords[[2]]), config)
  new_species[[2]]$traits[ , "dispersal"] <- 7.5
  new_species[[2]]$traits[ , "temp_niche_centre"] <- landscape$environment[rownames(species_coords[[2]]), "mean_temp"]
  new_species[[2]]$traits[ , "temp_niche_width"] <- 1
  new_species[[2]]$traits[ , "start_island"] <- unique(landscape$environment[rownames(species_coords[[2]]), "patch"])
  new_species[[2]]$traits[ , "body_size"] <- 1
  new_species[[2]]$traits[ , "abundance"] <- 1
  new_species[[2]]$traits[ , "metabolic_rate"] <- (new_species[[2]]$traits[ , "body_size"]^(3/4)) *(exp(-E/(k*mean((new_species[[2]]$traits[ , "temp_niche_centre"]+273.15)))))
  new_species[[2]]$traits[ , "mass_specific_metabolic_rate"] <- new_species[[2]]$traits[ , "metabolic_rate"]/new_species[[2]]$traits[ , "body_size"]
  
  # Species 3 - medium-low dispersal, medium-high niche breadth
  new_species[[3]] <- create_species(rownames(species_coords[[3]]), config)
  new_species[[3]]$traits[ , "dispersal"] <- 5
  new_species[[3]]$traits[ , "temp_niche_centre"] <- landscape$environment[rownames(species_coords[[3]]), "mean_temp"]
  new_species[[3]]$traits[ , "temp_niche_width"] <- 2
  new_species[[3]]$traits[ , "start_island"] <- unique(landscape$environment[rownames(species_coords[[3]]), "patch"])
  new_species[[3]]$traits[ , "body_size"] <- 1
  new_species[[3]]$traits[ , "abundance"] <- 1
  new_species[[3]]$traits[ , "metabolic_rate"] <- (new_species[[3]]$traits[ , "body_size"]^(3/4)) *(exp(-E/(k*mean((new_species[[3]]$traits[ , "temp_niche_centre"]+273.15)))))
  new_species[[3]]$traits[ , "mass_specific_metabolic_rate"] <- new_species[[3]]$traits[ , "metabolic_rate"]/new_species[[3]]$traits[ , "body_size"]
  
  # Species 4 -low dispersal,high niche breadth
  new_species[[4]] <- create_species(rownames(species_coords[[4]]), config)
  new_species[[4]]$traits[ , "dispersal"] <- 2.5
  new_species[[4]]$traits[ , "temp_niche_centre"] <- landscape$environment[rownames(species_coords[[4]]), "mean_temp"]
  new_species[[4]]$traits[ , "temp_niche_width"] <- 4
  new_species[[4]]$traits[ , "start_island"] <- unique(landscape$environment[rownames(species_coords[[4]]), "patch"])
  new_species[[4]]$traits[ , "body_size"] <- 1
  new_species[[4]]$traits[ , "abundance"] <- 1
  new_species[[4]]$traits[ , "metabolic_rate"] <- (new_species[[4]]$traits[ , "body_size"]^(3/4)) *(exp(-E/(k*mean((new_species[[4]]$traits[ , "temp_niche_centre"]+273.15)))))
  new_species[[4]]$traits[ , "mass_specific_metabolic_rate"] <- new_species[[4]]$traits[ , "metabolic_rate"]/new_species[[4]]$traits[ , "body_size"]
  
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
  #browser()
  trait_evolutionary_power_tempniche <-0.05
  trait_evolutionary_power_bodysize <-0.1
  
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
  traits[, "temp_niche_centre"] <- traits[, "temp_niche_centre"] + rnorm(length(traits[, ti]), mean=0, sd=trait_evolutionary_power_tempniche)
  
  # mutate body size
  traits[, "body_size"] <- traits[, "body_size"] + rnorm(length(traits[, ti]), mean=0, sd=trait_evolutionary_power_bodysize)
  
  # calculate metabolic rate
  B <- (traits[, "body_size"]^(3/4)) *(exp(-E/(k*(landscape$environment[rownames(traits),"mean_temp"]+273.15))))
  traits[, "metabolic_rate"] <- B
  traits[, "mass_specific_metabolic_rate"] <- B/traits[, "body_size"]
  
  return(traits)
}
#################################################
### Environmental and Ecological Interactions ###
#################################################

# M0 - No Ecological Limits - Null model
#apply_ecology <- function(abundance, traits, landscape, config) {
#  # get the difference between species and site mean temp
#  diff <- abs(traits[, "temp_niche_centre"]-landscape[,"mean_temp"])
#  # set the abundance of species with a difference in mean temp larger than the species temp width to zero 
#  # and the ones below to one
#  abundance[diff>traits[,"temp_niche_width"]] <- 0
#  abundance[diff<=traits[,"temp_niche_width"]] <- 1
#  return(abundance)
#}
#
## function to rescale temp so no negatives
#range01 <- function(x,low, high, ...){(x - low) / (high - low)}

# M1 - Hard Ecological Limit on Species - Random Extirpation
#apply_ecology <- function(abundance, traits, landscape, config) {
#  
#  #browser()
#  
#  # Set Species level carrying capacity based on temp and precipitation
#  K_opt <-  range01(landscape[, "mean_temp"], 6, 26)  * 100
#  
#  # get the difference between species and site mean temp
#  diff <- abs(traits[, "temp_niche_centre"]-landscape[,"mean_temp"])
#  
#  # set the abundance of species with a difference in mean temp larger than the species temp width to zero and the ones below to one
#  abundance[diff>traits[,"temp_niche_width"]] <- 0
#  abundance[diff<=traits[,"temp_niche_width"]] <- 1
#  
#  # if the number of species is higher than the carrying capacity, randomly select species to extirpate
#  if( sum(abundance) > K_opt){
#    k_diff <- k_opt - sum(abundance)
#    abundance[which(abundance==1)][sample(1:sum(abundance), k_diff, replace=F)] <- 0
#  }
#  
#  return(abundance)
#}

# M2 - More-Individuals Hypothesis
#apply_ecology <- function(abundance, traits, landscape, config) {
#  
#  # presence / absence is dependant on the niche width
#  #browser()
#  x0 <- 1000 # inflexion point of extinction probability curve
#  omega <- 0.02 # environmental filtering (higher values equal less restrictive niche)
#  decay <- 0.001
#
#  # Set Individual-level carrying capacity based on temp and precipitation
#  K_opt <-  range01(landscape[, "mean_temp"], 6, 26) * 10000
#  
#  # difference between landscape and species optimum temperature niche
#  diff_t <- abs(traits[, "temp_niche_centre"] - landscape[, "mean_temp"])
#  
#  # potential population size based on resource use efficieny 
#  # omega is strength of environmental filtering
#  # will equal K_opt when tthe species is perfectly adapted
#  # From McPeek 2008/2007
#  Nij <- K_opt * exp(-(diff_t/omega)^2)
#  
#  # realised abundance given a zero sum constraint
#  # if the RUE sum to less than 1 then each species gets that proportion of the NPP
#  # but if the RUE sums to more than one it needs to be rescaled so each species gets a proportion of the NPP based on its RUE
#  
#  Nij_sum <- sum(Nij)
#  Nij_hat <- Nij * (sort(c(Nij_sum, K_opt),partial=1)[1] / Nij_sum)
#  # now do an extinction filter based on population size
#  prob_extinction <- (1/(1+exp(-decay*(x0 -  Nij_hat))))
#  Nij_hat[which(sapply(prob_extinction, FUN=function(x){sample(c(0,1),1, prob= c(x, 1-(x)))})==0)] <- 0 
#  return(Nij_hat)
#}

# M3 - Competition defined soft limits
#apply_ecology <- function(abundance, traits, landscape, config) {
#  
#  # presence / absence is dependant on the niche width
#  browser()
#  x0 <- 1000 # inflexion point of extinction probability curve
#  omega <- 0.02 # environmental filtering (higher values equal less restrictive niche)
#  decay <- 0.001
#
#  # Set Individual-level carrying capacity based on temp and precipitation
#  K_opt <-  range01(landscape[, "mean_temp"], 6, 26) * 10000
#  
#  # difference between landscape and species optimum temperature niche
#  diff_t <- abs(traits[, "temp_niche_centre"] - landscape[, "mean_temp"])
#  
#  # get mean pairwise distances between species traits
#  diff_b <- rowMeans(as.matrix(dist(traits[, "body_size"])))
#  diff_b <- abs(diff_b - max(diff_b, na.rm=T))
#  
#  # potential population size based on resource use efficieny 
#  # omega is strength of environmental filtering
#  # will equal K_opt when tthe species is perfectly adapted
#  # Modified from McPeek 2008/2007
#  
#  # here we nerf the species with the most overlap with other species
#  Nij <- K_opt * exp(-(diff_t/omega)^2) * exp(-(diff_b/omega)^2)
#  
#  # realised abundance given a zero sum constraint
#  # if the RUE sum to less than 1 then each species gets that proportion of the NPP
#  # but if the RUE sums to more than one it needs to be rescaled so each species gets a proportion of the NPP based on its RUE
#  
#  Nij_sum <- sum(Nij)
#  Nij_hat <- Nij * (sort(c(Nij_sum, K_opt),partial=1)[1] / Nij_sum)
#  # now do an extinction filter based on population size
#  prob_extinction <- (1/(1+exp(-decay*(x0 -  Nij_hat))))
#  Nij_hat[which(sapply(prob_extinction, FUN=function(x){sample(c(0,1),1, prob= c(x, 1-(x)))})==0)] <- 0 
#  return(Nij_hat)
#}

# M4 - Metabolic Theory of Ecology
#apply_ecology <- function(abundance, traits, landscape, config) {
#  
#  # presence / absence is dependant on the niche width
#  browser()
#  x0 <- 1000 # inflexion point of extinction probability curve
#  omega <- 2 # environmental filtering (higher values equal less restrictive niche)
#  decay <- 0.001
#
#  # Set Individual-level carrying capacity based on temp and precipitation
#  K_opt <-  range01(landscape[, "mean_temp"], 6, 26) * 1000
#  
#  # difference between landscape and species optimum temperature niche
#  diff_t <- abs(traits[, "temp_niche_centre"] - landscape[, "mean_temp"])
#  
#  # potential population size based on resource use efficieny 
#  # omega is strength of environmental filtering
#  # will equal K_opt when the species is perfectly adapted
#  # modefied From McPeek 2008/2007
#  
#  # Modification
#  # here we also get a value for the metabolic rate and body size (resource use requirements of each species - log_BM)
#  # there is an exponentially negative effect of having high metabolic needs 
#  # e.g., coz your big and/or coz your metabolism is fast you need more food
#  log_BM <- log(traits[, "metabolic_rate"]^-1)*traits[, "body_size"]
#  Nij <- K_opt * exp(-(diff_t/omega)^2) * exp((log_BM)^-0.01) 
#
#  # realised abundance given a zero sum constraint
#  # if the RUE sum to less than 1 then each species gets that proportion of the NPP
#  # but if the RUE sums to more than one it needs to be rescaled so each species gets a proportion of the NPP based on its RUE
#  Nij_sum <- sum(Nij)
#  Nij_hat <- Nij * (sort(c(Nij_sum, K_opt),partial=1)[1] / Nij_sum)
#  
#  # now do an extinction filter based on population size
#  prob_extinction <- (1/(1+exp(-decay*(x0 -  Nij_hat))))
#  Nij_hat[which(sapply(prob_extinction, FUN=function(x){sample(c(0,1),1, prob= c(x, 1-(x)))})==0)] <- 0 
#  return(Nij_hat)
#}