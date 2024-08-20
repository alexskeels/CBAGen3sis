########################
### General settings ###
########################

random_seed =51801
start_time = 250
end_time = NA
max_number_of_species =20000
max_number_of_coexisting_species =20000
initial_abundance =  0.05

# custom function
{
  # ecological local equilibria variable J*
  get_J <- function(a_ff, a_fh, K_f){
    J <- sum((a_ff*K_f)/(a_ff-a_fh), na.rm=T)/(1+sum((a_fh/(a_ff-a_fh)), na.rm=T)) # new
    return(J)
  }
  # ecological environmental function
  fg <- function(x,a,b,c,ns=1){
    a <- a/c
    v <- a*exp(-((x-b)^2/(2*c^2)))
    return(v)
  }
}

apply_trs_tradeoff <- function(traits){
  #balance out traits in case they are above the traitoff surface
  cd <- traits[,c("dispersal", "competition")]
  # dispersal values on trade-off linear function. e.g. lines(x=c(0,1), y=c(0.9,1))
  dsurf <- -10*cd[,"competition"]+10
  csurf <- (10-cd[,"dispersal"])/10
  Dy <- cd[,"dispersal"]-dsurf
  # Dx <- cd[,"competition"]-csurf
  # get traits above and to the right of this linear function
  above_surface_mask <- Dy>0
  if (any(above_surface_mask)){ # then push them back to surface
    # print("traits were above the trait tradeoff surface... bringing them back to traide-off surface")
    # sample if either surface should be corresponding to d or c. 
    c_or_d <- sample(list(c(T,F), c(F,T)), size=sum(above_surface_mask), replace=TRUE)
    c_or_d <- do.call(rbind, c_or_d)
    cd[above_surface_mask,][c_or_d] <-   cbind(dsurf[above_surface_mask], csurf[above_surface_mask])[c_or_d]
    new_traits <- traits
    new_traits[, colnames(cd)] <- cd
    
    return(new_traits)
  }else{ #do nothing
    
    return(traits)
  }
}

# a list of traits to include with each species, traits with ss_eff_ are hack traits to extract in cell processes
trait_names = c("dispersal", "mean_temp", "temp_width", "competition")

# ranges to scale the input environemts with:
# not listed variable:         no scaling takes place
# listed, set to NA:           the environmental variable will be scaled from [min, max] to [0, 1]
# lsited with a given range r: the environmental variable will be scaled from [r1, r2] to [0, 1]
environmental_ranges = list("mean_temp"=c(9,26), "min_temp"=c(9,26),  "max_temp"=c(9,26))

#########################
### Observer Function ###
#########################

end_of_timestep_observer = function(data, vars, config){
  plot_richness(data$all_species, data$landscape)
}


######################
### Initialization ###
######################

create_ancestor_species <- function(landscape, config) {
  co <- landscape$coordinates
  new_species <- list()
  for(i in unique(landscape$environment[,"patch"])){
    initial_cells <- rownames(co)[landscape$environment[,"patch"]==i]
    # initial_cells <- sample(initial_cells, 1)
    new_species[[i]] <- create_species(initial_cells, config)
    #set local adaptation to max optimal temp equals local temp
    new_species[[i]]$traits[ , "dispersal"] <-0.65
    new_species[[i]]$traits[ , "competition"] <-0.2
    new_species[[i]]$traits[ , "mean_temp"] <- landscape$environment[initial_cells,"mean_temp"]
    new_species[[i]]$traits[ , "temp_width"] <- 0.2
    plot_species_presence(landscape, species=new_species[[i]])
  }
  return(new_species)
}


#################
### Dispersal ###
#################

# returns n dispersal values (proba distrib function)
get_dispersal_values <- function(n, species, landscape, config) {
  # hist(rweibull(117546,shape=2, scale=50)) # sample of # of drawn max for this simulation per time step.
  mean_abd <- mean(species$abundance)
  weight_abd <- species$abundance/mean_abd
  values <- rweibull(n, shape = 2, scale = 1+(mean(species$traits[,"dispersal"]*weight_abd)*49))  #from 1 to 50
  return(values)
}


##################
### Speciation ###
##################

# threshold for genetic distance after which a speciation event takes place
divergence_threshold =65 # between 10 and 50 ? as 0.1 to 0.5 Myrs or 100 - 500 kyrs

# adds a value of 1 to each geographic population cluster
get_divergence_factor <- function(species, cluster_indices, landscape, config) {
  return(1)
}

#######################
### Trait Evolution ###
#######################

apply_evolution <- function(species, cluster_indices, landscape, config) {
  trait_evolutionary_power <-0.01
  pw_tr_hom <- 0.5 # percentage of movement of local trait towards weighted trait mean within each population cluster 
  # pw_tr_hom = ZERO means no change while a value of ONE means that traits are equal within each populations cluster)
  traits <- species[["traits"]]
  # site names
  sites <- rownames(traits)
  #homogenize trait based on abundance
  # Homogenize all traits by weighted abundance, attention, exclude here any trait that should be more neutral
  # Bring traits closer to 50% of it's distance to the abundance weighted mean of the demographically connected 
  # population cluster
  trn <- config$gen3sis$general$trait_names
  for(cluster_index in unique(cluster_indices)){
    # cluster_index <- 1
    sites_cluster <- sites[which(cluster_indices == cluster_index)]
    B_t <- sum(species$abundance[sites_cluster])# total biomass
    B_r <- species$abundance[sites_cluster]/B_t# relative biomass
    for (ti in trn){
      # par(mfrow=c(1,2))
      # hist(traits[sites_cluster, ti], main="before")
      traits_ti <- traits[sites_cluster, ti]
      mean_trait <- sum(traits_ti*B_r)
      tr_hom <- mean_trait-traits_ti
      traits[sites_cluster, ti] <- traits_ti+(tr_hom*pw_tr_hom)
      # hist(traits[sites_cluster, ti], main="after")
    }
  }
  
  #mutate all traits except dispersal and competitive ability M0 and summary effective traits, which where excluded previously
  for (ti in c("mean_temp", "temp_width", "competition", "dispersal")){#trn[!trn%in%"dispersal"]){ # do not evolve dispersal
    if(ti=="competition"){ # check if competition and scale proportionally 0.1 variance /10 = 1, while 1 variance /1 = 1
      tep_sd=trait_evolutionary_power/10
    } else {
      tep_sd=trait_evolutionary_power
    }
    mutation_deltas <-rnorm(length(traits[, ti]), mean=0, sd=tep_sd)
    traits[, ti] <- traits[, ti] + mutation_deltas
  }
  # set bounds so that the species cant evolve a niche beyond that present in the data 
  # note that all temperature values are scaled between 0 and 1
  if(any(traits[, "temp_width"] > 1)){traits[which(traits[,"temp_width"]>1), "temp_width"] <- 1}
  if(any(traits[, "temp_width"] <= 0)){traits[which(traits[,"temp_width"]<=0), "temp_width"] <- 0.001} #limit of zero to avoid zero divisions on specialist/generalist trade-off
  if(any(traits[, "competition"] > 1)){traits[which(traits[,"competition"]>1), "competition"] <- 1}
  if(any(traits[, "competition"] < 0.9)){traits[which(traits[,"competition"]<0.9), "competition"] <- 0.9}
  if(any(traits[, "dispersal"] > 1)){traits[which(traits[,"dispersal"]>1), "dispersal"] <- 1}
  if(any(traits[, "dispersal"] < 0)){traits[which(traits[,"dispersal"]<0), "dispersal"] <- 0}

  # apply traits trade-off
  traits <- apply_trs_tradeoff(traits)
  
  return(traits)
}

#################################################
### Environmental and Ecological Interactions ###
#################################################

apply_ecology <- function(abundance, traits, landscape, config) {
  ns <- length(abundance)
  #### get rf, here r_f is the per capita growth rate of biomass that depends on the local site conditions 
  # set env niche
  env_min_fg <- fg(x=landscape[,"min_temp"], a=1, b=traits[, "mean_temp"], c=traits[, "temp_width"])
  env_max_fg <- fg(x=landscape[,"max_temp"], a=1, b=traits[, "mean_temp"], c=traits[, "temp_width"])
  # set growth rate
  g <- .1
  # abundance_tii first is only what the env. determines to be the new abundances
  r_f <- g*sqrt(env_min_fg*env_max_fg) # geometric mean
  ###### get (a_ff) = species interaction coefficient and (afh)= heterospecific interaction coefficient 
  # get traits Competition
  c_c <- rep(0.8,ns) # intra competition
  c_l <- traits[,"competition"]
  # set a_ff and a_fh
  a_ff <- 1-c_c
  a_fh <- 1-c_l
  # check if conditions are met in order to continue
  if (any(a_ff<=a_fh)){
    stop("a_ff has to be bigger than a_fh for this equilibrium function to be used! Check your intial and evolutionary conditions of the competition traits.")
  }
  ####### get K_f = the carrying capacity of species f (that is the equilibrium for the case without heterospecific biomass
  K_f <- r_f/a_ff
  ###### get J = the total biomass J* of the community in equilibrium
  J <- get_J(a_ff, a_fh, K_f)
  wistop <- FALSE
  keep_on_while <- rep(TRUE, ns)
  while(wistop==FALSE){
    shall_live <- (a_ff*K_f)>(a_fh*J)
    if (all(!shall_live)){ # in case all shall die
      B_f <- as.numeric(shall_live) # set all to zero, since we only have shall_live==FALSE
      wistop <- TRUE
    }
    if (all(shall_live[keep_on_while])){
      B_f <- ((a_ff*K_f)-(a_fh*J))/(a_ff-a_fh)
      B_f[!shall_live] <- 0
      wistop <- TRUE
    } else{
      a_ff[!shall_live] <- 0
      a_fh[!shall_live] <- 0
      K_f[!shall_live] <- 0
      keep_on_while <- shall_live
      if (sum(keep_on_while)==0){
        B_f[!shall_live] <- 0
        wistop <- TRUE
      }
      J <- get_J(a_ff, a_fh, K_f)
    }
  }
  names(B_f) <- names(abundance)
  B_f[B_f<0] <- 0 # set extinct species abundance to zero
  return(B_f)
}
