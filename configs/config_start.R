########################
### General settings ###
########################

random_seed =87527
start_time = NA
end_time = NA
max_number_of_species =20000
max_number_of_coexisting_species =20000
initial_abundance =  0.05

# serves as a template
#ss_eff_emp <- list("dispersal"=NA, "dispersal_success"=NA, "environmental_filter"=NA, "competition_c"=NA, "competition_l"=NA)
#ss_eff <- list()

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

# temperature performance curve
ftpc <- function(x, opt, breadth, rate, amplitude){
  # hard coded trade-off between amplitude and breadth
  amplitude <- amplitude/breadth
  p <- amplitude*exp(-(x-opt)^2/(amplitude*breadth^2))*pnorm(rate*(x-opt)/breadth)
  return(p)
}
# plot the temperature performance curve along a gradient
# plot(ftpc(seq(5,40,0.3), 26, 2, -8, 8), type="l")
# lines(ftpc(seq(5,40,0.3), 26, 5, -8, 8), type="l")


# set for first time
# assign("ss_eff", ss_eff_emp, envir = .GlobalEnv)
# a list of traits to include with each species, traits with ss_eff_ are hack traits to extract in site processes
trait_names = c("dispersal", "tpc_breadth", "tpc_opt", "tpc_rate", "plasticity", "sd_breath", "sd_opt", "sd_rate")

# ranges to scale the input environemts with:
# not listed variable:         no scaling takes place
# listed, set to NA:           the environmental variable will be scaled from [min, max] to [0, 1]
# lsited with a given range r: the environmental variable will be scaled from [r1, r2] to [0, 1]
# environmental_ranges = list("mean_temp"=c(9,26), "min_temp"=c(9,26),  "max_temp"=c(9,26))

#########################
### Observer Function ###
#########################

end_of_timestep_observer = function(data, vars, config){
  # browser()
  # jpeg(paste0(config$directories$output_plots, "/ranges_t", vars$ti, ".jpeg"))
  # {
  #   par(mfrow=c(1,1))
    plot_ranges(data$all_species, data$landscape)
  # }
  # dev.off()
  
  #  plot_richness(data$all_species, data$landscape)
  # save_species()
  # save_abundance()
  # save_divergence()
  # save_occupancy()
  # save_phylogeny()
  save_traits()
  
  # # make p/a matrices if necessary
  # if(!file.exists(file.path(config$directories$output, "occs"))){dir.create(file.path(config$directories$output, "occs"))}
  # # site names
  # all_sites <- rownames(data$landscape$coordinates)
  # # get 0 for absence and 1 for presence in each grid site
  # all_species_presence <- do.call( cbind, lapply(data$all_species, FUN = function(x) {ifelse(all_sites %in% names(x$abundance), 1, 0)}))
  # # colnames are species names
  # colnames(all_species_presence ) <- unlist(lapply(data$all_species, function(x){x$id}))
  # # column bind with x/y coordinates
  # presence_absence_matrix <- cbind(data$landscape$coordinates, all_species_presence)
  # saveRDS(presence_absence_matrix, file=file.path(config$directories$output,"occs",  paste0("pa_t_",vars$ti, ".rds")))
  
  
  # if(!file.exists(file.path(config$directories$output, "mean_traits_sp"))){dir.create(file.path(config$directories$output, "mean_traits_sp"))}
  # saveRDS(data$eco_by_sp, file=file.path(config$directories$output,"mean_traits_sp",  paste0("mean_trs_t_",vars$ti, ".rds")))
}


######################
### Initialization ###
######################

create_ancestor_species <- function(landscape, config) {
  co <- landscape$coordinates
  new_species <- list()
  for(i in unique(landscape$environment[,"island_id"])){
    # i <- 1
    initial_sites <- rownames(co)[landscape$environment[,"island_id"]==i]
    # initial_sites <- sample(initial_sites, 1)
    new_species[[i]] <- create_species(initial_sites, config)
    #set local adaptation to max optimal temp equals local temp
    new_species[[i]]$traits[ , "dispersal"] <- 1
    new_species[[i]]$traits[ ,  "tpc_breadth"] <- 3
    new_species[[i]]$traits[ ,  "tpc_opt"] <- landscape$environment[initial_sites,"mean_temp"]
    new_species[[i]]$traits[ ,  "tpc_rate"] <- -8
    new_species[[i]]$traits[ ,  "plasticity"] <- 2
    new_species[[i]]$traits[ ,  "sd_breath"] <- 0.8
    new_species[[i]]$traits[ ,  "sd_opt"] <- 2
    new_species[[i]]$traits[ ,  "sd_rate"] <- 0.01
    # plot_species_presence(landscape, species=new_species[[i]])
  }
  return(new_species)
}


#################
### Dispersal ###
#################

# returns n dispersal values (proba distrib function)
get_dispersal_values <- function(n, species, landscape, config){
  # hist(rweibull(117546,shape=2, scale=50)) # sample of # of drawn max for this simulation per time step.
  mean_abd <- mean(species$abundance)
  weight_abd <- species$abundance/mean_abd
  values <- rweibull(n, shape = 2, scale = 1+(mean(species$traits[,"dispersal"]*weight_abd)*10))  #from 5 to 50
  return(values)
}


##################
### Speciation ###
##################

# threshold for genetic distance after which a speciation event takes place
divergence_threshold = 4 # between 10 and 50 ? as 0.1 to 0.5 Myrs or 100 - 500 kyrs

# adds a value of 1 to each geographic population cluster
get_divergence_factor <- function(species, cluster_indices, landscape, config) {
  return(1)
}

#######################
### Trait Evolution ###
#######################

apply_evolution <- function(species, cluster_indices, landscape, config) {
  trait_evolutionary_power <-0.01 # variance 
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
    N_t <- sum(species$abundance[sites_cluster])# total abundance
    N_r <- species$abundance[sites_cluster]/N_t# relative abundance vector
    for (ti in trn){
      # ti <- trn[1]
      # par(mfrow=c(1,2))
      # hist(traits[sites_cluster, ti], main="before")
      traits_ti <- traits[sites_cluster, ti]
      mean_trait <- sum(traits_ti*N_r)
      tr_hom <- mean_trait-traits_ti
      traits[sites_cluster, ti] <- traits_ti+(tr_hom*pw_tr_hom)
    }
  }
  
  #mutate all traits except dispersal and competitive ability M0 and summary effective traits, which where excluded previously
  for (ti in c("dispersal", "tpc_breadth", "tpc_opt")){#trn[!trn%in%"dispersal"]){ # do not evolve dispersal
    # hsqr <- (traits[,"sd_breath"]^2+trait_evolutionary_power)
    # hsqr <- hsqr/(hsqr+landscape$environment[sites, "sd_breath"]^2)
    # traits[, ti] <- traits[, ti]*hsqr
    mutation_deltas <-rnorm(length(traits[, ti]), mean=0, sd=trait_evolutionary_power)
    traits[, ti] <- traits[, ti] + mutation_deltas
  }
  # set bounds so that the species cant evolve a niche beyond that present in the data 
  # note that all temperature values are scaled between 0 and 1
  if(any(traits[, "dispersal"] > 1)){traits[which(traits[,"dispersal"]>1), "dispersal"] <- 1}
  if(any(traits[, "dispersal"] < 0)){traits[which(traits[,"dispersal"]<0), "dispersal"] <- 0}
  if(any(traits[, "tpc_breadth"] < 0.1)){traits[which(traits[,"tpc_breadth"]>0.1), "tpc_breadth"] <- 0.1}

  
  return(traits)
}

#################################################
### Environmental and Ecological Interactions ###
#################################################

apply_ecology <- function(abundance, traits, landscape, config) {
  ns <- length(abundance)
  #### get rf, here r_f is the per capita growth rate of biomass that depends on the local site conditions 
  # set env niche
  
  random_local_temps <- rnorm(10, mean=landscape[,"mean_temp"], sd=landscape[,"sd_temp"])
  
  
  
  r_f <- rep(NA, ns)
  for (sp_i in 1:ns){
     env_performance <- ftpc(x=random_local_temps, 
                     opt=traits[sp_i, "tpc_opt"], 
                     breadth=traits[sp_i, "tpc_breadth"], 
                     rate=traits[sp_i, "tpc_rate"], 
                     amplitude=20)
     r_f[sp_i] <- mean(env_performance)
  }
  
  # intra and inter specific competition coeficient
  # set a_ff and a_fh
  a_ff <- rep(0.5, ns)
  a_fh <- rep(0.2, ns)
  
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
    # 
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
