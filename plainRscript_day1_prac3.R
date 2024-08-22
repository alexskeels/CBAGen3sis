# Dynamic Island Landscape

# libraries
library(gen3sis)
library(terra)
library(ape)

# load landscapes
lc <- readRDS(file.path("data", "landscapes", "islands","landscapes.rds"))

# get names of landscape variables
names(lc)

# get first time step
first_step_pos <- ncol(lc$mean_temp)

# get first 10 sites of mean temperature for the 2 last time steps and the first (oldest) time step
lc$mean_temp[100:110, c(1:4, first_step_pos)]

# plot mean_temp for first and last time step
plot(rast(lc$mean_temp[ ,c(1:3, first_step_pos)]))

# plot the change in mean tempereature through time
mean_temperatures <- colMeans(lc$mean_temp[, 3:ncol(lc$mean_temp)], na.rm=T)
plot(0:500,mean_temperatures ,  xlab="time", ylab="mean temperature")

# plot the change in mean elevation through time
mean_elevation <- colMeans(lc$elevation[, 3:ncol(lc$elevation)], na.rm=T)
plot(0:500,mean_elevation ,  xlab="time", ylab="mean elevation")

# load config
conf <- create_input_config("configs/config_islands_simple_Day1Prac3_M1.R")
# list all main elements of the config file
names(conf$gen3sis)
# list all elements of the general section, i.e. the main settings and not so much on the eco-evolutionary processes
names(conf$gen3sis$general)


# Run a basic simulation

# run simulation
sim <- run_simulation(config = "configs/config_islands_simple_Day1Prac3_M1.R", 
                      landscape = "data/landscapes/islands", 
                      output_directory = "output/islands")


sim <- readRDS("output/islands/config_islands_simple_Day1Prac3_M1/sgen3sis.rds")

#check elements inside the sim object
names(sim)

# visualize the outputs
plot_summary(sim)

# plot richness from summary in custom fashion
na_mask <- is.na(lc$elevation[,"0"])
rich <- sim$summary$`richness-final`
rich[na_mask,3] <- NA
plot(rast(rich, type="xyz"), main="Richness")


# plot richness at time step 32 using saved data
sps32 <- readRDS("output/islands/config_islands_simple_Day1Prac3_M1/species/species_t_32.rds")
lc32 <- readRDS("output/islands/config_islands_simple_Day1Prac3_M1/landscapes/landscape_t_32.rds")
plot_richness(sps32, lc32)

# plot richness at time step 12 using saved data
sps12 <- readRDS("output/islands/config_islands_simple_Day1Prac3_M1/species/species_t_12.rds")
lc12 <- readRDS("output/islands/config_islands_simple_Day1Prac3_M1/landscapes/landscape_t_12.rds")
plot_richness(sps12, lc12)

# Now lets look at the Phylogeny
library(ape)
phy <- read.nexus("output/islands/config_islands_simple_Day1Prac3_M1/phy.nex")
plot(phy)


# Customize simulations

#For an example of a more complex simulation, with species abundances, traits trading off and evolving check the config_M2_TH.R.

# load in config
conf_m2 <- create_input_config("configs/config_islands_simple_Day1Prac3_M2.R")

# currently populations need 10 time steps to complete speciation
conf$gen3sis$speciation$divergence_threshold

# currently populations need 10 time steps to complete speciation
conf_m2$gen3sis$speciation$divergence_threshold


#We have also change the observer function to save the presence/absence matrix for each time step.
# see original config
conf$gen3sis$general$end_of_timestep_observer 

# see modeified config
conf_m2$gen3sis$general$end_of_timestep_observer 

# run the new model
sim_m2 <- run_simulation(config="configs/config_islands_simple_Day1Prac3_M2.R", landscape="data/landscapes/islands", output_directory="output/islands")

# From the small tweak to the divergence factor we can see more species being generated

# original dynamics
plot_summary(readRDS("output/islands/config_islands_simple_Day1Prac3_M1/sgen3sis.rds"))

# modified dynamics
plot_summary(readRDS("output/islands/config_islands_simple_Day1Prac3_M2/sgen3sis.rds"))


#Exercise

#Review the M2 file and try to understand what it's doing. How are the dynamics different to M1, and how are they similar? Consider how you might modify the configuration or apply it to a specific research question.

# lets make some more drastic changes now to see if we the patterns can be more radically changed

#

# run the new model
sim_m3 <- run_simulation(config="configs/config_islands_simple_Day1Prac3_M3_2.R", landscape="data/landscapes/islands", output_directory="output/islands")

sim_m4 <- run_simulation(config="configs/config_islands_simple_Day1Prac3_M4.R", landscape="data/landscapes/islands", output_directory="output/islands")


## 5. Troubleshoot

#Creating or modifying a gen3sis configuration can definitely lead to some weird errors, especially since it's so flexible. That's the downside, along with the steep learning curve. But if you're not too overwhelmed, you're doing great!

#Here are some handy debugging tips for when you run into those pesky errors:

#browser(): This function lets you pause the execution and explore what's going on. It's like a pit stop where you can check out variables and step through the code.

#To make your R session enter browser mode whenever you hit an error, you can use: \*options(error = recover)\*. This will help you diagnose and fix issues more easily.

#You can also condition browser calls, which is very useful when you want to stop execution at a specific time step or when a certain condition is met. Here's an example:
  

# Use 'stop_time' to halt execution at a specific timestep in the landscape object:

stop_time <- "25"

get_dispersal_values <- function(n, species, landscape, config) {
  
  if (landscape$timestep == stop_time) {
    
    browser()
    
  }
  
  # You can also check the 'vars' object for the current timestep:
  
  vars <- dynGet("vars", inherits = TRUE)
  
  if (vars$ti == stop_time) {
    
    browser()
    
  }
  
}


