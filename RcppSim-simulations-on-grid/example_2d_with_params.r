library(pacman)
p_load(tidyverse)
library(MathBioSim)
library (MASS)
library (mnormt)
library(msm)

experiment = 'c'
switch(experiment, 
       a = {
         sw=0.15
         sm=0.15
         filename = "step_01_a.csv"
       },
       b = {
         sw=0.05
         sm=0.05
         filename = "step_01_b.csv"
       },
       c = {
         sw = 0.015
         sm = 0.15
         filename = "step_01_c.csv"
       }
)

death_radius_cutoff = sw*5
birth_radius_cutoff = sm*5

r <- seq(0, death_radius_cutoff, length.out = 1001)

sim<-initialize_simulator(area_length_x = 1, 
                          area_length_y = 1,
                          b=0.4, d=0.2, dd=0.001, 
                          initial_population_x = runif(200, 0, 1), 
                          initial_population_y =runif(200, 0, 1), 
                          death_r = death_radius_cutoff,
                          death_y = sapply(r, function(r_value) {
                            1 / (2 * sw * sw *pi) * exp((-1 / (2 * sw * sw)) * r_value * r_value)
                          }),
                          birth_ircdf_y = VGAM::qrayleigh(seq(0,0.9999, length.out = 1001), scale= sm),
                          realtime_limit = 100, 
                          ndim=2)

sim_results <- run_simulation(sim, 100, calculate.pcf = TRUE, pcf_grid=seq(0, pmax(sm,sw), length.out = 101))

ggplot(sim_results$population,aes(x=time,y=pop))+
  geom_line()
ggplot(head(sim_results$time_steps, n=1000),aes(x=time,y=pop))+
  geom_line()

write.csv(sim_results[["time_steps"]], file = filename)

#ggplot(sim_results$pattern,aes(x=x))+
  #geom_histogram(bins=10,breaks=seq(0,10,length.out = 11))
#ggplot(sim_results$pattern,aes(x=y))+
  #geom_histogram(bins=10,breaks=seq(0,10,length.out = 11))



