# https://github.com/Health-Economics-in-R/CEdecisiontree
# http://darthworkgroup.com/r-code/
# https://github.com/DARTH-git/cohort-modeling-tutorial-intro
# https://github.com/DARTH-git/Cohort-modeling-tutorial/blob/master/R/old/Functions%20STM_01.R
################################ Initial setup ############################### 
rm(list = ls())    # remove any variables in R's memory 


### Load packages
library(dplyr)    
library(tidyr)
library(reshape2) 
library(ggplot2)
library(gridExtra)
library(scales)    
library(boot)
library(dampack) 
library(darthtools)

# General setup for Markov model disease states 
cycle.length <- 0.5  # 60months (=5 years) for basic cycle length
n.age.init   <- 60  # age60 as a baseline
n.age.max    <- 100 # age100 as a max
n.cycles     <- (n.age.max - n.age.init)/cycle.length  # time horizon, number of cycles
state.names  <- c("<1cm","1-2cm","2-3cm","WF","Cancer","Death") # 3diff cyst sizes, worrisome feature, cancer, death
n.states     <- length(state.names) # number of states in disease progression
surgery.name <- c("No comp", "Comp", "Cancer", "Cancer", "Death") # states in surgery markov 
s.states     <- length(surgery.name)
d_e          <- 0.03 # annual discount rate for DALY 3%
d_c          <- 0.03 # annual discount rate for cost 3%
names.strat  <- c("surgery","surveillance")
n.strat      <- length(names.strat)


# Transition rates and hazard ratios (HRs)
r_HD         # rate of dying for <1cm cohort (all-cause mortality)
r_HWF        # rate of transitioning from <1cm to WF
r_HC         # rate of transitioning from <1cm to cancer
r_HS1        # rate of transitioning from <1cm to 1-2cm
r_S1S2       # rate of transitioning from 1-2cm to 2-3cm
hrd_S1       # hazard ratio of death in <1cm vs. 1-2cm
hrd_S2       # hazard ratio of death in <1cm vs. 2-3cm
hrw_S1       # hazard ratio of WF in <1cm vs. 1-2cm
hrw_S2       # hazard ratio of WF in <1cm vs. 2-3cm
hrc_S1       # hazard ratio of cancer in <1cm vs. 1-2cm
hrc_S2       # hazard ratio of cancer in <1cm vs. 2-3cm

#state rewards ---
###costs 


###utilities


#initial state vector
v_init       <- c(H1=0.36, H2=0.3, H3=0.1, H4=0.6*0.4, H5=0, H6=0) # initial state: <1cm w/o WF, <1cm+WF, 1-2cm, 2-3cm, cancer, death


#cohort trace matrix (No WF surveillance): Distribution of the cohorts
m_M          <- matrix(NA, nrow= (n.cycles+1), ncol=n.states, 
                       dimnames = list(0:n.cycles, state.names))
m_M[1,]      <- v_init # store initial state vector in the first row of the cohort trace

#transition probability matrix for disease progression: 6*6 matrix 
m_P          <- matrix(0, nrow     =n.states, ncol=n.states, 
                      dimnames    = list(state.names, state.names))

#transition probability matrix for surgery 
m_S          <- matrix(0, nrow = s.states, ncol=s.states,
                       dimnames = list(surgery.name, surgery.name))

#cohort trace matrix for those who undergo surgery



