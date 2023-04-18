# https://github.com/Health-Economics-in-R/CEdecisiontree
# http://darthworkgroup.com/r-code/
# https://github.com/DARTH-git/cohort-modeling-tutorial-intro
# https://github.com/DARTH-git/Cohort-modeling-tutorial/blob/master/R/old/Functions%20STM_01.R
################################ Initial setup ############################### 
rm(list = ls())    # remove any variables in R's memory 
setwd("D:/OneDrive - London School of Hygiene and Tropical Medicine/SNU")

### Load packages
library(dplyr)    
library(tidyr)
library(reshape2) 
library(ggplot2)
library(gridExtra)
library(scales)    
library(boot)
library(dampack) 
library(ggplot2)
library(ggforce)
library(data.table)
library(ggridges)
library(ggpubr)
options(scipen = 999)

#-------------------------------------------------------------------------------
# General setup for Markov model disease states 
#-------------------------------------------------------------------------------
cycle.length <- 0.5  # 6months (=0.5 year) for basic cycle length
n.age.init   <- 60   # age60 as a baseline
n.age.max    <- 100  # age100 as a max
n.cycles     <- (n.age.max - n.age.init)/cycle.length  # time horizon, number of cycles
state.names  <- c("<1cm","1-2cm","2-3cm","WF","Cancer","Death") # 3diff cyst sizes, worrisome feature, cancer, death
n.states     <- length(state.names) # number of states in disease progression
surgery.name <- c("Surg", "No comp", "Comp", "Cancer", "Death") # states in surgery markov 
s.states     <- length(surgery.name)
d_e          <- 0.03 # annual discount rate for DALY 3%
d_c          <- 0.03 # annual discount rate for cost 3%
names.strat  <- c("surgery","surveillance")
n.strat      <- length(names.strat)

#------------------------------------------------------------------------------- 
#Transition rates and hazard ratios for Disease Progression Markov (HRs)
#-------------------------------------------------------------------------------

rate1        <- fread("D:/OneDrive - London School of Hygiene and Tropical Medicine/SNU/rates.csv")
rate2        <- fread("D:/OneDrive - London School of Hygiene and Tropical Medicine/SNU/rate2.csv")
rate3        <- fread("D:/OneDrive - London School of Hygiene and Tropical Medicine/SNU/rate3.csv")
cancerdeath  <- fread("D:/OneDrive - London School of Hygiene and Tropical Medicine/SNU/cancerdeath.csv")


my_mod <- lm(y ~ x, rate1) 
plot(rate1$x,                       # Draw line plot with regression line
     rate1$y,
     type = "l")
lines(rate1$x,
      predict(my_mod),
      col = 2,
      lwd = 2)
my_coef <- coef(my_mod)            # Extract coefficients of model
my_coef 
my_equation <- paste("y =",        # Extract equation of model
                     coef(my_mod)[[1]],
                     "+",
                     coef(my_mod)[[2]],
                     "* x")
my_equation                        # Print equation of model

my_mod2 <- lm(y ~ x, rate2) 
plot(rate2$x,                       # Draw line plot with regression line
     rate2$y,
     type = "l")
lines(rate2$x,
      predict(my_mod2),
      col = 2,
      lwd = 2)
my_coef2 <- coef(my_mod2)            # Extract coefficients of model
my_coef2 

my_mod3 <- lm(y ~ x, cancerdeath) 
plot(cancerdeath$x,                       # Draw line plot with regression line
     cancerdeath$y,
     type = "l")
lines(cancerdeath$x,
      predict(my_mod3),
      col = 2,
      lwd = 2)
my_coef3 <- coef(my_mod3)            # Extract coefficients of model
my_coef3 


r_H1D        <- (0.042533081/2) # rate of dying for <1cm cohort (all-cause mortality)
r_H1WF       <- (0.0158/2) # rate of transitioning from <1cm to WF
r_H1C        <- (0.0016/2) # rate of transitioning from <1cm to cancer
r_H1H2       <- (0.000283/0.5)  # rate of transitioning from <1cm to 1-2cm 
r_H1H3       <- 0 # rate of transitioning form <1cm to 2-3cm
r_H2H3       <- (0.0004587781/0.5)  # rate of transitioning from 1-2cm to 2-3cm
r_H2WF       <- (0.025/2) # rate of transitioning from 1-2cm to WF
r_H3WF       <- (0.0608/2) # rate of transitioning from 2-3cm to WF
r_H2C        <- (0.0021/2) # rate of transitioning from 1-2cm to cancer
r_H3C        <- (0.0042/2) # rate of transitioning from 2-3cm to cancer
r_H2D        <- (0.0064/2) # rate of dying for 1-2cm cohort (all-cause mortality)
r_H3D        <- (0.0135/2) # rate of dying for 2-3cm cohort (all-cause mortality)
hrc_WF       <- 5.592157   # hazard ratio of cancer in <1cm vs. WF
hrd_WF       <- 1   # hazard ratio of death in <1cm vs. WF (same as benign death = natural death)
hrd_C        <- 2.515   # hazard ratio of death in <1cm vs. Cancer 

#------------------------------------------------------------------------------- 
#Transition rates and hazard ratios for Surgery Markov 
#-------------------------------------------------------------------------------
r_surgComp    <- 0.2098871
r_surgDeath   <- 0
r_CompCancer  <- 0.3784219
r_CompDeath   <- 0.0027369
r_CancerDeath <- 0.0002472906/0.5
r_surgNocomp  <- 1 - r_surgComp
r_NocompCancer<- 0.3784219
r_NocompDeath <- 11.1/100000 # https://seer.cancer.gov/statfacts/html/pancreas.html

#transform rates into transition probs
p_surgComp     <- 1 - exp(-r_surgComp*cycle.length)
p_surgDeath    <- 0
p_CompCancer   <- 1 - exp(-r_CompCancer*cycle.length)
p_CompDeath    <- 1 - exp(-r_CompDeath*cycle.length)
p_CancerDeath  <- 1 - exp(-r_CancerDeath*cycle.length)
p_surgNocomp   <- 1 - exp(-r_surgNocomp*cycle.length)
p_NocompCancer <- 1 - exp(-r_NocompCancer*cycle.length)
p_NocompDeath  <- 1 - exp(-r_NocompDeath*cycle.length)
  
#initial distribution of different cysts 
d_H1             <- 1058/(1058+1728+617)
d_H2             <- 1728/(1058+1728+617)
d_H3             <- 617/(1058+1728+617)
d_H1WF           <- 23/1058
d_H2WF           <- 69/1728
d_H3WF           <- 103/617
d_H1WFSurg       <- (7/23)*d_H1WF
d_H2WFSurg       <- (44/69)*d_H2WF
d_H3WFSurg       <- (94/103)*d_H3WF

#transform rates into transition probs
p_H1H2           <- 1 - exp(-r_H1H2*cycle.length)
p_H1H3           <- 1 - exp(-r_H1H3*cycle.length)
p_H2H3           <- 1 - exp(-r_H2H3*cycle.length)
p_H1WF           <- 1 - exp(-r_H1WF*cycle.length)
p_H1C            <- 1 - exp(-r_H1C*cycle.length)
p_H1D            <- 1 - exp(-r_H1D*cycle.length)
p_H2WF           <- 1 - exp(-r_H2WF*cycle.length)
p_H2C            <- 1 - exp(-r_H2C*cycle.length)
p_H2D            <- 1 - exp(-r_H2D*cycle.length)
p_H3WF           <- 1 - exp(-r_H3WF*cycle.length)
p_H3C            <- 1 - exp(-r_H3C*cycle.length)
p_H3D            <- 1 - exp(-r_H3D*cycle.length)
p_WFC            <- 1 - exp(-r_H1C*hrc_WF*cycle.length)
p_WFD            <- 1 - exp(-r_H1D*hrd_WF*cycle.length)
p_CD             <- 1 - exp(-r_H1D*hrd_C*cycle.length)



#-------------------------------------------------------------------------------
# START PROGRAM: <1CM
#-------------------------------------------------------------------------------
#initial state vector (<1cm): total sum should be d_H1
init.surgery     <- c("Surg" =d_H1WFSurg, "Nocomp"= 0, "Comp"=0, "Cancer"=0, "Death"=0) # apply surgery markov
init.surv        <- c("H1"= 0,"H2"= 0,"H3"= 0,"WF"= d_H1WF - d_H1WFSurg,"Cancer"= 0,"Death"= 0)# apply disease markov
init.noWF        <- c("H1" = d_H1-d_H1WF, "H2"= 0,"H3"= 0,"WF"= 0,"Cancer"= 0,"Death"= 0)# apply disease markov

#initial state vector (1-2cm): total sum should be d_H2
init.surgeryH2    <- c("Surg" =d_H2WFSurg, "Nocomp"= 0, "Comp"=0, "Cancer"=0, "Death"=0) # apply surgery markov
init.survH2       <- c("H1"= 0,"H2"= 0,"H3"= 0,"WF"= d_H2WF - d_H2WFSurg,"Cancer"= 0,"Death"= 0)# apply disease markov
init.noWFH2       <- c("H1" = 0, "H2"= d_H2-d_H2WF,"H3"= 0,"WF"= 0,"Cancer"= 0,"Death"= 0)# apply disease markov

#initial state vector (2-3cm): total sum should be d_H3
init.surgeryH3    <- c("Surg" =d_H3WFSurg, "Nocomp"= 0, "Comp"=0, "Cancer"=0, "Death"=0) # apply surgery markov
init.survH3       <- c("H1"= 0,"H2"= 0,"H3"= 0,"WF"= d_H3WF - d_H3WFSurg,"Cancer"= 0,"Death"= 0)# apply disease markov
init.noWFH3       <- c("H1" = 0, "H2"= 0,"H3"= d_H3-d_H3WF,"WF"= 0,"Cancer"= 0,"Death"= 0)# apply disease markov


#cohort trace matrix: Distribution of the cohorts (<1cm)
m_M_surgery1     <- matrix(NA, nrow = (n.cycles+1), ncol=s.states, 
                       dimnames = list(0:n.cycles, surgery.name))
m_M_surgery1[1,] <- init.surgery

m_M_surv         <- matrix(NA, nrow = (n.cycles+1), ncol=n.states, 
                           dimnames = list(0:n.cycles, state.names))
m_M_surv[1,]     <- init.surv

m_M_noWf         <- matrix(NA, nrow = (n.cycles+1), ncol=n.states, 
                           dimnames = list(0:n.cycles, state.names))
m_M_noWf[1,]     <- init.noWF

#cohort trace matrix: Distribution of the cohorts (1-2cm)
m_M_surgery2      <- matrix(NA, nrow = (n.cycles+1), ncol=s.states, 
                           dimnames = list(0:n.cycles, surgery.name))
m_M_surgery2[1,]  <- init.surgeryH2

m_M_surv2         <- matrix(NA, nrow = (n.cycles+1), ncol=n.states, 
                           dimnames = list(0:n.cycles, state.names))
m_M_surv2[1,]     <- init.survH2

m_M_noWf2         <- matrix(NA, nrow = (n.cycles+1), ncol=n.states, 
                           dimnames = list(0:n.cycles, state.names))
m_M_noWf2[1,]     <- init.noWFH2

#cohort trace matrix: Distribution of the cohorts (2-3cm)
m_M_surgery3      <- matrix(NA, nrow = (n.cycles+1), ncol=s.states, 
                            dimnames = list(0:n.cycles, surgery.name))
m_M_surgery3[1,]  <- init.surgeryH3

m_M_surv3         <- matrix(NA, nrow = (n.cycles+1), ncol=n.states, 
                            dimnames = list(0:n.cycles, state.names))
m_M_surv3[1,]     <- init.survH3

m_M_noWf3         <- matrix(NA, nrow = (n.cycles+1), ncol=n.states, 
                            dimnames = list(0:n.cycles, state.names))
m_M_noWf3[1,]     <- init.noWFH3


#transition probability matrix for disease progression: 6*6 matrix 
m_P                   <- matrix(0, nrow    = n.states, ncol=n.states, 
                      dimnames    = list(state.names, state.names))
m_P["<1cm","1-2cm"]   <- p_H1H2
m_P["<1cm","2-3cm"]   <- p_H1H3
m_P["<1cm","WF"]      <- p_H1WF
m_P["<1cm","Cancer"]  <- p_H1C
m_P["<1cm","Death"]   <- p_H1D
m_P["<1cm","<1cm"]    <- 1 - (p_H1H2 + p_H1H3 + p_H1WF + p_H1C + p_H1D)
m_P["1-2cm","<1cm"]   <- 0
m_P["1-2cm","2-3cm"]  <- p_H2H3
m_P["1-2cm","WF"]     <- p_H2WF
m_P["1-2cm","Cancer"] <- p_H2C
m_P["1-2cm","Death"]  <- p_H2D
m_P["1-2cm","1-2cm"]  <- 1 - (p_H2H3 + p_H2WF + p_H2C + p_H2D)
m_P["2-3cm","<1cm"]   <- 0
m_P["2-3cm","1-2cm"]  <- 0
m_P["2-3cm","WF"]     <- p_H3WF
m_P["2-3cm","Cancer"] <- p_H3C
m_P["2-3cm","Death"]  <- p_H3D
m_P["2-3cm","2-3cm"]  <- 1 - (p_H3WF + p_H3C + p_H3D)
m_P["WF","Cancer"]    <- p_WFC
m_P["WF","Death"]     <- p_WFD
m_P["WF","WF"]        <- 1 - (p_WFC + p_WFD)
m_P["Cancer","Death"] <- p_CD
m_P["Cancer","Cancer"]<- 1 - p_CD
m_P["Death", "Death"] <- 1

rowSums(m_P) == 1

#transition probability matrix for surgery 
m_S                        <- matrix(0, nrow  = s.states, ncol=s.states,
                              dimnames = list(surgery.name, surgery.name))
m_S["Surg", "No comp"]     <- p_surgNocomp
m_S["Surg", "Comp"]        <- p_surgComp
m_S["Surg", "Cancer"]      <- 0
m_S["Surg", "Death"]       <- p_surgDeath
m_S["Surg", "Surg"]        <- 1 - (p_surgNocomp + p_surgComp)
m_S["No comp", "Comp"]     <- 0
m_S["No comp", "Cancer"]   <- p_NocompCancer
m_S["No comp", "Death"]    <- p_NocompDeath
m_S["No comp", "No comp"]  <- 1 - (p_NocompCancer + p_NocompDeath)
m_S["Comp", "Cancer"]      <- p_CompCancer
m_S["Comp", "Death"]       <- p_CompDeath
m_S["Comp", "Comp"]        <- 1 - (p_CompCancer + p_CompDeath)
m_S["Cancer", "Death"]     <- p_CancerDeath
m_S["Cancer", "Cancer"]    <- 1 - p_CancerDeath
m_S["Death", "Death"]      <- 1

rowSums(m_S) == 1

#-------------------------------------------------------------------------------
# <1cm
#-------------------------------------------------------------------------------
#cohort trace matrix for those who undergo surgery 

#iterative solutions for disease progression
for(t in 1:n.cycles) {
  m_M_noWf[t+1, ] <- m_M_noWf[t, ] %*% m_P
}
# check
rowSums(m_M_noWf)

#cohort distribution graph
m_M_noWf <- as.data.frame(m_M_noWf)
noWfStack <- stack(m_M_noWf)
noWfStack$cycle <- rep(0:80, time(6))
colnames(noWfStack) <- c("Proportion", "State", "Cycle")


ggplot(noWfStack, aes(x = Cycle, y = Proportion))+
  geom_line()+
  facet_wrap(~`State`)+
  #labs(title = col)+
    theme_minimal()
  
noWf1 <- ggplot(noWfStack, aes(x = Cycle, y = Proportion, color = State)) +
            geom_line() +
             theme_minimal()+
               ylim(0,0.5)

#iterative solutions for surgery markov
for(t in 1:n.cycles) {
  m_M_surgery1[t+1, ] <- m_M_surgery1[t, ] %*% m_S
}
# check
rowSums(m_M_surgery1)

m_M_surgery1 <- as.data.frame(m_M_surgery1)
surg1Stack <- stack(m_M_surgery1)
surg1Stack$cycle <- rep(0:80, time(5))
colnames(surg1Stack) <- c("Proportion", "State", "Cycle")


ggplot(surg1Stack, aes(x = Cycle, y = Proportion))+
  geom_line()+
  facet_wrap(~`State`)+
  #labs(title = col)+
  theme_minimal()

surg1 <- ggplot(surg1Stack, aes(x = Cycle, y = Proportion, color = State)) +
            geom_line() +
             theme_minimal()+
               ylim(0,0.2)

#iterative solutions for surveillance + WF markov 
for(t in 1:n.cycles) {
  m_M_surv[t+1, ] <- m_M_surv[t, ] %*% m_P
}
# check
rowSums(m_M_surv)

m_M_surv <- as.data.frame(m_M_surv)
surv1Stack <- stack(m_M_surv)
surv1Stack$cycle <- rep(0:80, time(6))
colnames(surv1Stack) <- c("Proportion", "State", "Cycle")

ggplot(surv1Stack, aes(x = Cycle, y = Proportion))+
  geom_line()+
  facet_wrap(~`State`)+
  #labs(title = col)+
  theme_minimal()

wf1surv <- ggplot(surv1Stack, aes(x = Cycle, y = Proportion, color = State)) +
            geom_line() +
              theme_minimal()+
                ylim(0,0.02)
  

combinedH1 <- cbind(m_M_noWf, m_M_surv, m_M_surgery1)
write.csv(combinedH1, file="combined_H1.csv", row.names = F)

#-------------------------------------------------------------------------------
# 1-2cm
#-------------------------------------------------------------------------------

#cohort trace matrix for those who undergo surgery 

#iterative solutions for disease progression
for(t in 1:n.cycles) {
  m_M_noWf2[t+1, ] <- m_M_noWf2[t, ] %*% m_P
}
# check
rowSums(m_M_noWf2)

#cohort distribution graph
m_M_noWf2 <- as.data.frame(m_M_noWf2)
noWf2Stack <- stack(m_M_noWf2)
noWf2Stack$cycle <- rep(0:80, time(6))
colnames(noWf2Stack) <- c("Proportion", "State", "Cycle")


ggplot(noWf2Stack, aes(x = Cycle, y = Proportion))+
  geom_line()+
  facet_wrap(~`State`)+
  #labs(title = col)+
  theme_minimal()

noWf2 <- ggplot(noWf2Stack, aes(x = Cycle, y = Proportion, color = State)) +
           geom_line() +
             theme_minimal()+
              ylim(0,0.5)

#iterative solutions for surgery markov
for(t in 1:n.cycles) {
  m_M_surgery2[t+1, ] <- m_M_surgery2[t, ] %*% m_S
}
# check
rowSums(m_M_surgery2)

m_M_surgery2 <- as.data.frame(m_M_surgery2)
surg2Stack <- stack(m_M_surgery2)
surg2Stack$cycle <- rep(0:80, time(5))
colnames(surg2Stack) <- c("Proportion", "State", "Cycle")


ggplot(surg2Stack, aes(x = Cycle, y = Proportion))+
  geom_line()+
  facet_wrap(~`State`)+
  #labs(title = col)+
  theme_minimal()

surg2 <- ggplot(surg2Stack, aes(x = Cycle, y = Proportion, color = State)) +
          geom_line() +
           theme_minimal()+
             ylim(0,0.2)

#iterative solutions for surveillance + WF markov 
for(t in 1:n.cycles) {
  m_M_surv2[t+1, ] <- m_M_surv2[t, ] %*% m_P
}
# check
rowSums(m_M_surv2)

m_M_surv2 <- as.data.frame(m_M_surv2)
surv2Stack <- stack(m_M_surv2)
surv2Stack$cycle <- rep(0:80, time(6))
colnames(surv2Stack) <- c("Proportion", "State", "Cycle")

ggplot(surv2Stack, aes(x = Cycle, y = Proportion))+
  geom_line()+
  facet_wrap(~`State`)+
  #labs(title = col)+
  theme_minimal()

wf2surv <-  ggplot(surv2Stack, aes(x = Cycle, y = Proportion, color = State)) +
             geom_line() +
               theme_minimal()+
                 ylim(0,0.02)

combinedH2 <- cbind(m_M_noWf2, m_M_surv2, m_M_surgery2)
write.csv(combinedH2, file="combined_H2.csv", row.names = F)

figure1 <- ggarrange(noWf1, noWf2, noWf3,
                     labels = c("<1cm", "1-2cm", "2-3cm"),
                     ncol = 3, nrow = 1)
ggsave("nowf.pdf", figure1, dpi=1000, device= "pdf", height=8, width=15,units="in", bg=NULL, limitsize = F)

figure2 <- ggarrange(surg1, surg2, surg3,
                     labels = c("<1cm", "1-2cm", "2-3cm"),
                     ncol = 3, nrow = 1)
ggsave("surg.pdf", figure2, dpi=1000, device= "pdf", height=8, width=15,units="in", bg=NULL, limitsize = F)

figure2 <- ggarrange(wf1surv, wf2surv, wf3surv,
                     labels = c("<1cm", "1-2cm", "2-3cm"),
                     ncol = 3, nrow = 1)
ggsave("wfsurv.pdf", figure2, dpi=1000, device= "pdf", height=8, width=15,units="in", bg=NULL, limitsize = F)


#-------------------------------------------------------------------------------
# 2-3cm
#-------------------------------------------------------------------------------

#cohort trace matrix for those who undergo surgery 

#iterative solutions for disease progression
for(t in 1:n.cycles) {
  m_M_noWf3[t+1, ] <- m_M_noWf3[t, ] %*% m_P
}
# check
rowSums(m_M_noWf3)

#cohort distribution graph
m_M_noWf3 <- as.data.frame(m_M_noWf3)
noWf3Stack <- stack(m_M_noWf3)
noWf3Stack$cycle <- rep(0:80, time(6))
colnames(noWf3Stack) <- c("Proportion", "State", "Cycle")


ggplot(noWf3Stack, aes(x = Cycle, y = Proportion))+
  geom_line()+
  facet_wrap(~`State`)+
  #labs(title = col)+
  theme_minimal()

noWf3 <-ggplot(noWf3Stack, aes(x = Cycle, y = Proportion, color = State)) +
           geom_line() +
            theme_minimal()+
               ylim(0,0.5)

#iterative solutions for surgery markov
for(t in 1:n.cycles) {
  m_M_surgery3[t+1, ] <- m_M_surgery3[t, ] %*% m_S
}
# check
rowSums(m_M_surgery3)

m_M_surgery3 <- as.data.frame(m_M_surgery3)
surg3Stack <- stack(m_M_surgery3)
surg3Stack$cycle <- rep(0:80, time(5))
colnames(surg3Stack) <- c("Proportion", "State", "Cycle")


ggplot(surg3Stack, aes(x = Cycle, y = Proportion))+
  geom_line()+
  facet_wrap(~`State`)+
  #labs(title = col)+
  theme_minimal()

surg3 <- ggplot(surg3Stack, aes(x = Cycle, y = Proportion, color = State)) +
             geom_line() +
               theme_minimal()+
                ylim(0,0.2)

#iterative solutions for surveillance + WF markov 
for(t in 1:n.cycles) {
  m_M_surv3[t+1, ] <- m_M_surv3[t, ] %*% m_P
}
# check
rowSums(m_M_surv3)

m_M_surv3 <- as.data.frame(m_M_surv3)
surv3Stack <- stack(m_M_surv3)
surv3Stack$cycle <- rep(0:80, time(6))
colnames(surv3Stack) <- c("Proportion", "State", "Cycle")

ggplot(surv3Stack, aes(x = Cycle, y = Proportion))+
  geom_line()+
  facet_wrap(~`State`)+
  #labs(title = col)+
  theme_minimal()

wf3surv <- ggplot(surv3Stack, aes(x = Cycle, y = Proportion, color = State)) +
            geom_line() +
             theme_minimal()+
              ylim(0,0.02)

combinedH3 <- cbind(m_M_noWf3, m_M_surv3, m_M_surgery3)
write.csv(combinedH3, file="combined_H3.csv", row.names = F)

allcombined <- cbind(combinedH1,combinedH2,combinedH3)
write.csv(allcombined, file="allcombined.csv", row.names = F)

figure1 <-  ggarrange(noWf1, noWf2, noWf3, ncol = 1, nrow = 1)
ggarrange(noWf1, noWf2, noWf3,
          labels = c("A", "B", "C"),
          ncol = 1, nrow = 1)

#-------------------------------------------------------------------------------
# State Rewards
#-------------------------------------------------------------------------------

#costs (annual)
c_H1
c_H2
c_H3
c_WF
c_Cancer
c_Death    <-0

#utilities (annual)
u_H1
u_H2
u_H3
u_WF
u_Cancer
u_Death    <-0



u_SOC <- c(H1 = u_H, H2 =, H3 = , WF = , Cancer = , Death = )


#utilities

#-------------------------------------------------------------------------------
# PSA https://github.com/DARTH-git/cohort-modeling-tutorial-intro/blob/main/R/Functions_cSTM_time_indep.R
# https://github.com/DARTH-git/cohort-modeling-tutorial-intro/blob/main/R/Functions.R
#-------------------------------------------------------------------------------



