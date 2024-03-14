###########
# MATLAB parameters used for deb model
pars <- readMat(paste0('R/DEB/Lampropholis_delicata/results_Lampropholis_delicata.mat'))
par.names <- unlist(labels(pars$par))
for(i in 1:length(par.names)){
  assign(par.names[i], unlist(pars$par[i]))
}
stages <- 4
if(exists("T.L")==FALSE){T_L <- T.L <- 173.15}else{T_L <- T.L}
if(exists("T.H")==FALSE){T_H <- T.H <- 373.15}else{T_H <- T.H}
if(exists("T.AL")==FALSE){T_AL <- T.AL <- 5E04}else{T_AL <- T.AL}
if(exists("T.AH")==FALSE){T_AH <- T.AH <- 9E04}else{T_AH <- T.AH}
if(exists("T.REF")==FALSE){T_REF <- T.REF <- 20 + 273.15}else{T_REF <- T.REF}


######## Grab these from MATLAB OUTPUT
L.b <- 0.4213
E.0 <- 1145

# assign possible missing parameters
if(exists("E.Hj")==FALSE){E.Hj <- E.Hb}
if(exists("E.He")==FALSE){E.He <- E.Hb}
if(exists("L.j")==FALSE){L.j <- L.b}
p.Am <- p.M / kap*z
E.m <- p.Am / v
F.m <- p.Am / kap.X # redefining F.m to max possible value
z.mult <- 1         # DEB body size scaling parameter

# morph, behav and water loss
soilnode <- 7
rundeb <- 1 # run the DEB model?
reset <- 0
metab_mode <- 0

V_i <- z ^ 3
Ww_g <- (E.m * V_i) / mu.E * 23.9 / d.E + V_i
L_i <- z / del.M * 10

# energy content
#E_insect_g <- 23850 # J/g, insects, Buckley 2008 citing Reichle 1971; Andrews and Asato 1977
E_insect_g <- 25600 # J/g, insects, Waldschmidt 1986 cited by Angilletta 2001
V_prey <- 0.082 # cm3 from Pianka 1969, mean of largest 10 items
V_prey <- 0.2 # cm3 from Pianka 1969, largest item

# energy per insect
f_dry.insect <- 0.3 # g dry / g wet
rho.insect <- 1 # g / cm3
W_insect_dry <- V_prey * rho.insect * f_dry.insect # cm3 to g to dry
E_insect <- W_insect_dry * E_insect_g

## Stomach volume
V_stomach <- V_prey # largest prey item
V_stomach <- 0.33E-5 * L_i ^ 2.64 # Avery 1973 Lacerta vivipara 
V_stomach <- 0.048 * Ww_g  ^ 0.99 # Avery 1973 Lacerta vivipara
E_stomach <- E_insect_g * 0.3 * V_stomach
E_sm <- E_stomach / V_i # J / cm3 structure
E_sm <- 2000

d_Egg <- 0.5 # Dry mass fraction of egg (0-1)
pct_H_X <- 70 # Water content of food (%)
pct_H_P <- 30 # Water in faeces (product) (%)
pct_H_R <- 10 # Minimum tolerated dehydration (% of wet mass)
aestivate <- 0
aestdepth <- 5
depress <- 0.7
raindrink <- .1 # mm

K <- 0.05
X_max <- 3
pct_H_X <- 70
X <- 100000
F_m_ref <- p.Am/kap.X * 10000

# egg
V_init <- 3e-9
E_init <- E.0 / V_init
E_H_init <- 0
stage <- 0


# reproduction parameters
viviparous <- 0 # live bearing (1) or egg laying (0)
clutchsize <- 5 # how many eggs per clutch? 
#clutch_ab <- c(0.21354, -4.39046) # regression results from shine paper - in raw_data folder "SVL_Clutchsize.csv" 
photostart <- 3 
photofinish <- 2 
batch <- 1
#clutchsize = 6
#minclutch = 3
#maxclutch = 8

# overwrite nitrogenous waste indices with those of uric acid (currently ammonia by default)
n.NC <- 1
n.NH <- 4/5
n.NO <- 3/5
n.NN <- 4/5

