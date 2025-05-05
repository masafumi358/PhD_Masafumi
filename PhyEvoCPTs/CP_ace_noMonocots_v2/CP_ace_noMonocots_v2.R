#######################################################
# Example Ancestral Character Estimation ("ace")
# on Carnivorous Plant Trap traits
#######################################################

# Change nsim for stochastic mapping = 100 for getting more average (so you dont have to run many times)

# BioGeoBEARS installation guide > https://github.com/nmatzke/BioGeoBEARS

# install.packages('phytools')
# install.packages('maps')

library(ape)
library(maps)
library(phytools)

# set working directory
wd = "~/Github/PhD_Masafumi/PhyEvoCPTs/CP_ace_noMonocots_v2/"
setwd(wd)

# get tree file and character states file for each species on the tree
trfn = "gbotb_tr14_wSister_genera_edited2_minusMonocots.newick"
states_fn = "gbotb_tr14_wSister_genera_edited_states_noMonocots_v2.txt"

# read tree files
tr = read.tree(trfn)

# check the tree (how many tips?)
length(tr$tip.label)

# read discrete data frame
ddf = read.table(states_fn, header=TRUE, sep="\t")
head(ddf)
dim(ddf)

# ensure appropriate sorting
# all state numbers are -1 (0 in data is non-carnivorous, but 1 for non-carnivorous)
# +1 for all states
states = ddf$state
states = c(as.numeric(states) + 1)
names(states) = ddf$species
head(states)

#how many species for each state?
summary(as.factor(states))

# 1    2    3    4    5    6    7    8    9   10   11 
# 1447   36   94    3    1    6   29  144   38   23   58 

#check tip labels
tr$tip.label



#########################################
# Ancestral Character Estimation ("ace")#
#########################################

# see "11state_rate_matrix_noMonocots_v2" excel spreadsheet for how this works!



###################
#Equal rates model#
###################

# all transitions between different character states occur at the same rate
equal_rates_model = matrix(data=1, nrow=11, ncol=11, byrow=TRUE)

#set 0 for diagonal (no transition within the same character)
diag(equal_rates_model) = 0

#check the model
equal_rates_model

# get results, maximum likelihood to estimate the parameters of a discrete character evolution model
resER = fitMk.parallel(tree=tr, x=states, model=equal_rates_model, fixedQ=NULL, ncores=11)

# save results
save(resER, file="resER.Rdata")

# resER_noParallel = fitMk.parallel(tree=tr, x=states, model=equal_rates_model, fixedQ=NULL)
# save(resER_noParallel, file="resER_noParallel.Rdata")
load(file="resER.Rdata")
resER



#################
#Symmetric model#
#################

# all transitions between different character states can occur at different rates,
# but the rate from state i to j is the same as from j to i (i.e., symmetric rates)
symmetric_model = matrix(data=c(0,1,2,3,4,5,6,7,8,9,10,
1,0,11,12,13,14,15,16,17,18,19,
2,11,0,20,21,22,23,24,25,26,27,
3,12,20,0,28,29,30,31,32,33,34,
4,13,21,28,0,35,36,37,38,39,40,
5,14,22,29,35,0,41,42,43,44,45,
6,15,23,30,36,41,0,46,47,48,49,
7,16,24,31,37,42,46,0,50,51,52,
8,17,25,32,38,43,47,50,0,53,54,
9,18,26,33,39,44,48,51,53,0,55,
10,19,27,34,40,45,49,52,54,55,0), nrow=11, ncol=11, byrow=TRUE)

symmetric_model

resSYM = fitMk.parallel(tree=tr, x=states, model=symmetric_model, fixedQ=NULL, ncores=11)
save(resSYM, file="resSYM.Rdata")

#loads to resSYM
load(file="resSYM.Rdata")
resSYM



#################################
#ARD (all-rates-different) model#
#################################

# ARD (All-Rates-Different) model: all transitions between different character states 
# can occur at different rates, and the rate from state i to j can differ from j to i.
ARD_model = matrix(data=c(0,1,2,3,4,5,6,7,8,9,10,
56,0,11,12,13,14,15,16,17,18,19,
57,66,0,20,21,22,23,24,25,26,27,
58,67,75,0,28,29,30,31,32,33,34,
59,68,76,83,0,35,36,37,38,39,40,
60,69,77,84,90,0,41,42,43,44,45,
61,70,78,85,91,96,0,46,47,48,49,
62,71,79,86,92,97,101,0,50,51,52,
63,72,80,87,93,98,102,105,0,53,54,
64,73,81,88,94,99,103,106,108,0,55,
65,74,82,89,95,100,104,107,109,110,0), nrow=11, ncol=11, byrow=TRUE)

ARD_model

resARD = fitMk.parallel(tree=tr, x=states, model=ARD_model, fixedQ=NULL, ncores=11)
save(resARD, file="resARD.Rdata")

#loads to resARD
load(file="resARD.Rdata")
resARD



######################################
#Gain-loss-change (constrained) model#
######################################
GLCC_model = matrix(data=c(0,2,2,2,2,2,2,2,2,2,2,
1,0,3,0,0,0,0,0,0,0,0,
1,3,0,0,0,0,0,0,0,0,0,
1,0,0,0,3,0,0,0,0,0,0,
1,0,0,3,0,0,0,0,0,0,0,
1,0,0,0,0,0,3,0,0,0,0,
1,0,0,0,0,3,0,0,0,0,0,
1,0,0,0,0,0,0,0,3,3,3,
1,0,0,0,0,0,0,3,0,3,3,
1,0,0,0,0,0,0,3,3,0,3,
1,0,0,0,0,0,0,3,3,3,0), nrow=11, ncol=11, byrow=TRUE)

GLCC_model

resGLCC = fitMk.parallel(tree=tr, x=states, model=GLCC_model, fixedQ=NULL, ncores=11)
save(resGLCC, file="resGLCC.Rdata")

load(file="resGLCC.Rdata")
resGLCC



#Gain-loss-change (unconstrained) model
GLCU_model = matrix(data=c(0,2,2,2,2,2,2,2,2,2,2,
1,0,3,3,3,3,3,3,3,3,3,
1,3,0,3,3,3,3,3,3,3,3,
1,3,3,0,3,3,3,3,3,3,3,
1,3,3,3,0,3,3,3,3,3,3,
1,3,3,3,3,0,3,3,3,3,3,
1,3,3,3,3,3,0,3,3,3,3,
1,3,3,3,3,3,3,0,3,3,3,
1,3,3,3,3,3,3,3,0,3,3,
1,3,3,3,3,3,3,3,3,0,3,
1,3,3,3,3,3,3,3,3,3,0), nrow=11, ncol=11, byrow=TRUE)

GLCU_model

resGLCU = fitMk.parallel(tree=tr, x=states, model=GLCU_model, fixedQ=NULL, ncores=11)
save(resGLCU, file="resGLCU.Rdata")

load(file="resGLCU.Rdata")
resGLCU



#Asymmetric Rate Variation by Trapping Zone (ARVT) model
ARVT_model = matrix(data=c(0,2,2,2,2,2,2,2,2,2,2,
1,0,3,4,5,0,3,0,3,4,5,
1,6,0,7,8,6,0,6,0,7,8,
1,9,10,0,11,9,10,9,10,0,11,
1,12,13,14,0,12,13,12,13,14,0,
1,0,3,4,5,0,3,0,3,4,5,
1,6,0,7,8,6,0,6,0,7,8,
1,0,3,4,5,0,3,0,3,4,5,
1,6,0,7,8,6,0,6,0,7,8,
1,9,10,0,11,9,10,9,10,0,11,
1,12,13,14,0,12,13,12,13,14,0), nrow=11, ncol=11, byrow=TRUE)

ARVT_model

resARVT = fitMk.parallel(tree=tr, x=states, model=ARVT_model, fixedQ=NULL, ncores=11)
save(resARVT, file="resARVT.Rdata")

load(file="resARVT.Rdata")
resARVT



#Symmetric Rate Variation by Trapping zone (SRVT) model
SRVT_model = matrix(data=c(0,2,2,2,2,2,2,2,2,2,2,
1,0,3,4,5,0,3,0,3,4,5,
1,3,0,6,7,3,0,3,0,6,7,
1,4,6,0,8,4,6,4,6,0,8,
1,5,7,8,0,5,7,5,7,8,0,
1,0,3,4,5,0,3,0,3,4,5,
1,3,0,6,7,3,0,3,0,6,7,
1,0,3,4,5,0,3,0,3,4,5,
1,3,0,6,7,3,0,3,0,6,7,
1,4,6,0,8,4,6,4,6,0,8,
1,5,7,8,0,5,7,5,7,8,0), nrow=11, ncol=11, byrow=TRUE)

SRVT_model

resSRVT = fitMk.parallel(tree=tr, x=states, model=SRVT_model, fixedQ=NULL, ncores=11)
save(resSRVT, file="resSRVT.Rdata")

load(file="resSRVT.Rdata")
resSRVT



#Gain-loss-change-trapping-zone model
GLCTZ_model = matrix(data=c(0,2,2,2,2,2,2,2,2,2,2,
1,0,0,0,0,3,0,3,0,0,0,
1,0,0,0,0,0,3,0,3,0,0,
1,0,0,0,0,0,0,0,0,3,0,
1,0,0,0,0,0,0,0,0,0,3,
1,3,0,0,0,0,0,3,0,0,0,
1,0,3,0,0,0,0,0,3,0,0,
1,3,0,0,0,3,0,0,0,0,0,
1,0,3,0,0,0,3,0,0,0,0,
1,0,0,3,0,0,0,0,0,0,0,
1,0,0,0,3,0,0,0,0,0,0), nrow=11, ncol=11, byrow=TRUE)

GLCTZ_model

resGLCTZ = fitMk.parallel(tree=tr, x=states, model=GLCTZ_model, fixedQ=NULL, ncores=11)
save(resGLCTZ, file="resGLCTZ.Rdata")

load(file="resGLCTZ.Rdata")
resGLCTZ



#####################################
#Complex trap evolution - Jump model#
#####################################
CTE_Jump_model = matrix(data=c(0,2,2,0,0,0,0,0,0,0,0,
1,0,0,5,0,0,0,3,0,0,0,
1,0,0,0,6,0,0,0,4,0,0,
1,5,0,0,0,0,0,0,0,0,0,
1,0,6,0,0,0,0,0,0,0,0,
1,0,0,0,0,0,0,0,7,0,0,
1,0,0,0,0,0,0,0,0,8,0,
1,3,0,0,0,0,0,0,0,9,0,
1,0,4,0,0,7,0,0,0,0,10,
1,0,0,0,0,0,8,9,0,0,0,
1,0,0,0,0,0,0,0,10,0,0), nrow=11, ncol=11, byrow=TRUE)

CTE_Jump_model

resCTEJ = fitMk.parallel(tree=tr, x=states, model=CTE_Jump_model, fixedQ=NULL, ncores=11)
save(resCTEJ, file="resCTEJ.Rdata")

load(file="resCTEJ.Rdata")
resCTEJ



##############################################################
#Complex trap evolution model for the origin of bladder traps#
##############################################################
Complex_trap_evolution_model = matrix(data=c(0,2,2,0,0,0,0,0,0,0,0,
1,0,3,0,0,11,0,0,0,0,0,
1,3,0,4,0,0,9,0,0,0,0,
1,0,0,0,5,0,0,0,0,0,0,
1,0,0,0,0,0,0,0,0,0,0,
1,0,0,0,0,0,0,12,0,0,0,
1,0,0,0,0,0,0,0,10,0,0,
1,0,0,0,0,0,0,0,6,0,0,
1,0,0,0,0,0,0,6,0,7,0,
1,0,0,0,0,0,0,0,0,0,8,
1,0,0,0,0,0,0,0,0,0,0), nrow=11, ncol=11, byrow=TRUE)

Complex_trap_evolution_model

resCTE = fitMk.parallel(tree=tr, x=states, model=Complex_trap_evolution_model, fixedQ=NULL, ncores=11)
save(resCTE, file="resCTE.Rdata")

load(file="resCTE.Rdata")
resCTE



#########################################################################
#Reversible complex trap evolution model for the origin of bladder traps#
#########################################################################
Reversible_complex_trap_evolution_model = matrix(data=c(0,2,2,0,0,0,0,0,0,0,0,
1,0,3,0,0,11,0,0,0,0,0,
1,3,0,4,0,0,9,0,0,0,0,
1,0,0,0,5,0,0,0,0,0,0,
1,0,0,0,0,0,0,0,0,0,0,
1,11,0,0,0,0,0,12,0,0,0,
1,0,9,0,0,0,0,0,10,0,0,
1,0,0,0,0,12,0,0,6,0,0,
1,0,0,0,0,0,10,6,0,7,0,
1,0,0,0,0,0,0,0,0,0,8,
1,0,0,0,0,0,0,0,0,0,0), nrow=11, ncol=11, byrow=TRUE)

Reversible_complex_trap_evolution_model

resRCTE = fitMk.parallel(tree=tr, x=states, model=Reversible_complex_trap_evolution_model, fixedQ=NULL, ncores=11)
save(resRCTE, file="resRCTE.Rdata")

load(file="resRCTE.Rdata")
resRCTE



###################################################################################
#6 rates (one-way ground transition) complex trap evolution model#
###################################################################################
sixG_rates_complex_trap_evolution_model = matrix(data=c(0,2,2,0,0,0,0,0,0,0,0,
1,0,3,0,0,11,0,0,0,0,0,
1,3,0,4,0,0,9,0,0,0,0,
1,0,0,0,5,0,0,0,0,0,0,
1,0,0,0,0,0,0,0,0,0,0,
1,13,0,0,0,0,0,12,0,0,0,
1,0,0,0,0,0,0,0,10,0,0,
1,0,0,0,0,14,0,0,6,0,0,
1,0,0,0,0,0,0,6,0,7,0,
1,0,0,0,0,0,0,0,0,0,8,
1,0,0,0,0,0,0,0,0,0,0), nrow=11, ncol=11, byrow=TRUE)

sixG_rates_complex_trap_evolution_model

res6GCTE = fitMk.parallel(tree=tr, x=states, model=sixG_rates_complex_trap_evolution_model, fixedQ=NULL, ncores=11)
save(res6GCTE, file="res6GCTE.Rdata")

load(file="res6GCTE.Rdata")
res6GCTE



###################################################################################
#6 rates (one-way aerial transition) complex trap evolution model#
###################################################################################
sixA_rates_complex_trap_evolution_model = matrix(data=c(0,2,2,0,0,0,0,0,0,0,0,
                                                        1,0,3,0,0,11,0,0,0,0,0,
                                                        1,3,0,4,0,0,9,0,0,0,0,
                                                        1,0,0,0,5,0,0,0,0,0,0,
                                                        1,0,0,0,0,0,0,0,0,0,0,
                                                        1,0,0,0,0,0,0,12,0,0,0,
                                                        1,0,13,0,0,0,0,0,10,0,0,
                                                        1,0,0,0,0,0,0,0,6,0,0,
                                                        1,0,0,0,0,0,14,6,0,7,0,
                                                        1,0,0,0,0,0,0,0,0,0,8,
                                                        1,0,0,0,0,0,0,0,0,0,0), nrow=11, ncol=11, byrow=TRUE)

sixA_rates_complex_trap_evolution_model

res6ACTE = fitMk.parallel(tree=tr, x=states, model=sixA_rates_complex_trap_evolution_model, fixedQ=NULL, ncores=11)
save(res6ACTE, file="res6ACTE.Rdata")

load(file="res6ACTE.Rdata")
res6ACTE



###################################################################################
#8 rates complex trap evolution model#
###################################################################################
eight_rates_complex_trap_evolution_model = matrix(data=c(0,2,2,0,0,0,0,0,0,0,0,
                                                         1,0,3,0,0,11,0,0,0,0,0,
                                                         1,3,0,4,0,0,9,0,0,0,0,
                                                         1,0,0,0,5,0,0,0,0,0,0,
                                                         1,0,0,0,0,0,0,0,0,0,0,
                                                         1,15,0,0,0,0,0,12,0,0,0,
                                                         1,0,13,0,0,0,0,0,10,0,0,
                                                         1,0,0,0,0,16,0,0,6,0,0,
                                                         1,0,0,0,0,0,14,6,0,7,0,
                                                         1,0,0,0,0,0,0,0,0,0,8,
                                                         1,0,0,0,0,0,0,0,0,0,0), nrow=11, ncol=11, byrow=TRUE)

eight_rates_complex_trap_evolution_model

res8CTE = fitMk.parallel(tree=tr, x=states, model=eight_rates_complex_trap_evolution_model, fixedQ=NULL, ncores=11)
save(res8CTE, file="res8CTE.Rdata")

load(file="res8CTE.Rdata")
res8CTE



###################################################################################
#7 rates (aerial: both-pitcher or reverse) complex trap evolution model#
###################################################################################
seven_rates_abpr_complex_trap_evolution_model = matrix(data=c(0,2,2,0,0,0,0,0,0,0,0,
                                                              1,0,3,0,0,11,0,0,0,0,0,
                                                              1,3,0,4,0,0,9,0,0,0,0,
                                                              1,0,0,0,5,0,0,0,0,0,0,
                                                              1,0,0,0,0,0,0,0,0,0,0,
                                                              1,0,0,0,0,0,0,12,0,0,0,
                                                              1,0,13,0,0,0,0,0,10,0,0,
                                                              1,0,0,0,0,15,0,0,6,0,0,
                                                              1,0,0,0,0,0,14,6,0,7,0,
                                                              1,0,0,0,0,0,0,0,0,0,8,
                                                              1,0,0,0,0,0,0,0,0,0,0), nrow=11, ncol=11, byrow=TRUE)

seven_rates_abpr_complex_trap_evolution_model

res7abprCTE = fitMk.parallel(tree=tr, x=states, model=seven_rates_abpr_complex_trap_evolution_model, fixedQ=NULL, ncores=11)
save(res7abprCTE, file="res7abprCTE.Rdata")

load(file="res7abprCTE.Rdata")
res7abprCTE



###################################################################################
#7 rates (aerial: sticky-both or reverse) complex trap evolution model#
###################################################################################
seven_rates_asbr_complex_trap_evolution_model = matrix(data=c(0,2,2,0,0,0,0,0,0,0,0,
                                                              1,0,3,0,0,11,0,0,0,0,0,
                                                              1,3,0,4,0,0,9,0,0,0,0,
                                                              1,0,0,0,5,0,0,0,0,0,0,
                                                              1,0,0,0,0,0,0,0,0,0,0,
                                                              1,15,0,0,0,0,0,12,0,0,0,
                                                              1,0,13,0,0,0,0,0,10,0,0,
                                                              1,0,0,0,0,0,0,0,6,0,0,
                                                              1,0,0,0,0,0,14,6,0,7,0,
                                                              1,0,0,0,0,0,0,0,0,0,8,
                                                              1,0,0,0,0,0,0,0,0,0,0), nrow=11, ncol=11, byrow=TRUE)

seven_rates_asbr_complex_trap_evolution_model

res7asbrCTE = fitMk.parallel(tree=tr, x=states, model=seven_rates_asbr_complex_trap_evolution_model, fixedQ=NULL, ncores=11)
save(res7asbrCTE, file="res7asbrCTE.Rdata")

load(file="res7asbrCTE.Rdata")
res7asbrCTE



###################################################################################
#7 rates (ground: both-pitcher or reverse) complex trap evolution model#
###################################################################################
seven_rates_gbpr_complex_trap_evolution_model = matrix(data=c(0,2,2,0,0,0,0,0,0,0,0,
                                                              1,0,3,0,0,11,0,0,0,0,0,
                                                              1,3,0,4,0,0,9,0,0,0,0,
                                                              1,0,0,0,5,0,0,0,0,0,0,
                                                              1,0,0,0,0,0,0,0,0,0,0,
                                                              1,13,0,0,0,0,0,12,0,0,0,
                                                              1,0,0,0,0,0,0,0,10,0,0,
                                                              1,0,0,0,0,14,0,0,6,0,0,
                                                              1,0,0,0,0,0,15,6,0,7,0,
                                                              1,0,0,0,0,0,0,0,0,0,8,
                                                              1,0,0,0,0,0,0,0,0,0,0), nrow=11, ncol=11, byrow=TRUE)

seven_rates_gbpr_complex_trap_evolution_model

res7gbprCTE = fitMk.parallel(tree=tr, x=states, model=seven_rates_gbpr_complex_trap_evolution_model, fixedQ=NULL, ncores=11)
save(res7gbprCTE, file="res7gbprCTE.Rdata")

load(file="res7gbprCTE.Rdata")
res7gbprCTE



###################################################################################
#7 rates (ground: both-sticky or reverse) complex trap evolution model#
###################################################################################
seven_rates_gsbr_complex_trap_evolution_model = matrix(data=c(0,2,2,0,0,0,0,0,0,0,0,
                                                              1,0,3,0,0,11,0,0,0,0,0,
                                                              1,3,0,4,0,0,9,0,0,0,0,
                                                              1,0,0,0,5,0,0,0,0,0,0,
                                                              1,0,0,0,0,0,0,0,0,0,0,
                                                              1,13,0,0,0,0,0,12,0,0,0,
                                                              1,0,15,0,0,0,0,0,10,0,0,
                                                              1,0,0,0,0,14,0,0,6,0,0,
                                                              1,0,0,0,0,0,0,6,0,7,0,
                                                              1,0,0,0,0,0,0,0,0,0,8,
                                                              1,0,0,0,0,0,0,0,0,0,0), nrow=11, ncol=11, byrow=TRUE)

seven_rates_gsbr_complex_trap_evolution_model

res7gsbrCTE = fitMk.parallel(tree=tr, x=states, model=seven_rates_gsbr_complex_trap_evolution_model, fixedQ=NULL, ncores=11)
save(res7gsbrCTE, file="res7gsbrCTE.Rdata")

load(file="res7gsbrCTE.Rdata")
res7gsbrCTE

# see "Model_Selection_AIC_noMonocots.xlsx" for statistical model comparison (AIC)
# the best model is "res7abprCTE"
# compare the model with ER model

#################################
#Stochastic mapping for ER model#
#################################


set.seed(34321)
stochastic_maps_ER= make.simmap(tree=tr, x=states, model=equal_rates_model, nsim=100)
summary_stochastic_maps_ER = summary(stochastic_maps_ER)

save(summary_stochastic_maps_ER, file="stochasticER.Rdata")
load(file="stochasticER.Rdata")
summary_stochastic_maps_ER

names(summary_stochastic_maps_ER)

summary_stochastic_maps_ER$count

# ace = ancestral character estimation
ancstates = summary_stochastic_maps_ER$ace

pdf(file="plot_ACE_ER.pdf", width=12, height=40)

cols = c("white","lightblue","blue", "yellow","orange",
         "orange3","red","lightgrey", "darkgrey", "green3", "darkgreen")

plot.phylo(tr, show.tip.label=FALSE)

# get the internal node numbers
tipnode_nums = 1:length(tr$tip.label)
internal_nodenums = (length(tr$tip.label)+1):(length(tr$tip.label)+tr$Nnode)

# Cut to just carnivorous nodes
carnivorous_TF = ancstates[,"1"] < 0.50
internal_nodenums_CPs = internal_nodenums[carnivorous_TF]
ancstates_CPs = ancstates[carnivorous_TF,]

nodelabels(node=internal_nodenums_CPs, pie=ancstates_CPs, piecol=cols, cex=0.5)

dev.off()
system("open plot_ACE_ER.pdf")


#Getting mean counts on phylogenetic stochastic mapping
#means of each column
apply(X=summary_stochastic_maps_ER$count, MARGIN=2, FUN=mean)

mean_counts <- apply(X=summary_stochastic_maps_ER$count[,2:111], MARGIN=2, FUN=mean)

mean_counts

max(mean_counts)

sort_means <- sort(mean_counts, decreasing = TRUE)

top_10_means <- head(sort_means, 10)
top_10_means

# Standard deviations
apply(X=summary_stochastic_maps_ER$count, MARGIN=2, FUN=sd)

# 95% confidence intervals
1.96 * apply(X=summary_stochastic_maps_ER$count, MARGIN=2, FUN=sd)
-1.96 * apply(X=summary_stochastic_maps_ER$count, MARGIN=2, FUN=sd)

# 95% credible with quantiles
apply(X=summary_stochastic_maps_ER$count [,2:111], MARGIN=2, FUN=quantile, prob=c(0.025, 0.975))



#####################
###7abprCTE_v2 map###
#####################

runslow = FALSE
if (runslow == TRUE)
{
  set.seed(34321)
  
  trtable = BioGeoBEARS::prt(tr)
  save(trtable, file="gbotb_tr14_wSister_genera_edited2_minusMonocots_trtable.Rdata")
  
  # ~1 hour to run
  stochastic_maps_7abprCTE = make.simmap(tree=tr, x=states, model=seven_rates_abpr_complex_trap_evolution_model, nsim=100)
  # save the full stochastic maps
  save(stochastic_maps_7abprCTE, file="stochastic_maps_7abprCTE.Rdata")
  
  # make and save a summary as well 
  # several minutes
  summary_stochastic_maps_7abprCTE = summary(stochastic_maps_7abprCTE)
  save(summary_stochastic_maps_7abprCTE, file="summary_stochastic_maps_7abprCTE.Rdata")
} else {
  # loads to: trtable
  load(file="gbotb_tr14_wSister_genera_edited2_minusMonocots_trtable.Rdata")
  
  # loads to: stochastic_maps_7abprCTE
  load(file="stochastic_maps_7abprCTE.Rdata")	
  # loads to: stochastic_maps_7abprCTE
  load(file="summary_stochastic_maps_7abprCTE.Rdata")	
} # END runslow


summary_stochastic_maps_7abprCTE$count
names(summary_stochastic_maps_7abprCTE)


# ace = ancestral character estimation
ancstates = summary_stochastic_maps_7abprCTE$ace


# lca_of_Lentibulariaceae
2205

# ancnode_of_lca_Lentibulariaceae
2114

ancstates[2205,]
ancstates[2114,]

# Row of prt table
trtable[2205,]
#    node ord_ndname node_lvl node.type parent_br edge.length ancestor daughter_nds node_ht  time_bp fossils label
#2205 2205       2205       12  internal       640    10.53411     2114   2206, 2285 90.8475 32.88674      NA      

# edge/branch number:
640

stochastic_maps_7abprCTE[[1]]$maps[640]
stochastic_maps_7abprCTE[[2]]$maps[640]
stochastic_maps_7abprCTE[[3]]$maps[640]
stochastic_maps_7abprCTE[[4]]$maps[640]
stochastic_maps_7abprCTE[[5]]$maps[640]














pdf(file="plot_ACE_7abprCTE_v2.pdf", width=12, height=40)

cols = c("white","lightblue","blue", "yellow","orange",
         "orange3","red","lightgrey", "darkgrey", "green3", "darkgreen")

# get the internal node numbers
tipnode_nums = 1:length(tr$tip.label)
internal_nodenums = (length(tr$tip.label)+1):(length(tr$tip.label)+tr$Nnode)

# Cut to just carnivorous nodes
carnivorous_TF = ancstates[,"1"] < 0.50
CPnums = (1:length(carnivorous_TF))[carnivorous_TF]
CPnums
#internal_nodenums_CPs = internal_nodenums[carnivorous_TF]
#ancstates_CPs = ancstates[carnivorous_TF,]



# get the node ancestor of the LCA of each clade
list_of_ancnodes_of_lca = NULL
# Lentibulariaceae
lca_of_Lentibulariaceae = getMRCA(phy=tr, tip=c("Pinguicula_alpina","Utricularia_albiflora"))
lca_of_Lentibulariaceae
edgeTF = tr$edge[,2] == lca_of_Lentibulariaceae
edgenum = (1:(length(edgeTF)))[edgeTF]
edgenum
ancnode_of_lca_Lentibulariaceae = tr$edge[edgeTF,1]
ancnode_of_lca_Lentibulariaceae
tr$edge[edgeTF,]


# Nepentheaceae-Droseraceae (to see the transtion from non-carnivorous to sticky traps)
lca_of_NeDro = getMRCA(phy=tr, tip = c("Drosera_lanata", "Nepenthes_gracilis"))
lca_of_NeDro

edgeTF_NeDro = tr$edge[,2] == lca_of_NeDro
edgeTF_NeDro
edgenum_NeDro = (1:(length(edgeTF_NeDro)))[edgeTF_NeDro]
edgenum_NeDro
ancnode_of_lca_NeDro = tr$edge[edgeTF_NeDro,1]
ancnode_of_lca_NeDro
tr$edge[edgeTF_NeDro,]


# Nepentheaceae (see transition from sticky to pitcher throught ancestral transitional state)
#trtable$tipnames[1583:1732]
lca_of_Nepentheaceae = getMRCA(phy=tr, tip = c("Nepenthes_madagascariensis", "Nepenthes_pervillei"))
lca_of_Nepentheaceae

edgeTF_Nepentheaceae = tr$edge[,2] == lca_of_Nepentheaceae
edgeTF_Nepentheaceae
edgenum_Nepentheaceae = (1:(length(edgeTF_Nepentheaceae)))[edgeTF_Nepentheaceae]
edgenum_Nepentheaceae
ancnode_of_lca_Nepentheaceae = tr$edge[edgeTF_Nepentheaceae,1]
ancnode_of_lca_Nepentheaceae
tr$edge[edgeTF_Nepentheaceae,]


#Sarraceniaceae
lca_of_Sarrac = getMRCA(phy=tr, tip = c("Sarracenia_rosea", "Darlingtonia_californica"))
lca_of_Sarrac

edgeTF_Sarrac = tr$edge[,2] == lca_of_Sarrac
edgeTF_Sarrac
edgenum_Sarrac = (1:(length(edgeTF_Sarrac)))[edgeTF_Sarrac]
edgenum_Sarrac
ancnode_of_lca_Sarrac = tr$edge[edgeTF_Sarrac,1]
ancnode_of_lca_Sarrac
tr$edge[edgeTF_Sarrac,]


#Droseraceae
lca_of_Droseraceae = getMRCA(phy=tr, tip = c("Drosera_arcturi", "Aldrovanda_vesiculosa"))
lca_of_Droseraceae

edgeTF_Droseraceae = tr$edge[,2] == lca_of_Droseraceae
edgeTF_Droseraceae
edgenum_Droseraceae = (1:(length(edgeTF_Droseraceae)))[edgeTF_Droseraceae]
edgenum_Droseraceae
ancnode_of_lca_Droseraceae = tr$edge[edgeTF_Droseraceae,1]
ancnode_of_lca_Droseraceae
tr$edge[edgeTF_Droseraceae,]


# Lentibulariaceae-transitional-pitcher-eel
lca_of_Lentibulariaceae_2 = getMRCA(phy=tr, tip=c("Genlisea_uncinata","Utricularia_uniflora"))
lca_of_Lentibulariaceae_2
edgeTF_2 = tr$edge[,2] == lca_of_Lentibulariaceae_2
edgenum_2 = (1:(length(edgeTF_2)))[edgeTF_2]
edgenum_2
ancnode_of_lca_Lentibulariaceae_2 = tr$edge[edgeTF_2,1]
ancnode_of_lca_Lentibulariaceae_2
tr$edge[edgeTF_2,]

# make sure the secondarily noncarnivorous relatives
# of Droseraceae and Nepenthes are in there:
TF = grepl(tr$tip.label, pattern="Triphyophyllum")
num = (1:length(TF))[TF]
Triphyophyllum_num = num
tr$tip.label[1725:1775]
# Ancistrocladus_likoko
# Habropetalum_dawei
# Triphyophyllum_peltatum
getMRCA(phy=tr, tip=c("Ancistrocladus_likoko", "Habropetalum_dawei"))
getMRCA(phy=tr, tip=c("Ancistrocladus_likoko", "Triphyophyllum_peltatum"))
getMRCA(phy=tr, tip=c("Habropetalum_dawei", "Triphyophyllum_peltatum"))

nodes_of_nonCarnivorous_Droseraceae_relatives = BioGeoBEARS::get_daughter_nodes(nodenum=3615, tr=tr, nodes=NULL)
nodes_of_nonCarnivorous_Droseraceae_relatives

# Remove Triphyophyllum_num & put last, so it plots better
nodes_of_nonCarnivorous_Droseraceae_relatives = nodes_of_nonCarnivorous_Droseraceae_relatives[nodes_of_nonCarnivorous_Droseraceae_relatives != Triphyophyllum_num]
nodes_of_nonCarnivorous_Droseraceae_relatives = c(nodes_of_nonCarnivorous_Droseraceae_relatives, Triphyophyllum_num)


list_of_ancnodes_of_lca = c(ancnode_of_lca_Lentibulariaceae, ancnode_of_lca_NeDro, ancnode_of_lca_Nepentheaceae, ancnode_of_lca_Sarrac, ancnode_of_lca_Droseraceae, ancnode_of_lca_Lentibulariaceae_2, nodes_of_nonCarnivorous_Droseraceae_relatives)

length(internal_nodenums)
length(carnivorous_TF)

ancstates[1:50,]
ancstates[1850:1900,] # <- changes to tips here
ancstates[3700:3717,]

ancstates[carnivorous_TF,][1:50,]
ancstates[carnivorous_TF,][410:440,] # <- changes to tips here
ancstates[carnivorous_TF,][830:860,]


# these issues lead to this having NA at end:
internal_nodenums[carnivorous_TF]


# re-order ancstates, makes ancstates2 with
# standard APE node order
tip_nodes = 1:length(tr$tip.label)
internal_nodes = (length(tr$tip.label)+1):(length(tr$tip.label)+tr$Nnode)
all_nodes = c(tip_nodes, internal_nodes)
ancstates_internal = summary_stochastic_maps_7abprCTE$ace[1:length(internal_nodes),]
ancstates_tips = summary_stochastic_maps_7abprCTE$ace[(length(internal_nodes)+1):nrow(summary_stochastic_maps_7abprCTE$ace),]
# These should have APE order
ancstates2 = rbind(ancstates_tips, ancstates_internal)

# Check ancstates2
ancstates2[1:50,]
ancstates2[1850:1900,]
ancstates2[3700:3717,]

carnivorous_TF2 = ancstates2[,"1"] < 0.50

internal_nodenums_CPs = c(all_nodes[carnivorous_TF2], list_of_ancnodes_of_lca)

ancstates_CPs2 = rbind(ancstates2[carnivorous_TF2,], ancstates2[list_of_ancnodes_of_lca,])

ancstates2[carnivorous_TF,][1:100,]
ancstates2[carnivorous_TF,][101:200,]
ancstates2[carnivorous_TF,][201:300,]
ancstates2[carnivorous_TF,][301:400,]
ancstates2[carnivorous_TF,][401:500,]
ancstates2[carnivorous_TF,][501:600,]
ancstates2[carnivorous_TF,][601:700,]
ancstates2[carnivorous_TF,][701:800,]
ancstates2[carnivorous_TF,][750:860,]



#ancstates_CPs = ancstates[carnivorous_TF,]

plot.phylo(tr, show.tip.label=FALSE)

nodelabels(node=internal_nodenums_CPs, pie=ancstates_CPs2, piecol=cols, cex=0.5)

dev.off()
system("open plot_ACE_7abprCTE_v2.pdf")








###########################
##State Distribution Plot##
###########################

# look at Lentibulariaceae ancestral branch stochastic maps
dev.off()
edgenum

list_of_maps_for_Lentibulariaceae_ancestral_branch = list()
for (i in 1:length(stochastic_maps_7abprCTE))
  {
  list_of_maps_for_Lentibulariaceae_ancestral_branch[[i]] = stochastic_maps_7abprCTE[[i]]$maps[edgenum][[1]]
  }
list_of_maps_for_Lentibulariaceae_ancestral_branch

# number of steps in each of the 100 stochastic maps on this branch
sapply(X=list_of_maps_for_Lentibulariaceae_ancestral_branch, FUN=length)

# let's get percentages of each state at several timepoints along a branch/edge
branchlength = tr$edge.length[edgenum]
branchlength

# make a list of fractions along the branch
fractions_of_branch2 = seq(from=0, to=1, by=0.01)

# for a mapped branch, get the CUMULATIVE time when a change occurs
cumsum(list_of_maps_for_Lentibulariaceae_ancestral_branch[[1]])

# table of states for each timepoint along branch
tmpmat = matrix(data=NA, nrow=length(stochastic_maps_7abprCTE), ncol=length(fractions_of_branch2))
dim(tmpmat)

# fill out the table
for (i in 1:nrow(tmpmat))
  {
  for (j in 1:ncol(tmpmat))
    {
    fraction_time = fractions_of_branch2[j]
    time_on_branch = branchlength * fraction_time
    time_on_branch
    cumulative_times_mapped = cumsum(list_of_maps_for_Lentibulariaceae_ancestral_branch[[i]])
    
    # which time bin are we in:
    # time_on_branch = 10.5
    lessTF = time_on_branch <= cumulative_times_mapped
    lessTF
    
    # if all FALSE, state is the last state
    if (all(lessTF==FALSE))
    {
      lessTF[length(lessTF)] = TRUE
    }
    state_residences = 1:length(cumulative_times_mapped)
    state_residences
    state_residences[lessTF]
    the_state_at_the_time = min(state_residences[lessTF])
    the_state_at_the_time
    
    tmpmat[i,j] = names(cumulative_times_mapped)[the_state_at_the_time]
    }
  }
tmpmat

# Get list of unique states in stochastic maps
uniq_states = sort(unique(c(tmpmat)))
uniq_states

# Get percentages of each, for each timepoint
percentages_by_timepoint = matrix(data=NA, ncol=length(fractions_of_branch2), nrow=length(uniq_states))
percentages_by_timepoint

for (i in 1:nrow(percentages_by_timepoint))
{
  for (j in 1:ncol(percentages_by_timepoint))
  {
    TF = tmpmat[,j] == uniq_states[i]
    percentages_by_timepoint[i,j] = sum(TF)
  }
}

percentages_by_timepoint

# Barplot showing percentages of each state at 100 timepoints along the branch
# below Lentibulariaceae
statecols = c("grey50","orange","yellow2","grey90")
barplot(percentages_by_timepoint[4:1,], col=statecols,border=NA, space=0)

# Density plots for each state
timepoints = fractions_of_branch2*branchlength
plot(x=timepoints, y=1:101, main="State Distribution Over Time Along the Lentibulariaceae Ancestral Branch", xlab="Time along branch (units of length)", ylab="Percentage of each state", pch=".", col="white")
for (i in 1:nrow(percentages_by_timepoint))
{
  lines(x=timepoints,y=percentages_by_timepoint[i,], lwd = 3, col=rev(statecols)[i])
}

names(stochastic_maps_7abprCTE[[1]])






# Look at Lentibulariaceae ancestral branch stochastic maps
dev.off()
edgenum_2

list_of_maps_for_Lentibulariaceae_ancestral_branch = list()
for (i in 1:length(stochastic_maps_7abprCTE))
{
  list_of_maps_for_Lentibulariaceae_ancestral_branch[[i]] = stochastic_maps_7abprCTE[[i]]$maps[edgenum_2][[1]]
}
list_of_maps_for_Lentibulariaceae_ancestral_branch

# Number of steps in each of the 100 stochastic maps on this branch
sapply(X=list_of_maps_for_Lentibulariaceae_ancestral_branch, FUN=length)

# Let's get percentages of each state at several timepoints along a branch/edge
branchlength = tr$edge.length[edgenum_2]
branchlength

# Make a list of fractions along the branch
fractions_of_branch2 = seq(from=0, to=1, by=0.01)

# For a mapped branch, get the CUMULATIVE time when a change occurs
cumsum(list_of_maps_for_Lentibulariaceae_ancestral_branch[[1]])

# Table of states for each timepoint along branch
tmpmat = matrix(data=NA, nrow=length(stochastic_maps_7abprCTE), ncol=length(fractions_of_branch2))
dim(tmpmat)

# Fill out the table
for (i in 1:nrow(tmpmat))
{
  for (j in 1:ncol(tmpmat))
  {
    fraction_time = fractions_of_branch2[j]
    time_on_branch = branchlength * fraction_time
    time_on_branch
    cumulative_times_mapped = cumsum(list_of_maps_for_Lentibulariaceae_ancestral_branch[[i]])
    
    # Which time bin are we in:
    #time_on_branch = 10.5
    lessTF = time_on_branch <= cumulative_times_mapped
    lessTF
    
    # If all FALSE, state is the last state
    if (all(lessTF==FALSE))
    {
      lessTF[length(lessTF)] = TRUE
    }
    state_residences = 1:length(cumulative_times_mapped)
    state_residences
    state_residences[lessTF]
    the_state_at_the_time = min(state_residences[lessTF])
    the_state_at_the_time
    
    tmpmat[i,j] = names(cumulative_times_mapped)[the_state_at_the_time]
  }
}
tmpmat

# Get list of unique states in stochastic maps
uniq_states = sort(unique(c(tmpmat)))
uniq_states

# Get percentages of each, for each timepoint
percentages_by_timepoint = matrix(data=NA, ncol=length(fractions_of_branch2), nrow=length(uniq_states))
percentages_by_timepoint

for (i in 1:nrow(percentages_by_timepoint))
{
  for (j in 1:ncol(percentages_by_timepoint))
  {
    TF = tmpmat[,j] == uniq_states[i]
    percentages_by_timepoint[i,j] = sum(TF)
  }
}

percentages_by_timepoint

# Barplot showing percentages of each state at 100 timepoints along the branch
# below Lentibulariaceae
statecols = c("darkgreen","grey50","orange","lightblue")
barplot(percentages_by_timepoint[4:1,], col=statecols,border=NA, space=0)

# Density plots for each state
timepoints = fractions_of_branch2*branchlength
plot(x=timepoints, y=1:101, main="State Distribution Over Time Along the Lentibulariaceae Ancestral Branch (Transitional-Pitcher-Eel)", xlab="Time along branch (units of length)", ylab="Percentage of each state", pch=".", col="white")
for (i in 1:nrow(percentages_by_timepoint))
{
  lines(x=timepoints,y=percentages_by_timepoint[i,], lwd = 3, col=rev(statecols)[i])
}

names(stochastic_maps_7abprCTE[[1]])





#######################################################################
# Look at Droceraceae & Nepenthaceae ancestral branch stochastic maps #
#######################################################################
dev.off()
edgenum_NeDro

list_of_maps_for_Nepenthaceae_Droceraceae_ancestral_branch = list()
for (i in 1:length(stochastic_maps_7abprCTE))
{
  list_of_maps_for_Nepenthaceae_Droceraceae_ancestral_branch[[i]] = stochastic_maps_7abprCTE[[i]]$maps[edgenum_NeDro][[1]]
}
list_of_maps_for_Nepenthaceae_Droceraceae_ancestral_branch

# Number of steps in each of the 100 stochastic maps on this branch
sapply(X=list_of_maps_for_Nepenthaceae_Droceraceae_ancestral_branch, FUN=length)

# Let's get percentages of each state at several timepoints along a branch/edge
branchlength = tr$edge.length[edgenum_NeDro]
branchlength

# Make a list of fractions along the branch
fractions_of_branch2 = seq(from=0, to=1, by=0.01)

# For a mapped branch, get the CUMULATIVE time when a change occurs
cumsum(list_of_maps_for_Nepenthaceae_Droceraceae_ancestral_branch[[1]])

# Table of states for each timepoint along branch
tmpmat = matrix(data=NA, nrow=length(stochastic_maps_7abprCTE), ncol=length(fractions_of_branch2))
dim(tmpmat)

# Fill out the table
for (i in 1:nrow(tmpmat))
{
  for (j in 1:ncol(tmpmat))
  {
    fraction_time = fractions_of_branch2[j]
    time_on_branch = branchlength * fraction_time
    time_on_branch
    cumulative_times_mapped = cumsum(list_of_maps_for_Nepenthaceae_Droceraceae_ancestral_branch[[i]])
    
    # Which time bin are we in:
    #time_on_branch = 10.5
    lessTF = time_on_branch <= cumulative_times_mapped
    lessTF
    
    # If all FALSE, state is the last state
    if (all(lessTF==FALSE))
    {
      lessTF[length(lessTF)] = TRUE
    }
    state_residences = 1:length(cumulative_times_mapped)
    state_residences
    state_residences[lessTF]
    the_state_at_the_time = min(state_residences[lessTF])
    the_state_at_the_time
    
    tmpmat[i,j] = names(cumulative_times_mapped)[the_state_at_the_time]
  }
}
tmpmat

# Get list of unique states in stochastic maps
uniq_states = sort(unique(c(tmpmat)))
uniq_states

# Get percentages of each, for each timepoint
percentages_by_timepoint = matrix(data=NA, ncol=length(fractions_of_branch2), nrow=length(uniq_states))
percentages_by_timepoint

for (i in 1:nrow(percentages_by_timepoint))
{
  for (j in 1:ncol(percentages_by_timepoint))
  {
    TF = tmpmat[,j] == uniq_states[i]
    percentages_by_timepoint[i,j] = sum(TF)
  }
}

percentages_by_timepoint

# Barplot showing percentages of each state at 100 timepoints along the branch
# below Lentibulariaceae
statecols = c("lightgrey","orange","yellow3","grey90")
barplot(percentages_by_timepoint[4:1,], col=statecols,border=NA, space=0)

# Density plots for each state
timepoints = fractions_of_branch2*branchlength
plot(x=timepoints, y=1:101, main="State Distribution Over Time Along the Nepenthaceae+Droseraceae Ancestral Branch", xlab="Time along branch (units of length)", ylab="Percentage of each state", pch=".", col="white")
for (i in 1:nrow(percentages_by_timepoint))
{
  lines(x=timepoints,y=percentages_by_timepoint[i,], lwd = 3, col=rev(statecols)[i])
}

names(stochastic_maps_7abprCTE[[1]])




#######################################################################
# Look at Sarraceniaceae ancestral branch stochastic maps #
#######################################################################
dev.off()
edgenum_Sarrac

list_of_maps_for_Sarraceniaceae_ancestral_branch = list()
for (i in 1:length(stochastic_maps_7abprCTE))
{
  list_of_maps_for_Sarraceniaceae_ancestral_branch[[i]] = stochastic_maps_7abprCTE[[i]]$maps[edgenum_Sarrac][[1]]
}
list_of_maps_for_Sarraceniaceae_ancestral_branch

# Number of steps in each of the 100 stochastic maps on this branch
sapply(X=list_of_maps_for_Sarraceniaceae_ancestral_branch, FUN=length)

# Let's get percentages of each state at several timepoints along a branch/edge
branchlength = tr$edge.length[edgenum_Sarrac]
branchlength

# Make a list of fractions along the branch
fractions_of_branch2 = seq(from=0, to=1, by=0.01)

# For a mapped branch, get the CUMULATIVE time when a change occurs
cumsum(list_of_maps_for_Sarraceniaceae_ancestral_branch[[1]])

# Table of states for each timepoint along branch
tmpmat = matrix(data=NA, nrow=length(stochastic_maps_7abprCTE), ncol=length(fractions_of_branch2))
dim(tmpmat)

# Fill out the table
for (i in 1:nrow(tmpmat))
{
  for (j in 1:ncol(tmpmat))
  {
    fraction_time = fractions_of_branch2[j]
    time_on_branch = branchlength * fraction_time
    time_on_branch
    cumulative_times_mapped = cumsum(list_of_maps_for_Sarraceniaceae_ancestral_branch[[i]])
    
    # Which time bin are we in:
    #time_on_branch = 10.5
    lessTF = time_on_branch <= cumulative_times_mapped
    lessTF
    
    # If all FALSE, state is the last state
    if (all(lessTF==FALSE))
    {
      lessTF[length(lessTF)] = TRUE
    }
    state_residences = 1:length(cumulative_times_mapped)
    state_residences
    state_residences[lessTF]
    the_state_at_the_time = min(state_residences[lessTF])
    the_state_at_the_time
    
    tmpmat[i,j] = names(cumulative_times_mapped)[the_state_at_the_time]
  }
}
tmpmat

# Get list of unique states in stochastic maps
uniq_states = sort(unique(c(tmpmat)))
uniq_states

# Get percentages of each, for each timepoint
percentages_by_timepoint = matrix(data=NA, ncol=length(fractions_of_branch2), nrow=length(uniq_states))
percentages_by_timepoint

for (i in 1:nrow(percentages_by_timepoint))
{
  for (j in 1:ncol(percentages_by_timepoint))
  {
    TF = tmpmat[,j] == uniq_states[i]
    percentages_by_timepoint[i,j] = sum(TF)
  }
}

percentages_by_timepoint

# Barplot showing percentages of each state at 100 timepoints along the branch
# below Lentibulariaceae
statecols = c("darkgreen","green3","grey50","lightgrey","orange","yellow2","grey90")
barplot(percentages_by_timepoint[7:1,], col=statecols,border=NA, space=0)

# Density plots for each state
timepoints = fractions_of_branch2*branchlength
plot(x=timepoints, y=1:101, main="State Distribution Over Time Along the Sarraceniaeceae Ancestral Branch", xlab="Time along branch (units of length)", ylab="Percentage of each state", pch=".", col="white")
for (i in 1:nrow(percentages_by_timepoint))
{
  lines(x=timepoints,y=percentages_by_timepoint[i,], lwd = 3, col=rev(statecols)[i])
}

names(stochastic_maps_7abprCTE[[1]])





#######################################################################
# Look at Nepenthaceae ancestral branch (thtough transitional state) stochastic maps #
#######################################################################
dev.off()
edgenum_Nepentheaceae

list_of_maps_for_Nepenthaceae_ancestral_branch = list()
for (i in 1:length(stochastic_maps_7abprCTE))
{
  list_of_maps_for_Nepenthaceae_ancestral_branch[[i]] = stochastic_maps_7abprCTE[[i]]$maps[edgenum_Nepentheaceae][[1]]
}
list_of_maps_for_Nepenthaceae_ancestral_branch

# Number of steps in each of the 100 stochastic maps on this branch
sapply(X=list_of_maps_for_Nepenthaceae_ancestral_branch, FUN=length)

# Let's get percentages of each state at several timepoints along a branch/edge
branchlength = tr$edge.length[edgenum_Nepentheaceae]
branchlength

# Make a list of fractions along the branch
fractions_of_branch2 = seq(from=0, to=1, by=0.01)

# For a mapped branch, get the CUMULATIVE time when a change occurs
cumsum(list_of_maps_for_Nepenthaceae_ancestral_branch[[1]])

# Table of states for each timepoint along branch
tmpmat = matrix(data=NA, nrow=length(stochastic_maps_7abprCTE), ncol=length(fractions_of_branch2))
dim(tmpmat)

# Fill out the table
for (i in 1:nrow(tmpmat))
{
  for (j in 1:ncol(tmpmat))
  {
    fraction_time = fractions_of_branch2[j]
    time_on_branch = branchlength * fraction_time
    time_on_branch
    cumulative_times_mapped = cumsum(list_of_maps_for_Nepenthaceae_ancestral_branch[[i]])
    
    # Which time bin are we in:
    #time_on_branch = 10.5
    lessTF = time_on_branch <= cumulative_times_mapped
    lessTF
    
    # If all FALSE, state is the last state
    if (all(lessTF==FALSE))
    {
      lessTF[length(lessTF)] = TRUE
    }
    state_residences = 1:length(cumulative_times_mapped)
    state_residences
    state_residences[lessTF]
    the_state_at_the_time = min(state_residences[lessTF])
    the_state_at_the_time
    
    tmpmat[i,j] = names(cumulative_times_mapped)[the_state_at_the_time]
  }
}
tmpmat

# Get list of unique states in stochastic maps
uniq_states = sort(unique(c(tmpmat)))
uniq_states

# Get percentages of each, for each timepoint
percentages_by_timepoint = matrix(data=NA, ncol=length(fractions_of_branch2), nrow=length(uniq_states))
percentages_by_timepoint

for (i in 1:nrow(percentages_by_timepoint))
{
  for (j in 1:ncol(percentages_by_timepoint))
  {
    TF = tmpmat[,j] == uniq_states[i]
    percentages_by_timepoint[i,j] = sum(TF)
  }
}

percentages_by_timepoint

# Barplot showing percentages of each state at 100 timepoints along the branch
# below Lentibulariaceae
dev.off()
statecols = c("darkgreen","green3","grey50","lightgrey","orange","yellow3")
barplot(percentages_by_timepoint[3:1,], col=statecols,border=NA, space=0)

# Density plots for each state
timepoints = fractions_of_branch2*branchlength
plot(x=timepoints, y=1:101, main="State Distribution Over Time Along the Nepenthaceae Ancestral Branch", xlab="Time along branch (units of length)", ylab="Percentage of each state", pch=".", col="white")
for (i in 1:nrow(percentages_by_timepoint))
{
  lines(x=timepoints,y=percentages_by_timepoint[i,], lwd = 3, col=rev(statecols)[i])
}

names(stochastic_maps_7abprCTE[[1]])







#####################################################################################
# Look at Droseraceae ancestral branch (thtough transitional state) stochastic maps #
#####################################################################################
dev.off()
edgenum_Droseraceae

list_of_maps_for_Droseraceae_ancestral_branch = list()
for (i in 1:length(stochastic_maps_7abprCTE))
{
  list_of_maps_for_Droseraceae_ancestral_branch[[i]] = stochastic_maps_7abprCTE[[i]]$maps[edgenum_Droseraceae][[1]]
}
list_of_maps_for_Droseraceae_ancestral_branch

# Number of steps in each of the 100 stochastic maps on this branch
sapply(X=list_of_maps_for_Droseraceae_ancestral_branch, FUN=length)

# Let's get percentages of each state at several timepoints along a branch/edge
branchlength = tr$edge.length[edgenum_Droseraceae]
branchlength

# Make a list of fractions along the branch
fractions_of_branch2 = seq(from=0, to=1, by=0.01)

# For a mapped branch, get the CUMULATIVE time when a change occurs
cumsum(list_of_maps_for_Droseraceae_ancestral_branch[[1]])

# Table of states for each timepoint along branch
tmpmat = matrix(data=NA, nrow=length(stochastic_maps_7abprCTE), ncol=length(fractions_of_branch2))
dim(tmpmat)

# Fill out the table
for (i in 1:nrow(tmpmat))
{
  for (j in 1:ncol(tmpmat))
  {
    fraction_time = fractions_of_branch2[j]
    time_on_branch = branchlength * fraction_time
    time_on_branch
    cumulative_times_mapped = cumsum(list_of_maps_for_Droseraceae_ancestral_branch[[i]])
    
    # Which time bin are we in:
    #time_on_branch = 10.5
    lessTF = time_on_branch <= cumulative_times_mapped
    lessTF
    
    # If all FALSE, state is the last state
    if (all(lessTF==FALSE))
    {
      lessTF[length(lessTF)] = TRUE
    }
    state_residences = 1:length(cumulative_times_mapped)
    state_residences
    state_residences[lessTF]
    the_state_at_the_time = min(state_residences[lessTF])
    the_state_at_the_time
    
    tmpmat[i,j] = names(cumulative_times_mapped)[the_state_at_the_time]
  }
}
tmpmat

# Get list of unique states in stochastic maps
uniq_states = sort(unique(c(tmpmat)))
uniq_states

# Get percentages of each, for each timepoint
percentages_by_timepoint = matrix(data=NA, ncol=length(fractions_of_branch2), nrow=length(uniq_states))
percentages_by_timepoint

for (i in 1:nrow(percentages_by_timepoint))
{
  for (j in 1:ncol(percentages_by_timepoint))
  {
    TF = tmpmat[,j] == uniq_states[i]
    percentages_by_timepoint[i,j] = sum(TF)
  }
}

percentages_by_timepoint

# Barplot showing percentages of each state at 100 timepoints along the branch
# below Lentibulariaceae
statecols = c("grey50","orange","yellow3","grey90")
barplot(percentages_by_timepoint[3:1,], col=statecols,border=NA, space=0)

# Density plots for each state
timepoints = fractions_of_branch2*branchlength
plot(x=timepoints, y=1:101, main="State Distribution Over Time Along the Droseraceae Ancestral Branch", xlab="Time along branch (units of length)", ylab="Percentage of each state", pch=".", col="white")
for (i in 1:nrow(percentages_by_timepoint))
{
  lines(x=timepoints,y=percentages_by_timepoint[i,], lwd = 3, col=rev(statecols)[i])
}

names(stochastic_maps_7abprCTE[[1]])





#Getting mean counts on phylogenetic stochastic mapping
#means of each column
apply(X=summary_stochastic_maps_7abprCTE$count, MARGIN=2, FUN=mean)

mean_counts <- apply(X=summary_stochastic_maps_7abprCTE$count[,2:111], MARGIN=2, FUN=mean)

mean_counts

max(mean_counts)

sort_means <- sort(mean_counts, decreasing = TRUE)

top_10_means <- head(sort_means, 10)
top_10_means

# Standard deviations
apply(X=summary_stochastic_maps_7abprCTE$count, MARGIN=2, FUN=sd)

# 95% confidence intervals
1.96 * apply(X=summary_stochastic_maps_7abprCTE$count, MARGIN=2, FUN=sd)
-1.96 * apply(X=summary_stochastic_maps_7abprCTE$count, MARGIN=2, FUN=sd)

# 95% credible with quantiles
apply(X=summary_stochastic_maps_7abprCTE$count [,2:111], MARGIN=2, FUN=quantile, prob=c(0.025, 0.975))



?getMRCA
tr
tr$edge
