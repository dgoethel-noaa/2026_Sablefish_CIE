########################################################################################################################

# Purpose: To create an RTMB model for the sablefish assessment using SPoRC
# Creator: Matthew LH. Cheng (UAF CFOS)
# Date 5/28/25
# NOTES: This uses SPoRC release v1.0.0


########################################################################################################################
# Setup -------------------------------------------------------------------

#install.packages("devtools") # install dev tools
#install.packages("TMB") # install TMB
#install.packages("RTMB") # install RTMB
#TMB:::install.contrib("https://github.com/vtrijoulet/OSA_multivariate_dists/archive/main.zip") # get multivariate OSA distributions

# optional packages to install
#devtools::install_github("fishfollower/compResidual/compResidual")
#devtools::install_github("noaa-afsc/afscOSA", dependencies = TRUE)

# install SPoRC
# devtools::install_github("chengmatt/SPoRC", dependencies = TRUE)
# devtools::install_github("chengmatt/SPoRC", dependencies = TRUE, lib = Sys.getenv("R_LIBS_USER"))

library(here)
library(SPoRC)
library(ggplot2)
library(RTMB)
library(tidyr)
library(tidyverse)

#############################################################################################################################################################################################################
#---------USER INPUTS--NEED TO UPDATE

mod_name <- '26.2b_FAA_3Flt_LLS_Only_No_Francis'                                                    # Model name to use for saving files and labels
root_folder <- here("Mods","__FAA_mod_build","Sensitivity Runs")                                            # The root folder where runs will be stored (sub folders named with mod_name)

do_francis <- 0                                                                # whether or not to do francis reweighting for this run, ==0 NO, ==1 YES
do_diags <- 0                                                                  # whether or not to do extra diagnostics (retrospective, jitter, profile likelihoods, etc) for this run; ==0 NO, ==1, YES

do_high_res_comps <- 0                                                         # whether to make comp graphs as tiffs (==1) or not (==0); these take up a lot of space ~1gb
do_OSA <- 1                                                                    # whether to do One Step Ahead (OSA) residual calcs and plots; Yes==1, NO==0
do_pearson <- 1                                                                # whether to calculate and print out Pearson residual bubble plots ; Yes==1, No==0
do_bias_ramp <- 0                                                              # whether to update the recruitment bias ramp parameters; NOTE: this is done after FRANCIS reweighting
##############################################################################################################################################################################################################

# Load in assessment data
SPoRC_one_reg <- readRDS(here(root_folder,mod_name,"SPoRC_FAA_3LLS_1LLF.rds"))

# Make plotting directory
dir.create(here(root_folder,mod_name,"plots"))


#######################################################################################################################
#------------------ Setup model weighting
# NOTE: These are replaced by the Francis Reweighting Procedure as long as do_francis <- 1
#######################################################################################################################

comps_ESS <- SPoRC_one_reg$comps_ESS                                                               # Input ESS for composition data
comps_ESS_FAA <- SPoRC_one_reg$comps_ESS_FAA                                                       # Input ESS for composition data by fleet

Wt_FishAgeComps <- array(1, dim = c(SPoRC_one_reg$nreg,
                                     SPoRC_one_reg$nyrs,
                                     SPoRC_one_reg$nsex,
                                     SPoRC_one_reg$nflt))

Wt_FishLenComps <- array(1, dim = c(SPoRC_one_reg$nreg,
                                     SPoRC_one_reg$nyrs,
                                     SPoRC_one_reg$nsex,
                                     SPoRC_one_reg$nflt))

Wt_SrvAgeComps <- array(1, dim = c(SPoRC_one_reg$nreg,
                                    SPoRC_one_reg$nyrs,
                                    SPoRC_one_reg$nsex,
                                    SPoRC_one_reg$nsrv_flt))

Wt_SrvLenComps <- array(1, dim = c(SPoRC_one_reg$nreg,
                                   SPoRC_one_reg$nyrs,
                                   SPoRC_one_reg$nsex,
                                   SPoRC_one_reg$nsrv_flt))

##################################################################################################################################################################################
# Prepare Model Inputs ----------------------------------------------------
# NOTE: this is all generalized so nothing needs to be updated unless you are changing the model structure; fully generalized based on updates to the data pull code
##################################################################################################################################################################################

#------------------ Initialize model dimensions and data list
input_list <- Setup_Mod_Dim(years = 1:SPoRC_one_reg$nyrs,                      # vector of years
                            ages = 1:SPoRC_one_reg$nages,                      # vector of ages
                            lens = SPoRC_one_reg$length_bins,                  # length bins
                            n_regions = SPoRC_one_reg$nreg,                    # number of regions, if nreg==1 then AK-wide, if n==5 then 1==BS, 2==AI, 3==WG, 4==CG, 5==EG
                            n_sexes = SPoRC_one_reg$nsex,                      # number of sexes, 1==female, 2==male
                            n_fish_fleets = SPoRC_one_reg$nflt,                # number of fishery fleet, if FAA default is 1=BS_FG, 2=AI_FG, 3=GOA_FG, 4=TF (aggregated); alt FAA models have 2 or 3 TF fleets
                            n_srv_fleets = SPoRC_one_reg$nsrv_flt,             # number of survey fleets, if FAA default is 1=BS_LLS, 2=AI_LLS, 3=GOA_LLS, 4=TS (aggregated, not fit), 5=JPN_LLS (aggregated)
                            verbose = TRUE
)

#------------------- Setup recruitment stuff (using defaults for other stuff)
input_list <- Setup_Mod_Rec(input_list = input_list,                           # input data list from above
                            do_rec_bias_ramp = 1,                              # do bias ramp (0 == don't do bias ramp, 1 == do bias ramp)
                            bias_year = c(8,     # breakpoints for bias ramp--NOTE these are flexibly coded and don't need to be changed as add more years to model (first entry == years no bias ramp - 1960 to 1980--years 1 to 21;
                                          24,     # second entry == years ascending limb of bias ramp - 1981 to 1990--years 21 to 30;
                                          63,
                                          65) ,
                            max_bias_ramp_fct = 0.9,
                            sigmaR_switch = as.integer(length(                 # when to switch from early to late sigmaR
                              SPoRC_one_reg$styr:1975)),
                            dont_est_recdev_last = 1,                          # don't estimate last recruitment deviate (set equal to mean recruitment)
                            ln_sigmaR = log(c(0.4, 0.9)),                      # early and late sigma_R starting values
                            rec_model = "mean_rec",                            # recruitment model, use mean recruitment (not BH SRR)
                            sigmaR_spec = "fix",                               # fix early sigmaR, estimate late sigmaR
                            InitDevs_spec = NULL,                              # estimate all initial deviations
                            RecDevs_spec = NULL,                               # estimate all recruitment deviations
                            init_age_strc = 1,                                 # Initial age structure is derived by assuming a geometric series (init_age_strc = 1; the alternative is iterating age structure to some equilibrium init_age_strc = 0),
                            init_F_prop = 0.1                                  # A 10% fraction of the mean fishing mortality rate is assumed for initializing the age structure.
)

#bias_year = c(12, 19, 63, 65),
#max_bias_ramp_fct = 0.9

#------------------- Setup Biological Inputs
M_prior <- data.frame(
  regionblk = 1,
  yearblk = 1,
  ageblk = 1,
  sexblk = 1:2,
  mu = 0.085,
  sd = 0.1
)

input_list <- Setup_Mod_Biologicals(input_list = input_list,
                                    WAA = SPoRC_one_reg$WAA,                   # weight-at-age
                                    MatAA = SPoRC_one_reg$MatAA,               # maturity at age
                                    AgeingError = as.matrix(                   # ageing error
                                      SPoRC_one_reg$AgeError),
                                    SizeAgeTrans = SPoRC_one_reg$SizeAgeTrans, # size age transition matrix

                                    # Model options
                                    fit_lengths = 1,                           # fitting length compositions
                                    Use_M_prior = 1,                           # use natural mortality prior
                                    M_prior = M_prior,                         # mean and sd for M prior
                                    M_spec = "est_ln_M",                       # estimate M
                                    M_sexblk_spec = list(1,2)                  # estimate sex specific M (one par per sex)
)

#------------------- Setup Movement (No Movement in Operational Assessment)
input_list <- Setup_Mod_Movement(input_list = input_list,
                                 use_fixed_movement = 1,                       # don't est move
                                 Fixed_Movement = NA,                          # use identity move
                                 do_recruits_move = 0                          # recruits don't move
)

#------------------- Setup Tagging (No Tagging Fit in Operational Assessment)
input_list <- Setup_Mod_Tagging(input_list = input_list,
                                UseTagging = 0                                 # don't fit/use tagging data
)

#------------------- Setup Catch Specifications
input_list <- Setup_Mod_Catch_and_F(input_list = input_list,
                                    ObsCatch = SPoRC_one_reg$ObsCatch,         # Observed Catch in mt
                                    Catch_Type =                               # just a switch to determine whether to fit region-specific catch, ==1 aggregate across regions (does nothing for single region model)
                                      array(1, dim = c(length(input_list$data$years),
                                                       input_list$data$n_fish_fleets)),
                                    UseCatch = SPoRC_one_reg$UseCatch,         # Years in which catch data is available and should be fit in the model

                                    # Model options
                                    Use_F_pen = 1,                             # whether to use f penalty, == 0 don't use, == 1 use
                                    ln_sigmaF = array(log(0.75),
                                                      dim = c(SPoRC_one_reg$nreg,
                                                              SPoRC_one_reg$nflt)),          # the sd for the F prior (high value gives low 'weight' to the F prior...1.0 is the default value)
                                    sigmaC_spec = 'fix',                       # Whether the variance term for fitting catch data is estimated or fixed
                                    ln_sigmaC = array(log(.05),
                                                      dim = c(SPoRC_one_reg$nreg,
                                                              SPoRC_one_reg$nyrs,
                                                              SPoRC_one_reg$nflt)),        # 0.05 is common value to fit catch relatively closely, but allow for some error
                                    Catch_Constant = array(0,
                                                           dim = c(SPoRC_one_reg$nreg,
                                                                   SPoRC_one_reg$nflt)),              # The small constant to add to the catch term to avoid log(0)--Fournier robust likelihood approach; fixed gear then trawl fishery
)

#------------------- Setup fishery indices and compositions
input_list <- Setup_Mod_FishIdx_and_Comps(input_list = input_list,
                                          ObsFishIdx =                         # Observed fishery indices (CPUE)
                                            SPoRC_one_reg$ObsFishIdx,
                                          ObsFishIdx_SE =                      # Fishery index standard errors
                                            (SPoRC_one_reg$ObsFishIdx_SE) * 2,  # converting to true variance for use in the lognormal likelihood; *2 is scaling fishery CPUE to mean CV of 0.2 (calc historically as mean of 0.1 which is too small)
                                          UseFishIdx =                         # Switch to identify years in which fishery indices exist and whether to fit them in the model
                                            SPoRC_one_reg$UseFishIdx,
                                          ObsFishAgeComps =                    # Observed fishery age comps
                                            SPoRC_one_reg$ObsFishAgeComps,
                                          UseFishAgeComps =                    # Switch identifying whether fishery age comps exist and whether to fit them
                                            SPoRC_one_reg$UseFishAgeComps,
                                          ISS_FishAgeComps =                   # The input sample sizes for fishery age comps
                                            SPoRC_one_reg$ISS_FishAgeComps,
                                          ObsFishLenComps =                    # Observed fishery length comps
                                            SPoRC_one_reg$ObsFishLenComps,
                                          UseFishLenComps =                    # Whether to fit fishery length comps in a given year
                                            SPoRC_one_reg$UseFishLenComps,
                                          ISS_FishLenComps =                   # Input sample size for fishery length comps
                                            SPoRC_one_reg$ISS_FishLenComps,

                                          # Model options
                                          fish_idx_type =                      # Define units (e.g. biomass vs. abundance/numbers) for fishery indices (none is used if no index for a given fleet)
                                            c("biom","none"),
                                          FishAgeComps_LikeType =              # age comp likelihoods by fleet (none if not fit or not available for given fleet)
                                            c("Multinomial","none"),
                                          FishLenComps_LikeType =              # length comp likelihoods for fishery fleet 1 and 2
                                            c("Multinomial","Multinomial"),
                                          FishAgeComps_Type =                  # age comp structure for fishery fleet 1 and 2 (agg_Year_1-terminal_Fleet_1 indicates that age comps are aggregated across sex for all years for fleet 1--fixed gear)
                                            c("spltRjntS_Year_1-terminal_Fleet_1",   # NOTE: these entries are needed for all years and all fleets, so if want to agg for some years then disagg for remaining years, would need two entries for that fleet
                                              "none_Year_1-terminal_Fleet_2"), # NOTE: if want joint by sex would use spltRjntS_Year_1-terminal_Fleet_1
                                          FishLenComps_Type =                  # length comp structure for fishery fleet 1 and 2 (spltRspltS_Year_1-terminal_Fleet_1 indicates that length comps are disaggregated across sex and fit using 'split' approach for all years for fleet 1--fixed gear)
                                            c("spltRjntS_Year_1-terminal_Fleet_1",
                                              "spltRjntS_Year_1-terminal_Fleet_2"),
                                          FishAge_comp_agg_type = c(1,NA),     # for age comps, identifies the sequence of the normalizing, aggregation, and aging error process (matches what is done in ADMB sablefish currently)
                                          FishLen_comp_agg_type = c(0,0)       # for length comps, identifies the sequence of the normalizing, aggregation, and aging error process (matches what is done in ADMB sablefish currently)
)

#------------------- Setup survey indices and compositions

input_list <- Setup_Mod_SrvIdx_and_Comps(input_list = input_list,
                                         ObsSrvIdx =                           # Observed fishery independent survey indices
                                           SPoRC_one_reg$ObsSrvIdx,
                                         ObsSrvIdx_SE =                        # Observed fishery independent survey standard errors
                                           (SPoRC_one_reg$ObsSrvIdx_SE)*1.5,     # generally LLS (JPN and DOM) CV ~ 0.05, TS ~ 0.1; the *2 is to give more reasonable uncertainty and weights to survey data
                                         UseSrvIdx =                           # Whether to fit fishery independent surveys in each year
                                           SPoRC_one_reg$UseSrvIdx,
                                         ObsSrvAgeComps =                      # Observed age comps for fishery independent surveys
                                           SPoRC_one_reg$ObsSrvAgeComps,
                                         ISS_SrvAgeComps =                     # Input sample size for fishery independent survey age comps
                                           SPoRC_one_reg$ISS_SrvAgeComps,
                                         UseSrvAgeComps =                      # Whether to fit fishery independent survey age comps in a given year
                                           SPoRC_one_reg$UseSrvAgeComps,
                                         ObsSrvLenComps =                      # Observed fishery independent survey length comps
                                           SPoRC_one_reg$ObsSrvLenComps,
                                         UseSrvLenComps =                      # Whether to fit fishery independent survey length comps in a given year
                                           SPoRC_one_reg$UseSrvLenComps,
                                         ISS_SrvLenComps =                     # Input sample size for fishery independent survey length comps
                                           SPoRC_one_reg$ISS_SrvLenComps,

                                         # Model options
                                         srv_idx_type =                        # Define units (e.g. biomass vs. abundance/numbers) for fishery independent survey indices (none is used if no index for a given fleet)
                                           c("abd", "abd","abd","biom", "abd"),            # NOAA LLS uses RPNs, Trawl Survey uses biomass, and JPN LLS uses RPNs
                                         SrvAgeComps_LikeType =                # age comp likelihoods by survey fleet (none if not fit or not available for given fleet)
                                           c("Multinomial", "Multinomial","Multinomial","none", "Multinomial"),
                                         SrvLenComps_LikeType =                # length comp likelihoods for fishery survey fleets
                                           c("Multinomial","Multinomial","Multinomial", "Multinomial", "Multinomial"),
                                         SrvAgeComps_Type =                    # age comp structure for survey fleets (agg_Year_1-terminal_Fleet_1 indicates that age comps are aggregated across sex for all years for fleet 1--NOAA LLS)
                                           c("spltRjntS_Year_1-terminal_Fleet_1",
                                             "spltRjntS_Year_1-terminal_Fleet_2",
                                             "spltRjntS_Year_1-terminal_Fleet_3",
                                             "none_Year_1-terminal_Fleet_4",
                                             "spltRjntS_Year_1-terminal_Fleet_5"),
                                         SrvLenComps_Type =                    # length comp structure for survey fleets (spltRspltS_Year_1-terminal_Fleet_1 indicates that length comps are disaggregated across sex and fit using 'split' approach for all years for fleet 1--NOAA LLS)
                                           c("spltRjntS_Year_1-terminal_Fleet_1",
                                             "spltRjntS_Year_1-terminal_Fleet_2",
                                             "spltRjntS_Year_1-terminal_Fleet_3",
                                             "spltRjntS_Year_1-terminal_Fleet_4",
                                             "spltRjntS_Year_1-terminal_Fleet_5"),
                                         SrvAge_comp_agg_type = c(1,1,1,NA,1),     # for age comps, identifies the sequence of the normalizing, aggregation, and aging error process (matches what is done in ADMB sablefish currently)
                                         SrvLen_comp_agg_type = c(0,0,0,0,0)       # for length comps, identifies the sequence of the normalizing, aggregation, and aging error process (matches what is done in ADMB sablefish currently)
)


#------------------- Setup data weighting
input_list <- Setup_Mod_Weighting(input_list = input_list,
                                  Wt_Catch = 1,                                # Weight (lambda) for catch data
                                  Wt_FishIdx = 1,                              # Weight (lambda) for fishery indices data
                                  Wt_SrvIdx = 1,                               # Weight (lambda) for survey indices data
                                  Wt_Rec = 1,                                  # Weight (lambda) for stock-recruit penalty term
                                  Wt_F = 1,                                    # Weight (lambda) for fishing mortality penalty term
                                  Wt_Tagging = 0,                              # Weight (lambda) for tagging data (not fit in single region model)
                                  Wt_FishAgeComps = Wt_FishAgeComps,           # Weights (lambda) for fishery age comp data (note these are updated by Francis reweighting and input at start of this script)
                                  Wt_FishLenComps = Wt_FishLenComps,           # Weights (lambda) for fishery length comp data (note these are updated by Francis reweighting and input at start of this script)
                                  Wt_SrvAgeComps = Wt_SrvAgeComps,             # Weights (lambda) for survey age comp data (note these are updated by Francis reweighting and input at start of this script)
                                  Wt_SrvLenComps = Wt_SrvLenComps              # Weights (lambda) for survey length comp data (note these are updated by Francis reweighting and input at start of this script)
)

#------------------- Setup up fishery selectivity
# Setup priors for selectivity
# Define valid fleet-block combinations
fleet_blocks <- data.frame(
  fleet = c(1,1, 2),                                        # fleets corresponding to time blocks (only 2 time blocks for GOA LLF, none for other fleets)
  block = c(1,2, 1)                                       # corresponding time blocks
)

# Define sex and parameter combinations
sex_par <- expand.grid(sex = 1:2, par = 1:2)

# Merge to get all valid combinations
fish_selex_structure <- merge(fleet_blocks, sex_par) 
# dplyr::filter(!(fleet == 1 & block == 1 & sex == 2 & par == 2)) %>%          # remove priors for any unestimated pars -- par1=a50, par2=delta; NEEDS TO MATCH PARAMETER MAPPING
#  dplyr::filter(!(fleet == 2 & block == 1 & sex == 2 & par == 2)) %>%          # remove priors for any unestimated pars -- par1=a50, par2=delta; NEEDS TO MATCH PARAMETER MAPPING
# dplyr::filter(!(fleet == 3 & block == 1 & sex == 2 & par == 2)) %>%          # remove priors for any unestimated pars -- par1=a50, par2=delta; NEEDS TO MATCH PARAMETER MAPPING
# dplyr::filter(!(fleet == 4 & block == 1 & sex == 2 & par == 1))              # remove priors for any unestimated pars -- par1=a50, par2=delta; NEEDS TO MATCH PARAMETER MAPPING

# Add the lognormal prior values - creates a dataframe, each row is a unique parameter combination to apply the prior to
fish_selex_prior <- cbind(
  region = 1,
  fish_selex_structure,
  mu = 4,                                                                      # All selex means = 1 (means should be defined in normal space)
  sd = 1.5                                                                       # All selex sd = 5
)


#########################################################################################################################################
# Quick plot of gamma to help determine appropriate priors...
# Extract out and exponentiate the parameters here
#bmax = 5 # Bin at max selex
#delta = 3 # slope parameter

#Bin = 1:30
# Now, calculate/derive power parameter + selex values
#p = 0.5 * (sqrt( bmax^2 + (4 * delta^2)) - bmax)
#selex = (Bin / bmax)^(bmax/p) * exp( (bmax - Bin) / p ) # return parametric form

#plot(selex, type = 'l')
##########################################################################################################################################                     

input_list <- Setup_Mod_Fishsel_and_Q(input_list = input_list,
                                      cont_tv_fish_sel =                       # whether continuous time-varying selectivity is used for either fishery fleet
                                        c("none_Fleet_1", "none_Fleet_2"),
                                      fish_sel_blocks =                        # fishery selectivity time blocks if not TV specified above for a given fleet
                                        c("Block_1_Year_1-56_Fleet_1",          # Early time block LLS-- 1960 to 2016
                                          "Block_2_Year_57-terminal_Fleet_1",   # Recent time block for LLS--2017 to terminal year                      # Recent time block for fixed gear fishery--2016 to terminal year
                                          "none_Fleet_2"),                     # no blocks for trawl fishery
                                      fish_sel_model =                         # fishery selectivity form
                                        c("logist1_Block_1_Fleet_1",                   # logistic selectivity for fixed gear fleet
                                          "logist1_Block_2_Fleet_1",                   # logistic selectivity for fixed gear fleet
                                          "logist1_Fleet_2"),                    # gamma function selectivity for trawl gear fleet
                                      fish_q_blocks =                          # fishery indices catchability time blocks
                                        c(
                                          "none_Fleet_1",         # pre-IFQ time block for fixed gear fishery 1994 and before
                                          "none_Fleet_2"),                     # no blocks for trawl fishery since no CPUE index used
                                      fish_fixed_sel_pars =                    # Whether to estimate all fixed effects for fishery selectivity, but can later modify to fix and share parameters via the parameter mapping
                                        c("est_all", "est_all"),               # Estimate all selectivity parameters for both fishery fleets
                                      fish_q_spec =                            # Whether to estimate all fixed effects for fishery catchability
                                        c("fix","fix"),           # DO NOT Estimate fishery q for fixed gear or for trawl (no CPUE index)
                                      Use_fish_selex_prior = 1,                # Using selex priors
                                      fish_selex_prior = fish_selex_prior
)

#------------------- Setup survey selectivity
# Setup priors for selectivity
# Define valid fleet-block combinations
fleet_blocks <- data.frame(
  fleet = c(1, 1,2,2, 3,3, 4, 5),                                             # fleets for each time block (LLS has 2 timeb blocks TS and JPN have 1 each)
  block = c(1, 2,1,2, 1,2, 1, 1)                                              # corresponding blocks
)

# Define sex and parameter combinations
sex_par <- expand.grid(sex = 1:2, par = 1:2)

# Merge to get all valid combinations
srv_selex_structure <- merge(fleet_blocks, sex_par) %>%
  dplyr::filter(!(fleet == 4))                                                 # remove priors for any unestimated pars -- par1=a50, par2=delta; NEEDS TO MATCH PARAMETER MAPPING

# Add the lognormal prior values - creates a dataframe, each row is a unique parameter combination to apply the prior to
srv_selex_prior <- cbind(
  region = 1,
  srv_selex_structure,
  mu = 3,                                                                      # All selex means = 1 (means should be defined in normal space)
  sd = 1.5                                                                       # All selex sd = 5
)


input_list <- Setup_Mod_Srvsel_and_Q(input_list = input_list,
                                     cont_tv_srv_sel =                         # whether continuous time-varying selectivity is used for either survey fleet
                                       c("none_Fleet_1",
                                         "none_Fleet_2",
                                         "none_Fleet_3",
                                         "none_Fleet_4",
                                         "none_Fleet_5"),
                                     srv_sel_blocks =                          # survey selectivity time blocks if not TV specified above for a given fleet
                                       c("Block_1_Year_1-56_Fleet_1",          # Early time block LLS-- 1960 to 2016
                                          "Block_2_Year_57-terminal_Fleet_1",   # Recent time block for LLS--2017 to terminal year
                                          "Block_1_Year_1-56_Fleet_2",          # Early time block LLS-- 1960 to 2016
                                          "Block_2_Year_57-terminal_Fleet_2",   # Recent time block for LLS--2017 to terminal year
                                          "Block_1_Year_1-56_Fleet_3",          # Early time block LLS-- 1960 to 2016
                                          "Block_2_Year_57-terminal_Fleet_3",   # Recent time block for LLS--2017 to terminal year
                                          "none_Fleet_4",                       # No blocks for trawl survey
                                          "none_Fleet_5"                        # No blocks for JPN LLS
                                       ),
                                     srv_sel_model =                           # Survey selectivity form
                                       c("logist1_Block_1_Fleet_1",
                                         "logist1_Block_2_Fleet_1",                      # gamma for BS
                                         "logist1_Block_1_Fleet_2",                    # logistic selectivity for LLS
                                         "logist1_Block_2_Fleet_2",                    # logistic selectivity for LLS
                                         "logist1_Block_1_Fleet_3",
                                         "logist1_Block_2_Fleet_3",
                                          "exponential_Fleet_4",                # exponential declining selectivity for trawl survey
                                         "logist1_Fleet_5"),                   # logistic selectivity for JPN LLS
                                     srv_q_blocks =                            # survey indices catchability time blocks (no q time blocks for any fleets)
                                       c("none_Fleet_1",
                                         "none_Fleet_2",
                                         "none_Fleet_3",
                                         "none_Fleet_4",
                                         "none_Fleet_5"),
                                     srv_fixed_sel_pars_spec =                 # Whether to estimate all fixed effects for survey selectivity, but can later modify to fix and share parameters via the parameter mapping
                                       c("est_all",                            # Estimate all selectivity parameters for all survey fleets
                                         "est_all",
                                         "est_all",
                                         "fix",
                                         "est_all"),
                                     srv_q_spec =                              # Whether to estimate all fixed effects for survey catchability
                                       c("est_all",                            # Estimate q for all survey fleets
                                         "est_all",
                                         "est_all",
                                         "fix",
                                         "est_all"),
                                     Use_srv_selex_prior = 1,                 # Using selex priors
                                     srv_selex_prior = srv_selex_prior
)


##################################################################################################################################################################################
# Do some parameter mapping for sharing of selectivity pars
# NOTE: the parameter mapping is a bit complex, but more digestable when look at how the parameter array is filled
#          The dimensions of the selectivity array are n_regions, max_fish_pars (i.e., the max number of parameters in any given selectivity function being used--2 pars for logistic),
#           max_fishsel_blks (i.e., the maximum number of time blocks in any given fishery--3 in the case of fixed gear), n_sexes, n_fish_fleets
##################################################################################################################################################################################

#-------------------------------- Fishery Selectivity
# To fill the fishery selectivity par array, for each gear/sex/region there needs to be a matrix with dims of max_fishsel_blks X max_fish_pars (for sablefish this is 3 rows X 2 cols)
# When entering factors, use a unique number for each parameter being estimated, and repeat a number to share the current par to an earlier par (i.e., if 2 is repeated then the current par will use the same value as par 2)
# Since LLF (gear 1) has estimated pars for all female pars across all time blocks (2 selex pars X 3 time blocks ==6 pars total), factors 1:6 are input for gear 1 sex 1 selectivities
# Since LLF (gear 1) males (sex 2) share the logistic selex delta par in the first time block, factor 2 is repeated after factor 7 (i.e., LLF male first time block delta sel par is the same as the corresponding par for LLF females)
# For gear/time block combos where pars are the same, then just repeat the values (e.g., no time blocks for trawl fleet--gear2--selex, so factors/pars are repeated for each of the 3 time blocks)
# For selex functions with < max_fish_pars, just repeat the par number in the place of the par that isn't being used (e.g., for trawl survey--nsrvy_flt2--use a single par exponential function, so just repeat 1st par # for that survey in the second par spot)
# If a parameter is set as a fixed input, then just put NA to indicate that it is not estimated (e.g., survey fleet 3--JPN LLS uses fixed input for a50, so NA put in that par spot)

# Can use this code to see how the parameter array is constructed
# array(c(1:7,2,8:11, rep(12:13,3),rep(c(14,13),3)), dim=dim(input_list$par$ln_fish_fixed_sel_pars))

# CONT SETUP
input_list$map$ln_fish_fixed_sel_pars <- factor(c(rep(c(1:4),1),               # estimate a50 and a95 for females (2pars), but repeat 2xs for the 2 time blocks in the LLF fleet 
                                                  rep(c(5:8),1),               # estimate a50 and a95 for males (2pars), but repeat 2xs for the 2 time blocks in the LLF fleet 
                                                  rep(9:10,2),                # for trawl fishery estimate both female gamma pars, but no time blocks (repeat 2 times to fill the 2 time blocks associated with the LLF GOA)
                                                  rep(c(11,12),2)))            # for trawl fishery males, share the bmax par of the gamma function with the female estimate par (repeat 2 times to fill the 2 time blocks associated with the LLF)

#-------------------------------- Survey Selectivity
# For survey selectivity max_srvsel_blks==2 and max_srv_pars=2, so need 2 X 2 matrix for each reg/fleet/sex
# Trawl survey (flt2) assumes a 1 par exponential with no time blocks, so repeat the factor/par number 4 times for each sex
# JPN LLS (flt3) uses a fixed input a50, so put NA for that value, then it shares the delta par with the Dom LLS (fleet 1), so use the factor associated with sex-specfic delta from the first LLS time block

# Can use this code to see how the parameter array is constructed
#array(c(1:3,2,4:6,5,rep(7,4),rep(8, 4),rep(c(NA,2), 2),rep(c(NA, 5), 2)),dim=dim(input_list$par$ln_srv_fixed_sel_pars))

# Share JPN sex-delta
input_list$map$ln_srv_fixed_sel_pars <-
  factor(c(rep(c(1:4),1),              # females BS 2 pars, 2 blocks                                                         
           rep(c(5:8),1),              # males BS 2 pars, 2 blocks      
           rep(c(9:12),1),             # females AI 2 pars, 2 blocks                                                          
           rep(c(13:16),1),            # males AI 2 pars, 2 blocks      
           rep(c(17:20),1),            # females GOA 2 pars, 2 blocks                                                                  
           rep(c(21:24),1),            # males GOA 2 pars, 2 blocks      
           rep(NA,4),                                                          # domestic trawl survey (single time block) and only one parameter (NOT ESTIMATED) (exponential selex fxn), so only one parameter estimated per sex across blocks (so repeat value 2 times for each sex)
           rep(NA, 4),                                                         # domestic trawl survey (single time block) and only one parameter (NOT ESTIMATED) (exponential selex fxn), so only one parameter estimated per sex across blocks (so repeat value 2 times for each sex)
           rep(c(25:26),2),                                                    # coop JPN LL survey (single time block, females)
           rep(c(27:28),2)))                                                   # coop JPN LL survey (single time block, males), share delta

# Selex starting values
input_list$par$ln_srv_fixed_sel_pars[1,,,,] <- log(3)
input_list$par$ln_fish_fixed_sel_pars[1,,,,] <- log(3)

##################################################################################################################################################################################
# Extract out lists updated with helper functions ---------------------------------------------------------------
# NOTE: these three assignments must be done right before fitting the model, otherwise the updated mapping, etc. won't be recorded
##################################################################################################################################################################################

data <- input_list$data
parameters <- input_list$par
mapping <- input_list$map



if(do_bias_ramp == 1){
  
  sabie_rtmb_model <- fit_model(data,
                                parameters,
                                mapping,
                                random = NULL,
                                newton_loops = 3,
                                silent = TRUE
  )
  
  SPoRC:::SPoRC_rtmb(parameters, data)
  
  # Get standard error report
  sabie_rtmb_model$sd_rep <- RTMB::sdreport(sabie_rtmb_model)
  
  source(here("R","mt_rec_bias_correction.r"))
  sigr2 <- c(rep(exp(parameters$ln_sigmaR[1])^2, 15), rep(exp(parameters$ln_sigmaR[2])^2, length(16:(SPoRC_one_reg$nyrs-1))))
  dont_est_recdev_last <- 1
  bias_ramp <- mt_rec_bias_correction(
    data, 
    sabie_rtmb_model$sd_rep,
    sigr2,
    dont_est_recdev_last = 1
  )
  
  bias_ramp$best_params
  write.csv(bias_ramp$best_params,here(root_folder,mod_name,"plots",paste0(mod_name,"_recr_bias_corr.csv")))
  
  # Saves the current screen plot to a PNG file
  dev.copy(png, here(root_folder,mod_name,"plots",paste0(mod_name,"_recr_bias_corr.png")))
  dev.off()  # Mandatory to close the file and finish the save
  
  
  data$bias_year <- c(as.numeric(bias_ramp$best_params[1]),     # breakpoints for bias ramp--NOTE these are flexibly coded and don't need to be changed as add more years to model (first entry == years no bias ramp - 1960 to 1980--years 1 to 21;
                      as.numeric(bias_ramp$best_params[2]),     # second entry == years ascending limb of bias ramp - 1981 to 1990--years 21 to 30;
                      as.numeric(bias_ramp$best_params[3]),
                      as.numeric(bias_ramp$best_params[4]))
  #data$max_bias_ramp_fct <- as.numeric(bias_ramp$best_params[5])   # just leave at 0.9 
  
}


##################################################################################################################################################################################
# Run Francis  ------------------------------------------------------------

if(do_francis == 1){

  n_francis_iter <- 10                                                           # number of francis iterations to run

  st_wts <- data.frame(Comp_Type = c("GOA_LLF_Age_F", "GOA_LLF_Age_M",
                                     "BS_LLS_Age_F","BS_LLS_Age_M", "AI_LLS_Age_F","AI_LLS_Age_M", "GOA_LLS_Age_F","GOA_LLS_Age_M", 
                                     "JPN_LLS_Age_F", "JPN_LLS_Age_M",
                                     "TF_Len_F", "TF_Len_M",
                                     #"BS_LLS_Len_F", "BS_LLS_Len_M","AI_LLS_Len_F", "AI_LLS_Len_M",
                                      "JPN_LLS_Len_F", "JPN_LLS_Len_M"),
                       Init_Weight = c(unique(data$Wt_FishAgeComps[1,which(data$UseFishAgeComps[,,1] == 1),1,1]), 
                                       unique(data$Wt_FishAgeComps[1,which(data$UseFishAgeComps[,,1] == 1),2,1]),
                                       unique(data$Wt_SrvAgeComps[1,which(data$UseSrvAgeComps[,,1] == 1),1,1]), 
                                       unique(data$Wt_SrvAgeComps[1,which(data$UseSrvAgeComps[,,1] == 1),2,1]),
                                       unique(data$Wt_SrvAgeComps[1,which(data$UseSrvAgeComps[,,2] == 1),1,2]), 
                                       unique(data$Wt_SrvAgeComps[1,which(data$UseSrvAgeComps[,,2] == 1),2,2]),
                                       unique(data$Wt_SrvAgeComps[1,which(data$UseSrvAgeComps[,,3] == 1),1,3]), 
                                       unique(data$Wt_SrvAgeComps[1,which(data$UseSrvAgeComps[,,3] == 1),2,3]),
                                       unique(data$Wt_SrvAgeComps[1,which(data$UseSrvAgeComps[,,5] == 1),1,5]), 
                                       unique(data$Wt_SrvAgeComps[1,which(data$UseSrvAgeComps[,,5] == 1),2,5]),
                                      # unique(data$Wt_FishLenComps[1,which(data$UseFishLenComps[,,1] == 1),1,1]), 
                                      # unique(data$Wt_FishLenComps[1,which(data$UseFishLenComps[,,1] == 1),2,1]),
                                      # unique(data$Wt_FishLenComps[1,which(data$UseFishLenComps[,,2] == 1),1,2]), 
                                      # unique(data$Wt_FishLenComps[1,which(data$UseFishLenComps[,,2] == 1),2,2]),
                                      # unique(data$Wt_FishLenComps[1,which(data$UseFishLenComps[,,3] == 1),1,3]), 
                                      # unique(data$Wt_FishLenComps[1,which(data$UseFishLenComps[,,3] == 1),2,3]),
                                       unique(data$Wt_FishLenComps[1,which(data$UseFishLenComps[,,2] == 1),1,2]), 
                                       unique(data$Wt_FishLenComps[1,which(data$UseFishLenComps[,,2] == 1),2,2]),
                                      # unique(data$Wt_SrvLenComps[1,which(data$UseSrvLenComps[,,1] == 1),1,1]),      # no length comps used for BS and AI because always had ages for domestic survey since started sampling BSAI
                                      # unique(data$Wt_SrvLenComps[1,which(data$UseSrvLenComps[,,1] == 1),2,1]),
                                      # unique(data$Wt_SrvLenComps[1,which(data$UseSrvLenComps[,,2] == 1),1,2]), 
                                      # unique(data$Wt_SrvLenComps[1,which(data$UseSrvLenComps[,,2] == 1),2,2]),
                                      # unique(data$Wt_SrvLenComps[1,which(data$UseSrvLenComps[,,3] == 1),1,3]), 
                                      # unique(data$Wt_SrvLenComps[1,which(data$UseSrvLenComps[,,3] == 1),2,3]),
                                       unique(data$Wt_SrvLenComps[1,which(data$UseSrvLenComps[,,5] == 1),1,5]), 
                                       unique(data$Wt_SrvLenComps[1,which(data$UseSrvLenComps[,,5] == 1),2,5])))

  # get francis weights
  francis_run <- SPoRC::run_francis(data,
                                    parameters,
                                    mapping,
                                    random = NULL,
                                    n_francis_iter = n_francis_iter,
                                    newton_loops = 3
  )

  # get objects
  sabie_rtmb_model <- francis_run$obj # get RTMB object
  data <- sabie_rtmb_model$data # get data list with updated francis weights (NOTE! Outputting the new data list back out)
  rep <- sabie_rtmb_model$rep # get report file

  # get final weights and compare to starting weights
  end_wts <- data.frame(Comp_Type = c("GOA_LLF_Age_F", "GOA_LLF_Age_M",
                                      "BS_LLS_Age_F","BS_LLS_Age_M", "AI_LLS_Age_F","AI_LLS_Age_M", "GOA_LLS_Age_F","GOA_LLS_Age_M", 
                                      "JPN_LLS_Age_F", "JPN_LLS_Age_M",
                                      "TF_Len_F", "TF_Len_M",
                                      #"BS_LLS_Len_F", "BS_LLS_Len_M","AI_LLS_Len_F", "AI_LLS_Len_M",
                                      "JPN_LLS_Len_F", "JPN_LLS_Len_M"),
                        Final_Weight = c(unique(data$Wt_FishAgeComps[1,which(data$UseFishAgeComps[,,1] == 1),1,1]), 
                                        unique(data$Wt_FishAgeComps[1,which(data$UseFishAgeComps[,,1] == 1),2,1]),
                                        unique(data$Wt_SrvAgeComps[1,which(data$UseSrvAgeComps[,,1] == 1),1,1]), 
                                        unique(data$Wt_SrvAgeComps[1,which(data$UseSrvAgeComps[,,1] == 1),2,1]),
                                        unique(data$Wt_SrvAgeComps[1,which(data$UseSrvAgeComps[,,2] == 1),1,2]), 
                                        unique(data$Wt_SrvAgeComps[1,which(data$UseSrvAgeComps[,,2] == 1),2,2]),
                                        unique(data$Wt_SrvAgeComps[1,which(data$UseSrvAgeComps[,,3] == 1),1,3]), 
                                        unique(data$Wt_SrvAgeComps[1,which(data$UseSrvAgeComps[,,3] == 1),2,3]),
                                        unique(data$Wt_SrvAgeComps[1,which(data$UseSrvAgeComps[,,5] == 1),1,5]), 
                                        unique(data$Wt_SrvAgeComps[1,which(data$UseSrvAgeComps[,,5] == 1),2,5]),
                                        # unique(data$Wt_FishLenComps[1,which(data$UseFishLenComps[,,1] == 1),1,1]), 
                                        # unique(data$Wt_FishLenComps[1,which(data$UseFishLenComps[,,1] == 1),2,1]),
                                        # unique(data$Wt_FishLenComps[1,which(data$UseFishLenComps[,,2] == 1),1,2]), 
                                        # unique(data$Wt_FishLenComps[1,which(data$UseFishLenComps[,,2] == 1),2,2]),
                                        # unique(data$Wt_FishLenComps[1,which(data$UseFishLenComps[,,3] == 1),1,3]), 
                                        # unique(data$Wt_FishLenComps[1,which(data$UseFishLenComps[,,3] == 1),2,3]),
                                        unique(data$Wt_FishLenComps[1,which(data$UseFishLenComps[,,2] == 1),1,2]), 
                                        unique(data$Wt_FishLenComps[1,which(data$UseFishLenComps[,,2] == 1),2,2]),
                                        # unique(data$Wt_SrvLenComps[1,which(data$UseSrvLenComps[,,1] == 1),1,1]),      # no length comps used for BS and AI because always had ages for domestic survey since started sampling BSAI
                                        # unique(data$Wt_SrvLenComps[1,which(data$UseSrvLenComps[,,1] == 1),2,1]),
                                        # unique(data$Wt_SrvLenComps[1,which(data$UseSrvLenComps[,,2] == 1),1,2]), 
                                        # unique(data$Wt_SrvLenComps[1,which(data$UseSrvLenComps[,,2] == 1),2,2]),
                                        # unique(data$Wt_SrvLenComps[1,which(data$UseSrvLenComps[,,3] == 1),1,3]), 
                                        # unique(data$Wt_SrvLenComps[1,which(data$UseSrvLenComps[,,3] == 1),2,3]),
                                        unique(data$Wt_SrvLenComps[1,which(data$UseSrvLenComps[,,5] == 1),1,5]), 
                                        unique(data$Wt_SrvLenComps[1,which(data$UseSrvLenComps[,,5] == 1),2,5])))

  comp_wts <- st_wts %>% full_join(end_wts)

  Francis_plot <- grid::grid.grabExpr(gridExtra::grid.table(comp_wts))
  Francis_plot1 <- cowplot::ggdraw() + cowplot::draw_grob(Francis_plot)

  write.csv(comp_wts,file = here(root_folder,mod_name,"plots",paste0(mod_name,"_Francis_Comp_Wts.csv")))

  # Do francis mean plot
  if(unique(data$FishAgeComps_Type[,1]) != 2){        # joint sex comps
    plot_francis <- francis_run$end_mean_francis %>%
      drop_na() %>%
      ggplot() +
      geom_point(francis_run$end_mean_francis %>% drop_na(), mapping = aes(x = Comp_Year, y = obs)) +
      geom_line(francis_run$end_mean_francis %>% drop_na(), mapping = aes(x = Comp_Year, y = pred, color = 'Last Iter')) +
      geom_line(francis_run$start_mean_francis %>% drop_na(), mapping = aes(x = Comp_Year, y = pred, color = 'First Iter')) +
      facet_grid(Type~Fleet+Sex, labeller = labeller(
        Fleet = function(x) paste("Fleet", x, sep = " ")
      ), scales = 'free_y') +
      theme_sablefish() +
      labs(x = 'Year', y = 'Mean Composition Bin')
  }
  if(unique(data$FishAgeComps_Type[,1]) == 2){      # split sex comps
    plot_francis <- francis_run$end_mean_francis %>%
      drop_na() %>%
      ggplot() +
      geom_point(francis_run$end_mean_francis %>% drop_na(), mapping = aes(x = Comp_Year, y = obs)) +
      geom_line(francis_run$end_mean_francis %>% drop_na(), mapping = aes(x = Comp_Year, y = pred, color = 'Last Iter')) +
      geom_line(francis_run$start_mean_francis %>% drop_na(), mapping = aes(x = Comp_Year, y = pred, color = 'First Iter')) +
      facet_grid(Type~Fleet, labeller = labeller(
        Fleet = function(x) paste("Fleet", x, sep = " ")
      ), scales = 'free_y') +
      theme_sablefish() +
      labs(x = 'Year', y = 'Mean Composition Bin')
  }

  tiff(filename=here(root_folder,mod_name,"plots",paste0(mod_name,"_Francis_Wts.png")),units = "in", width=10,height=8, res = 300)  # these tiff commands weren't working correctly within a loop, so had to pull them out
  print(Francis_plot1)
  dev.off()
  tiff(filename=here(root_folder,mod_name,"plots",paste0(mod_name,"_Francis_mean.png")),units = "in", width=16,height=14, res = 300)
  print(plot_francis)
  dev.off()
}

##################################################################################################################################################################################
# Fit Model ---------------------------------------------------------------
##################################################################################################################################################################################

sabie_rtmb_model <- fit_model(data,
                              parameters,
                              mapping,
                              random = NULL,
                              newton_loops = 3,
                              silent = TRUE
)

SPoRC:::SPoRC_rtmb(parameters, data)

# Get standard error report
sabie_rtmb_model$sd_rep <- RTMB::sdreport(sabie_rtmb_model)

#write.csv(sabie_rtmb_model$sd_rep, here(root_folder,mod_name,"plots",paste0(mod_name,"_sd_rep.csv")))

# Save Model --------------------------------------------------------------

# Save data, parameters, and mapping in model object
sabie_rtmb_model$data <- data
sabie_rtmb_model$parameters <- parameters
sabie_rtmb_model$mapping <- mapping

# Write out RDS file
saveRDS(sabie_rtmb_model, file = here(root_folder,mod_name,paste0(mod_name,"_model_results.RDS")))
saveRDS(input_list, file = here(root_folder,mod_name, paste0(mod_name,"_input_list.RDS")))

# Check convergence, will let you know if high correlation, large gradient, or high SE for any pars
messages <- capture.output({
  post_optim_sanity_checks(sabie_rtmb_model$sd_rep,sabie_rtmb_model$rep,gradient_tol = 0.00001,se_tol = 5,corr_tol = 0.95)
}, type = "message")

writeLines(messages,  here(root_folder,mod_name,"plots", paste0(mod_name,"_convergence.txt")))

##################################################################################################################################################################################
# Reference Points --------------------------------------------------------
##################################################################################################################################################################################

ref_pts <- SPoRC::Get_Reference_Points(data = data,
                                       rep = sabie_rtmb_model$rep,
                                       SPR_x = 0.4,
                                       t_spwn = 0,
                                       sex_ratio_f = 0.5,
                                       calc_rec_st_yr = 20,
                                       rec_age = 2, type = "single_region",
                                       what = "SPR")

b40 <- ref_pts$b_ref_pt
f40 <- ref_pts$f_ref_pt
sabie_rtmb_model$b40 <- ref_pts$b_ref_pt
sabie_rtmb_model$f40 <- ref_pts$f_ref_pt

##################################################################################################################################################################################
# Catch Advice ------------------------------------------------------------
# NOTE: that these projections are just rough approximations for comparison purposes since don't use stochastic (inverse guassian simulated) recruitment--here just use mean recruitment
#       For official projections (and 7 scenarios needed for MSRA), then use the full Do_Population_Projection function (https://chengmatt.github.io/SPoRC/reference/Do_Population_Projection.html)
##################################################################################################################################################################################


# Define HCR to use (NPFMC F40)
HCR_function <- function(x, frp, brp, alpha = 0.05) {
  stock_status <- x / brp                                                      # define stock status
  if(stock_status >= 1) f <- frp                                               # If stock status is > 1, use full Fref
  if(stock_status > alpha && stock_status < 1){                                # If stock status is between brp and alpha (limit/lower BRP), then linearly scale down Fref
    f <- frp * (stock_status - alpha) / (1 - alpha)
  }
  if(stock_status < alpha) f <- 0                                              # If stock status is less than alpha (limit/lower BRP), then no fishing
  return(f)
}

# Pull necessary inputs from RTMB assessment run
t_spawn <- 0                                                                   # time of spawning for discounting mortality for SSB calcs (assume winter, Jan 1, spawning for sablefish)
n_proj_yrs <- 15                                                               # number of projection years to model
terminal_NAA <- array(sabie_rtmb_model$rep$NAA[,length(data$years),,],         # Pull the terminal year numbers at age by sex
                      dim = c(data$n_regions, length(data$ages),
                              data$n_sexes))
terminal_NAA0 <- array(sabie_rtmb_model$rep$NAA0[,length(data$years),,],         # Pull the terminal year numbers at age unfished by sex
                       dim = c(data$n_regions, length(data$ages),
                               data$n_sexes))
WAA <- array(rep(data$WAA[,length(data$years),,], each = n_proj_yrs),          # Pull WAA
             dim = c(data$n_regions, n_proj_yrs, length(data$ages),
                     data$n_sexes))
WAA_fish <- array(rep(data$WAA[,length(data$years),,], each = n_proj_yrs),          # Pull WAA for fisheries
                  dim = c(data$n_regions, n_proj_yrs, length(data$ages),
                          data$n_sexes, data$n_fish_fleets))
MatAA <- array(rep(data$MatAA[,length(data$years),,], each = n_proj_yrs),      # Pull maturity
               dim = c(data$n_regions, n_proj_yrs, length(data$ages),
                       data$n_sexes))
fish_sel <- array(rep(sabie_rtmb_model$rep$fish_sel[,length(data$years),,,],   # Pull estimated fishery selectivity
                      each = n_proj_yrs), dim = c(data$n_regions, n_proj_yrs,
                                                  length(data$ages), data$n_sexes, data$n_fish_fleets))
Movement <- array(rep(sabie_rtmb_model$rep$Movement[,,length(data$years),,],   # Pull estimated movement (not used for single region models)
                      each = n_proj_yrs), dim = c(data$n_regions,
                                                  data$n_regions, n_proj_yrs, length(data$ages), data$n_sexes))
terminal_F <- array(sabie_rtmb_model$rep$Fmort[,length(data$years),],          # Pull terminal year Fishing mortality
                    dim = c(data$n_regions, data$n_fish_fleets))
natmort <- array(rep(sabie_rtmb_model$rep$natmort[,length(data$years),,],      # Pull estimated natural mortality
                     each = n_proj_yrs), dim = c(data$n_regions, n_proj_yrs,
                                                 length(data$ages), data$n_sexes))
recruitment <- array(sabie_rtmb_model$rep$Rec[,20:(length(data$years) - 2)],   # Pull estimated recruitment time series from years 20 - terminal year-1 (corresponds to 1979 to 2023--don't use terminal year because don't estimate terminal recruitment)
                     dim = c(data$n_regions, length(20:(length(data$years)-2))))

# Sex ratio
sexratio <- array(dim = c(data$n_regions, n_proj_yrs, data$n_sexes)) # empty array
for(yr in 1:n_proj_yrs) sexratio[, yr, ] <- sabie_rtmb_model$rep$sexratio[,length(data$years),] # populate empty array

# Use FABC to project
f_ref_pt <- array(f40, dim = c(data$n_regions, n_proj_yrs))
b_ref_pt <- array(b40, dim = c(data$n_regions, n_proj_yrs))

# Do projection to get ABC
get_abc_vals <- Do_Population_Projection(n_proj_yrs = n_proj_yrs,              # projection years
                                         n_regions = data$n_regions,           # number of regions
                                         n_ages = length(data$ages),           # number of ages
                                         n_sexes = data$n_sexes,               # number of sexes
                                         sexratio = sexratio,                  # sex ratio
                                         n_fish_fleets = data$n_fish_fleets,   # number of fishery fleets
                                         do_recruits_move = 0,                 # whether recruits move, ==0 no movement
                                         recruitment = recruitment,            # recruitment vector
                                         terminal_NAA = terminal_NAA,          # terminal numbers at age
                                         terminal_NAA0 = terminal_NAA0,
                                         terminal_F = terminal_F,              # terminal F array
                                         natmort = natmort,                    # natural mortality array
                                         WAA = WAA,                            # weight at age array
                                         WAA_fish = WAA_fish,
                                         MatAA = MatAA,                        # maturity at age array
                                         fish_sel = fish_sel,                  # fishery selex array
                                         Movement = Movement,                  # movement array
                                         f_ref_pt = f_ref_pt,                  # fishery reference point
                                         b_ref_pt = b_ref_pt,                  # biological reference point
                                         HCR_function = HCR_function,          # HCR function that takes x, frp, and brp
                                         recruitment_opt = "mean_rec",         # recruitment assumption
                                         t_spawn = t_spawn                     # spawn fraction
)

# ABC for subsequent years (dimensions, region, year, fleet)
sum(get_abc_vals$proj_Catch[1,2,])                                             # year 2 = the actual first projected year

saveRDS(get_abc_vals, file = here(root_folder,mod_name,paste0(mod_name,"_ABC_proj.RDS")))


##################################################################################################################################################################################
# Some Basic Plots and Diagnostics ------------------------------------------------------------
##################################################################################################################################################################################

bad_pars <- which(!is.finite(sabie_rtmb_model$gr()))                           # check which pars have NA gradients

data <- sabie_rtmb_model$data
rep <- sabie_rtmb_model$rep
sd_rep <- sabie_rtmb_model$sd_rep

# Define dimensions
years <- SPoRC_one_reg$yrs
ages <- SPoRC_one_reg$ages
lens <- SPoRC_one_reg$length_bins


# get key quantities
reference_points_opt <- list(SPR_x = 0.4,
                             t_spwn = 0,
                             sex_ratio_f = 0.5,
                             calc_rec_st_yr = 20,
                             rec_age = 2,
                             type = "single_region",
                             what = "SPR"
)
proj_model_opt <- list(n_proj_yrs = 2,
                       n_avg_yrs = 1,
                       HCR_function = HCR_function,
                       recruitment_opt = 'mean_rec',
                       fmort_opt = 'HCR')
out <- get_key_quants(list(data), list(rep), reference_points_opt, proj_model_opt, 1)

out[[1]]  # key quantities data.frame

tiff(filename=here(root_folder,mod_name,"plots",paste0(mod_name,"_ABC.png")),units = "in", width=16,height=4, res = 300)
print(out[[2]]+ labs(title="Management Quantities"))  # table plot as cowplot object
dev.off()


# plot all basic info (basically just a wrapper function for a subset of functions defined above) ------------------------------------------------------------------
all_plots <- plot_all_basic(list(data), list(rep), list(sd_rep), mod_name, out_path = here(root_folder,mod_name,"plots"))

bios <- get_biological_plot(list(data), list(rep), mod_name)
tiff(filename=here(root_folder,mod_name,"plots",paste0(mod_name,"_M.png")),units = "in", width=12,height=8, res = 300)
print(bios[[2]]+
        ggthemes::scale_color_colorblind() +
        ggthemes::scale_fill_colorblind()+ labs(title="Natural Mortality"))
dev.off()


# SSB and Recruitment -----------------------------------------------------

# ssb estimates
ssb_plot_df <- reshape2::melt(rep$SSB) %>%
  dplyr::rename(Region = Var1, Year = Var2) %>%
  dplyr::bind_cols(se = sd_rep$sd[names(sd_rep$value) ==  "log(SSB)"]) %>%
  dplyr::mutate(lwr = exp(log(value) - 1.96 * se),
                upr = exp(log(value) + 1.96 * se), Region = paste("Region", Region), Type = "SSB (kt)")

# bio estimates
bio_plot_df <- reshape2::melt(rep$Total_Biom) %>%
  dplyr::rename(Region = Var1, Year = Var2) %>%
  dplyr::bind_cols(se = sd_rep$sd[names(sd_rep$value) ==  "log(Total_Biom)"]) %>%
  dplyr::mutate(lwr = exp(log(value) - 1.96 * se),
                upr = exp(log(value) + 1.96 * se), Region = paste("Region", Region), Type = "Biomass (kt)") %>%
  mutate(Year = Year + 1959)

# recruitment estimates
rec_plot_df <- reshape2::melt(rep$Rec) %>%
  dplyr::rename(Region = Var1,  Year = Var2) %>%
  dplyr::bind_cols(se = sd_rep$sd[names(sd_rep$value) == "log(Rec)"]) %>%
  dplyr::mutate(lwr = exp(log(value) - 1.96 * se),
                upr = exp(log(value) + 1.96 * se),
                Region = paste("Region", Region), Type = "Recruitment (mil. fish)")

# bind
biom_rec_df <- rbind(ssb_plot_df, rec_plot_df) %>%
  mutate(Year = Year + 1959)

# plot
ssb_recr <- ggplot(biom_rec_df, aes(x = Year, y = value, ymin = lwr, ymax = upr)) +
  geom_line(lwd = 0.9) +
  geom_ribbon(alpha = 0.3) +
  facet_wrap(~Type, scales = "free") + labs(x = "Year", y = "Value") +
  coord_cartesian(ylim = c(0,NA)) +
  theme_sablefish()

tiff(filename=here(root_folder,mod_name,"plots",paste0(mod_name,"_ssb_recr.png")),units = "in", width=10,height=6, res = 300)
print(ssb_recr)
dev.off()

bio <- ggplot(bio_plot_df, aes(x = Year, y = value, ymin = lwr, ymax = upr)) +
  geom_line(lwd = 0.9) +
  geom_ribbon(alpha = 0.3) +
  facet_wrap(~Type, scales = "free") + labs(x = "Year", y = "Value") +
  coord_cartesian(ylim = c(0,NA)) +
  labs(x = 'Year', y = 'Biomass (kt)') +
  theme_sablefish()

tiff(filename=here(root_folder,mod_name,"plots",paste0(mod_name,"_Biomass.png")),units = "in", width=10,height=6, res = 300)
print(bio)
dev.off()

# Numbers at Age ----------------------------------------------------------

naa_df <- reshape2::melt(rep$NAA) %>%
  rename(Region = Var1, Year = Var2, Age = Var3, Sex = Var4) %>%
  mutate(Year = Year + 1959,
         Sex = ifelse(Sex == 1, "Female", "Male"))

naa <- ggplot(naa_df, aes(x = Year, y = value, color = Sex)) +
  geom_line() +
  facet_wrap(~Age, scales = 'free') +
  labs(x = 'Year', y = 'Numbers-at-age (millions)') +
  theme_sablefish()+ labs(title="Numbers-at-Age")

tiff(filename=here(root_folder,mod_name,"plots",paste0(mod_name,"_naa.png")),units = "in", width=14,height=12, res = 300)
print(naa)
dev.off()

# Selex --------------------------------------------------------------

selex_plot <- get_selex_plot(list(sabie_rtmb_model$rep),  Selex_Type = "age",1, year_indx = c(1:length(years)))

tiff(filename=here(root_folder,mod_name,"plots",paste0(mod_name,"_selex_fsh.png")),units = "in", width=14,height=12, res = 300)
print(selex_plot[[1]]+ labs(title="Fishery Selectivity"))
dev.off()
tiff(filename=here(root_folder,mod_name,"plots",paste0(mod_name,"_selex_srvy.png")),units = "in", width=14,height=12, res = 300)
print(selex_plot[[2]]+ labs(title="Survey Selectivity"))
dev.off()

# fits to catch --------------------------------------------------------------

catch_plot <- get_catch_fits_plot(list(data), list(sabie_rtmb_model$rep), 1)

tiff(filename=here(root_folder,mod_name,"plots",paste0(mod_name,"_llf_catch_fit.png")),units = "in", width=14,height=12, res = 300)
print(catch_plot+labs(title="Catch"))
dev.off()


# Index Fits --------------------------------------------------------------

idx_fits <- get_idx_fits(data, rep, year_labs = years) %>%
  # Renaming stuff
  mutate(Index_Type =
           case_when(
             Category == "Fishery1, Q1" ~ "Japanese Fishery CPUE Index",
             Category %in% c("Fishery1, Q2", "Fishery1, Q3") ~ "Domestic Fishery CPUE Index",
             Category == "Survey1, Q1" ~ "BS LL Survey Relative Population Numbers",
             Category == "Survey2, Q1" ~ "AI LL Survey Relative Population Numbers",
             Category == "Survey3, Q1" ~ "GOA LL Survey Relative Population Numbers",
             Category == "Survey4, Q1" ~ "GOA Trawl Survey Biomass (kt)",
             Category == "Survey5, Q1" ~ "Japanese LL Survey Relative Population Numbers"
           ))

# index plot
idx_plot <- ggplot() +
  geom_line(idx_fits %>% filter(obs != 0),
            mapping = aes(x = Year, y = value), lwd = 1.3, col = 'red') +
  geom_pointrange(idx_fits %>% filter(obs != 0),
                  mapping = aes(x = Year, y = obs, ymin = lci, ymax= uci), color = 'blue') +
  labs(x = 'Year', y = 'Index') +
  facet_wrap(~Index_Type, scales = "free") +
  theme_sablefish() +
  theme(strip.text = element_text(
    size = 12))+labs(title="Index Fits")

tiff(filename=here(root_folder,mod_name,"plots",paste0(mod_name,"_idx_fit.png")),units = "in", width=14,height=12, res = 300)
print(idx_plot)
dev.off()

# Plot residuals and runs test
unique_idx <- unique(idx_fits$Index_Type)
runs_all <- data.frame()
for(i in 1:length(unique(idx_fits$Index_Type))) {
  tmp <- idx_fits %>% filter(Index_Type == unique_idx[i])
  runstest <- SPoRC::do_runs_test(x=as.numeric(tmp$resid),type="resid", mixing = "less")
  tmp_runs <- data.frame(p = runstest$p.runs, lwr = runstest$sig3lim[1], upr = runstest$sig3lim[2], Index_Type = unique_idx[i])
  runs_all <- rbind(runs_all, tmp_runs)
} # end i

runs_tst <- ggplot() +
  geom_point(idx_fits, mapping = aes(x = Year, y = resid)) +
  geom_segment(idx_fits, mapping = aes(x = Year, xend = Year, y = 0, yend = resid)) +
  geom_smooth(idx_fits, mapping = aes(x = Year, y = resid), se = F) +
  geom_hline(yintercept = 0, lty = 2) +
  geom_hline(runs_all, mapping = aes(yintercept = upr), lty = 2) +
  geom_hline(runs_all, mapping = aes(yintercept = lwr), lty = 2) +
  geom_text(data = runs_all, aes(x = -Inf, y = Inf, label = paste("p = ", round(p, 3))), hjust = -0.5, vjust = 2, size = 5)+
  labs(x = "Year", y = 'Residuals') +
  theme_bw(base_size = 20) +
  facet_wrap(~Index_Type, scales = 'free', ncol = 2)+
  theme(strip.text = element_text(
    size = 12))+labs(title="Index Residuals Runs Test")

tiff(filename=here(root_folder,mod_name,"plots",paste0(mod_name,"_idx_runs_test.png")),units = "in", width=14,height=12, res = 300)
print(runs_tst)
dev.off()

# if want to graph RPW
#  srv_baa <- rep$srv_q[1,,1] * rep$NAA[,-66,,] * exp(rep$ZAA[,,,] * 0.5) * rep$srv_sel[1,,,,1] * data$WAA[1,,,]
# srv_biom <- apply(srv_baa, 1, sum)


# Composition Fits --------------------------------------------------------

# Get composition fits
comp_fits <- get_comp_prop(data, rep, ages, lens, years)

### Fishery Ages ------------------------------------------------------------
# Female
# Fits by Year
llf_age_fem <- comp_fits$Fishery_Ages %>%
  filter(Fleet == 1, Sex == 1) %>%
  mutate(Cohort = Year - Age) %>%
  ggplot() +
  geom_col(aes(x = Age, y = obs), alpha = 0.5) +
  geom_line(aes(x = Age, y = pred)) +
  facet_wrap(~Year) +
  theme_sablefish() +
  labs(x = 'Age', y = 'Proportion')+labs(title="LLF Age Comps Female")

# Fits by Cohort
llf_age_coh_fem <- comp_fits$Fishery_Ages %>%
  filter(Fleet == 1, Sex == 1) %>%
  mutate(Cohort = Year - Age) %>%
  ggplot() +
  geom_col(aes(x = Age, y = obs), alpha = 0.5) +
  geom_line(aes(x = Age, y = pred)) +
  facet_wrap(~Cohort) +
  theme_sablefish() +
  labs(x = 'Age', y = 'Proportion')+labs(title="LLF Age Comps Cohort Fits Female")

# Male
# Fits by Year
if(unique(data$FishAgeComps_Type[,1]) > 0){
  llf_age_mal <- comp_fits$Fishery_Ages %>%
    filter(Fleet == 1, Sex == 2) %>%
    mutate(Cohort = Year - Age) %>%
    ggplot() +
    geom_col(aes(x = Age, y = obs), alpha = 0.5) +
    geom_line(aes(x = Age, y = pred)) +
    facet_wrap(~Year) +
    theme_sablefish() +
    labs(x = 'Age', y = 'Proportion')+labs(title="LLF Age Comps Male")
  
  # Fits by Cohort
  llf_age_coh_mal <- comp_fits$Fishery_Ages %>%
    filter(Fleet == 1, Sex == 2) %>%
    mutate(Cohort = Year - Age) %>%
    ggplot() +
    geom_col(aes(x = Age, y = obs), alpha = 0.5) +
    geom_line(aes(x = Age, y = pred)) +
    facet_wrap(~Cohort) +
    theme_sablefish() +
    labs(x = 'Age', y = 'Proportion')+labs(title="LLF Age Comps Cohort Fits Male")
  
}
### Fishery Lengths ----------------------------------------------------------
if(any(data$UseFishLenComps[,,1] == 1)){
  
  # Fits by Year - Fixed-Gear, Females
  llf_len_fem <- comp_fits$Fishery_Lens %>%
    filter(Fleet == 1, Sex == 1) %>%
    ggplot() +
    geom_col(aes(x = Len, y = obs), alpha = 0.5) +
    geom_line(aes(x = Len, y = pred)) +
    facet_wrap(~Year) +
    theme_sablefish() +
    labs(x = 'Length', y = 'Proportion')+labs(title="LLF Length Comps Female")
  
  
  # Fits by Year - Fixed-Gear, Males
  llf_len_male <- comp_fits$Fishery_Lens %>%
    filter(Fleet == 1, Sex == 2) %>%
    ggplot() +
    geom_col(aes(x = Len, y = obs), alpha = 0.5) +
    geom_line(aes(x = Len, y = pred)) +
    facet_wrap(~Year) +
    theme_sablefish() +
    labs(x = 'Length', y = 'Proportion')+labs(title="LLF Length Comps Male")
}

# Fits by Year - Trawl-Gear, Females
tf_len_fem <- comp_fits$Fishery_Lens %>%
  filter(Fleet == 2, Sex == 1) %>%
  ggplot() +
  geom_col(aes(x = Len, y = obs), alpha = 0.5) +
  geom_line(aes(x = Len, y = pred)) +
  facet_wrap(~Year) +
  theme_sablefish() +
  labs(x = 'Length', y = 'Proportion')+labs(title="TF Length Comps Female")


# Fits by Year - Trawl-Gear, Males
tf_len_male <- comp_fits$Fishery_Lens %>%
  filter(Fleet == 2, Sex == 2) %>%
  ggplot() +
  geom_col(aes(x = Len, y = obs), alpha = 0.5) +
  geom_line(aes(x = Len, y = pred)) +
  facet_wrap(~Year) +
  theme_sablefish() +
  labs(x = 'Length', y = 'Proportion')+labs(title="TF Length Comps Male")




### Survey Ages -------------------------------------------------------------

# Fits by Year - Domestic LL Survey
# Female
# BS
BS_lls_age_fem <- comp_fits$Survey_Ages %>%
  filter(Fleet == 1, Sex ==1) %>%
  mutate(Cohort = Year - Age) %>%
  ggplot() +
  geom_col(aes(x = Age, y = obs), alpha = 0.5) +
  geom_line(aes(x = Age, y = pred)) +
  facet_wrap(~Year) +
  theme_sablefish() +
  labs(x = 'Age', y = 'Proportion')+labs(title="BS LLS Age Comps Female")


# Fits by Cohort - Domestic LL Survey
BS_lls_age_coh_fem <- comp_fits$Survey_Ages %>%
  filter(Fleet == 1, Sex == 1) %>%
  mutate(Cohort = Year - Age) %>%
  ggplot() +
  geom_col(aes(x = Age, y = obs), alpha = 0.5) +
  geom_line(aes(x = Age, y = pred)) +
  facet_wrap(~Cohort) +
  theme_sablefish() +
  labs(x = 'Age', y = 'Proportion')+labs(title="BS LLS Age Comps Cohort Fit Female")


# Male
if(unique(data$SrvAgeComps_Type[,1]) > 0){
  BS_lls_age_mal <- comp_fits$Survey_Ages %>%
    filter(Fleet == 1, Sex ==2) %>%
    mutate(Cohort = Year - Age) %>%
    ggplot() +
    geom_col(aes(x = Age, y = obs), alpha = 0.5) +
    geom_line(aes(x = Age, y = pred)) +
    facet_wrap(~Year) +
    theme_sablefish() +
    labs(x = 'Age', y = 'Proportion')+labs(title="BS LLS Age Comps Male")
  
  
  # Fits by Cohort - Domestic LL Survey
  BS_lls_age_coh_mal <- comp_fits$Survey_Ages %>%
    filter(Fleet == 1, Sex == 2) %>%
    mutate(Cohort = Year - Age) %>%
    ggplot() +
    geom_col(aes(x = Age, y = obs), alpha = 0.5) +
    geom_line(aes(x = Age, y = pred)) +
    facet_wrap(~Cohort) +
    theme_sablefish() +
    labs(x = 'Age', y = 'Proportion')+labs(title="BS LLS Age Comps Cohort Fit Male")
  
}

# Fits by Year - Domestic LL Survey
# Female
# AI
AI_lls_age_fem <- comp_fits$Survey_Ages %>%
  filter(Fleet == 2, Sex ==1) %>%
  mutate(Cohort = Year - Age) %>%
  ggplot() +
  geom_col(aes(x = Age, y = obs), alpha = 0.5) +
  geom_line(aes(x = Age, y = pred)) +
  facet_wrap(~Year) +
  theme_sablefish() +
  labs(x = 'Age', y = 'Proportion')+labs(title="AI LLS Age Comps Female")


# Fits by Cohort - Domestic LL Survey
AI_lls_age_coh_fem <- comp_fits$Survey_Ages %>%
  filter(Fleet == 2, Sex == 1) %>%
  mutate(Cohort = Year - Age) %>%
  ggplot() +
  geom_col(aes(x = Age, y = obs), alpha = 0.5) +
  geom_line(aes(x = Age, y = pred)) +
  facet_wrap(~Cohort) +
  theme_sablefish() +
  labs(x = 'Age', y = 'Proportion')+labs(title="AI LLS Age Comps Cohort Fit Female")


# Male
if(unique(data$SrvAgeComps_Type[,2]) > 0){
  AI_lls_age_mal <- comp_fits$Survey_Ages %>%
    filter(Fleet == 2, Sex ==2) %>%
    mutate(Cohort = Year - Age) %>%
    ggplot() +
    geom_col(aes(x = Age, y = obs), alpha = 0.5) +
    geom_line(aes(x = Age, y = pred)) +
    facet_wrap(~Year) +
    theme_sablefish() +
    labs(x = 'Age', y = 'Proportion')+labs(title="AI LLS Age Comps Male")
  
  
  # Fits by Cohort - Domestic LL Survey
  AI_lls_age_coh_mal <- comp_fits$Survey_Ages %>%
    filter(Fleet == 2, Sex == 2) %>%
    mutate(Cohort = Year - Age) %>%
    ggplot() +
    geom_col(aes(x = Age, y = obs), alpha = 0.5) +
    geom_line(aes(x = Age, y = pred)) +
    facet_wrap(~Cohort) +
    theme_sablefish() +
    labs(x = 'Age', y = 'Proportion')+labs(title="AI LLS Age Comps Cohort Fit Male")
  
}

# Fits by Year - Domestic LL Survey
# Female
# GOA
GOA_lls_age_fem <- comp_fits$Survey_Ages %>%
  filter(Fleet == 3, Sex ==1) %>%
  mutate(Cohort = Year - Age) %>%
  ggplot() +
  geom_col(aes(x = Age, y = obs), alpha = 0.5) +
  geom_line(aes(x = Age, y = pred)) +
  facet_wrap(~Year) +
  theme_sablefish() +
  labs(x = 'Age', y = 'Proportion')+labs(title="GOA LLS Age Comps Female")


# Fits by Cohort - Domestic LL Survey
GOA_lls_age_coh_fem <- comp_fits$Survey_Ages %>%
  filter(Fleet == 3, Sex == 1) %>%
  mutate(Cohort = Year - Age) %>%
  ggplot() +
  geom_col(aes(x = Age, y = obs), alpha = 0.5) +
  geom_line(aes(x = Age, y = pred)) +
  facet_wrap(~Cohort) +
  theme_sablefish() +
  labs(x = 'Age', y = 'Proportion')+labs(title="GOA LLS Age Comps Cohort Fit Female")


# Male
if(unique(data$SrvAgeComps_Type[,3]) > 0){
  GOA_lls_age_mal <- comp_fits$Survey_Ages %>%
    filter(Fleet == 3, Sex ==2) %>%
    mutate(Cohort = Year - Age) %>%
    ggplot() +
    geom_col(aes(x = Age, y = obs), alpha = 0.5) +
    geom_line(aes(x = Age, y = pred)) +
    facet_wrap(~Year) +
    theme_sablefish() +
    labs(x = 'Age', y = 'Proportion')+labs(title="GOA LLS Age Comps Male")
  
  
  # Fits by Cohort - Domestic LL Survey
  GOA_lls_age_coh_mal <- comp_fits$Survey_Ages %>%
    filter(Fleet == 3, Sex == 2) %>%
    mutate(Cohort = Year - Age) %>%
    ggplot() +
    geom_col(aes(x = Age, y = obs), alpha = 0.5) +
    geom_line(aes(x = Age, y = pred)) +
    facet_wrap(~Cohort) +
    theme_sablefish() +
    labs(x = 'Age', y = 'Proportion')+labs(title="GOA LLS Age Comps Cohort Fit Male")
  
}


# Fits by Year - Japanese LL Survey
# Female
jpn_lls_age_fem <- comp_fits$Survey_Ages %>%
  filter(Fleet == 5, Sex == 1) %>%
  mutate(Cohort = Year - Age) %>%
  ggplot() +
  geom_col(aes(x = Age, y = obs), alpha = 0.5) +
  geom_line(aes(x = Age, y = pred)) +
  facet_wrap(~Year) +
  theme_sablefish() +
  labs(x = 'Age', y = 'Proportion')+labs(title="JPN LLS Age Comps Female")


# Fits by Cohort - Japanese LL Survey
jpn_lls_age_coh_fem <- comp_fits$Survey_Ages %>%
  filter(Fleet == 5, Sex == 1) %>%
  mutate(Cohort = Year - Age) %>%
  ggplot() +
  geom_col(aes(x = Age, y = obs), alpha = 0.5) +
  geom_line(aes(x = Age, y = pred)) +
  facet_wrap(~Cohort) +
  theme_sablefish() +
  labs(x = 'Age', y = 'Proportion')+labs(title="JPN LLS Age Comps Cohort Fit Female")


# Male
if(unique(data$SrvAgeComps_Type[,5]) > 0){
  jpn_lls_age_mal <- comp_fits$Survey_Ages %>%
    filter(Fleet == 5, Sex == 2) %>%
    mutate(Cohort = Year - Age) %>%
    ggplot() +
    geom_col(aes(x = Age, y = obs), alpha = 0.5) +
    geom_line(aes(x = Age, y = pred)) +
    facet_wrap(~Year) +
    theme_sablefish() +
    labs(x = 'Age', y = 'Proportion')+labs(title="JPN LLS Age Comps Male")
  
  
  # Fits by Cohort - Japanese LL Survey
  jpn_lls_age_coh_mal <- comp_fits$Survey_Ages %>%
    filter(Fleet == 5, Sex == 2) %>%
    mutate(Cohort = Year - Age) %>%
    ggplot() +
    geom_col(aes(x = Age, y = obs), alpha = 0.5) +
    geom_line(aes(x = Age, y = pred)) +
    facet_wrap(~Cohort) +
    theme_sablefish() +
    labs(x = 'Age', y = 'Proportion')+labs(title="JPN LLS Age Comps Cohort Fit Male")
  
}

### Survey Lengths ----------------------------------------------------------
# Fits by Year - Domestic LL Survey, 
# Females
# BS
if(any(data$UseSrvLenComps[,,1] == 1)){
  BS_lls_len_fem <- comp_fits$Survey_Lens %>%
    filter(Fleet == 1, Sex == 1) %>%
    ggplot() +
    geom_col(aes(x = Len, y = obs), alpha = 0.5) +
    geom_line(aes(x = Len, y = pred)) +
    facet_wrap(~Year) +
    theme_sablefish() +
    labs(x = 'Length', y = 'Proportion')+labs(title="BS LLS Length Comps Female")
  
  
  # Fits by Year - Domestic LL Survey, Males
  BS_lls_len_male <- comp_fits$Survey_Lens %>%
    filter(Fleet == 1, Sex == 2) %>%
    ggplot() +
    geom_col(aes(x = Len, y = obs), alpha = 0.5) +
    geom_line(aes(x = Len, y = pred)) +
    facet_wrap(~Year) +
    theme_sablefish() +
    labs(x = 'Length', y = 'Proportion')+labs(title="BS LLS Length Comps Male")
}

# Females
# AI
if(any(data$UseSrvLenComps[,,2] == 1)){
  AI_lls_len_fem <- comp_fits$Survey_Lens %>%
    filter(Fleet == 2, Sex == 1) %>%
    ggplot() +
    geom_col(aes(x = Len, y = obs), alpha = 0.5) +
    geom_line(aes(x = Len, y = pred)) +
    facet_wrap(~Year) +
    theme_sablefish() +
    labs(x = 'Length', y = 'Proportion')+labs(title="AI LLS Length Comps Female")
  
  
  # Fits by Year - Domestic LL Survey, Males
  AI_lls_len_male <- comp_fits$Survey_Lens %>%
    filter(Fleet == 2, Sex == 2) %>%
    ggplot() +
    geom_col(aes(x = Len, y = obs), alpha = 0.5) +
    geom_line(aes(x = Len, y = pred)) +
    facet_wrap(~Year) +
    theme_sablefish() +
    labs(x = 'Length', y = 'Proportion')+labs(title="AI LLS Length Comps Male")
}

# Females
# GOA
if(any(data$UseSrvLenComps[,,3] == 1)){
  GOA_lls_len_fem <- comp_fits$Survey_Lens %>%
    filter(Fleet == 3, Sex == 1) %>%
    ggplot() +
    geom_col(aes(x = Len, y = obs), alpha = 0.5) +
    geom_line(aes(x = Len, y = pred)) +
    facet_wrap(~Year) +
    theme_sablefish() +
    labs(x = 'Length', y = 'Proportion')+labs(title="GOA LLS Length Comps Female")
  
  
  # Fits by Year - Domestic LL Survey, Males
  GOA_lls_len_male <- comp_fits$Survey_Lens %>%
    filter(Fleet == 3, Sex == 2) %>%
    ggplot() +
    geom_col(aes(x = Len, y = obs), alpha = 0.5) +
    geom_line(aes(x = Len, y = pred)) +
    facet_wrap(~Year) +
    theme_sablefish() +
    labs(x = 'Length', y = 'Proportion')+labs(title="GOA LLS Length Comps Male")
}

if(any(data$UseSrvLenComps[,,4] == 1)){
  # Fits by Year - GOA Trawl Survey, Females
  ts_len_fem <- comp_fits$Survey_Lens %>%
    filter(Fleet == 4, Sex == 1) %>%
    ggplot() +
    geom_col(aes(x = Len, y = obs), alpha = 0.5) +
    geom_line(aes(x = Len, y = pred)) +
    facet_wrap(~Year) +
    theme_sablefish() +
    labs(x = 'Length', y = 'Proportion')+labs(title="GOA TS Length Comps Female")
  
  
  # Fits by Year - GOA Trawl Survey, Males
  ts_len_male <- comp_fits$Survey_Lens %>%
    filter(Fleet == 4, Sex == 2) %>%
    ggplot() +
    geom_col(aes(x = Len, y = obs), alpha = 0.5) +
    geom_line(aes(x = Len, y = pred)) +
    facet_wrap(~Year) +
    theme_sablefish() +
    labs(x = 'Length', y = 'Proportion')+labs(title="GOA TS LEngth Comps Male")
}


# Fits by Year - Japanese LL Survey, Females
if(any(data$UseSrvLenComps[,,5] == 1)){
  jpn_lls_len_fem <- comp_fits$Survey_Lens %>%
    filter(Fleet == 5, Sex == 1) %>%
    ggplot() +
    geom_col(aes(x = Len, y = obs), alpha = 0.5) +
    geom_line(aes(x = Len, y = pred)) +
    facet_wrap(~Year) +
    theme_sablefish() +
    labs(x = 'Length', y = 'Proportion')+labs(title="JPN LLS Length Comps Female")
  
  # Fits by Year - Japanese LL Survey, Males
  jpn_lls_len_male <- comp_fits$Survey_Lens %>%
    filter(Fleet == 5, Sex == 2) %>%
    ggplot() +
    geom_col(aes(x = Len, y = obs), alpha = 0.5) +
    geom_line(aes(x = Len, y = pred)) +
    facet_wrap(~Year) +
    theme_sablefish() +
    labs(x = 'Length', y = 'Proportion')+labs(title="JPN LLS Length Comps Male")
}


# Get OSA residuals ---------------------------------------------------------------------------------------------------------------------------------------

# NOTE: If do Francis reweighting need to be careful to remove any NA values introduced by using the which == 1 argument in teh N array

if(do_OSA == 1){
  ### Fishery Ages ------------------------------------------------------------
  #### Domestic Fishery LL Ages -------------------------------------------------
  
  if(unique(data$FishAgeComps_Type[,1]) == 0){
    osa_results_llfa <- SPoRC::get_osa(obs_mat = comp_fits$Obs_FishAge_mat + 1e-4,
                                       exp_mat = comp_fits$Pred_FishAge_mat + 1e-4,
                                       N = array(c(na.omit(data$ISS_FishAgeComps[1,which(data$UseFishAgeComps[,,1] == 1),1,1] *
                                                             unique(data$Wt_FishAgeComps[1,which(data$UseFishAgeComps[,,1] == 1),1,1])))),
                                       years = which(data$UseFishAgeComps[,,1] == 1),
                                       fleet = 1,
                                       bins = ages,
                                       comp_type = 0,
                                       bin_label = "Ages")
  }
  
  if(unique(data$FishAgeComps_Type[,1]) == 1){
    osa_results_llfa <- SPoRC::get_osa(obs_mat = comp_fits$Obs_FishAge_mat + 1e-4,
                                       exp_mat = comp_fits$Pred_FishAge_mat + 1e-4,
                                       N = array(c(na.omit(data$ISS_FishAgeComps[1,which(data$UseFishAgeComps[,,1] == 1),1,1] *
                                                             unique(data$Wt_FishAgeComps[1,which(data$UseFishAgeComps[,,1] == 1),1,1])),
                                                   na.omit(data$ISS_FishAgeComps[1,which(data$UseFishAgeComps[,,1] == 1),2,1] *
                                                             unique(data$Wt_FishAgeComps[1,which(data$UseFishAgeComps[,,1] == 1),2,1]))),
                                                 dim = c(data$n_regions, length(which(data$UseFishAgeComps[,,1] == 1)),
                                                         data$n_sexes)),
                                       years = which(data$UseFishAgeComps[,,1] == 1),
                                       fleet = 1,
                                       bins = ages,
                                       comp_type = 1,                              # ==0 is aggregated comps, ==1 is split comps, ==2 is joint comps
                                       bin_label = "Ages")
  }
  if(unique(data$FishAgeComps_Type[,1]) == 2){
    osa_results_llfa <- SPoRC::get_osa(obs_mat = comp_fits$Obs_FishAge_mat + 1e-4,
                                       exp_mat = comp_fits$Pred_FishAge_mat + 1e-4,
                                       N = array(c(na.omit(data$ISS_FishAgeComps[1,which(data$UseFishAgeComps[,,1] == 1),1,1] * 
                                                             unique(data$Wt_FishAgeComps[1,which(data$UseFishAgeComps[,,1] == 1),1,1]))),
                                                 dim = c(data$n_regions, length(which(data$UseFishAgeComps[,,1] == 1)))),
                                       years = which(data$UseFishAgeComps[,,1] == 1),
                                       fleet = 1,
                                       bins = ages,
                                       comp_type = 2,                              # ==0 is aggregated comps, ==1 is split comps, ==2 is joint comps
                                       bin_label = "Ages")
  }
  # Get SDNR QQ plot and Bubble Plot
  resid_plot_llfa <- SPoRC::plot_resids(osa_results = osa_results_llfa)
  
  # Get Aggregated Plot
  fishage_agg <- comp_fits$Fishery_Ages %>%
    group_by(Region, Age, Sex, Fleet) %>%
    summarize(obs = mean(obs), pred = mean(pred)) %>%
    filter(Fleet == 1) %>%
    ggplot() +
    geom_col(aes(x = Age, y = obs), fill = 'darkgreen', alpha = 0.8) +
    geom_line(aes(x = Age, y = pred), col = 'black', lwd = 1.3) +
    theme_bw(base_size = 15) +
    facet_grid( Region~Sex, labeller = labeller(
      Region = function(x) paste0("Region ", x),
      Sex = function(x) paste0("Sex ", x)
    )) +
    labs(x = 'Age', y = 'Proportions')
  
  osa_llf_age <- cowplot::plot_grid(resid_plot_llfa[[2]], # Bubble
                                    cowplot::plot_grid(resid_plot_llfa[[1]], fishage_agg, nrow = 1), # QQ and Aggregated
                                    ncol = 1, rel_heights = c(0.6, 0.4), labels = "LLF OSA Age Comps")
   
   
  #### Domestic Fixed-Gear Fishery Length --------------------------------------------------------
  # Get OSA residuals - need to add small constant to allow OSA calculation for multinomial
  if(any(data$UseFishLenComps[,,1] == 1)){
    
    if(unique(data$FishLenComps_Type[,1]) == 1){
      osa_results_llfl <- SPoRC::get_osa(obs_mat = comp_fits$Obs_FishLen_mat + 1e-4,
                                         exp_mat = comp_fits$Pred_FishLen_mat + 1e-4,
                                         N = array(c(na.omit(data$ISS_FishLenComps[1,which(data$UseFishLenComps[,,1] == 1),1,1] *
                                                               unique(data$Wt_FishLenComps[1,which(data$UseFishLenComps[,,1] == 1),1,1])),
                                                     na.omit(data$ISS_FishLenComps[1,which(data$UseFishLenComps[,,1] == 1),2,1] *
                                                               unique(data$Wt_FishLenComps[1,which(data$UseFishLenComps[,,1] == 1),2,1]))),
                                                   dim = c(data$n_regions, length(which(data$UseFishLenComps[,,1] == 1)),
                                                           data$n_sexes)),
                                         years = which(data$UseFishLenComps[,,1] == 1),
                                         fleet = 1,
                                         bins = lens,
                                         comp_type = 1,
                                         bin_label = "Lengths (cm)")
    }
    if(unique(data$FishLenComps_Type[,1]) == 2){
      osa_results_llfl <- SPoRC::get_osa(obs_mat = comp_fits$Obs_FishLen_mat + 1e-4,
                                         exp_mat = comp_fits$Pred_FishLen_mat + 1e-4,
                                         N = array(c(na.omit(data$ISS_FishLenComps[1,which(data$UseFishLenComps[,,1] == 1),1,1] *
                                                               unique(data$Wt_FishLenComps[1,which(data$UseFishLenComps[,,1] == 1),1,1]))),
                                                   dim = c(data$n_regions, length(which(data$UseFishLenComps[,,1] == 1)))),
                                         years = which(data$UseFishLenComps[,,1] == 1),
                                         fleet = 1,
                                         bins = lens,
                                         comp_type = 2,
                                         bin_label = "Lengths (cm)")
    }
    
    
    # Get SDNR QQ plot and Bubble Plot
    resid_plot_llfl <- SPoRC::plot_resids(osa_results = osa_results_llfl)
    
    # Get Aggregated Plot
    fishlens_agg <- comp_fits$Fishery_Lens %>%
      group_by(Region, Len, Sex, Fleet) %>%
      summarize(obs = mean(obs), pred = mean(pred)) %>%
      filter(Fleet == 1) %>%
      ggplot() +
      geom_col(aes(x = Len, y = obs), fill = 'darkgreen', alpha = 0.8) +
      geom_line(aes(x = Len, y = pred), col = 'black', lwd = 1.3) +
      theme_bw(base_size = 15) +
      facet_grid( Region~Sex, labeller = labeller(
        Region = function(x) paste0("Region ", x),
        Sex = function(x) paste0("Sex ", x)
      )) +
      theme_bw(base_size = 20) +
      labs(x = 'Lengths (cm)', y = 'Proportions')
    
    osa_llf_len <- cowplot::plot_grid(resid_plot_llfl[[2]], # Bubble
                                      cowplot::plot_grid(resid_plot_llfl[[1]], fishlens_agg, nrow = 1), # QQ and Aggregated
                                      ncol = 1, rel_heights = c(0.6, 0.4), labels = "LLF OSA Length Comps")
    
  }
   
  #### Domestic Trawl Fishery Lengths --------------------------------------------------
  # Get OSA residuals - need to add small constant to allow OSA calculation for multinomial
  
  if(unique(data$FishLenComps_Type[,2]) == 1){
    osa_results_tf <- SPoRC::get_osa(obs_mat = comp_fits$Obs_FishLen_mat + 1e-4,
                                     exp_mat = comp_fits$Pred_FishLen_mat + 1e-4,
                                     N = array(c(na.omit(data$ISS_FishLenComps[1,which(data$UseFishLenComps[,,2] == 1),1,2] *
                                                           unique(data$Wt_FishLenComps[1,which(data$UseFishLenComps[,,2] == 1),1,2])),
                                                 na.omit(data$ISS_FishLenComps[1,which(data$UseFishLenComps[,,2] == 1),2,2] *
                                                           unique(data$Wt_FishLenComps[1,which(data$UseFishLenComps[,,2] == 1),2,2]))),
                                               dim = c(data$n_regions, length(which(data$UseFishLenComps[,,2] == 1)),
                                                       data$n_sexes)),
                                     years = which(data$UseFishLenComps[,,2] == 1),
                                     fleet = 2,
                                     bins = lens,
                                     comp_type = 1,
                                     bin_label = "Lengths (cm)")
  }
  if(unique(data$FishLenComps_Type[,2]) == 2){
    osa_results_tf <- SPoRC::get_osa(obs_mat = comp_fits$Obs_FishLen_mat + 1e-4,
                                     exp_mat = comp_fits$Pred_FishLen_mat + 1e-4,
                                     N = array(c(na.omit(data$ISS_FishLenComps[1,which(data$UseFishLenComps[,,2] == 1),1,2] *
                                                           unique(data$Wt_FishLenComps[1,which(data$UseFishLenComps[,,2] == 1),1,2]))),
                                               dim = c(data$n_regions, length(which(data$UseFishLenComps[,,2] == 1)))),
                                     years = which(data$UseFishLenComps[,,2] == 1),
                                     fleet = 2,
                                     bins = lens,
                                     comp_type = 2,
                                     bin_label = "Lengths (cm)")
  }
  # Get SDNR QQ plot and Bubble Plot
  resid_plot_tf <- SPoRC::plot_resids(osa_results = osa_results_tf)
  
  # Get Aggregated Plot
  trwlens_agg <- comp_fits$Fishery_Lens %>%
    group_by(Region, Len, Sex, Fleet) %>%
    summarize(obs = mean(obs), pred = mean(pred)) %>%
    filter(Fleet == 2) %>%
    ggplot() +
    geom_col(aes(x = Len, y = obs), fill = 'darkgreen', alpha = 0.8) +
    geom_line(aes(x = Len, y = pred), col = 'black', lwd = 1.3) +
    theme_bw(base_size = 15) +
    facet_grid( Region~Sex, labeller = labeller(
      Region = function(x) paste0("Region ", x),
      Sex = function(x) paste0("Sex ", x)
    )) +
    theme_bw(base_size = 20) +
    labs(x = 'Lengths (cm)', y = 'Proportions')
  
  osa_tf <- cowplot::plot_grid(resid_plot_tf[[2]], # Bubble
                               cowplot::plot_grid(resid_plot_tf[[1]], trwlens_agg, nrow = 1), # QQ and Aggregated
                               ncol = 1, rel_heights = c(0.6, 0.4), labels = "TF OSA Length Comps")
  
  

  ### Survey Ages -------------------------------------------------------------
  #### Domestic Survey LL Ages -------------------------------------------------
  # Get OSA residuals - need to add small constant to allow OSA calculation for multinomial
  
  # BS
  
  if(unique(data$SrvAgeComps_Type[,1]) == 2){
    osa_results_llsa_BS <- SPoRC::get_osa(obs_mat = comp_fits$Obs_SrvAge_mat + 1e-4,
                                          exp_mat = comp_fits$Pred_SrvAge_mat + 1e-4,
                                          N = array(c(na.omit(data$ISS_SrvAgeComps[1,which(data$UseSrvAgeComps[,,1] == 1),1,1] *
                                                                unique(data$Wt_SrvAgeComps[1,which(data$UseSrvAgeComps[,,1] == 1),1,1]))),
                                                    dim = c(data$n_regions, length(which(data$UseSrvAgeComps[,,1] == 1)))),
                                          years = which(data$UseSrvAgeComps[,,1] == 1),
                                          fleet = 1,
                                          bins = ages,
                                          comp_type = 2,                              # ==0 is aggregated comps, ==1 is split comps, ==2 is joint comps
                                          bin_label = "Ages")
  }
  # Get SDNR QQ plot and Bubble Plot
  resid_plot_llsa_BS <- SPoRC::plot_resids(osa_results = osa_results_llsa_BS)
  
  # Get Aggregated Plot
  srvages_agg_BS <- comp_fits$Survey_Ages %>%
    group_by(Region, Age, Sex, Fleet) %>%
    summarize(obs = mean(obs), pred = mean(pred)) %>%
    filter(Fleet == 1) %>%
    ggplot() +
    geom_col(aes(x = Age, y = obs), fill = 'darkgreen', alpha = 0.8) +
    geom_line(aes(x = Age, y = pred), col = 'black', lwd = 1.3) +
    theme_bw(base_size = 15) +
    facet_grid( Region~Sex, labeller = labeller(
      Region = function(x) paste0("Region ", x),
      Sex = function(x) paste0("Sex ", x)
    )) +
    theme_bw(base_size = 20) +
    labs(x = 'Ages', y = 'Proportions')
  
  osa_lls_age_BS <- cowplot::plot_grid(resid_plot_llsa_BS[[2]], # Bubble
                                       cowplot::plot_grid(resid_plot_llsa_BS[[1]], srvages_agg_BS, nrow = 1), # QQ and Aggregated
                                       ncol = 1, rel_heights = c(0.6, 0.4), labels = "BS LLS OSA Age Comps")
  

  # AI
  
  if(unique(data$SrvAgeComps_Type[,2]) == 2){
    osa_results_llsa_AI <- SPoRC::get_osa(obs_mat = comp_fits$Obs_SrvAge_mat + 1e-4,
                                          exp_mat = comp_fits$Pred_SrvAge_mat + 1e-4,
                                          N = array(c(na.omit(data$ISS_SrvAgeComps[1,which(data$UseSrvAgeComps[,,2] == 1),1,2] *
                                                                unique(data$Wt_SrvAgeComps[1,which(data$UseSrvAgeComps[,,2] == 1),1,2]))),
                                                    dim = c(data$n_regions, length(which(data$UseSrvAgeComps[,,2] == 1)))),
                                          years = which(data$UseSrvAgeComps[,,2] == 1),
                                          fleet = 2,
                                          bins = ages,
                                          comp_type = 2,                              # ==0 is aggregated comps, ==1 is split comps, ==2 is joint comps
                                          bin_label = "Ages")
  }
  # Get SDNR QQ plot and Bubble Plot
  resid_plot_llsa_AI <- SPoRC::plot_resids(osa_results = osa_results_llsa_AI)
  
  # Get Aggregated Plot
  srvages_agg_AI <- comp_fits$Survey_Ages %>%
    group_by(Region, Age, Sex, Fleet) %>%
    summarize(obs = mean(obs), pred = mean(pred)) %>%
    filter(Fleet == 2) %>%
    ggplot() +
    geom_col(aes(x = Age, y = obs), fill = 'darkgreen', alpha = 0.8) +
    geom_line(aes(x = Age, y = pred), col = 'black', lwd = 1.3) +
    theme_bw(base_size = 15) +
    facet_grid( Region~Sex, labeller = labeller(
      Region = function(x) paste0("Region ", x),
      Sex = function(x) paste0("Sex ", x)
    )) +
    theme_bw(base_size = 20) +
    labs(x = 'Ages', y = 'Proportions')
  
  osa_lls_age_AI <- cowplot::plot_grid(resid_plot_llsa_AI[[2]], # Bubble
                                       cowplot::plot_grid(resid_plot_llsa_AI[[1]], srvages_agg_AI, nrow = 1), # QQ and Aggregated
                                       ncol = 1, rel_heights = c(0.6, 0.4), labels = "AI LLS OSA Age Comps")
  

  # GOA
  
  if(unique(data$SrvAgeComps_Type[,3]) == 2){
    osa_results_llsa_GOA <- SPoRC::get_osa(obs_mat = comp_fits$Obs_SrvAge_mat + 1e-4,
                                           exp_mat = comp_fits$Pred_SrvAge_mat + 1e-4,
                                           N = array(c(na.omit(data$ISS_SrvAgeComps[1,which(data$UseSrvAgeComps[,,3] == 1),1,3] *
                                                                 unique(data$Wt_SrvAgeComps[1,which(data$UseSrvAgeComps[,,3] == 1),1,3]))),
                                                     dim = c(data$n_regions, length(which(data$UseSrvAgeComps[,,3] == 1)))),
                                           years = which(data$UseSrvAgeComps[,,3] == 1),
                                           fleet = 3,
                                           bins = ages,
                                           comp_type = 2,                              # ==0 is aggregated comps, ==1 is split comps, ==2 is joint comps
                                           bin_label = "Ages")
  }
  # Get SDNR QQ plot and Bubble Plot
  resid_plot_llsa_GOA <- SPoRC::plot_resids(osa_results = osa_results_llsa_GOA)
  
  # Get Aggregated Plot
  srvages_agg_GOA <- comp_fits$Survey_Ages %>%
    group_by(Region, Age, Sex, Fleet) %>%
    summarize(obs = mean(obs), pred = mean(pred)) %>%
    filter(Fleet == 3) %>%
    ggplot() +
    geom_col(aes(x = Age, y = obs), fill = 'darkgreen', alpha = 0.8) +
    geom_line(aes(x = Age, y = pred), col = 'black', lwd = 1.3) +
    theme_bw(base_size = 15) +
    facet_grid( Region~Sex, labeller = labeller(
      Region = function(x) paste0("Region ", x),
      Sex = function(x) paste0("Sex ", x)
    )) +
    theme_bw(base_size = 20) +
    labs(x = 'Ages', y = 'Proportions')
  
  osa_lls_age_GOA <- cowplot::plot_grid(resid_plot_llsa_GOA[[2]], # Bubble
                                        cowplot::plot_grid(resid_plot_llsa_GOA[[1]], srvages_agg_GOA, nrow = 1), # QQ and Aggregated
                                        ncol = 1, rel_heights = c(0.6, 0.4), labels = "GOA LLS OSA Age Comps")
  
  
  ### Survey Lengths ----------------------------------------------------------
  #### Domestic LL Survey Lengths ----------------------------------------------
  # Get OSA residuals - need to add small constant to allow OSA calculation for multinomial
  
  # BS
  
  if(any(data$UseSrvLenComps[,,1] == 1)){
    if(unique(data$SrvLenComps_Type[,1]) == 2){
      osa_results_llsl_BS <- SPoRC::get_osa(obs_mat = comp_fits$Obs_SrvLen_mat + 1e-4,
                                            exp_mat = comp_fits$Pred_SrvLen_mat + 1e-4,
                                            N = array(c(na.omit(data$ISS_SrvLenComps[1,which(data$UseSrvLenComps[,,1] == 1),1,1] *
                                                                  unique(data$Wt_SrvLenComps[1,which(data$UseSrvLenComps[,,1] == 1),1,1]))),
                                                      dim = c(data$n_regions, length(which(data$UseSrvLenComps[,,1] == 1)))),
                                            years = which(data$UseSrvLenComps[,,1] == 1),
                                            fleet = 1,
                                            bins = lens,
                                            comp_type = 2,
                                            bin_label = "Lengths (cm)")
    }
    
    
    # Get SDNR QQ plot and Bubble Plot
    resid_plot_llsl_BS <- SPoRC::plot_resids(osa_results = osa_results_llsl_BS)
    
    # Get Aggregated Plot
    srvdom_lens_agg_BS <- comp_fits$Survey_Lens %>%
      group_by(Region, Len, Sex, Fleet) %>%
      summarize(obs = mean(obs), pred = mean(pred)) %>%
      filter(Fleet == 1) %>%
      ggplot() +
      geom_col(aes(x = Len, y = obs), fill = 'darkgreen', alpha = 0.8) +
      geom_line(aes(x = Len, y = pred), col = 'black', lwd = 1.3) +
      theme_bw(base_size = 15) +
      facet_grid( Region~Sex, labeller = labeller(
        Region = function(x) paste0("Region ", x),
        Sex = function(x) paste0("Sex ", x)
      )) +
      theme_bw(base_size = 20) +
      labs(x = 'Lengths (cm)', y = 'Proportions')
    
    osa_lls_len_BS <- cowplot::plot_grid(resid_plot_llsl_BS[[2]], # Bubble
                                         cowplot::plot_grid(resid_plot_llsl_BS[[1]], srvdom_lens_agg_BS, nrow = 1), # QQ and Aggregated
                                         ncol = 1, rel_heights = c(0.6, 0.4), labels = "BS LLS OSA Length Comps")
    
  }
  # AI
  
  if(any(data$UseSrvLenComps[,,2] == 1)){
    if(unique(data$SrvLenComps_Type[,2]) == 2){
      osa_results_llsl_AI <- SPoRC::get_osa(obs_mat = comp_fits$Obs_SrvLen_mat + 1e-4,
                                            exp_mat = comp_fits$Pred_SrvLen_mat + 1e-4,
                                            N = array(c(na.omit(data$ISS_SrvLenComps[1,which(data$UseSrvLenComps[,,2] == 1),1,2] *
                                                                  unique(data$Wt_SrvLenComps[1,which(data$UseSrvLenComps[,,2] == 1),1,2]))),
                                                      dim = c(data$n_regions, length(which(data$UseSrvLenComps[,,2] == 1)))),
                                            years = which(data$UseSrvLenComps[,,2] == 1),
                                            fleet = 2,
                                            bins = lens,
                                            comp_type = 2,
                                            bin_label = "Lengths (cm)")
    }
    
    
    # Get SDNR QQ plot and Bubble Plot
    resid_plot_llsl_AI <- SPoRC::plot_resids(osa_results = osa_results_llsl_AI)
    
    # Get Aggregated Plot
    srvdom_lens_agg_AI <- comp_fits$Survey_Lens %>%
      group_by(Region, Len, Sex, Fleet) %>%
      summarize(obs = mean(obs), pred = mean(pred)) %>%
      filter(Fleet == 2) %>%
      ggplot() +
      geom_col(aes(x = Len, y = obs), fill = 'darkgreen', alpha = 0.8) +
      geom_line(aes(x = Len, y = pred), col = 'black', lwd = 1.3) +
      theme_bw(base_size = 15) +
      facet_grid( Region~Sex, labeller = labeller(
        Region = function(x) paste0("Region ", x),
        Sex = function(x) paste0("Sex ", x)
      )) +
      theme_bw(base_size = 20) +
      labs(x = 'Lengths (cm)', y = 'Proportions')
    
    osa_lls_len_AI <- cowplot::plot_grid(resid_plot_llsl_AI[[2]], # Bubble
                                         cowplot::plot_grid(resid_plot_llsl_AI[[1]], srvdom_lens_agg_AI, nrow = 1), # QQ and Aggregated
                                         ncol = 1, rel_heights = c(0.6, 0.4), labels = "AI LLS OSA Length Comps")
    
  }
  
  
  # GOA
  
  if(any(data$UseSrvLenComps[,,3] == 1)){
    if(unique(data$SrvLenComps_Type[,3]) == 2){
      osa_results_llsl_GOA <- SPoRC::get_osa(obs_mat = comp_fits$Obs_SrvLen_mat + 1e-4,
                                             exp_mat = comp_fits$Pred_SrvLen_mat + 1e-4,
                                             N = array(c(na.omit(data$ISS_SrvLenComps[1,which(data$UseSrvLenComps[,,3] == 1),1,3] *
                                                                   unique(data$Wt_SrvLenComps[1,which(data$UseSrvLenComps[,,3] == 1),1,3]))),
                                                       dim = c(data$n_regions, length(which(data$UseSrvLenComps[,,3] == 1)))),
                                             years = which(data$UseSrvLenComps[,,3] == 1),
                                             fleet = 3,
                                             bins = lens,
                                             comp_type = 2,
                                             bin_label = "Lengths (cm)")
    }
    
    
    # Get SDNR QQ plot and Bubble Plot
    resid_plot_llsl_GOA <- SPoRC::plot_resids(osa_results = osa_results_llsl_GOA)
    
    # Get Aggregated Plot
    srvdom_lens_agg_GOA <- comp_fits$Survey_Lens %>%
      group_by(Region, Len, Sex, Fleet) %>%
      summarize(obs = mean(obs), pred = mean(pred)) %>%
      filter(Fleet == 3) %>%
      ggplot() +
      geom_col(aes(x = Len, y = obs), fill = 'darkgreen', alpha = 0.8) +
      geom_line(aes(x = Len, y = pred), col = 'black', lwd = 1.3) +
      theme_bw(base_size = 15) +
      facet_grid( Region~Sex, labeller = labeller(
        Region = function(x) paste0("Region ", x),
        Sex = function(x) paste0("Sex ", x)
      )) +
      theme_bw(base_size = 20) +
      labs(x = 'Lengths (cm)', y = 'Proportions')
    
    osa_lls_len_GOA <- cowplot::plot_grid(resid_plot_llsl_GOA[[2]], # Bubble
                                          cowplot::plot_grid(resid_plot_llsl_GOA[[1]], srvdom_lens_agg_GOA, nrow = 1), # QQ and Aggregated
                                          ncol = 1, rel_heights = c(0.6, 0.4), labels = "GOA LLS OSA Length Comps")
    
  }
  
  #### Domestic Trawl Survey ---------------------------------------------------
  # Get OSA residuals - need to add small constant to allow OSA calculation for multinomial
  
  if(any(data$UseSrvLenComps[,,4] == 1)){
    if(unique(data$SrvLenComps_Type[,4]) == 2){
      osa_results_tsl <- SPoRC::get_osa(obs_mat = comp_fits$Obs_SrvLen_mat + 1e-4,
                                        exp_mat = comp_fits$Pred_SrvLen_mat + 1e-4,
                                        N = array(c(na.omit(data$ISS_SrvLenComps[1,which(data$UseSrvLenComps[,,4] == 1),1,4] *
                                                              unique(data$Wt_SrvLenComps[1,which(data$UseSrvLenComps[,,4] == 1),1,4]))),
                                                  dim = c(data$n_regions, length(which(data$UseSrvLenComps[,,4] == 1)))),
                                        years = which(data$UseSrvLenComps[,,4] == 1),
                                        fleet = 4,
                                        bins = lens,
                                        comp_type = 2,
                                        bin_label = "Lengths (cm)")
    }
    
    # Get SDNR QQ plot and Bubble Plot
    resid_plot_tsl <- SPoRC::plot_resids(osa_results = osa_results_tsl)
    
    # Get Aggregated Plot
    srvtrw_lens_agg <- comp_fits$Survey_Lens %>%
      group_by(Region, Len, Sex, Fleet) %>%
      summarize(obs = mean(obs), pred = mean(pred)) %>%
      filter(Fleet == 4) %>%
      ggplot() +
      geom_col(aes(x = Len, y = obs), fill = 'darkgreen', alpha = 0.8) +
      geom_line(aes(x = Len, y = pred), col = 'black', lwd = 1.3) +
      theme_bw(base_size = 15) +
      facet_grid( Region~Sex, labeller = labeller(
        Region = function(x) paste0("Region ", x),
        Sex = function(x) paste0("Sex ", x)
      )) +
      theme_bw(base_size = 20) +
      labs(x = 'Lengths (cm)', y = 'Proportions')
    
    osa_ts_len <- cowplot::plot_grid(resid_plot_tsl[[2]], # Bubble
                                     cowplot::plot_grid(resid_plot_tsl[[1]], srvtrw_lens_agg, nrow = 1), # QQ and Aggregated
                                     ncol = 1, rel_heights = c(0.6, 0.4), labels = "TS OSA Length Comps")
  }
  
  
  #### Japanese LL Survey ---------------------------------------------------
  # Ages--------------------------------
  
  if(unique(data$SrvAgeComps_Type[,5]) == 2){
    osa_results_jpn_llsa <- SPoRC::get_osa(obs_mat = comp_fits$Obs_SrvAge_mat + 1e-4,
                                           exp_mat = comp_fits$Pred_SrvAge_mat + 1e-4,
                                           N = array(c(na.omit(data$ISS_SrvAgeComps[1,which(data$UseSrvAgeComps[,,5] == 1),1,5] *
                                                                 unique(data$Wt_SrvAgeComps[1,which(data$UseSrvAgeComps[,,5] == 1),1,5]))),
                                                     dim = c(data$n_regions, length(which(data$UseSrvAgeComps[,,5] == 1)))),
                                           years = which(data$UseSrvAgeComps[,,5] == 1),
                                           fleet = 5,
                                           bins = ages,
                                           comp_type = 2,                          # ==0 is aggregated comps, ==1 is split comps, ==2 is joint comps
                                           bin_label = "Ages")
  }
  
  # Get SDNR QQ plot and Bubble Plot
  resid_plot_jpn_llsa <- SPoRC::plot_resids(osa_results = osa_results_jpn_llsa)
  
  # Get Aggregated Plot
  jp_ll_age_agg <- comp_fits$Survey_Ages %>%
    group_by(Region,Age,Sex, Fleet) %>%
    summarize(obs = mean(obs), pred = mean(pred)) %>%
    filter(Fleet == 5) %>%
    ggplot() +
    geom_col(aes(x = Age, y = obs), fill = 'darkgreen', alpha = 0.8) +
    geom_line(aes(x = Age, y = pred), col = 'black', lwd = 1.3) +
    theme_bw(base_size = 15) +
    facet_grid( Region~Sex, labeller = labeller(
      Region = function(x) paste0("Region ", x),
      Sex = function(x) paste0("Sex ", x)
    )) +
    labs(x = 'Age', y = 'Proportions')
  
  osa_jpn_lls_age <- cowplot::plot_grid(resid_plot_jpn_llsa[[2]], # Bubble
                                        cowplot::plot_grid(resid_plot_jpn_llsa[[1]], jp_ll_age_agg, nrow = 1), # QQ and Aggregated
                                        ncol = 1, rel_heights = c(0.6, 0.4), labels = "JPN LLS OSA Age Comps")
  

  # lengths---------------------------------------------------
  
  if(any(data$UseSrvLenComps[,,5] == 1)){
    if(unique(data$SrvLenComps_Type[,5]) == 2){
      osa_results_jpn_llsl <- SPoRC::get_osa(obs_mat = comp_fits$Obs_SrvLen_mat + 1e-4,
                                             exp_mat = comp_fits$Pred_SrvLen_mat + 1e-4,
                                             N = array(c(na.omit(data$ISS_SrvLenComps[1,which(data$UseSrvLenComps[,,5] == 1),1,5] *
                                                                   unique(data$Wt_SrvLenComps[1,which(data$UseSrvLenComps[,,5] == 1),1,5]))),
                                                       dim = c(data$n_regions, length(which(data$UseSrvLenComps[,,5] == 1)))),
                                             years = which(data$UseSrvLenComps[,,5] == 1),
                                             fleet = 5,
                                             bins = lens,
                                             comp_type = 2,
                                             bin_label = "Lengths (cm)")
    }
    
    # Get SDNR QQ plot and Bubble Plot
    resid_plot_jpn_llsl <- SPoRC::plot_resids(osa_results = osa_results_jpn_llsl)
    
    # Get Aggregated Plot
    jp_ll_lens_agg <- comp_fits$Survey_Lens %>%
      group_by(Region, Len, Sex, Fleet) %>%
      summarize(obs = mean(obs), pred = mean(pred)) %>%
      filter(Fleet == 5) %>%
      ggplot() +
      geom_col(aes(x = Len, y = obs), fill = 'darkgreen', alpha = 0.8) +
      geom_line(aes(x = Len, y = pred), col = 'black', lwd = 1.3) +
      theme_bw(base_size = 15) +
      facet_grid( Region~Sex, labeller = labeller(
        Region = function(x) paste0("Region ", x),
        Sex = function(x) paste0("Sex ", x)
      )) +
      theme_bw(base_size = 20) +
      labs(x = 'Lengths (cm)', y = 'Proportions')
    
    osa_jpn_lls_len <- cowplot::plot_grid(resid_plot_jpn_llsl[[2]], # Bubble
                                          cowplot::plot_grid(resid_plot_jpn_llsl[[1]], jp_ll_lens_agg, nrow = 1), # QQ and Aggregated
                                          ncol = 1, rel_heights = c(0.6, 0.4), labels = "JPN LLS OSA Length Comps")
    
  }
  
}




########################################################################################################################################################
# Diagnostic Runs -------------------------------------------------------------------------------------------------------------------------------------
#    Don't need to run for every model
########################################################################################################################################################

if(do_diags == 1){
  
  # get retrospective ------------------------------------------------------------------
  retro <- SPoRC::do_retrospective(n_retro = 7,                                  # number of retro peels to run
                                   data = sabie_rtmb_model$data,                 # rtmb data
                                   parameters = sabie_rtmb_model$parameters,     # rtmb parameters
                                   mapping = sabie_rtmb_model$mapping,           # rtmb mapping
                                   random = NULL,                                # if random effects are used
                                   do_par = TRUE,                                # whether or not to parralleize
                                   n_cores = 7,                                  # if parallel, number of cores to use
                                   do_francis = TRUE,                            # if we want tod o Francis
                                   n_francis_iter = 10                           # Number of francis iterations to do
  )
  
  retro_plot <- get_retrospective_plot(retro, Rec_Age = SPoRC_one_reg$rec_age)   # get retrospective plot
  
  tiff(filename=here(root_folder,mod_name,"plots",paste0(mod_name,"_retro_relative.png")),units = "in", width=14,height=8, res = 300)
  print(retro_plot[[1]]+labs(title = "Retrospective Relative Difference"))
  dev.off()
  tiff(filename=here(root_folder,mod_name,"plots",paste0(mod_name,"_retro_abs.png")),units = "in", width=14,height=8, res = 300)
  print(retro_plot[[2]]+labs(title = "Retrospective Absolute Difference"))
  dev.off()
  tiff(filename=here(root_folder,mod_name,"plots",paste0(mod_name,"_retro_squid.png")),units = "in", width=8,height=8, res = 300)
  print(retro_plot[[3]]+labs(title = "Retrospective Recruitment Squid Plot"))
  dev.off()
  
  
  
  # Jitter Analysis ---------------------------------------------------------
  # get jitter results
  jitter_res <- SPoRC::do_jitter(data = data,
                                 parameters = parameters,
                                 mapping = mapping,
                                 random = NULL,
                                 sd = 0.1,
                                 n_jitter = 50,
                                 n_newton_loops = 3,
                                 do_par = TRUE,
                                 n_cores = 8
  )
  
  # get proportion converged
  prop_converged <- jitter_res %>%
    filter(Year == 1, Type == 'Recruitment') %>%
    summarize(prop_conv = sum(Hessian) / length(Hessian))
  
  # get final model results
  final_mod <- reshape2::melt(rep$SSB) %>% rename(Region = Var1, Year = Var2) %>%
    mutate(Type = 'SSB') %>%
    bind_rows(reshape2::melt(rep$Rec) %>%
                rename(Region = Var1, Year = Var2) %>% mutate(Type = 'Recruitment'))
  
  # comparison of SSB and recruitment
  jitter_ssb <- ggplot() +
    geom_line(jitter_res, mapping = aes(x = Year + 1959, y = value, group = jitter, color = Hessian), lwd = 1) +
    geom_line(final_mod, mapping = aes(x = Year + 1959, y = value), color = "black", lwd = 1.3 , lty = 2) +
    facet_grid(Type~Region, scales = 'free',
               labeller = labeller(Region = function(x) paste0("Region ", x),
                                   Type = c("Recruitment" = "Age 2 Recruitment (millions)", "SSB" = 'SSB (kt)'))) +
    labs(x = "Year", y = "Value") +
    theme_bw(base_size = 20) +
    scale_color_manual(values = c("red", 'grey')) +
    geom_text(data = jitter_res %>% filter(Type == 'SSB', Year == 1, jitter == 1),
              aes(x = Inf, y = Inf, label = paste("Proportion Converged: ", round(prop_converged$prop_conv, 3))),
              hjust = 1.1, vjust = 1.9, size = 6, color = "black")+labs(title = "Jitter SSB and Recruitment")
  
  tiff(filename=here(root_folder,mod_name,"plots",paste0(mod_name,"_jitter_SSB_Recr.png")),units = "in", width=8,height=8, res = 300)
  print(jitter_ssb)
  dev.off()
  
  # compare jitter of max gradient and hessian PD
  jitter_grad <- ggplot(jitter_res, aes(x = jitter, y = jnLL, color = Max_Gradient, shape = Hessian)) +
    geom_point(size = 5, alpha = 0.3) +
    geom_hline(yintercept = min(rep$jnLL), lty = 2, linewidth = 2, color = "blue") +
    facet_wrap(~Hessian, labeller = labeller(
      Hessian = c("FALSE" = "non-PD Hessian", "TRUE" = 'PD Hessian')
    )) +
    scale_color_viridis_c() +
    theme_bw(base_size = 20) +
    theme(legend.position = "bottom") +
    guides(color = guide_colorbar(barwidth = 15, barheight = 0.5)) +
    labs(x = 'Jitter') +
    geom_text(data = jitter_res %>% filter(Hessian == TRUE, Year == 1, jitter == 1),
              aes(x = Inf, y = Inf, label = paste("Proportion Converged: ", round(prop_converged$prop_conv, 3))),
              hjust = 1.1, vjust = 1.9, size = 6, color = "black")+labs(title = "Jitter nLL and Convergence")
  
  tiff(filename=here(root_folder,mod_name,"plots",paste0(mod_name,"_jitter_grad.png")),units = "in", width=8,height=8, res = 300)
  print(jitter_grad)
  dev.off()
  
  
  # Likelihood Profiles ------------------------------------------------------
  
  # extract out mean recruitment value
  mean_rec_val <- sabie_rtmb_model$env$last.par.best[names(sabie_rtmb_model$env$last.par.best) == 'ln_global_R0']
  
  # do mean recruitment profile
  mean_rec_profile <- SPoRC::do_likelihood_profile(data = data,
                                                   parameters = parameters,
                                                   mapping = mapping,
                                                   random = NULL,
                                                   what = "ln_global_R0",
                                                   min_val = 0.15 * mean_rec_val,
                                                   max_val = 1.5 * mean_rec_val,
                                                   inc = 0.05,
                                                   do_par = TRUE,
                                                   n_cores = 8
  )
  
  # Plot mean recruitment profile
  mean_rec_unsmoothed_prof <- mean_rec_profile$agg_nLL %>%
    mutate(Summarized_Type = case_when(
      str_detect(type, "Pen|Prior") ~ "Other",
      str_detect(type, "Len") ~ "Length Comps",
      str_detect(type, "Age") ~ "Age Comps",
      str_detect(type, "Idx") ~ "Indices",
      str_detect(type, "Catch") ~ "Catch",
      str_detect(type, "jnLL") ~ "jnLL",
    )) %>%
    filter(value != 0) %>%
    group_by(Summarized_Type, prof_val) %>%
    summarize(value = sum(value), .groups = "drop") %>%
    group_by(Summarized_Type) %>%
    mutate(value = value - min(value))
  
  # Plot!
  unsmoothed_prof <- ggplot(mean_rec_unsmoothed_prof, aes(x = prof_val, y = value, color = Summarized_Type)) +
    geom_line(lwd = 1.3) +
    geom_vline(xintercept = mean_rec_val, lty = 2) +
    labs(x = 'Log Mean Recruitment', y = "Scaled nLL", color = "Type") +
    theme_bw(base_size = 15)+labs(title = "Profile Likelihood R0")
  
  tiff(filename=here(root_folder,mod_name,"plots",paste0(mod_name,"_likelihood_prof_unsmooth.png")),units = "in", width=8,height=8, res = 300)
  print(unsmoothed_prof)
  dev.off()
  
  ll_prof_unsmooth_plot_type <- ggplot(mean_rec_unsmoothed_prof, aes(x = prof_val, y = value, color = Summarized_Type)) +
    geom_line(lwd = 1.3) +
    facet_wrap(~Summarized_Type, scales = 'free') +
    geom_vline(xintercept = mean_rec_val, lty = 2) +
    labs(x = 'Log Mean Recruitment', y = "Scaled nLL", color = "Type") +
    theme_bw(base_size = 15)+labs(title = "Profile Likelihood R0")
  
  tiff(filename=here(root_folder,mod_name,"plots",paste0(mod_name,"_likelihood_prof_unsmooth_type.png")),units = "in", width=8,height=8, res = 300)
  print(ll_prof_unsmooth_plot_type)
  dev.off()
  
  # Get rescaled smooths for local minima
  mean_rec_prof <- mean_rec_unsmoothed_prof %>%
    group_by(Summarized_Type) %>%                                                # basically, predicting with a smooth, and then rescaling the smooth to have min of 0 - doing this because the smooths sometimes do not predict a min of 0
    group_modify(~ {
      fit <- loess(value ~ prof_val, data = .x, span = 3)                        # Fit smoother on grouped likelihoods
      pred_grid <- seq(min(.x$prof_val), max(.x$prof_val), length.out = 100)     # prediction grdi
      preds <- predict(fit, newdata = pred_grid)                                 # smoothed likelihood profile predictions
      tibble(prof_val = pred_grid, smooth_value = preds)}) %>%                   # return data
    group_by(Summarized_Type) %>%                                                # Now rescale smooth values so minimum is 0 for each group
    mutate(smooth_value = smooth_value - min(smooth_value, na.rm = TRUE))
  
  # Plot!
  ll_prof_plot <- ggplot(mean_rec_prof, aes(x = prof_val, y = smooth_value, color = Summarized_Type)) +
    geom_line(lwd = 1.3) +
    geom_vline(xintercept = mean_rec_val, lty = 2) +
    labs(x = 'Log Mean Recruitment', y = "Scaled nLL", color = "Type") +
    theme_bw(base_size = 15)+labs(title = "Profile Likelihood R0")
  
  tiff(filename=here(root_folder,mod_name,"plots",paste0(mod_name,"_likelihood_prof.png")),units = "in", width=8,height=8, res = 300)
  print(ll_prof_plot)
  dev.off()
  
  ll_prof_plot_type <- ggplot(mean_rec_prof, aes(x = prof_val, y = smooth_value, color = Summarized_Type)) +
    geom_line(lwd = 1.3) +
    facet_wrap(~Summarized_Type, scales = 'free') +
    geom_vline(xintercept = mean_rec_val, lty = 2) +
    labs(x = 'Log Mean Recruitment', y = "Scaled nLL", color = "Type") +
    theme_bw(base_size = 15)+labs(title = "Profile Likelihood R0")
  
  tiff(filename=here(root_folder,mod_name,"plots",paste0(mod_name,"_likelihood_prof_type.png")),units = "in", width=8,height=8, res = 300)
  print(ll_prof_plot_type)
  dev.off()
  
  
  # extract out M value
  M_val_fem <- sabie_rtmb_model$env$last.par.best[names(sabie_rtmb_model$env$last.par.best) == 'ln_M'][1]
  M_val_mal <- sabie_rtmb_model$env$last.par.best[names(sabie_rtmb_model$env$last.par.best) == 'ln_M'][2]
  
  # do mean recruitment profile
  M_profile_fem <- SPoRC::do_likelihood_profile(data = data,
                                                parameters = parameters,
                                                mapping = mapping,
                                                random = NULL,
                                                idx = 1,
                                                what = "ln_M",
                                                min_val = 1.55 * M_val_fem,
                                                max_val = 0.8 * M_val_fem,
                                                inc = 0.05,
                                                do_par = TRUE,
                                                n_cores = 8
  )
  
  M_profile_mal <- SPoRC::do_likelihood_profile(data = data,
                                                parameters = parameters,
                                                mapping = mapping,
                                                random = NULL,
                                                idx = 2,
                                                what = "ln_M",
                                                min_val = 1.55 * M_val_mal,
                                                max_val = 0.8 * M_val_mal,
                                                inc = 0.05,
                                                do_par = TRUE,
                                                n_cores = 8
  )
  
  # Plot mean recruitment profile
  M_unsmoothed_prof_fem <- M_profile_fem$agg_nLL %>%
    mutate(Summarized_Type = case_when(
      str_detect(type, "Pen|Prior") ~ "Other",
      str_detect(type, "Len") ~ "Length Comps",
      str_detect(type, "Age") ~ "Age Comps",
      str_detect(type, "Idx") ~ "Indices",
      str_detect(type, "Catch") ~ "Catch",
      str_detect(type, "jnLL") ~ "jnLL",
    )) %>%
    filter(value != 0) %>%
    group_by(Summarized_Type, prof_val) %>%
    summarize(value = sum(value), .groups = "drop") %>%
    group_by(Summarized_Type) %>%
    mutate(value = value - min(value))
  
  # Plot!
  unsmoothed_prof_M_fem <- ggplot(M_unsmoothed_prof_fem, aes(x = prof_val, y = value, color = Summarized_Type)) +
    geom_line(lwd = 1.3) +
    geom_vline(xintercept = M_val_fem, lty = 2) +
    labs(x = 'Log_M', y = "Scaled nLL", color = "Type") +
    theme_bw(base_size = 15)+labs(title = "Profile Likelihood Female M")
  
  tiff(filename=here(root_folder,mod_name,"plots",paste0(mod_name,"_likelihood_prof_unsmooth_M_fem.png")),units = "in", width=8,height=8, res = 300)
  print(unsmoothed_prof_M_fem)
  dev.off()
  
  ll_prof_unsmooth_plot_type_M_fem <- ggplot(M_unsmoothed_prof_fem, aes(x = prof_val, y = value, color = Summarized_Type)) +
    geom_line(lwd = 1.3) +
    facet_wrap(~Summarized_Type, scales = 'free') +
    geom_vline(xintercept = M_val_fem, lty = 2) +
    labs(x = 'Log_M', y = "Scaled nLL", color = "Type") +
    theme_bw(base_size = 15)+labs(title = "Profile Likelihood Female M")
  
  tiff(filename=here(root_folder,mod_name,"plots",paste0(mod_name,"_likelihood_prof_unsmooth_type_M_Fem.png")),units = "in", width=8,height=8, res = 300)
  print(ll_prof_unsmooth_plot_type_M_fem)
  dev.off()
  
  # Get rescaled smooths for local minima
  M_prof_fem <- M_unsmoothed_prof_fem %>%
    group_by(Summarized_Type) %>%                                                # basically, predicting with a smooth, and then rescaling the smooth to have min of 0 - doing this because the smooths sometimes do not predict a min of 0
    group_modify(~ {
      fit <- loess(value ~ prof_val, data = .x, span = 3)                        # Fit smoother on grouped likelihoods
      pred_grid <- seq(min(.x$prof_val), max(.x$prof_val), length.out = 100)     # prediction grdi
      preds <- predict(fit, newdata = pred_grid)                                 # smoothed likelihood profile predictions
      tibble(prof_val = pred_grid, smooth_value = preds)}) %>%                   # return data
    group_by(Summarized_Type) %>%                                                # Now rescale smooth values so minimum is 0 for each group
    mutate(smooth_value = smooth_value - min(smooth_value, na.rm = TRUE))
  
  # Plot!
  ll_prof_plot_M_fem <- ggplot(M_prof_fem, aes(x = prof_val, y = smooth_value, color = Summarized_Type)) +
    geom_line(lwd = 1.3) +
    geom_vline(xintercept = M_val_fem, lty = 2) +
    labs(x = 'Log_M', y = "Scaled nLL", color = "Type") +
    theme_bw(base_size = 15)+labs(title = "Profile Likelihood Female M")
  
  tiff(filename=here(root_folder,mod_name,"plots",paste0(mod_name,"_likelihood_prof_M_fem.png")),units = "in", width=8,height=8, res = 300)
  print(ll_prof_plot_M_fem)
  dev.off()
  
  ll_prof_plot_type_M_fem <- ggplot(M_prof_fem, aes(x = prof_val, y = smooth_value, color = Summarized_Type)) +
    geom_line(lwd = 1.3) +
    facet_wrap(~Summarized_Type, scales = 'free') +
    geom_vline(xintercept = M_val_fem, lty = 2) +
    labs(x = 'Log_M', y = "Scaled nLL", color = "Type") +
    theme_bw(base_size = 15)+labs(title = "Profile Likelihood Female M")
  
  tiff(filename=here(root_folder,mod_name,"plots",paste0(mod_name,"_likelihood_prof_type_M_fem.png")),units = "in", width=8,height=8, res = 300)
  print(ll_prof_plot_type_M_fem)
  dev.off()
  
  
  # Plot mean recruitment profile
  M_unsmoothed_prof_mal <- M_profile_mal$agg_nLL %>%
    mutate(Summarized_Type = case_when(
      str_detect(type, "Pen|Prior") ~ "Other",
      str_detect(type, "Len") ~ "Length Comps",
      str_detect(type, "Age") ~ "Age Comps",
      str_detect(type, "Idx") ~ "Indices",
      str_detect(type, "Catch") ~ "Catch",
      str_detect(type, "jnLL") ~ "jnLL",
    )) %>%
    filter(value != 0) %>%
    group_by(Summarized_Type, prof_val) %>%
    summarize(value = sum(value), .groups = "drop") %>%
    group_by(Summarized_Type) %>%
    mutate(value = value - min(value))
  
  # Plot!
  unsmoothed_prof_M_mal <- ggplot(M_unsmoothed_prof_mal, aes(x = prof_val, y = value, color = Summarized_Type)) +
    geom_line(lwd = 1.3) +
    geom_vline(xintercept = M_val_mal, lty = 2) +
    labs(x = 'Log_M', y = "Scaled nLL", color = "Type") +
    theme_bw(base_size = 15)+labs(title = "Profile Likelihood Male M")
  
  tiff(filename=here(root_folder,mod_name,"plots",paste0(mod_name,"_likelihood_prof_unsmooth_M_mal.png")),units = "in", width=8,height=8, res = 300)
  print(unsmoothed_prof_M_mal)
  dev.off()
  
  ll_prof_unsmooth_plot_type_M_mal <- ggplot(M_unsmoothed_prof_mal, aes(x = prof_val, y = value, color = Summarized_Type)) +
    geom_line(lwd = 1.3) +
    facet_wrap(~Summarized_Type, scales = 'free') +
    geom_vline(xintercept = M_val_mal, lty = 2) +
    labs(x = 'Log_M', y = "Scaled nLL", color = "Type") +
    theme_bw(base_size = 15)+labs(title = "Profile Likelihood Male M")
  
  tiff(filename=here(root_folder,mod_name,"plots",paste0(mod_name,"_likelihood_prof_unsmooth_type_M_mal.png")),units = "in", width=8,height=8, res = 300)
  print(ll_prof_unsmooth_plot_type_M_mal)
  dev.off()
  
  # Get rescaled smooths for local minima
  M_prof_mal <- M_unsmoothed_prof_mal %>%
    group_by(Summarized_Type) %>%                                                # basically, predicting with a smooth, and then rescaling the smooth to have min of 0 - doing this because the smooths sometimes do not predict a min of 0
    group_modify(~ {
      fit <- loess(value ~ prof_val, data = .x, span = 3)                        # Fit smoother on grouped likelihoods
      pred_grid <- seq(min(.x$prof_val), max(.x$prof_val), length.out = 100)     # prediction grdi
      preds <- predict(fit, newdata = pred_grid)                                 # smoothed likelihood profile predictions
      tibble(prof_val = pred_grid, smooth_value = preds)}) %>%                   # return data
    group_by(Summarized_Type) %>%                                                # Now rescale smooth values so minimum is 0 for each group
    mutate(smooth_value = smooth_value - min(smooth_value, na.rm = TRUE))
  
  # Plot!
  ll_prof_plot_M_mal <- ggplot(M_prof_mal, aes(x = prof_val, y = smooth_value, color = Summarized_Type)) +
    geom_line(lwd = 1.3) +
    geom_vline(xintercept = M_val_mal, lty = 2) +
    labs(x = 'Log_M', y = "Scaled nLL", color = "Type") +
    theme_bw(base_size = 15)+labs(title = "Profile Likelihood Male M")
  
  tiff(filename=here(root_folder,mod_name,"plots",paste0(mod_name,"_likelihood_prof_M_mal.png")),units = "in", width=8,height=8, res = 300)
  print(ll_prof_plot_M_mal)
  dev.off()
  
  ll_prof_plot_type_M_mal <- ggplot(M_prof_mal, aes(x = prof_val, y = smooth_value, color = Summarized_Type)) +
    geom_line(lwd = 1.3) +
    facet_wrap(~Summarized_Type, scales = 'free') +
    geom_vline(xintercept = M_val_mal, lty = 2) +
    labs(x = 'Log_M', y = "Scaled nLL", color = "Type") +
    theme_bw(base_size = 15)+labs(title = "Profile Likelihood Male M")
  
  tiff(filename=here(root_folder,mod_name,"plots",paste0(mod_name,"_likelihood_prof_type_M_mal.png")),units = "in", width=8,height=8, res = 300)
  print(ll_prof_plot_type_M_mal)
  dev.off()
}  # do_diags loop




if(do_high_res_comps == 1){
  
  tiff(filename=here(root_folder,mod_name,"plots",paste0(mod_name,"_llf_age_fem.png")),units = "in", width=12,height=12, res = 300)
  print(llf_age_fem)
  dev.off()
  
  tiff(filename=here(root_folder,mod_name,"plots",paste0(mod_name,"_llf_age_cohort_fem.png")),units = "in", width=12,height=12, res = 300)
  print(llf_age_coh_fem)
  dev.off()
  
  # Male
  # Fits by Year
  if(unique(data$FishAgeComps_Type[,1]) > 0){
    
    tiff(filename=here(root_folder,mod_name,"plots",paste0(mod_name,"_llf_age_male.png")),units = "in", width=12,height=12, res = 300)
    print(llf_age_mal)
    dev.off()
    
    tiff(filename=here(root_folder,mod_name,"plots",paste0(mod_name,"_llf_age_cohort_male.png")),units = "in", width=12,height=12, res = 300)
    print(llf_age_coh_mal)
    dev.off()
  }
  
  if(any(data$UseFishLenComps[,,1] == 1)){
    
    tiff(filename=here(root_folder,mod_name,"plots",paste0(mod_name,"_llf_len_fem.png")),units = "in", width=12,height=12, res = 300)
    print(llf_len_fem)
    dev.off()
    
    tiff(filename=here(root_folder,mod_name,"plots",paste0(mod_name,"_llf_len_male.png")),units = "in", width=12,height=12, res = 300)
    print(llf_len_male)
    dev.off()
  }
  
  tiff(filename=here(root_folder,mod_name,"plots",paste0(mod_name,"_tf_len_fem.png")),units = "in", width=12,height=12, res = 300)
  print(tf_len_fem)
  dev.off()
  
  tiff(filename=here(root_folder,mod_name,"plots",paste0(mod_name,"_tf_len_male.png")),units = "in", width=12,height=12, res = 300)
  print(tf_len_male)
  dev.off()
  
  
  
  tiff(filename=here(root_folder,mod_name,"plots",paste0(mod_name,"_BS_lls_age_fem.png")),units = "in", width=12,height=12, res = 300)
  print(BS_lls_age_fem)
  dev.off()
  
  tiff(filename=here(root_folder,mod_name,"plots",paste0(mod_name,"_BS_lls_age_cohort_fem.png")),units = "in", width=12,height=12, res = 300)
  print(BS_lls_age_coh_fem)
  dev.off()
  
  # Male
  if(unique(data$SrvAgeComps_Type[,1]) > 0){
    
    tiff(filename=here(root_folder,mod_name,"plots",paste0(mod_name,"_BS_lls_age_mal.png")),units = "in", width=12,height=12, res = 300)
    print(BS_lls_age_mal)
    dev.off()
    
    tiff(filename=here(root_folder,mod_name,"plots",paste0(mod_name,"_BS_lls_age_cohort_mal.png")),units = "in", width=12,height=12, res = 300)
    print(BS_lls_age_coh_mal)
    dev.off()
  }
  
  tiff(filename=here(root_folder,mod_name,"plots",paste0(mod_name,"_AI_lls_age_fem.png")),units = "in", width=12,height=12, res = 300)
  print(AI_lls_age_fem)
  dev.off()
  
  tiff(filename=here(root_folder,mod_name,"plots",paste0(mod_name,"_AI_lls_age_cohort_fem.png")),units = "in", width=12,height=12, res = 300)
  print(AI_lls_age_coh_fem)
  dev.off()
  
  # Male
  if(unique(data$SrvAgeComps_Type[,2]) > 0){
    
    tiff(filename=here(root_folder,mod_name,"plots",paste0(mod_name,"_AI_lls_age_mal.png")),units = "in", width=12,height=12, res = 300)
    print(AI_lls_age_mal)
    dev.off()
    
    tiff(filename=here(root_folder,mod_name,"plots",paste0(mod_name,"_AI_lls_age_cohort_mal.png")),units = "in", width=12,height=12, res = 300)
    print(AI_lls_age_coh_mal)
    dev.off()
  }
  
  tiff(filename=here(root_folder,mod_name,"plots",paste0(mod_name,"_GOA_lls_age_fem.png")),units = "in", width=12,height=12, res = 300)
  print(GOA_lls_age_fem)
  dev.off()
  
  tiff(filename=here(root_folder,mod_name,"plots",paste0(mod_name,"_GOA_lls_age_cohort_fem.png")),units = "in", width=12,height=12, res = 300)
  print(GOA_lls_age_coh_fem)
  dev.off()
  
  # Male
  if(unique(data$SrvAgeComps_Type[,3]) > 0){
    
    tiff(filename=here(root_folder,mod_name,"plots",paste0(mod_name,"_GOA_lls_age_mal.png")),units = "in", width=12,height=12, res = 300)
    print(GOA_lls_age_mal)
    dev.off()
    
    tiff(filename=here(root_folder,mod_name,"plots",paste0(mod_name,"_GOA_lls_age_cohort_mal.png")),units = "in", width=12,height=12, res = 300)
    print(GOA_lls_age_coh_mal)
    dev.off()
  }
  
  tiff(filename=here(root_folder,mod_name,"plots",paste0(mod_name,"_jpn_lls_age_fem.png")),units = "in", width=12,height=12, res = 300)
  print(jpn_lls_age_fem)
  dev.off()
  
  tiff(filename=here(root_folder,mod_name,"plots",paste0(mod_name,"_jpn_lls_age_cohort_fem.png")),units = "in", width=12,height=12, res = 300)
  print(jpn_lls_age_coh_fem)
  dev.off()
  
  # Male
  if(unique(data$SrvAgeComps_Type[,5]) > 0){
    
    tiff(filename=here(root_folder,mod_name,"plots",paste0(mod_name,"_jpn_lls_age_mal.png")),units = "in", width=12,height=12, res = 300)
    print(jpn_lls_age_mal)
    dev.off()
    
    tiff(filename=here(root_folder,mod_name,"plots",paste0(mod_name,"_jpn_lls_age_cohort_mal.png")),units = "in", width=12,height=12, res = 300)
    print(jpn_lls_age_coh_mal)
    dev.off()
  }
  
  if(any(data$UseSrvLenComps[,,1] == 1)){
    tiff(filename=here(root_folder,mod_name,"plots",paste0(mod_name,"_BS_lls_len_fem.png")),units = "in", width=12,height=12, res = 300)
    print(BS_lls_len_fem)
    dev.off()
    
    
    tiff(filename=here(root_folder,mod_name,"plots",paste0(mod_name,"_BS_lls_len_male.png")),units = "in", width=12,height=12, res = 300)
    print(BS_lls_len_male)
    dev.off()
  }
  
  if(any(data$UseSrvLenComps[,,2] == 1)){
    tiff(filename=here(root_folder,mod_name,"plots",paste0(mod_name,"_AI_lls_len_fem.png")),units = "in", width=12,height=12, res = 300)
    print(AI_lls_len_fem)
    dev.off()
    
    
    tiff(filename=here(root_folder,mod_name,"plots",paste0(mod_name,"_AI_lls_len_male.png")),units = "in", width=12,height=12, res = 300)
    print(AI_lls_len_male)
    dev.off()
  }
  
  if(any(data$UseSrvLenComps[,,3] == 1)){
    tiff(filename=here(root_folder,mod_name,"plots",paste0(mod_name,"_GOA_lls_len_fem.png")),units = "in", width=12,height=12, res = 300)
    print(GOA_lls_len_fem)
    dev.off()
    
    
    tiff(filename=here(root_folder,mod_name,"plots",paste0(mod_name,"_GOA_lls_len_male.png")),units = "in", width=12,height=12, res = 300)
    print(GOA_lls_len_male)
    dev.off()
  }
  
  if(any(data$UseSrvLenComps[,,4] == 1)){
    
    tiff(filename=here(root_folder,mod_name,"plots",paste0(mod_name,"_ts_len_fem.png")),units = "in", width=12,height=12, res = 300)
    print(ts_len_fem)
    dev.off()
    
    tiff(filename=here(root_folder,mod_name,"plots",paste0(mod_name,"_ts_len_male.png")),units = "in", width=12,height=12, res = 300)
    print(ts_len_male)
    dev.off()
  }
  
  if(any(data$UseSrvLenComps[,,5] == 1)){
    tiff(filename=here(root_folder,mod_name,"plots",paste0(mod_name,"_jpn_lls_len_fem.png")),units = "in", width=12,height=12, res = 300)
    print(jpn_lls_len_fem)
    dev.off()
    
    tiff(filename=here(root_folder,mod_name,"plots",paste0(mod_name,"_jpn_lls_len_male.png")),units = "in", width=12,height=12, res = 300)
    print(jpn_lls_len_male)
    dev.off()
  }
  
  
  # OSA COMPS
  if(do_OSA == 1){
    tiff(filename=here(root_folder,mod_name,"plots",paste0(mod_name,"_osa_llf_age.png")),units = "in", width=12,height=12, res = 300)
    print(osa_llf_age)
    dev.off()
    
    if(any(data$UseFishLenComps[,,1] == 1)){
      
      tiff(filename=here(root_folder,mod_name,"plots",paste0(mod_name,"_osa_llf_len.png")),units = "in", width=12,height=12, res = 300)
      print(osa_llf_len)
      dev.off()
    }
    
    tiff(filename=here(root_folder,mod_name,"plots",paste0(mod_name,"_osa_tf_len.png")),units = "in", width=12,height=12, res = 300)
    print(osa_tf)
    dev.off()
    
    
    
    tiff(filename=here(root_folder,mod_name,"plots",paste0(mod_name,"_osa_lls_age_BS.png")),units = "in", width=12,height=12, res = 300)
    print(osa_lls_age_BS)
    dev.off()
    
    tiff(filename=here(root_folder,mod_name,"plots",paste0(mod_name,"_osa_lls_age_AI.png")),units = "in", width=12,height=12, res = 300)
    print(osa_lls_age_AI)
    dev.off()
    
    tiff(filename=here(root_folder,mod_name,"plots",paste0(mod_name,"_osa_lls_age_GOA.png")),units = "in", width=12,height=12, res = 300)
    print(osa_lls_age_GOA)
    dev.off()
    
    if(any(data$UseSrvLenComps[,,1] == 1)){
      tiff(filename=here(root_folder,mod_name,"plots",paste0(mod_name,"_osa_lls_len_BS.png")),units = "in", width=12,height=12, res = 300)
      print(osa_lls_len_BS)
      dev.off()
    }
    
    if(any(data$UseSrvLenComps[,,2] == 1)){
      tiff(filename=here(root_folder,mod_name,"plots",paste0(mod_name,"_osa_lls_len_AI.png")),units = "in", width=12,height=12, res = 300)
      print(osa_lls_len_AI)
      dev.off()  
    }
    
    if(any(data$UseSrvLenComps[,,3] == 1)){
      tiff(filename=here(root_folder,mod_name,"plots",paste0(mod_name,"_osa_lls_len_GOA.png")),units = "in", width=12,height=12, res = 300)
      print(osa_lls_len_GOA)
      dev.off()
    }
    
    if(any(data$UseSrvLenComps[,,4] == 1)){
      tiff(filename=here(root_folder,mod_name,"plots",paste0(mod_name,"_osa_ts_len.png")),units = "in", width=12,height=12, res = 300)
      print(osa_ts_len)
      dev.off()
    }
    
    tiff(filename=here(root_folder,mod_name,"plots",paste0(mod_name,"_osa_jpn_lls_age.png")),units = "in", width=12,height=12, res = 300)
    print(osa_jpn_lls_age)
    dev.off()
    
    if(any(data$UseSrvLenComps[,,5] == 1)){
      tiff(filename=here(root_folder,mod_name,"plots",paste0(mod_name,"_osa_jpn_lls_len.png")),units = "in", width=12,height=12, res = 300)
      print(osa_jpn_lls_len)
      dev.off()  
    }
    
  }
  
  
}


# Pearson Resids
if(do_pearson == 1){

source(here("R","pearson_residuals.r"))
  
obs_mat = comp_fits$Obs_FishAge_mat + 1e-4
exp_mat = comp_fits$Pred_FishAge_mat + 1e-4
N = array(c(na.omit(data$ISS_FishAgeComps[1,which(data$UseFishAgeComps[,,1] == 1),1,1] * 
                      unique(data$Wt_FishAgeComps[1,which(data$UseFishAgeComps[,,1] == 1),1,1]))),
          dim = c(data$n_regions, length(data$years)))
years = which(data$UseFishAgeComps[,,1] == 1)
fleet = 1
bins = ages
comp_type = 2
bin_label = "Ages"
data_source = "Fixed Gear Fishery Age Compositions"

pearson_stuff <- get_pearson(
  obs_mat,
  exp_mat,
  N,
  years,
  fleet,
  bins,
  comp_type,
  bin_label,
  data_source
)

llf_age_pear <- pearson_stuff[[2]]

tiff(filename=here(root_folder,mod_name,"plots",paste0(mod_name,"_Pearson_Resids_LLF_Age.png")),units = "in", width=12,height=8, res = 300)
print(llf_age_pear)
dev.off()

if(any(data$UseFishLenComps[,,1] == 1)){
  
  obs_mat = comp_fits$Obs_FishLen_mat + 1e-4
  exp_mat = comp_fits$Pred_FishLen_mat + 1e-4
  N = array(c(na.omit(data$ISS_FishLenComps[1,which(data$UseFishLenComps[,,1] == 1),1,1] * 
                        unique(data$Wt_FishLenComps[1,which(data$UseFishLenComps[,,1] == 1),1,1]))),
            dim = c(data$n_regions, length(data$years)))
  years = which(data$UseFishLenComps[,,1] == 1)
  fleet = 1
  bins = lens
  comp_type = 2
  bin_label = "Lengths"
  data_source = "Fixed Gear Fishery Length Compositions"
  
  pearson_stuff <- get_pearson(
    obs_mat,
    exp_mat,
    N,
    years,
    fleet,
    bins,
    comp_type,
    bin_label,
    data_source
  )
  
  llf_len_pear <- pearson_stuff[[2]]
  
  tiff(filename=here(root_folder,mod_name,"plots",paste0(mod_name,"_Pearson_Resids_LLF_Length.png")),units = "in", width=12,height=8, res = 300)
  print(llf_len_pear)
  dev.off()
}


obs_mat = comp_fits$Obs_FishLen_mat + 1e-4
exp_mat = comp_fits$Pred_FishLen_mat + 1e-4
N = array(c(na.omit(data$ISS_FishLenComps[1,which(data$UseFishLenComps[,,2] == 1),1,2] * 
                      unique(data$Wt_FishLenComps[1,which(data$UseFishLenComps[,,2] == 1),1,2]))),
          dim = c(data$n_regions, length(data$years)))
years = which(data$UseFishLenComps[,,2] == 1)
fleet = 2
bins = lens
comp_type = 2
bin_label = "Lengths"
data_source = "Trawl Fishery Length Compositions"

pearson_stuff <- get_pearson(
  obs_mat,
  exp_mat,
  N,
  years,
  fleet,
  bins,
  comp_type,
  bin_label,
  data_source
)

tf_len_pear <- pearson_stuff[[2]]

tiff(filename=here(root_folder,mod_name,"plots",paste0(mod_name,"_Pearson_Resids_TF_Length.png")),units = "in", width=12,height=8, res = 300)
print(tf_len_pear)
dev.off()


obs_mat = comp_fits$Obs_SrvAge_mat + 1e-4
exp_mat = comp_fits$Pred_SrvAge_mat + 1e-4
N = array(c(na.omit(data$ISS_SrvAgeComps[1,which(data$UseSrvAgeComps[,,1] == 1),1,1] * 
                      unique(data$Wt_SrvAgeComps[1,which(data$UseSrvAgeComps[,,1] == 1),1,1]))),
          dim = c(data$n_regions, length(data$years)))
years = which(data$UseSrvAgeComps[,,1] == 1)
fleet = 1
bins = ages
comp_type = 2
bin_label = "Ages"
data_source = "BS LLS Age Compositions"

pearson_stuff <- get_pearson(
  obs_mat,
  exp_mat,
  N,
  years,
  fleet,
  bins,
  comp_type,
  bin_label,
  data_source
)

bs_lls_age_pear <- pearson_stuff[[2]]

tiff(filename=here(root_folder,mod_name,"plots",paste0(mod_name,"_Pearson_Resids_BS_LLS_Age.png")),units = "in", width=12,height=8, res = 300)
print(bs_lls_age_pear)
dev.off()

obs_mat = comp_fits$Obs_SrvAge_mat + 1e-4
exp_mat = comp_fits$Pred_SrvAge_mat + 1e-4
N = array(c(na.omit(data$ISS_SrvAgeComps[1,which(data$UseSrvAgeComps[,,2] == 1),1,2] * 
                      unique(data$Wt_SrvAgeComps[1,which(data$UseSrvAgeComps[,,2] == 1),1,2]))),
          dim = c(data$n_regions, length(data$years)))
years = which(data$UseSrvAgeComps[,,2] == 1)
fleet = 2
bins = ages
comp_type = 2
bin_label = "Ages"
data_source = "AI LLS Age Compositions"

pearson_stuff <- get_pearson(
  obs_mat,
  exp_mat,
  N,
  years,
  fleet,
  bins,
  comp_type,
  bin_label,
  data_source
)

ai_lls_age_pear <- pearson_stuff[[2]]

tiff(filename=here(root_folder,mod_name,"plots",paste0(mod_name,"_Pearson_Resids_AI_LLS_Age.png")),units = "in", width=12,height=8, res = 300)
print(ai_lls_age_pear)
dev.off()

obs_mat = comp_fits$Obs_SrvAge_mat + 1e-4
exp_mat = comp_fits$Pred_SrvAge_mat + 1e-4
N = array(c(na.omit(data$ISS_SrvAgeComps[1,which(data$UseSrvAgeComps[,,3] == 1),1,3] * 
                      unique(data$Wt_SrvAgeComps[1,which(data$UseSrvAgeComps[,,3] == 1),1,3]))),
          dim = c(data$n_regions, length(data$years)))
years = which(data$UseSrvAgeComps[,,3] == 1)
fleet = 3
bins = ages
comp_type = 2
bin_label = "Ages"
data_source = "GOA LLS Age Compositions"

pearson_stuff <- get_pearson(
  obs_mat,
  exp_mat,
  N,
  years,
  fleet,
  bins,
  comp_type,
  bin_label,
  data_source
)

goa_lls_age_pear <- pearson_stuff[[2]]

tiff(filename=here(root_folder,mod_name,"plots",paste0(mod_name,"_Pearson_Resids_GOA_LLS_Age.png")),units = "in", width=12,height=8, res = 300)
print(goa_lls_age_pear)
dev.off()

if(any(data$UseSrvLenComps[,,1] == 1)){
  obs_mat = comp_fits$Obs_SrvLen_mat + 1e-4
  exp_mat = comp_fits$Pred_SrvLen_mat + 1e-4
  N = array(c(na.omit(data$ISS_SrvLenComps[1,which(data$UseSrvLenComps[,,1] == 1),1,1] * 
                        unique(data$Wt_SrvLenComps[1,which(data$UseSrvLenComps[,,1] == 1),1,1]))),
            dim = c(data$n_regions, length(data$years)))
  years = which(data$UseSrvLenComps[,,1] == 1)
  fleet = 1
  bins = lens
  comp_type = 2
  bin_label = "Lengths"
  data_source = "BS LLS Length Compositions"
  
  pearson_stuff <- get_pearson(
    obs_mat,
    exp_mat,
    N,
    years,
    fleet,
    bins,
    comp_type,
    bin_label,
    data_source
  )
  
  bs_lls_len_pear <- pearson_stuff[[2]]
  
  tiff(filename=here(root_folder,mod_name,"plots",paste0(mod_name,"_Pearson_Resids_BS_LLS_Length.png")),units = "in", width=12,height=8, res = 300)
  print(bs_lls_len_pear)
  dev.off()
}

if(any(data$UseSrvLenComps[,,2] == 1)){
  obs_mat = comp_fits$Obs_SrvLen_mat + 1e-4
  exp_mat = comp_fits$Pred_SrvLen_mat + 1e-4
  N = array(c(na.omit(data$ISS_SrvLenComps[1,which(data$UseSrvLenComps[,,2] == 1),1,2] * 
                        unique(data$Wt_SrvLenComps[1,which(data$UseSrvLenComps[,,2] == 1),1,2]))),
            dim = c(data$n_regions, length(data$years)))
  years = which(data$UseSrvLenComps[,,2] == 1)
  fleet = 2
  bins = lens
  comp_type = 2
  bin_label = "Lengths"
  data_source = "AI LLS Length Compositions"
  
  pearson_stuff <- get_pearson(
    obs_mat,
    exp_mat,
    N,
    years,
    fleet,
    bins,
    comp_type,
    bin_label,
    data_source
  )
  
  ai_lls_len_pear <- pearson_stuff[[2]]
  
  tiff(filename=here(root_folder,mod_name,"plots",paste0(mod_name,"_Pearson_Resids_AI_LLS_Length.png")),units = "in", width=12,height=8, res = 300)
  print(ai_lls_len_pear)
  dev.off()
}

if(any(data$UseSrvLenComps[,,3] == 1)){
  obs_mat = comp_fits$Obs_SrvLen_mat + 1e-4
  exp_mat = comp_fits$Pred_SrvLen_mat + 1e-4
  N = array(c(na.omit(data$ISS_SrvLenComps[1,which(data$UseSrvLenComps[,,3] == 1),1,3] * 
                        unique(data$Wt_SrvLenComps[1,which(data$UseSrvLenComps[,,3] == 1),1,3]))),
            dim = c(data$n_regions, length(data$years)))
  years = which(data$UseSrvLenComps[,,3] == 1)
  fleet = 3
  bins = lens
  comp_type = 2
  bin_label = "Lengths"
  data_source = "GOA LLS Length Compositions"
  
  pearson_stuff <- get_pearson(
    obs_mat,
    exp_mat,
    N,
    years,
    fleet,
    bins,
    comp_type,
    bin_label,
    data_source
  )
  
  goa_lls_len_pear <- pearson_stuff[[2]]
  
  tiff(filename=here(root_folder,mod_name,"plots",paste0(mod_name,"_Pearson_Resids_GOA_LLS_Length.png")),units = "in", width=12,height=8, res = 300)
  print(goa_lls_len_pear)
  dev.off()
}

obs_mat = comp_fits$Obs_SrvAge_mat + 1e-4
exp_mat = comp_fits$Pred_SrvAge_mat + 1e-4
N = array(c(na.omit(data$ISS_SrvAgeComps[1,which(data$UseSrvAgeComps[,,5] == 1),1,5] * 
                      unique(data$Wt_SrvAgeComps[1,which(data$UseSrvAgeComps[,,5] == 1),1,5]))),
          dim = c(data$n_regions, length(data$years)))
years = which(data$UseSrvAgeComps[,,5] == 1)
fleet = 5
bins = ages
comp_type = 2
bin_label = "Ages"
data_source = "JPN LLS Age Compositions"

pearson_stuff <- get_pearson(
  obs_mat,
  exp_mat,
  N,
  years,
  fleet,
  bins,
  comp_type,
  bin_label,
  data_source
)

jpn_lls_age_pear <- pearson_stuff[[2]]

tiff(filename=here(root_folder,mod_name,"plots",paste0(mod_name,"_Pearson_Resids_JPN_LLS_Age.png")),units = "in", width=12,height=8, res = 300)
print(jpn_lls_age_pear)
dev.off()

if(any(data$UseSrvLenComps[,,5] == 1)){
  obs_mat = comp_fits$Obs_SrvLen_mat + 1e-4
  exp_mat = comp_fits$Pred_SrvLen_mat + 1e-4
  N = array(c(na.omit(data$ISS_SrvLenComps[1,which(data$UseSrvLenComps[,,5] == 1),1,5] * 
                        unique(data$Wt_SrvLenComps[1,which(data$UseSrvLenComps[,,5] == 1),1,5]))),
            dim = c(data$n_regions, length(data$years)))
  years = which(data$UseSrvLenComps[,,5] == 1)
  fleet = 5
  bins = lens
  comp_type = 2
  bin_label = "Lengths"
  data_source = "JPN LLS Length Compositions"
  
  pearson_stuff <- get_pearson(
    obs_mat,
    exp_mat,
    N,
    years,
    fleet,
    bins,
    comp_type,
    bin_label,
    data_source
  )
  
  jpn_lls_len_pear <- pearson_stuff[[2]]
  
  tiff(filename=here(root_folder,mod_name,"plots",paste0(mod_name,"_Pearson_Resids_JPN_LLS_Length.png")),units = "in", width=12,height=8, res = 300)
  print(jpn_lls_len_pear)
  dev.off()
}

}



# ISS, ESS, SE ---------------------------------------------------------------

ISS_LLF_age = c(na.omit(data$ISS_FishAgeComps[1,which(data$UseFishAgeComps[,,1] == 1),1,1]))

wt_LLF_age = c(na.omit(unique(data$Wt_FishAgeComps[1,which(data$UseFishAgeComps[,,1] == 1),1,1])))

ESS_LLF_age = c(na.omit(data$ISS_FishAgeComps[1,which(data$UseFishAgeComps[,,1] == 1),1,1] * 
                          unique(data$Wt_FishAgeComps[1,which(data$UseFishAgeComps[,,1] == 1),1,1])))

ESS_tab_LLF_Age <- data.frame(ISS = ISS_LLF_age, Data_Wt = rep(wt_LLF_age,length(ISS_LLF_age)), ESS = ESS_LLF_age) 
write.csv(ESS_tab_LLF_Age, here(root_folder,mod_name,"plots",paste0(mod_name,"_LLF_ESS_Age.csv")))

ESS_tab_LLF_Age_mean <- ESS_tab_LLF_Age %>% summarize(ISS = mean(ISS), Data_Wt = mean(Data_Wt), ESS = mean(ESS)) %>% 
  mutate(Source = "Fixed Gear Fish Ages") %>% dplyr::relocate(Source) %>%
  mutate(dplyr::across(where(is.numeric), round, 2))

ISS_TF_len = c(na.omit(data$ISS_FishLenComps[1,which(data$UseFishLenComps[,,2] == 1),1,2]))

wt_TF_len = c(na.omit(unique(data$Wt_FishLenComps[1,which(data$UseFishLenComps[,,2] == 1),1,2])))

ESS_TF_len = c(na.omit(data$ISS_FishLenComps[1,which(data$UseFishLenComps[,,2] == 1),1,2] * 
                         unique(data$Wt_FishLenComps[1,which(data$UseFishLenComps[,,2] == 1),1,2])))

ESS_tab_TF_len <- data.frame(ISS = ISS_TF_len, Data_Wt = rep(wt_TF_len,length(ISS_TF_len)), ESS = ESS_TF_len) 
write.csv(ESS_tab_TF_len, here(root_folder,mod_name,"plots",paste0(mod_name,"_TF_ESS_Len.csv")))

ESS_tab_TF_len_mean <- ESS_tab_TF_len %>% summarize(ISS = mean(ISS), Data_Wt = mean(Data_Wt), ESS = mean(ESS)) %>% 
  mutate(Source = "Trawl Fish Lengths") %>% dplyr::relocate(Source) %>%
  mutate(dplyr::across(where(is.numeric), round, 2))

ISS_BS_LLS_age = c(na.omit(data$ISS_SrvAgeComps[1,which(data$UseSrvAgeComps[,,1] == 1),1,1]))

wt_BS_LLS_age = c(na.omit(unique(data$Wt_SrvAgeComps[1,which(data$UseSrvAgeComps[,,1] == 1),1,1])))

ESS_BS_LLS_age = c(na.omit(data$ISS_SrvAgeComps[1,which(data$UseSrvAgeComps[,,1] == 1),1,1] * 
                             unique(data$Wt_SrvAgeComps[1,which(data$UseSrvAgeComps[,,1] == 1),1,1])))

ESS_tab_BS_LLS_Age <- data.frame(ISS = ISS_BS_LLS_age, Data_Wt = rep(wt_BS_LLS_age,length(ISS_BS_LLS_age)), ESS = ESS_BS_LLS_age) 
write.csv(ESS_tab_BS_LLS_Age, here(root_folder,mod_name,"plots",paste0(mod_name,"_BS_LLS_ESS_Age.csv")))

ESS_tab_BS_LLS_Age_mean <- ESS_tab_BS_LLS_Age %>% summarize(ISS = mean(ISS), Data_Wt = mean(Data_Wt), ESS = mean(ESS)) %>% 
  mutate(Source = "BS LL Survey Ages") %>% dplyr::relocate(Source) %>%
  mutate(dplyr::across(where(is.numeric), round, 2))

ISS_AI_LLS_age = c(na.omit(data$ISS_SrvAgeComps[1,which(data$UseSrvAgeComps[,,2] == 1),1,2]))

wt_AI_LLS_age = c(na.omit(unique(data$Wt_SrvAgeComps[1,which(data$UseSrvAgeComps[,,2] == 1),1,2])))

ESS_AI_LLS_age = c(na.omit(data$ISS_SrvAgeComps[1,which(data$UseSrvAgeComps[,,2] == 1),1,2] * 
                             unique(data$Wt_SrvAgeComps[1,which(data$UseSrvAgeComps[,,2] == 1),1,2])))

ESS_tab_AI_LLS_Age <- data.frame(ISS = ISS_AI_LLS_age, Data_Wt = rep(wt_AI_LLS_age,length(ISS_AI_LLS_age)), ESS = ESS_AI_LLS_age) 
write.csv(ESS_tab_AI_LLS_Age, here(root_folder,mod_name,"plots",paste0(mod_name,"_AI_LLS_ESS_Age.csv")))

ESS_tab_AI_LLS_Age_mean <- ESS_tab_AI_LLS_Age %>% summarize(ISS = mean(ISS), Data_Wt = mean(Data_Wt), ESS = mean(ESS)) %>% 
  mutate(Source = "AI LL Survey Ages") %>% dplyr::relocate(Source) %>%
  mutate(dplyr::across(where(is.numeric), round, 2))

ISS_GOA_LLS_age = c(na.omit(data$ISS_SrvAgeComps[1,which(data$UseSrvAgeComps[,,3] == 1),1,3]))

wt_GOA_LLS_age = c(na.omit(unique(data$Wt_SrvAgeComps[1,which(data$UseSrvAgeComps[,,3] == 1),1,3])))

ESS_GOA_LLS_age = c(na.omit(data$ISS_SrvAgeComps[1,which(data$UseSrvAgeComps[,,3] == 1),1,3] * 
                              unique(data$Wt_SrvAgeComps[1,which(data$UseSrvAgeComps[,,3] == 1),1,3])))

ESS_tab_GOA_LLS_Age <- data.frame(ISS = ISS_GOA_LLS_age, Data_Wt = rep(wt_GOA_LLS_age,length(ISS_GOA_LLS_age)), ESS = ESS_GOA_LLS_age) 
write.csv(ESS_tab_GOA_LLS_Age, here(root_folder,mod_name,"plots",paste0(mod_name,"_GOA_LLS_ESS_Age.csv")))

ESS_tab_GOA_LLS_Age_mean <- ESS_tab_GOA_LLS_Age %>% summarize(ISS = mean(ISS), Data_Wt = mean(Data_Wt), ESS = mean(ESS)) %>% 
  mutate(Source = "GOA LL Survey Ages") %>% dplyr::relocate(Source) %>%
  mutate(dplyr::across(where(is.numeric), round, 2))

ISS_JPN_LLS_age = c(na.omit(data$ISS_SrvAgeComps[1,which(data$UseSrvAgeComps[,,5] == 1),1,5]))

wt_JPN_LLS_age = c(na.omit(unique(data$Wt_SrvAgeComps[1,which(data$UseSrvAgeComps[,,5] == 1),1,5])))

ESS_JPN_LLS_age = c(na.omit(data$ISS_SrvAgeComps[1,which(data$UseSrvAgeComps[,,5] == 1),1,5] * 
                              unique(data$Wt_SrvAgeComps[1,which(data$UseSrvAgeComps[,,5] == 1),1,5])))

ESS_tab_JPN_LLS_Age <- data.frame(ISS = ISS_JPN_LLS_age, Data_Wt = rep(wt_JPN_LLS_age,length(ISS_JPN_LLS_age)), ESS = ESS_JPN_LLS_age) 
write.csv(ESS_tab_JPN_LLS_Age, here(root_folder,mod_name,"plots",paste0(mod_name,"_JPN_LLS_ESS_Age.csv")))

ESS_tab_JPN_LLS_Age_mean <- ESS_tab_JPN_LLS_Age %>% summarize(ISS = mean(ISS), Data_Wt = mean(Data_Wt), ESS = mean(ESS)) %>% 
  mutate(Source = "JPN LL Survey Ages") %>% dplyr::relocate(Source) %>%
  mutate(dplyr::across(where(is.numeric), round, 2))

ISS_JPN_LLS_len = c(na.omit(data$ISS_SrvLenComps[1,which(data$UseSrvLenComps[,,5] == 1),1,5]))

wt_JPN_LLS_len = c(na.omit(unique(data$Wt_SrvLenComps[1,which(data$UseSrvLenComps[,,5] == 1),1,5])))

ESS_JPN_LLS_len = c(na.omit(data$ISS_SrvLenComps[1,which(data$UseSrvLenComps[,,5] == 1),1,5] * 
                              unique(data$Wt_SrvLenComps[1,which(data$UseSrvLenComps[,,5] == 1),1,5])))

ESS_tab_JPN_LLS_len <- data.frame(ISS = ISS_JPN_LLS_len, Data_Wt = rep(wt_JPN_LLS_len,length(ISS_JPN_LLS_len)), ESS = ESS_JPN_LLS_len) 
write.csv(ESS_tab_JPN_LLS_len, here(root_folder,mod_name,"plots",paste0(mod_name,"_JPN_LLS_ESS_Age.csv")))

ESS_tab_JPN_LLS_Len_mean <- ESS_tab_JPN_LLS_len %>% summarize(ISS = mean(ISS), Data_Wt = mean(Data_Wt), ESS = mean(ESS)) %>% 
  mutate(Source = "JPN LL Survey Lengths") %>% dplyr::relocate(Source) %>%
  mutate(dplyr::across(where(is.numeric), round, 2))

ESS_tab <- rbind(ESS_tab_LLF_Age_mean, ESS_tab_BS_LLS_Age_mean, ESS_tab_AI_LLS_Age_mean, ESS_tab_GOA_LLS_Age_mean, 
                 ESS_tab_JPN_LLS_Age_mean, ESS_tab_TF_len_mean, ESS_tab_JPN_LLS_Len_mean)



SE_BS_LLS_idx = c(na.omit(data$ObsSrvIdx_SE[1,which(data$UseSrvIdx[,,1] == 1),1]))

wt_BS_LLS_idx = data$Wt_SrvIdx

ESE_BS_LLS_idx = c(na.omit(data$ObsSrvIdx_SE[1,which(data$UseSrvIdx[,,1] == 1),1] * 
                             data$Wt_SrvIdx))

SE_tab_BS_LLS_idx <- data.frame(SE = SE_BS_LLS_idx, Data_Wt = rep(wt_BS_LLS_idx,length(SE_BS_LLS_idx)), Eff_SE = ESE_BS_LLS_idx) 
write.csv(SE_tab_BS_LLS_idx, here(root_folder,mod_name,"plots",paste0(mod_name,"_BS_LLS_Idx_SE.csv")))

SE_tab_BS_LLS_mean <- SE_tab_BS_LLS_idx %>% summarize(SE = mean(SE), Data_Wt = mean(Data_Wt), Eff_SE = mean(Eff_SE)) %>% 
  mutate(Source = "BS LL Survey Index") %>% dplyr::relocate(Source) %>%
  mutate(dplyr::across(where(is.numeric), round, 2))

SE_AI_LLS_idx = c(na.omit(data$ObsSrvIdx_SE[1,which(data$UseSrvIdx[,,2] == 1),2]))

wt_AI_LLS_idx = data$Wt_SrvIdx

ESE_AI_LLS_idx = c(na.omit(data$ObsSrvIdx_SE[1,which(data$UseSrvIdx[,,2] == 1),2] * 
                             data$Wt_SrvIdx))

SE_tab_AI_LLS_idx <- data.frame(SE = SE_AI_LLS_idx, Data_Wt = rep(wt_AI_LLS_idx,length(SE_AI_LLS_idx)), Eff_SE = ESE_AI_LLS_idx) 
write.csv(SE_tab_AI_LLS_idx, here(root_folder,mod_name,"plots",paste0(mod_name,"_AI_LLS_Idx_SE.csv")))

SE_tab_AI_LLS_mean <- SE_tab_AI_LLS_idx %>% summarize(SE = mean(SE), Data_Wt = mean(Data_Wt), Eff_SE = mean(Eff_SE)) %>% 
  mutate(Source = "AI LL Survey Index") %>% dplyr::relocate(Source) %>%
  mutate(dplyr::across(where(is.numeric), round, 2))

SE_GOA_LLS_idx = c(na.omit(data$ObsSrvIdx_SE[1,which(data$UseSrvIdx[,,3] == 1),3]))

wt_GOA_LLS_idx = data$Wt_SrvIdx

ESE_GOA_LLS_idx = c(na.omit(data$ObsSrvIdx_SE[1,which(data$UseSrvIdx[,,3] == 1),3] * 
                              data$Wt_SrvIdx))

SE_tab_GOA_LLS_idx <- data.frame(SE = SE_GOA_LLS_idx, Data_Wt = rep(wt_GOA_LLS_idx,length(SE_GOA_LLS_idx)), Eff_SE = ESE_GOA_LLS_idx) 
write.csv(SE_tab_GOA_LLS_idx, here(root_folder,mod_name,"plots",paste0(mod_name,"_GOA_LLS_Idx_SE.csv")))

SE_tab_GOA_LLS_mean <- SE_tab_GOA_LLS_idx %>% summarize(SE = mean(SE), Data_Wt = mean(Data_Wt), Eff_SE = mean(Eff_SE)) %>% 
  mutate(Source = "GOA LL Survey Index") %>% dplyr::relocate(Source) %>%
  mutate(dplyr::across(where(is.numeric), round, 2))

SE_JPN_LLS_idx = c(na.omit(data$ObsSrvIdx_SE[1,which(data$UseSrvIdx[,,5] == 1),5]))

wt_JPN_LLS_idx = data$Wt_SrvIdx

ESE_JPN_LLS_idx = c(na.omit(data$ObsSrvIdx_SE[1,which(data$UseSrvIdx[,,5] == 1),5] * 
                              data$Wt_SrvIdx))

SE_tab_JPN_LLS_idx <- data.frame(SE = SE_JPN_LLS_idx, Data_Wt = rep(wt_JPN_LLS_idx,length(SE_JPN_LLS_idx)), Eff_SE = ESE_JPN_LLS_idx) 
write.csv(SE_tab_JPN_LLS_idx, here(root_folder,mod_name,"plots",paste0(mod_name,"_JPN_LLS_Idx_SE.csv")))

SE_tab_JPN_LLS_mean <- SE_tab_JPN_LLS_idx %>% summarize(SE = mean(SE), Data_Wt = mean(Data_Wt), Eff_SE = mean(Eff_SE)) %>% 
  mutate(Source = "JPN LL Survey Index") %>% dplyr::relocate(Source) %>%
  mutate(dplyr::across(where(is.numeric), round, 2))


SE_tab <- rbind(SE_tab_BS_LLS_mean, SE_tab_AI_LLS_mean, SE_tab_GOA_LLS_mean, SE_tab_JPN_LLS_mean)

ESS_SE_tab <- full_join(SE_tab,ESS_tab) %>% dplyr::relocate(Eff_SE, .after = SE)
write.csv(ESS_SE_tab,here(root_folder,mod_name,"plots",paste0(mod_name,"_ESS_SE.csv")))

tiff(filename=here(root_folder,mod_name,"plots",paste0(mod_name,"_ESS_SE_Tab.png")),units = "in", width=10,height=6, res = 300)
ESS_SE_tab_plot <- gridExtra::grid.table(ESS_SE_tab, rows=NULL)
dev.off()

##############################################################

pdf(file = here(root_folder,mod_name,"plots",paste0(mod_name,"_all_results.pdf")),  width = 25, height = 13)

gridExtra::grid.table(ESS_SE_tab, rows=NULL)

plot.new()
text(0.1, 0.5, messages, cex = 1.5)

if(do_francis ==1){
  #print(Francis_plot1)
  print(plot_francis)
}
print(messages)
print(out[[2]])
print(ssb_recr)
print(bio)
print(naa)
print(selex_plot[[1]])
print(selex_plot[[2]])
print(catch_plot)
print(idx_plot)
print(runs_tst)
print(llf_age_fem)
print(llf_age_coh_fem)

if(unique(data$FishAgeComps_Type[,1]) > 0){
  print(llf_age_mal)
  print(llf_age_coh_mal)
}
if(any(data$UseFishLenComps[,,1] == 1)){
  print(llf_len_fem)
  print(llf_len_male)
}
print(tf_len_fem)
print(tf_len_male)


print(BS_lls_age_fem)
print(BS_lls_age_coh_fem)

if(unique(data$SrvAgeComps_Type[,1]) > 0){
  print(BS_lls_age_mal)
  print(BS_lls_age_coh_mal)
}

print(AI_lls_age_fem)
print(AI_lls_age_coh_fem)
if(unique(data$SrvAgeComps_Type[,2]) > 0){
  print(AI_lls_age_mal)
  print(AI_lls_age_coh_mal)
}
print(GOA_lls_age_fem)
print(GOA_lls_age_coh_fem)
if(unique(data$SrvAgeComps_Type[,3]) > 0){
  print(GOA_lls_age_mal)
  print(GOA_lls_age_coh_mal)
}

print(jpn_lls_age_fem)
print(jpn_lls_age_coh_fem)

if(unique(data$SrvAgeComps_Type[,5]) > 0){
  print(jpn_lls_age_mal)
  print(jpn_lls_age_coh_mal)
}

if(any(data$UseSrvLenComps[,,1] == 1)){
  print(BS_lls_len_fem)
  print(BS_lls_len_male)
}
if(any(data$UseSrvLenComps[,,2] == 1)){
  print(AI_lls_len_fem)
  print(AI_lls_len_male)
}
if(any(data$UseSrvLenComps[,,3] == 1)){
  print(GOA_lls_len_fem)
  print(GOA_lls_len_male)
}

if(any(data$UseSrvLenComps[,,4] == 1)){
  print(ts_len_fem)
  print(ts_len_male)
}

if(any(data$UseSrvLenComps[,,5] == 1)){
  print(jpn_lls_len_fem)
  print(jpn_lls_len_male)
}

if(do_OSA == 1){
  
  print(osa_llf_age)
  if(any(data$UseFishLenComps[,,1] == 1)){
    print(osa_llf_len)
  }
  print(osa_tf)
  print(osa_lls_age_BS)
  print(osa_lls_age_AI)
  print(osa_lls_age_GOA)
  
  if(any(data$UseSrvLenComps[,,1] == 1)){
    print(osa_lls_len_BS)
  }
  if(any(data$UseSrvLenComps[,,2] == 1)){
    print(osa_lls_len_AI)
  }
  if(any(data$UseSrvLenComps[,,3] == 1)){
    print(osa_lls_len_GOA)
  }
  
  
  if(any(data$UseSrvLenComps[,,4] == 1)){
    print(osa_ts_len)
  }
  print(osa_jpn_lls_age)
  
  if(any(data$UseSrvLenComps[,,5] == 1)){
    print(osa_jpn_lls_len)
  }
}

if(do_pearson == 1){
  
  print(llf_age_pear)
  if(any(data$UseFishLenComps[,,1] == 1)){
    print(llf_len_pear)
  }
  print(tf_len_pear)
  print(bs_lls_age_pear)
  print(ai_lls_age_pear)
  print(goa_lls_age_pear)
  
  if(any(data$UseSrvLenComps[,,1] == 1)){
    print(bs_lls_len_pear)
  }
  if(any(data$UseSrvLenComps[,,2] == 1)){
    print(ai_lls_len_pear)
  }
  if(any(data$UseSrvLenComps[,,3] == 1)){
    print(goa_lls_len_pear)
  }
  
  print(jpn_lls_age_pear)
  
  if(any(data$UseSrvLenComps[,,5] == 1)){
    print(jpn_lls_len_pear)
  }
}

if(do_diags == 1){
  print(retro_plot[[1]])
  print(retro_plot[[2]])
  print(retro_plot[[3]])
  print(jitter_ssb)
  print(jitter_grad)
  print(unsmoothed_prof)
  print(ll_prof_unsmooth_plot_type)
  print(ll_prof_plot)
  print(ll_prof_plot_type)
  print(unsmoothed_prof_M_fem)
  print(ll_prof_unsmooth_plot_type_M_fem)
  print(ll_prof_plot_M_fem)
  print(ll_prof_plot_type_M_fem)
  print(unsmoothed_prof_M_mal)
  print(ll_prof_unsmooth_plot_type_M_mal)
  print(ll_prof_plot_M_mal)
  print(ll_prof_plot_type_M_mal)
}  # do_diags loop

dev.off()

save.image(file=here(root_folder,mod_name,paste0(mod_name,".RData")))        # Save all the data frames in case want to reload data without doing full data pulls

