# Purpose: Spatial Sablefish Sensitivity Run - Multinomial Release Conditioned (3 region model)
# Creator: Matthew LH. Cheng
# Date Created: 3/31/26

# Set up ------------------------------------------------------------------
library(here)
library(tidyverse)
library(RTMB)
library(SPoRC)

# Setup folder
root_folder <- here("Sept PT Model Runs", "25_15_Spatial")
mod_name <- "Spatial_MltRel"

# read in 2025 assessment data
dat_new_2025 <- readRDS(here("Mods", "Spat_Cont", "Spat_3_Reg", "SPoRC_Spatial_3_Regions.rds"))

  # Setup Model -------------------------------------------------------------
# Initialize model dimensions and data list
input_list <- Setup_Mod_Dim(years = 1:length(dat_new_2025$yrs),
                            # vector of years (1 - 62)
                            ages = 1:length(dat_new_2025$ages),
                            # vector of ages (1 - 30)
                            lens = dat_new_2025$length_bins,
                            # number of lengths (41 - 99)
                            n_regions = dat_new_2025$nreg,
                            # number of regions (5)
                            n_sexes = dat_new_2025$nsex,
                            # number of sexes (2)
                            n_fish_fleets = dat_new_2025$nflt,
                            # number of fishery fleet (2)
                            n_srv_fleets = dat_new_2025$nsrv_flt,
                            # number of survey fleets (3)
                            verbose = TRUE
)

# Setup recruitment stuff (using defaults for other stuff)
input_list <- Setup_Mod_Rec(input_list = input_list, # input data list from above
                            do_rec_bias_ramp = 0,
                            sigmaR_switch = as.integer(length(dat_new_2025$styr:1975)),
                            
                            # Model options
                            rec_model = "mean_rec", # recruitment model
                            sigmaR_spec = "fix", # fixing
                            ln_sigmaR = log(c(0.4, 1.2)),
                            # values to fix sigmaR at, or starting values
                            ln_global_R0 = log(25), 
                            Use_Rec_prop_Prior = 1, 
                            Rec_prop_prior = 3.5
)

# Setup biological stuff (using defaults for other stuff)
input_list <- Setup_Mod_Biologicals(input_list = input_list,
                                    WAA = dat_new_2025$WAA, # weight at age
                                    MatAA = dat_new_2025$MatAA, # maturity at age
                                    AgeingError = as.matrix(dat_new_2025$AgeError),
                                    # ageing error matrix
                                    fit_lengths = 1, # fitting lengths
                                    SizeAgeTrans = dat_new_2025$SizeAgeTrans,
                                    # size age transition matrix
                                    M_spec = "fix", # fix natural mortality
                                    # values to fix natural mortality at
                                    Fixed_natmort = {
                                      Fixed_natmort <- array(NA,
                                                             dim = c(input_list$data$n_regions,
                                                                     length(input_list$data$years),
                                                                     length(input_list$data$ages),
                                                                     input_list$data$n_sexes))
                                      
                                      Fixed_natmort[,,,1] <- 0.09212236  # Females 
                                      Fixed_natmort[,,,2] <- 0.09676652  # Males 
                                      Fixed_natmort
                                    }
)

# setting up movement parameterization
input_list <- Setup_Mod_Movement(input_list = input_list,
                                 # Model options
                                 Movement_ageblk_spec = list(c(1:6), c(7:15), c(16:30)),
                                 # estimating movement in 3 age blocks
                                 # (ages 1-6, ages 7-15, ages 16-30)
                                 Movement_yearblk_spec = "constant", # time-invariant movement
                                 Movement_sexblk_spec = "constant", # sex-invariant movement
                                 do_recruits_move = 0, # recruits do not move
                                 use_fixed_movement = 0, # estimating movement
                                 Use_Movement_Prior = 0 # priors used for movement
)


# setting up tagging parameterization

# setup tagging priors
tag_prior <- data.frame(
  region = 1,
  block = 1:3,
  mu = NA, # no mean, since symmetric beta
  sd = 5, # sd = 5
  type = 0 # symmetric beta
)

input_list <- Setup_Mod_Tagging(input_list = input_list,
                                UseTagging = 1, # using tagging data
                                max_tag_liberty = 15, # maximum number of years to track a cohort
                                
                                # Data Inputs
                                tag_release_indicator = dat_new_2025$tag_release_indicator,
                                # tag release indicator (first col = tag region,
                                # second col = tag year),
                                # total number of rows = number of tagged cohorts
                                Tagged_Fish = dat_new_2025$Tagged_Fish, # Released fish
                                # dimensioned by total number of tagged cohorts, (implicitly
                                # tracks the release year and region), age, and sex
                                Obs_Tag_Recap = dat_new_2025$Obs_Tag_Recap,
                                # dimensioned by max tag liberty, tagged cohorts, regions,
                                # ages, and sexes
                                
                                # Model options
                                Tag_LikeType = "Multinomial_Release", # Multinomial Release Conditioned
                                mixing_period = 2, # Don't fit tagging until release year + 1
                                t_tagging = 0.5, # tagging happens midway through the year,
                                # movement does not occur within that year
                                tag_selex = "SexSp_DomFleet", # tagging recapture selectivity is the dominant fleet (fixed-gear)
                                tag_natmort = "AgeSp_SexSp", # tagging natural mortality is age and sex-specific
                                Use_TagRep_Prior = 1, # tag reporting rate priors are used
                                TagRep_Prior = tag_prior,
                                move_age_tag_pool = as.list(1:30), # whether or not to pool tagging data when fitting (for computational cost)
                                move_sex_tag_pool = list(c(1:2)), # whether or not to pool sex-specific data when fitting
                                Init_Tag_Mort_spec = "fix", # fixing initial tag mortality
                                Tag_Shed_spec = "fix", # fixing chronic shedding
                                TagRep_spec = "est_shared_r", # tag reporting rates are not region specific
                                Tag_Reporting_blocks = c(
                                  paste("Block_1_Year_1-35_Region_", 
                                        c(1:input_list$data$n_regions), sep = ''),
                                  paste("Block_2_Year_36-terminal_Region_", 
                                        c(1:input_list$data$n_regions), sep = ''),
                                  paste("Block_3_Year_56-terminal_Region_", 
                                        c(1:input_list$data$n_regions), sep = '')
                                ),
                                # Specify starting values or fixing values
                                ln_Init_Tag_Mort = log(0.1), # fixing initial tag mortality
                                ln_Tag_Shed = log(0.02),  # fixing tag shedding
                                ln_tag_theta = log(0.5), # starting value for tagging overdispersion
                                Tag_Reporting_Pars = array(log(0.2 / (1-0.2)), dim = c(input_list$data$n_regions, 3))  # starting values for tag reporting pars
)

# setting up catch data
dat_new_2025$UseCatch[is.na(dat_new_2025$ObsCatch)] <- 0 # replace w/ 0s

input_list <- Setup_Mod_Catch_and_F(input_list = input_list,
                                    # Data inputs
                                    ObsCatch = dat_new_2025$ObsCatch,
                                    Catch_Type = array(1, dim = c(length(input_list$data$years), input_list$data$n_fish_fleets)),
                                    UseCatch = dat_new_2025$UseCatch,
                                    # Model options
                                    Use_F_pen = 1,
                                    # whether to use f penalty, == 0 don't use, == 1 use
                                    sigmaC_spec = 'fix',
                                    ln_sigmaC =
                                      array(log(0.05), dim = c(input_list$data$n_regions, 
                                                               length(input_list$data$years), 
                                                               input_list$data$n_fish_fleets)),
                                    # fixing catch sd at small value
                                    ln_F_mean = array(-2, dim = c(input_list$data$n_regions,
                                                                  input_list$data$n_fish_fleets))
                                    # some starting values for fishing mortality
)

# Fishery Indices and Compositions

# Make sure the use indicators for compositions are 0s when NAs
dat_new_2025$UseFishAgeComps[is.na(apply(dat_new_2025$ObsFishAgeComps, c(1,2,5), sum))] <- 0
dat_new_2025$UseFishLenComps[is.na(apply(dat_new_2025$ObsFishLenComps, c(1,2,5), sum))] <- 0

input_list <- Setup_Mod_FishIdx_and_Comps(input_list = input_list,
                                          # data inputs
                                          ObsFishIdx = dat_new_2025$ObsFishIdx,
                                          ObsFishIdx_SE = dat_new_2025$ObsFishIdx_SE,
                                          UseFishIdx =  dat_new_2025$UseFishIdx,
                                          ObsFishAgeComps = dat_new_2025$ObsFishAgeComps,
                                          UseFishAgeComps = dat_new_2025$UseFishAgeComps,
                                          ISS_FishAgeComps = dat_new_2025$ISS_FishAgeComps,
                                          ObsFishLenComps = dat_new_2025$ObsFishLenComps,
                                          UseFishLenComps = dat_new_2025$UseFishLenComps,
                                          ISS_FishLenComps = dat_new_2025$ISS_FishLenComps,
                                          
                                          # Model options
                                          fish_idx_type = c("none", "none"),
                                          # fishery indices not used
                                          FishAgeComps_LikeType =
                                            c("Multinomial", "none"),
                                          # age comp likelihoods for fishery fleet 1 and 2
                                          FishLenComps_LikeType =
                                            c("Multinomial", "Multinomial"),
                                          # length comp likelihoods for fishery fleet 1 and 2
                                          FishAgeComps_Type =
                                            c("spltRjntS_Year_1-terminal_Fleet_1",
                                              "none_Year_1-terminal_Fleet_2"),
                                          # age comp structure for fishery fleet 1 and 2
                                          FishLenComps_Type =
                                            c("spltRjntS_Year_1-terminal_Fleet_1",
                                              "spltRjntS_Year_1-terminal_Fleet_2")
)

# Survey Indices and Compositions
# Note: These are now sds (i.e., the tmb likelihoods) for the domestic and JP LLS, but not for the trawl survey (TS is still CV * index).
# Not a big issue because the trawl survey is not used here.

# Make sure the use indicators for compositions are 0s when NAs
dat_new_2025$UseSrvAgeComps[is.na(apply(dat_new_2025$ObsSrvAgeComps, c(1,2,5), sum))] <- 0
dat_new_2025$UseSrvLenComps[is.na(apply(dat_new_2025$ObsSrvLenComps, c(1,2,5), sum))] <- 0

# Don't fit trawl survey
dat_new_2025$UseSrvLenComps[,,2] <- 0
dat_new_2025$UseSrvIdx[,,2] <- 0

input_list <- Setup_Mod_SrvIdx_and_Comps(input_list = input_list,
                                         # data inputs
                                         ObsSrvIdx = dat_new_2025$ObsSrvIdx,
                                         ObsSrvIdx_SE = dat_new_2025$ObsSrvIdx_SE * 1.5,
                                         UseSrvIdx =  dat_new_2025$UseSrvIdx,
                                         ObsSrvAgeComps = dat_new_2025$ObsSrvAgeComps,
                                         ISS_SrvAgeComps = dat_new_2025$ISS_SrvAgeComps,
                                         UseSrvAgeComps = dat_new_2025$UseSrvAgeComps,
                                         ObsSrvLenComps = dat_new_2025$ObsSrvLenComps,
                                         UseSrvLenComps = dat_new_2025$UseSrvLenComps,
                                         ISS_SrvLenComps = dat_new_2025$ISS_SrvLenComps,
                                         
                                         # Model options
                                         srv_idx_type = c("abd", 'biom', "abd"),
                                         # abundance and biomass for survey fleet 1 and 2
                                         SrvAgeComps_LikeType =
                                           c("Multinomial", 'none', "Multinomial"),
                                         # survey age composition likelihood for survey fleet
                                         # 1, and 2
                                         SrvLenComps_LikeType =
                                           c("none", "none", "none"),
                                         #  no length compositions used for survey
                                         SrvAgeComps_Type = c("spltRjntS_Year_1-terminal_Fleet_1",
                                                              "none_Year_1-terminal_Fleet_2",
                                                              "spltRjntS_Year_1-terminal_Fleet_3"),
                                         # survey age comp type
                                         SrvLenComps_Type = c("none_Year_1-terminal_Fleet_1",
                                                              "none_Year_1-terminal_Fleet_2",
                                                              "none_Year_1-terminal_Fleet_3")
)

# Fishery Selectivity and Catchability

# defining priors
sex_par <- expand.grid(sex = 1:2, par = 1:2)
fleet_blocks <- data.frame(
  fleet = c(1,1, 2),                                       
  block = c(1,2, 1)
)

# merge together (note that unlike the operational assessment, selectivity
# blocks are reduced from 3 to 2)
fish_selex_structure <- merge(fleet_blocks, sex_par) 

# Add the lognormal prior values - creates a dataframe, each row is a unique parameter combination to apply the prior to
fish_selex_prior <- cbind(
  region = 1,
  fish_selex_structure,
  mu = 4,                                                                      # All selex means = 1 (means should be defined in normal space)
  sd = 1.5                                                                       # All selex sd = 5
)

input_list <- Setup_Mod_Fishsel_and_Q(input_list = input_list,
                                      
                                      # Model options
                                      cont_tv_fish_sel = c("none_Fleet_1", "none_Fleet_2"),
                                      # fishery selectivity, whether continuous time-varying
                                      
                                      # fishery selectivity blocks
                                      fish_sel_blocks =
                                        c("Block_1_Year_1-56_Fleet_1",      
                                          "Block_2_Year_57-terminal_Fleet_1",
                                          "none_Fleet_2"),
                                      # no blocks for trawl fishery
                                      
                                      # fishery selectivity form
                                      fish_sel_model =
                                        c("logist1_Fleet_1", "logist1_Fleet_2"),
                                      
                                      # fishery catchability blocks
                                      fish_q_blocks =
                                        c("none_Fleet_1", "none_Fleet_2"),
                                      # no blocks since q is not estimated
                                      
                                      # sharing fishery selex parameters
                                      fish_fixed_sel_pars =
                                        c("est_shared_r", "est_shared_r"),
                                      
                                      # whether to estimate all fixed effects
                                      # for fishery catchability
                                      fish_q_spec =
                                        c("fix", "fix"),
                                      Use_fish_selex_prior = 1,
                                      fish_selex_prior = fish_selex_prior
)


# setup survey selectivity
# Define sex and parameter combinations
sex_par <- expand.grid(sex = 1:2, par = 1:2)

# Define valid fleet-block combinations (only estimating domestic and jp LLS)
fleet_blocks <- data.frame(
  fleet = c(1, 3),
  block = c(1, 1)
)

# Merge to get all valid combinations
srv_selex_structure <- merge(fleet_blocks, sex_par)

# Add the lognormal prior values - creates a dataframe, each row is a unique parameter combination to apply the prior to
srv_selex_prior <- cbind(
  region = 1,
  srv_selex_structure,
  mu = 3,
  sd = 1.5
) %>%
  filter(!(fleet == 3 & par == 2 & sex == 2))

input_list <- Setup_Mod_Srvsel_and_Q(input_list = input_list,
                                     
                                     # Model options
                                     # survey selectivity, whether continuous time-varying
                                     cont_tv_srv_sel =
                                       c("none_Fleet_1",
                                         "none_Fleet_2",
                                         "none_Fleet_3"),
                                     
                                     # survey selectivity blocks
                                     srv_sel_blocks =                          # survey selectivity time blocks if not TV specified above for a given fleet
                                       c("none_Fleet_1",
                                         "none_Fleet_2",                       # No blocks for trawl survey
                                         "none_Fleet_3"                        # No blocks for JPN LLS
                                       ),
                                     
                                     # survey selectivity form
                                     srv_sel_model =
                                       c("logist1_Fleet_1",
                                         "exponential_Fleet_2",
                                         "logist1_Fleet_3"
                                       ),
                                     
                                     # survey catchability blocks
                                     srv_q_blocks =
                                       c("none_Fleet_1",
                                         "none_Fleet_2",
                                         "none_Fleet_3"),
                                     
                                     # whether to estiamte all fixed effects
                                     # for survey selectivity and later
                                     # modify to fix/share parameters
                                     srv_fixed_sel_pars_spec =
                                       c("est_shared_r",
                                         "fix",
                                         "est_shared_r"),
                                     
                                     # whether to estiamte all
                                     # fixed effects for survey catchability
                                     # spatially-invariant q
                                     srv_q_spec =
                                       c("est_shared_r",
                                         "fix",
                                         "est_shared_r"),
                                     Use_srv_selex_prior = 1,
                                     srv_selex_prior = srv_selex_prior
)

# set up model weighting stuff
input_list <- Setup_Mod_Weighting(input_list = input_list,
                                  Wt_Catch = 1,
                                  Wt_FishIdx = 1,
                                  Wt_SrvIdx = 1,
                                  Wt_Rec = 1,
                                  Wt_F = 1,
                                  Wt_Tagging = 0.1,
                                  # Composition model weighting
                                  Wt_FishAgeComps =
                                    array(1, dim = c(input_list$data$n_regions,
                                                     length(input_list$data$years),
                                                     input_list$data$n_sexes,
                                                     input_list$data$n_fish_fleets)),
                                  Wt_FishLenComps =
                                    array(1, dim = c(input_list$data$n_regions,
                                                     length(input_list$data$years),
                                                     input_list$data$n_sexes,
                                                     input_list$data$n_fish_fleets)),
                                  Wt_SrvAgeComps =
                                    array(1, dim = c(input_list$data$n_regions,
                                                     length(input_list$data$years),
                                                     input_list$data$n_sexes,
                                                     input_list$data$n_srv_fleets)),
                                  Wt_SrvLenComps =
                                    array(1, dim = c(input_list$data$n_regions,
                                                     length(input_list$data$years),
                                                     input_list$data$n_sexes,
                                                     input_list$data$n_srv_fleets))
)

# extract out lists updated with helper functions
data <- input_list$data
parameters <- input_list$par
mapping <- input_list$map


# Additional Model Specifications -----------------------------------------

# Survey Ages (~100 total across all regions)
data$ISS_SrvAgeComps[] <- 33  

# Fishery Ages 
data$ISS_FishAgeComps[1,,,] <- 25  # BS
data$ISS_FishAgeComps[2,,,] <- 25  # AI
data$ISS_FishAgeComps[3,,,] <- 50  # GOA

# Fishery Lengths - Fixed Gear
data$ISS_FishLenComps[1,,,1] <- 13  # BS
data$ISS_FishLenComps[2,,,1] <- 13  # AI
data$ISS_FishLenComps[3,,,1] <- 18  # GOA

# Fishery Lengths - Trawl Gear
data$ISS_FishLenComps[1,,,1] <- 25  # BS
data$ISS_FishLenComps[2,,,1] <- 10  # AI
data$ISS_FishLenComps[3,,,1] <- 10  # GOA

# Map off delta for JP LLS
map_srv_fixed <- array(mapping$ln_srv_fixed_sel_pars, dim = dim(parameters$ln_srv_fixed_sel_pars))
map_srv_fixed[,2,1,2,3]  <- map_srv_fixed[,2,1,1,3] # share deltas
mapping$ln_srv_fixed_sel_pars <- factor(map_srv_fixed)

# Some starting values to help out the model
parameters$ln_srv_fixed_sel_pars[] <- log(1)
parameters$ln_fish_fixed_sel_pars[,,,,1] <- log(3) # fixed gear
parameters$ln_fish_fixed_sel_pars[,,,,2] <- log(1) # trawl gear

# Fit Model ---------------------------------------------------------------
sabie_rtmb_model <- fit_model(data,
                              parameters,
                              mapping,
                              random = NULL,
                              newton_loops = 3,
                              silent = F
)

sabie_rtmb_model$sd_rep <- sdreport(sabie_rtmb_model)
sabie_rtmb_model$data <- data
sabie_rtmb_model$parameters <- parameters
sabie_rtmb_model$mapping <- mapping

saveRDS(sabie_rtmb_model, here("Mods", "Spat_Cont","Spat_3_Reg", "3_rg_MltRel.RDS"))

sabie_rtmb_model <- readRDS(here("Mods", "Spat_Cont", "Spat_3_Reg", "3_rg_MltRel.RDS"))
data <- sabie_rtmb_model$data
pdf(here("Mods", "Spat_Cont", "Spat_3_Reg", "figs", "main_res.pdf"), width = 19, height = 10)
get_ts_plot(list(sabie_rtmb_model$rep), list(sabie_rtmb_model$sd_rep), 1)
get_catch_fits_plot(list(data), list(sabie_rtmb_model$rep), 1)
get_idx_fits_plot(list(data), list(sabie_rtmb_model$rep), 1)
get_selex_plot(list(sabie_rtmb_model$rep), 1, year_indx = c(1,60))
dev.off()


# Composition Fits --------------------------------------------------------
rep <- sabie_rtmb_model$rep
ages <- data$ages
lens <- data$lens
years <- data$years
comp_regions <- data.frame(Region = 1:3, Region_Name = c("BS", 'AI', "GOA"))
sex_names <- data.frame(Sex = 1:2, Sex_Name = c("Female", 'Male'))

# Get composition fits
comp_fits <- get_comp_prop(data, rep, ages, lens, years)

### Fishery Ages ------------------------------------------------------------

pdf(here("Mods", "Spat_Cont", "Spat_3_Reg", "figs", 'llf_acomp_fits.pdf'), width = 19, height = 10) 
for(r in 1:nrow(comp_regions)) {
  for(s in 1:nrow(sex_names)) {
    
    # yearly comp fits
    print(
      comp_fits$Fishery_Ages %>%
        filter(Fleet == 1, Sex == s, Region == r) %>%
        mutate(Cohort = Year - Age) %>%
        ggplot() +
        geom_col(aes(x = Age, y = obs), alpha = 0.5) +
        geom_line(aes(x = Age, y = pred)) +
        facet_wrap(~Year) +
        theme_sablefish() +
        labs(x = 'Age', y = 'Proportion')+
        labs(title=paste("LLF Age Comps", sex_names$Sex_Name[s], comp_regions$Region_Name[r]))
      
    )
  } # end s loop
  
  osas <-  SPoRC::get_osa(obs_mat = array(comp_fits$Obs_FishAge_mat[r,,,,] + 1e-4, dim = c(1, dim(comp_fits$Obs_FishAge_mat)[-1])),
                          exp_mat = array(comp_fits$Pred_FishAge_mat[r,,,,] + 1e-4, dim = c(1, dim(comp_fits$Obs_FishAge_mat)[-1])),
                          N = array(c(na.omit(data$ISS_FishAgeComps[r,which(data$UseFishAgeComps[r,,1] == 1),1,1] * 
                                                unique(data$Wt_FishAgeComps[r,which(data$UseFishAgeComps[r,,1] == 1),1,1]))),
                                    dim = c(1, length(which(data$UseFishAgeComps[r,,1] == 1)))),
                          years = which(data$UseFishAgeComps[r,,1] == 1),
                          fleet = 1,
                          bins = ages,
                          comp_type = 2,                
                          bin_label = "Ages"
  )
  
  # Get SDNR QQ plot and Bubble Plot
  osas$res$region <- r # rename region for plotting
  osa_res <- SPoRC::plot_resids(osa_results = osas)
  
  # Get Aggregated Plot
  fishage_agg <- comp_fits$Fishery_Ages %>%
    group_by(Region, Age, Sex, Fleet) %>%
    summarize(obs = mean(obs), pred = mean(pred)) %>%
    filter(Fleet == 1, Region == r) %>%
    ggplot() +
    geom_col(aes(x = Age, y = obs), fill = 'darkgreen', alpha = 0.8) +
    geom_line(aes(x = Age, y = pred), col = 'black', lwd = 1.3) +
    theme_bw(base_size = 15) +
    facet_grid( Region~Sex, labeller = labeller(
      Region = function(x) paste0("Region ", x),
      Sex = function(x) paste0("Sex ", x)
    )) +
    labs(x = 'Age', y = 'Proportions')
  
  print(
    cowplot::plot_grid(osa_res[[2]], # Bubble
                       cowplot::plot_grid(osa_res[[1]], fishage_agg, nrow = 1), # QQ and Aggregated
                       ncol = 1, rel_heights = c(0.6, 0.4), labels = paste("LLF OSA Age Comps", comp_regions$Region_Name[r]))
  )
} # end r loop
dev.off()

### Fishery Lengths ------------------------------------------------------------
#### LLF ---------------------------------------------------------------------
pdf(here("Mods", "Spat_Cont", "Spat_3_Reg", "figs", 'llf_lcomp_fits.pdf'), width = 19, height = 10) 
for(r in 1:nrow(comp_regions)) {
  for(s in 1:nrow(sex_names)) {
    
    # yearly comp fits
    print(
      comp_fits$Fishery_Lens %>%
        filter(Fleet == 1, Sex == s, Region == r) %>%
        ggplot() +
        geom_col(aes(x = Len, y = obs), alpha = 0.5) +
        geom_line(aes(x = Len, y = pred)) +
        facet_wrap(~Year) +
        theme_sablefish() +
        labs(x = 'Age', y = 'Proportion')+
        labs(title=paste("LLF Length Comps", sex_names$Sex_Name[s], comp_regions$Region_Name[r]))
      
    )
  } # end s loop
  
  osas <-  SPoRC::get_osa(obs_mat = array(comp_fits$Obs_FishLen_mat[r,,,,] + 1e-4, dim = c(1, dim(comp_fits$Obs_FishLen_mat)[-1])),
                          exp_mat = array(comp_fits$Obs_FishLen_mat[r,,,,] + 1e-4, dim = c(1, dim(comp_fits$Obs_FishLen_mat)[-1])),
                          N = array(c(na.omit(data$ISS_FishLenComps[r,which(data$UseFishLenComps[r,,1] == 1),1,1] * 
                                                unique(data$Wt_FishLenComps[r,which(data$UseFishLenComps[r,,1] == 1),1,1]))),
                                    dim = c(1, length(which(data$UseFishLenComps[r,,1] == 1)))),
                          years = which(data$UseFishLenComps[r,,1] == 1),
                          fleet = 1,
                          bins = lens,
                          comp_type = 2,                
                          bin_label = "Lengths"
  )
  
  # Get SDNR QQ plot and Bubble Plot
  osas$res$region <- r # rename region for plotting
  osa_res <- SPoRC::plot_resids(osa_results = osas)
  
  # Get Aggregated Plot
  fishlen_agg <- comp_fits$Fishery_Lens %>%
    group_by(Region, Len, Sex, Fleet) %>%
    summarize(obs = mean(obs), pred = mean(pred)) %>%
    filter(Fleet == 1, Region == r) %>%
    ggplot() +
    geom_col(aes(x = Len, y = obs), fill = 'darkgreen', alpha = 0.8) +
    geom_line(aes(x = Len, y = pred), col = 'black', lwd = 1.3) +
    theme_bw(base_size = 15) +
    facet_grid( Region~Sex, labeller = labeller(
      Region = function(x) paste0("Region ", x),
      Sex = function(x) paste0("Sex ", x)
    )) +
    labs(x = 'Length', y = 'Proportions')
  
  print(
    cowplot::plot_grid(osa_res[[2]], # Bubble
                       cowplot::plot_grid(osa_res[[1]], fishlen_agg, nrow = 1), # QQ and Aggregated
                       ncol = 1, rel_heights = c(0.6, 0.4), labels = paste("LLF OSA Length Comps", comp_regions$Region_Name[r]))
  )
} # end r loop
dev.off()

#### Trawl ---------------------------------------------------------------------
pdf(here("Mods", "Spat_Cont", "Spat_3_Reg", "figs", 'twlf_lcomp_fits.pdf'), width = 19, height = 10) 
for(r in 1:nrow(comp_regions)) {
  for(s in 1:nrow(sex_names)) {
    
    # yearly comp fits
    print(
      comp_fits$Fishery_Lens %>%
        filter(Fleet == 2, Sex == s, Region == r) %>%
        ggplot() +
        geom_col(aes(x = Len, y = obs), alpha = 0.5) +
        geom_line(aes(x = Len, y = pred)) +
        facet_wrap(~Year) +
        theme_sablefish() +
        labs(x = 'Age', y = 'Proportion')+
        labs(title=paste("Trawl Length Comps", sex_names$Sex_Name[s], comp_regions$Region_Name[r]))
      
    )
  } # end s loop
  
  osas <-  SPoRC::get_osa(obs_mat = array(comp_fits$Obs_FishLen_mat[r,,,,] + 1e-4, dim = c(1, dim(comp_fits$Obs_FishLen_mat)[-1])),
                          exp_mat = array(comp_fits$Obs_FishLen_mat[r,,,,] + 1e-4, dim = c(1, dim(comp_fits$Obs_FishLen_mat)[-1])),
                          N = array(c(na.omit(data$ISS_FishLenComps[r,which(data$UseFishLenComps[r,,2] == 1),1,1] * 
                                                unique(data$Wt_FishLenComps[r,which(data$UseFishLenComps[r,,2] == 1),1,2]))),
                                    dim = c(1, length(which(data$UseFishLenComps[r,,2] == 1)))),
                          years = which(data$UseFishLenComps[r,,2] == 1),
                          fleet = 2,
                          bins = lens,
                          comp_type = 2,                
                          bin_label = "Lengths"
  )
  
  # Get SDNR QQ plot and Bubble Plot
  osas$res$region <- r # rename region for plotting
  osa_res <- SPoRC::plot_resids(osa_results = osas)
  
  # Get Aggregated Plot
  fishlen_agg <- comp_fits$Fishery_Lens %>%
    group_by(Region, Len, Sex, Fleet) %>%
    summarize(obs = mean(obs), pred = mean(pred)) %>%
    filter(Fleet == 2, Region == r) %>%
    ggplot() +
    geom_col(aes(x = Len, y = obs), fill = 'darkgreen', alpha = 0.8) +
    geom_line(aes(x = Len, y = pred), col = 'black', lwd = 1.3) +
    theme_bw(base_size = 15) +
    facet_grid( Region~Sex, labeller = labeller(
      Region = function(x) paste0("Region ", x),
      Sex = function(x) paste0("Sex ", x)
    )) +
    labs(x = 'Length', y = 'Proportions')
  
  print(
    cowplot::plot_grid(osa_res[[2]], # Bubble
                       cowplot::plot_grid(osa_res[[1]], fishlen_agg, nrow = 1), # QQ and Aggregated
                       ncol = 1, rel_heights = c(0.6, 0.4), labels = paste("Trawl OSA Length Comps", comp_regions$Region_Name[r]))
  )
} # end r loop
dev.off()

pdf(here("Mods", "Spat_Cont", "Spat_3_Reg", "figs", 'domll_acomp_fits.pdf'), width = 19, height = 10) 
for(r in 1:nrow(comp_regions)) {
  for(s in 1:nrow(sex_names)) {
    
    # yearly comp fits
    print(
      comp_fits$Survey_Ages %>%
        filter(Fleet == 1, Sex == s, Region == r) %>%
        mutate(Cohort = Year - Age) %>%
        ggplot() +
        geom_col(aes(x = Age, y = obs), alpha = 0.5) +
        geom_line(aes(x = Age, y = pred)) +
        facet_wrap(~Year) +
        theme_sablefish() +
        labs(x = 'Age', y = 'Proportion')+
        labs(title=paste("domll Age Comps", sex_names$Sex_Name[s], comp_regions$Region_Name[r]))
      
    )
  } # end s loop
  
  osas <-  SPoRC::get_osa(obs_mat = array(comp_fits$Obs_SrvAge_mat[r,,,,] + 1e-4, dim = c(1, dim(comp_fits$Obs_SrvAge_mat)[-1])),
                          exp_mat = array(comp_fits$Pred_SrvAge_mat[r,,,,] + 1e-4, dim = c(1, dim(comp_fits$Obs_SrvAge_mat)[-1])),
                          N = array(c(na.omit(data$ISS_SrvAgeComps[r,which(data$UseSrvAgeComps[r,,1] == 1),1,1] * 
                                                unique(data$Wt_SrvAgeComps[r,which(data$UseSrvAgeComps[r,,1] == 1),1,1]))),
                                    dim = c(1, length(which(data$UseSrvAgeComps[r,,1] == 1)))),
                          years = which(data$UseSrvAgeComps[r,,1] == 1),
                          fleet = 1,
                          bins = ages,
                          comp_type = 2,                
                          bin_label = "Ages"
  )
  
  # Get SDNR QQ plot and Bubble Plot
  osas$res$region <- r # rename region for plotting
  osa_res <- SPoRC::plot_resids(osa_results = osas)
  
  # Get Aggregated Plot
  srvage_agg <- comp_fits$Survey_Ages %>%
    group_by(Region, Age, Sex, Fleet) %>%
    summarize(obs = mean(obs), pred = mean(pred)) %>%
    filter(Fleet == 1, Region == r) %>%
    ggplot() +
    geom_col(aes(x = Age, y = obs), fill = 'darkgreen', alpha = 0.8) +
    geom_line(aes(x = Age, y = pred), col = 'black', lwd = 1.3) +
    theme_bw(base_size = 15) +
    facet_grid( Region~Sex, labeller = labeller(
      Region = function(x) paste0("Region ", x),
      Sex = function(x) paste0("Sex ", x)
    )) +
    labs(x = 'Age', y = 'Proportions')
  
  print(
    cowplot::plot_grid(osa_res[[2]], # Bubble
                       cowplot::plot_grid(osa_res[[1]], srvage_agg, nrow = 1), # QQ and Aggregated
                       ncol = 1, rel_heights = c(0.6, 0.4), labels = paste("domll OSA Age Comps", comp_regions$Region_Name[r]))
  )
} # end r loop
dev.off()

pdf(here("Mods", "Spat_Cont", "Spat_3_Reg", "figs", 'jpll_acomp_fits.pdf'), width = 19, height = 10) 
for(r in 1:nrow(comp_regions)) {
  for(s in 1:nrow(sex_names)) {
    
    # yearly comp fits
    print(
      comp_fits$Survey_Ages %>%
        filter(Fleet == 3, Sex == s, Region == r) %>%
        mutate(Cohort = Year - Age) %>%
        ggplot() +
        geom_col(aes(x = Age, y = obs), alpha = 0.5) +
        geom_line(aes(x = Age, y = pred)) +
        facet_wrap(~Year) +
        theme_sablefish() +
        labs(x = 'Age', y = 'Proportion')+
        labs(title=paste("jpll Age Comps", sex_names$Sex_Name[s], comp_regions$Region_Name[r]))
      
    )
  } # end s loop
  
  osas <-  SPoRC::get_osa(obs_mat = array(comp_fits$Obs_SrvAge_mat[r,,,,] + 1e-4, dim = c(1, dim(comp_fits$Obs_SrvAge_mat)[-1])),
                          exp_mat = array(comp_fits$Pred_SrvAge_mat[r,,,,] + 1e-4, dim = c(1, dim(comp_fits$Obs_SrvAge_mat)[-1])),
                          N = array(c(na.omit(data$ISS_SrvAgeComps[r,which(data$UseSrvAgeComps[r,,3] == 1),1,3] * 
                                                unique(data$Wt_SrvAgeComps[r,which(data$UseSrvAgeComps[r,,3] == 1),1,3]))),
                                    dim = c(1, length(which(data$UseSrvAgeComps[r,,3] == 1)))),
                          years = which(data$UseSrvAgeComps[r,,3] == 1),
                          fleet = 3,
                          bins = ages,
                          comp_type = 2,                
                          bin_label = "Ages"
  )
  
  # Get SDNR QQ plot and Bubble Plot
  osas$res$region <- r # rename region for plotting
  osa_res <- SPoRC::plot_resids(osa_results = osas)
  
  # Get Aggregated Plot
  srvage_agg <- comp_fits$Survey_Ages %>%
    group_by(Region, Age, Sex, Fleet) %>%
    summarize(obs = mean(obs), pred = mean(pred)) %>%
    filter(Fleet == 3, Region == r) %>%
    ggplot() +
    geom_col(aes(x = Age, y = obs), fill = 'darkgreen', alpha = 0.8) +
    geom_line(aes(x = Age, y = pred), col = 'black', lwd = 1.3) +
    theme_bw(base_size = 15) +
    facet_grid( Region~Sex, labeller = labeller(
      Region = function(x) paste0("Region ", x),
      Sex = function(x) paste0("Sex ", x)
    )) +
    labs(x = 'Age', y = 'Proportions')
  
  print(
    cowplot::plot_grid(osa_res[[2]], # Bubble
                       cowplot::plot_grid(osa_res[[1]], srvage_agg, nrow = 1), # QQ and Aggregated
                       ncol = 1, rel_heights = c(0.6, 0.4), labels = paste("jpll OSA Age Comps", comp_regions$Region_Name[r]))
  )
} # end r loop
dev.off()

# Retrospective -----------------------------------------------------------
sabie_rtmb_model <- readRDS(here("Mods", "Spat_Cont", "Spat_3_Reg", "3_rg_MltRel.RDS"))

retro <- SPoRC::do_retrospective(n_retro = 7,                                 
                                 data = sabie_rtmb_model$data,                
                                 parameters = sabie_rtmb_model$parameters,    
                                 mapping = sabie_rtmb_model$mapping,          
                                 random = NULL,                               
                                 do_par = TRUE,                               
                                 n_cores = 7,                                 
                                 do_francis = FALSE                           
)

retro_plot <- get_retrospective_plot(retro, Rec_Age = 2)   # get retrospective plot

pdf(here("Mods", "Spat_Cont", "Spat_3_Reg", "figs", 'retro.pdf'), width = 13, height = 10) 
retro_plot
dev.off()
