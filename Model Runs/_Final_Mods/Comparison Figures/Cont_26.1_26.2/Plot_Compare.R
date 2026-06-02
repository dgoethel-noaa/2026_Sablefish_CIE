# Purpose: To demonstrate the use of different plotting functions from SPoRC
# Creator: Matthew LH. Cheng (UAF-CFOS)
# Date: 6/9/25


# Set up ------------------------------------------------------------------

library(here)
library(ggplot2)
library(SPoRC)

mod_name2 <- c(#"25.12_1Reg_Cont",
               "26.1i_1Reg_Cont_Upd_No_wAI_LLS", 
               #"26.2a_FAA_3Flt_LLS_Idx_Only",
               "26.2b_FAA_3Flt_LLS_Only"
               #"26.2c_FAA_3Flt",
               #"26.2d_FAA_3Flt_LLS_Only_log2domeS",
               #"26.2e_FAA_3Flt_log2domeS",
               #"26.2f_FAA_3Flt_No_Blks_F"
               #"26.2g_FAA_3Flt_log2domeS_No_Blks_F"
               #"26.2h_FAA_3Flt_LLS_Only_2x_CV_LLS"
               )                         # names of model want to compare, these must match the names used for the SPoRC outputs results rdata files

drop_cont <- c(-1)
mods_no_cont <- mod_name2[-c(1)]
comp_name <- "Final_Mods"                                                       # name for saving files

path <- here("Mods","___Final_Mods")                   # Assumes that comparison folder is inside the first model folder, but can update this to be elsewhere

do_high_res_plots <- 1

# extract out stuff------------------------------------------------------------------

info <- readRDS(here(path,"SPoRC_FAA_3FG_1TF.rds"))                            # read in data file to get years
mod_start_year <- info$yrs[1]
mod_years <- info$yrs


rtmb_out <- list()
data <- list()
rep <- list()
sd_rep <- list()

for(i in 1:length(mod_name2)){
  rtmb_out[[i]] <- readRDS(here(path, paste0(mod_name2[i],"_model_results.RDS")))     # read in base assessment created from SPoRC model run R file
  data[[i]] <- rtmb_out[[i]]$data
  rep[[i]] <- rtmb_out[[i]]$rep
  sd_rep[[i]] <- rtmb_out[[i]]$sd_rep
}

mod_name2 <- c("26.1_1Reg","26.2_FAA")

# Make Plots ------------------------------------------------------------------

# All the plots below are a list of lists - where each element of the list is a different model (allows for model comparison)

# get index fits ------------------------------------------------------------------
# get index fits ------------------------------------------------------------------
idx_fit <- get_idx_fits_plot(data , rep, mod_name2)

idx_plot <- idx_fit+
  ggthemes::scale_color_colorblind() +
  ggthemes::scale_fill_colorblind()+
  ggplot2::scale_x_continuous(breaks = seq(1, length(mod_years), 5),
                              labels = (mod_start_year  - 1)+ seq(1, length(mod_years), 5))+
  facet_wrap(~Category, scales = "free_y", labeller =
               labeller(Category = c(`Fishery1, Q1` = "JPN Fishery CPUE",
                                     `Fishery1, Q2` = "Fishery CPUE",
                                     `Fishery1, Q3` = "Fishery CPUE",
                                     `Survey1, Q1` = "BS LL Survey RPNs",
                                     `Survey2, Q1` = "AI LL Survey RPNs",
                                     `Survey3, Q1` = "GOA LL Survey RPNs",
                                     `Survey4, Q1` = "GOA Trawl Survey (kt)",
                                     `Survey5, Q1` = "JPN LL Survey RPNs")))+
  labs(title = "Index Fits")

idx_fit2 <- get_idx_fits_plot(data[drop_cont] , rep[drop_cont], mods_no_cont)

idx_plot2 <- idx_fit2+
  ggthemes::scale_color_colorblind() +
  ggthemes::scale_fill_colorblind()+
  ggplot2::scale_x_continuous(breaks = seq(1, length(mod_years), 5),
                              labels = (mod_start_year  - 1)+ seq(1, length(mod_years), 5))+
  facet_wrap(~Category, scales = "free_y", labeller =
               labeller(Category = c(`Fishery1, Q1` = "JPN Fishery CPUE",
                                     `Fishery1, Q2` = "Fishery CPUE",
                                     `Fishery1, Q3` = "Fishery CPUE",
                                     `Survey1, Q1` = "BS LL Survey RPNs",
                                     `Survey2, Q1` = "AI LL Survey RPNs",
                                     `Survey3, Q1` = "GOA LL Survey RPNs",
                                     `Survey4, Q1` = "GOA Trawl Survey (kt)",
                                     `Survey5, Q1` = "JPN LL Survey RPNs")))+
  labs(title = "Index Fits")

# get biologicals ------------------------------------------------------------------
bios <- get_biological_plot(data, rep, mod_name2)

M_plot <- bios[[2]]+
  ggthemes::scale_color_colorblind() +
  ggthemes::scale_fill_colorblind()+
  labs(title = "Natural Mortality")

WAA_plot <- bios[[3]]+
  ggthemes::scale_color_colorblind() +
  ggthemes::scale_fill_colorblind()+
  labs(title = "Weight-at-Age")

MAA_plot <- bios[[4]]+
  ggthemes::scale_color_colorblind() +
  ggthemes::scale_fill_colorblind()+
  labs(title = "Maturity")

# get selectivity------------------------------------------------------------------
selex <- get_selex_plot(rep, mod_name2,Selex_Type = "age",year_indx = c(1:length(mod_years)))

f_selex <- selex[[1]]+
  labs(title = "Fishery Selectivity")

s_selex <- selex[[2]]+
  labs(title = "Survey Selectivity")


selex2 <- get_selex_plot(rep[drop_cont], mods_no_cont,Selex_Type = "age",year_indx = c(1:length(mod_years)))

f_selex2 <- selex2[[1]]+
  labs(title = "Fishery Selectivity")

s_selex2 <- selex2[[2]]+
  labs(title = "Survey Selectivity")

# fits to catch --------------------------------------------------------------

catch <- get_catch_fits_plot(data, rep, mod_name2)

catch_plot <- catch + labs(title = "Catch")

catch2 <- get_catch_fits_plot(data[drop_cont], rep[drop_cont],mods_no_cont )

catch_plot2 <- catch2 + labs(title = "Catch")


# get time series plot ------------------------------------------------------------------
time_series <- get_ts_plot(rep, sd_rep, mod_name2,do_ci = TRUE)

ts_plot_CI <- time_series[[1]] +
  ggthemes::scale_color_colorblind() +
  ggthemes::scale_fill_colorblind()+
  ggplot2::scale_x_continuous(breaks = seq(1, length(mod_years), 5),
                              labels = (mod_start_year  - 1)+ seq(1, length(mod_years), 5))+ labs(title = "Time Series Plots")

F_plot_CI <- time_series[[2]]+  facet_wrap(~Type) +
  ggthemes::scale_color_colorblind() +
  ggthemes::scale_fill_colorblind()+
  ggplot2::scale_x_continuous(breaks = seq(1, length(mod_years), 5),
                              labels = (mod_start_year  - 1)+ seq(1, length(mod_years), 5))+ labs(title = "Fishing Mortality")

recr_CI <-      time_series[[3]]+  facet_wrap(~Type) +
  ggthemes::scale_color_colorblind() +
  ggthemes::scale_fill_colorblind()+
  ggplot2::scale_x_continuous(breaks = seq(1, length(mod_years), 5),
                              labels = (mod_start_year  - 1)+ seq(1, length(mod_years), 5))+ labs(title = "Recruitment")

ssb_CI <-      time_series[[4]]+  facet_wrap(~Type) +
  ggthemes::scale_color_colorblind() +
  ggthemes::scale_fill_colorblind()+
  ggplot2::scale_x_continuous(breaks = seq(1, length(mod_years), 5),
                              labels = (mod_start_year  - 1)+ seq(1, length(mod_years), 5))+ labs(title = "SSB")

bio_CI <-      time_series[[5]]+  facet_wrap(~Type) +
  ggthemes::scale_color_colorblind() +
  ggthemes::scale_fill_colorblind()+
  ggplot2::scale_x_continuous(breaks = seq(1, length(mod_years), 5),
                              labels = (mod_start_year  - 1)+ seq(1, length(mod_years), 5))+ labs(title = "Biomass")


# REPEAT PLOTS WITHOUT CIs
time_series <- get_ts_plot(rep, sd_rep, mod_name2,do_ci = FALSE)

ts_plot <-      time_series[[1]] +
  ggthemes::scale_color_colorblind() +
  ggthemes::scale_fill_colorblind()+
  ggplot2::scale_x_continuous(breaks = seq(1, length(mod_years), 5),
                              labels = (mod_start_year  - 1)+ seq(1, length(mod_years), 5)) + labs(title = "Time Series Plots")

F_plot <- time_series[[2]]+  facet_wrap(~Type) +
  ggthemes::scale_color_colorblind() +
  ggthemes::scale_fill_colorblind()+
  ggplot2::scale_x_continuous(breaks = seq(1, length(mod_years), 5),
                              labels = (mod_start_year  - 1)+ seq(1, length(mod_years), 5))+ labs(title = "Fishing Mortality")

recr <-      time_series[[3]]+  facet_wrap(~Type) +
  ggthemes::scale_color_colorblind() +
  ggthemes::scale_fill_colorblind()+
  ggplot2::scale_x_continuous(breaks = seq(1, length(mod_years), 5),
                              labels = (mod_start_year  - 1)+ seq(1, length(mod_years), 5))+ labs(title = "Recruitment")

ssb <-     time_series[[4]]+  facet_wrap(~Type) +
  ggthemes::scale_color_colorblind() +
  ggthemes::scale_fill_colorblind()+
  ggplot2::scale_x_continuous(breaks = seq(1, length(mod_years), 5),
                              labels = (mod_start_year  - 1)+ seq(1, length(mod_years), 5))+ labs(title = "SSB")

bio <-      time_series[[5]]+  facet_wrap(~Type) +
  ggthemes::scale_color_colorblind() +
  ggthemes::scale_fill_colorblind()+
  ggplot2::scale_x_continuous(breaks = seq(1, length(mod_years), 5),
                              labels = (mod_start_year  - 1)+ seq(1, length(mod_years), 5))+ labs(title = "Biomass")


# get data fitted plot ------------------------------------------------------------------
data_used <- get_data_fitted_plot(data, mod_name2)

data_plot <- data_used+
  ggplot2::scale_x_continuous(breaks = seq(1, length(mod_years), 5),
                              labels = (mod_start_year  - 1)+ seq(1, length(mod_years), 5))+ labs(title = "Data")

# get nLL plot ------------------------------------------------------------------
nll <- get_nLL_plot(data, rep, mod_name2)

nll_plot <- nll[[1]]+ labs(title = "Negative Log-Likelihood")


# Add table of terminal year values (SSB, stock status, ABC) and table of likelihood components, convergence status, max grad, etc.

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
                       HCR_function = HCR_function <- function(x, frp, brp, alpha = 0.05) {
                         stock_status <- x / brp # define stock status
                         # If stock status is > 1
                         if(stock_status >= 1) f <- frp
                         # If stock status is between brp and alpha
                         if(stock_status > alpha && stock_status < 1) f <- frp * (stock_status - alpha) / (1 - alpha)
                         # If stock status is less than alpha
                         if(stock_status < alpha) f <- 0
                         return(f)
                       },
                       recruitment_opt = 'mean_rec',
                       fmort_opt = 'HCR')
out <- get_key_quants(data, rep, reference_points_opt, proj_model_opt, mod_name2)

print(out[[1]])  # key quantities data.frame

tiff(filename=here(path,paste0(comp_name,"_status+BRPs.png")),units = "in", width=16,height=4, res = 300)
print(out[[2]])  # table plot as cowplot object
dev.off()


# plot all basic info (basically just a wrapper function for a subset of functions defined above) ------------------------------------------------------------------
all_plots <- plot_all_basic(data, rep, sd_rep, mod_name2, out_path = here(path))


pdf(file = here(path,paste0(comp_name,"_Full Comparison.pdf")),  width = 25, height = 13)

print(out[[2]])
print(nll[[2]])
print(nll_plot)

print(data_plot)

print(M_plot)
print(WAA_plot)
print(MAA_plot)

print(idx_plot)
print(idx_plot2)
print(catch_plot)
print(catch_plot2)
print(f_selex)
print(f_selex2)
print(s_selex)
print(s_selex2)
print(ts_plot)
print(F_plot)
print(recr)
print(ssb)
print(bio)
print(ts_plot_CI)
print(F_plot_CI)
print(recr_CI)
print(ssb_CI)
print(bio_CI)

dev.off()


if(do_high_res_plots == 1){
  tiff(filename=here(path,paste0(comp_name,"_nll.png")),units = "in", width=12,height=10, res = 300)
  print(nll_plot)
  dev.off()
  tiff(filename=here(path,paste0(comp_name,"_nll_tab.png")),units = "in", width=24,height=2, res = 300)
  print(nll[[2]])
  dev.off()
  tiff(filename=here(path,paste0(comp_name,"_data.png")),units = "in", width=14,height=12, res = 300)
  print(data_plot)
  dev.off()


  tiff(filename=here(path,paste0(comp_name,"_M.png")),units = "in", width=12,height=8, res = 300)
  print(M_plot)
  dev.off()
  tiff(filename=here(path,paste0(comp_name,"_WAA.png")),units = "in", width=12,height=8, res = 300)
  print(WAA_plot)
  dev.off()
  tiff(filename=here(path,paste0(comp_name,"_Mat.png")),units = "in", width=12,height=8, res = 300)
  print(MAA_plot)
  dev.off()

  tiff(filename=here(path,paste0(comp_name,"_idx_fit.png")),units = "in", width=24,height=14, res = 300)
  print(idx_plot)
  dev.off()
  tiff(filename=here(path,paste0(comp_name,"_selex_fish.png")),units = "in", width=12,height=10, res = 300)
  print(f_selex)
  dev.off()
  tiff(filename=here(path,paste0(comp_name,"_selex_srvy.png")),units = "in", width=12,height=10, res = 300)
  print(s_selex)
  dev.off()
  tiff(filename=here(path,paste0(comp_name,"_llf_catch_fit.png")),units = "in", width=14,height=12, res = 300)
  print(catch_plot)
  dev.off()


  tiff(filename=here(path,paste0(comp_name,"_ts_noCI.png")),units = "in", width=14,height=12, res = 300)
  print(ts_plot)
  dev.off()
  tiff(filename=here(path,paste0(comp_name,"_ts_F_noCI.png")),units = "in", width=24,height=6, res = 300)
  print(F_plot)
  dev.off()
  tiff(filename=here(path,paste0(comp_name,"_ts_Recr_noCI.png")),units = "in", width=22,height=10, res = 300)
  print(recr)
  dev.off()
  tiff(filename=here(path,paste0(comp_name,"_ts_SSB_noCI.png")),units = "in", width=22,height=10, res = 300)
  print(ssb)
  dev.off()
  tiff(filename=here(path,paste0(comp_name,"_ts_Bio_noCI.png")),units = "in", width=22,height=10, res = 300)
  print(bio)
  dev.off()

  tiff(filename=here(path,paste0(comp_name,"_ts.png")),units = "in", width=14,height=12, res = 300)
  print(ts_plot_CI)
  dev.off()
  tiff(filename=here(path,paste0(comp_name,"_ts_F.png")),units = "in", width=24,height=6, res = 300)
  print(F_plot_CI)
  dev.off()
  tiff(filename=here(path,paste0(comp_name,"_ts_Recr.png")),units = "in", width=22,height=10, res = 300)
  print(recr_CI)
  dev.off()
  tiff(filename=here(path,paste0(comp_name,"_ts_SSB.png")),units = "in", width=22,height=10, res = 300)
  print(ssb_CI)
  dev.off()
  tiff(filename=here(path,paste0(comp_name,"_ts_Bio.png")),units = "in", width=22,height=10, res = 300)
  print(bio_CI)
  dev.off()
}

