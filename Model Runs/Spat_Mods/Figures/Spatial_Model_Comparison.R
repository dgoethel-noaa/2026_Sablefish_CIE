# Purpose: To compare 2 region, 3 reigon, and 5 region spatial models
# Creator: Matthew LH. Cheng
# 4/2/26


# Setup -------------------------------------------------------------------

library(here)
library(tidyverse)
library(SPoRC)

# read in models
# two_rg <- readRDS( here("Mods", "Spat_Cont", "Spat_2_Reg", "2_rg_MltRel.RDS"))
three_rg <- readRDS( here("Mods", "Spat_Cont", "Spat_3_Reg", "3_rg_MltRel.RDS"))
five_rg <- readRDS( here("Mods", "Spat_Cont", "Spat_5_Reg", "5_rg_MltRel.RDS"))


# Helpers -----------------------------------------------------------------

# Helper to melt one array and tag it
melt_recap <- function(arr, type_label) {
  expand.grid(
    liberty      = seq_len(dim(arr)[1]),
    cohort       = seq_len(dim(arr)[2]),
    recap_region = seq_len(dim(arr)[3])
  )  %>% 
    mutate(value = as.vector(arr), type = type_label)
}

# build_ts_df <- function(metric) {
#   
#   n_yr <- ncol(two_rg$rep[[metric]])
#   yrs  <- seq_len(n_yr)
#   
#   data.frame(
#     year   = rep(yrs, 3),
#     model  = rep(c("2-region", "3-region", "5-region"), each = n_yr),
#     bsai   = c(
#       two_rg$rep[[metric]][1, ],
#       colSums(three_rg$rep[[metric]][1:2, ]),
#       colSums(five_rg$rep[[metric]][1:2, ])
#     ),
#     goa    = c(
#       two_rg$rep[[metric]][2, ],
#       three_rg$rep[[metric]][3, ],
#       colSums(five_rg$rep[[metric]][3:5, ])
#     ),
#     total  = c(
#       colSums(two_rg$rep[[metric]]),
#       colSums(three_rg$rep[[metric]]),
#       colSums(five_rg$rep[[metric]])
#     )
#   )  %>%
#     pivot_longer(c(bsai, goa, total),
#                  names_to  = "region",
#                  values_to = "value")  %>%
#     mutate(
#       region = factor(region,
#                       levels = c("bsai", "goa", "total"),
#                       labels = c("BS + AI", "GOA", "Aggregated")),
#       model  = factor(model, levels = c("2-region", "3-region", "5-region"))
#     )
# }

build_ts_df_3 <- function(metric) {
  
  n_yr <- ncol(three_rg$rep[[metric]])
  yrs  <- seq_len(n_yr)
  
  data.frame(
    year   = rep(yrs, 2),
    model  = rep(c("3-region", "5-region"), each = n_yr),
    bs   = c(
      three_rg$rep[[metric]][1, ],
      five_rg$rep[[metric]][1, ]
    ),
    ai   = c(
      three_rg$rep[[metric]][2, ],
      five_rg$rep[[metric]][2, ]
    ),
    goa    = c(
      three_rg$rep[[metric]][3, ],
      colSums(five_rg$rep[[metric]][3:5, ])
    ),
    total  = c(
      colSums(three_rg$rep[[metric]]),
      colSums(five_rg$rep[[metric]])
    )
  )  %>%
    pivot_longer(c(bs, ai, goa, total),
                 names_to  = "region",
                 values_to = "value")  %>%
    mutate(
      region = factor(region,
                      levels = c("bs", 'ai', "goa", "total"),
                      labels = c("BS", 'AI', "GOA", "Aggregated")),
      model  = factor(model, levels = c("3-region", "5-region"))
    )
}

build_regional_ts <- function(model_rep, region_labels, model_name, metric) {
  n_yr <- ncol(model_rep[[metric]])
  n_rg <- length(region_labels)
  
  data.frame(
    year   = rep(seq_len(n_yr), times = n_rg),
    model  = model_name,
    region = rep(region_labels, each = n_yr),
    value  = as.vector(t(model_rep[[metric]][seq_len(n_rg), ]))
  ) %>% 
    mutate(region = factor(region, levels = region_labels))
}


# Time Series -------------------------------------------------------------

# Compare aggregated time series (to the lowest resolution spatial model)
# ssb_df <- build_ts_df("SSB")  %>% mutate(metric = "SSB")
# rec_df <- build_ts_df("Rec")  %>% mutate(metric = "Recruitment")
# ts_df     <- bind_rows(ssb_df, rec_df) %>% 
#   dplyr::mutate(year = year + 1959)
# 
# ggplot(ts_df, 
#        aes(x = year, y = value, color = model, fill = model, lty = model)) +
#   geom_col(data = ~filter(.x, metric == "Recruitment"), alpha = 0.7,
#            position = position_dodge(), color = NA) +
#   geom_line(data = ~filter(.x, metric == "SSB"), linewidth = 1.3) +
#   facet_grid(metric ~ region, scales = "free_y") +
#   scale_x_continuous(breaks = scales::pretty_breaks(5)) +
#   scale_y_continuous(breaks = scales::pretty_breaks(4)) +
#   labs(x = 'Year', y = 'Value', lty = 'Model', color = 'Model', fill = 'Model') +
#   theme_bw(base_size = 20)
# 
# ggplot(ts_df, 
#        aes(x = year, y = value, color = model, lty = model)) +
#   geom_line() +
#   facet_grid(metric ~ region, scales = "free_y") +
#   scale_x_continuous(breaks = scales::pretty_breaks(5)) +
#   scale_y_continuous(breaks = scales::pretty_breaks(4)) +
#   labs(x = 'Year', y = 'Value', lty = 'Model', color = 'Model') +
#   theme_bw(base_size = 20)

# Compare aggregated time series (to the 3 region spatial model)
ssb_df <- build_ts_df_3("SSB")  %>% mutate(metric = "SSB")
rec_df <- build_ts_df_3("Rec")  %>% mutate(metric = "Recruitment")
ts_df     <- bind_rows(ssb_df, rec_df) %>% 
  dplyr::mutate(year = year + 1959)

ggplot(ts_df, 
       aes(x = year, y = value, color = model, fill = model, lty = model)) +
  geom_col(data = ~filter(.x, metric == "Recruitment"), alpha = 0.7,
           position = position_dodge(), color = NA) +
  geom_line(data = ~filter(.x, metric == "SSB"), linewidth = 1.3) +
  facet_grid(metric ~ region, scales = "free_y") +
  scale_x_continuous(breaks = scales::pretty_breaks(5)) +
  scale_y_continuous(breaks = scales::pretty_breaks(4)) +
  labs(x = 'Year', y = 'Value', lty = 'Model', color = 'Model', fill = 'Model') +
  theme_bw(base_size = 20)

pdf(here("Mods", "Spat_Cont", "ts_all.pdf"), width = 15, height = 10)
ggplot(ts_df, 
       aes(x = year, y = value, color = model, lty = model)) +
  geom_line() +
  facet_grid(metric ~ region, scales = "free_y") +
  scale_x_continuous(breaks = scales::pretty_breaks(5)) +
  scale_y_continuous(breaks = scales::pretty_breaks(4)) +
  labs(x = 'Year', y = 'Value', lty = 'Model', color = 'Model') +
  theme_bw(base_size = 20)
dev.off()

# look at three region model estiamtes
three_rg_labels <- c("BS", "AI", "GOA")
three_rg_df <- bind_rows(
  build_regional_ts(three_rg$rep, three_rg_labels, "3-region", "SSB"),
  build_regional_ts(three_rg$rep, three_rg_labels, "3-region", "Rec") %>%  mutate(metric = "Recruitment")
) %>%  
  mutate(metric = if_else(is.na(metric), "SSB", metric),
         year = year + 1959)

ggplot(three_rg_df, aes(x = year, y = value)) +
  geom_col(data = ~filter(.x, metric == "Recruitment"), color = 'black', 
           fill = "steelblue", alpha = 0.7, lwd = 0.5) +
  geom_line(data = ~filter(.x, metric == "SSB"), linewidth = 1.3) +
  facet_grid(metric ~ region, scales = "free_y") +
  scale_x_continuous(breaks = scales::pretty_breaks(5)) +
  scale_y_continuous(breaks = scales::pretty_breaks(4)) +
  labs(x = 'Year', y = 'Value') +
  theme_bw(base_size = 20)

pdf(here("Mods", "Spat_Cont", "ts_3.pdf"), width = 15, height = 10)
ggplot(three_rg_df, aes(x = year, y = value)) +
  geom_line(linewidth = 1.3) +
  facet_grid(metric ~ region, scales = "free_y") +
  ggthemes::scale_color_colorblind() +
  scale_x_continuous(breaks = scales::pretty_breaks(5)) +
  scale_y_continuous(breaks = scales::pretty_breaks(4)) +
  labs(x = 'Year', y = 'Value') +
  theme_bw(base_size = 20) 
dev.off()

# look at five region models
five_rg_labels <- c("BS", "AI", "WGOA", 'CGOA', 'EGOA')
five_rg_df <- bind_rows(
  build_regional_ts(five_rg$rep, five_rg_labels, "5-region", "SSB"),
  build_regional_ts(five_rg$rep, five_rg_labels, "5-region", "Rec") %>%  mutate(metric = "Recruitment")
) %>% 
  mutate(metric = if_else(is.na(metric), "SSB", metric),
         year = year + 1959)

ggplot(five_rg_df, aes(x = year, y = value)) +
  geom_col(data = ~filter(.x, metric == "Recruitment"), color = 'black', 
           fill = "steelblue", alpha = 0.7, lwd = 0.5) +
  geom_line(data = ~filter(.x, metric == "SSB"), linewidth = 1.3) +
  facet_grid(metric ~ region, scales = "free_y") +
  scale_x_continuous(breaks = scales::pretty_breaks(5)) +
  scale_y_continuous(breaks = scales::pretty_breaks(4)) +
  labs(x = 'Year', y = 'Value') +
  theme_bw(base_size = 20)

pdf(here("Mods", "Spat_Cont", "ts_5.pdf"), width = 15, height = 10)
ggplot(five_rg_df, aes(x = year, y = value)) +
  geom_line(linewidth = 1.3) +
  facet_grid(metric ~ region, scales = "free_y") +
  ggthemes::scale_color_colorblind() +
  scale_x_continuous(breaks = scales::pretty_breaks(5)) +
  scale_y_continuous(breaks = scales::pretty_breaks(4)) +
  labs(x = 'Year', y = 'Value') +
  theme_bw(base_size = 20) 
dev.off()

# Index Fits --------------------------------------------------------------

# Two region model
# two_rg$data$ObsSrvIdx[,,2] <- NA # remove trawl survey
# two_rg$data$ObsFishIdx[] <- NA # remove fishery index
# get_idx_fits_plot(list(two_rg$data), list(two_rg$rep), "Two Region") 

# Three region model
three_rg$data$ObsSrvIdx[,,2] <- NA # remove trawl survey
three_rg$data$ObsFishIdx[] <- NA # remove fishery index

pdf(here("Mods", "Spat_Cont", "idxfits_3.pdf"), width = 15, height = 10)
get_idx_fits_plot(list(three_rg$data), list(three_rg$rep), "Three Region") +
  facet_grid(Category ~ Region, scales = "free_y", labeller = labeller(
    Region = c("Region 1" = "BS", "Region 2" = "AI", "Region 3" = "GOA"),
    Category = c("Survey1, Q1" = "US Longline Survey", "Survey3, Q1" = "JP Longline Survey")
  )) +
  scale_x_continuous(labels = seq(1960, 2025, by = 10), breaks = seq(1, 65, 10)) +
  theme(legend.position = 'none')
dev.off()

# Five region model
five_rg$data$ObsSrvIdx[,,2] <- NA # remove trawl survey
five_rg$data$ObsFishIdx[] <- NA # remove fishery index

pdf(here("Mods", "Spat_Cont", "idxfits_5.pdf"), width = 15, height = 10)
get_idx_fits_plot(list(five_rg$data), list(five_rg$rep), "Five Region") +
  facet_grid(Category ~ Region, scales = "free_y", labeller = labeller(
    Region = c("Region 1" = "BS", "Region 2" = "AI", "Region 3" = "WGOA", "Region 4" = "EGOA", "Region 5" = "CGOA"),
    Category = c("Survey1, Q1" = "US Longline Survey", "Survey3, Q1" = "JP Longline Survey")
  )) +
  scale_x_continuous(labels = seq(1960, 2025, by = 10), breaks = seq(1, 65, 10)) +
  theme(legend.position = 'none')
dev.off()

# Catch Fits --------------------------------------------------------------

# Two region model
# get_catch_fits_plot(list(two_rg$data), list(two_rg$rep), "Two Region") 

# Three region model
pdf(here("Mods", "Spat_Cont", "catfits_3.pdf"), width = 15, height = 10)
get_catch_fits_plot(list(three_rg$data), list(three_rg$rep), "Three Region") +
  facet_grid(Fleet ~ Region, scales = "free_y", labeller = labeller(
    Region = c("Region 1" = "BS", "Region 2" = "AI", "Region 3" = "GOA"),
    Category = c("Fleet 1" = "Fixed Gear", "Fleet 2" = "Trawl Gear")
  )) +
  scale_x_continuous(labels = seq(1960, 2025, by = 10), breaks = seq(1, 65, 10)) +
  theme(legend.position = 'none')
dev.off()

# Five region model
pdf(here("Mods", "Spat_Cont", "catfits_5.pdf"), width = 15, height = 10)
get_catch_fits_plot(list(five_rg$data), list(five_rg$rep), "Five Region") +
  facet_grid(Fleet ~ Region, scales = "free_y", labeller = labeller(
    Region = c("Region 1" = "BS", "Region 2" = "AI", "Region 3" = "WGOA", "Region 4" = "EGOA", "Region 5" = "CGOA"),
    Category = c("Fleet 1" = "Fixed Gear", "Fleet 2" = "Trawl Gear")
  )) +
  scale_x_continuous(labels = seq(1960, 2025, by = 10), breaks = seq(1, 65, 10)) +
  theme(legend.position = 'none')
dev.off()


# Tag Fits ----------------------------------------------------------------

# two region model
# obs_collapsed  <- apply(two_rg$data$Obs_Tag_Recap,  1:3, sum)
# pred_collapsed <- apply(two_rg$rep$Pred_Tag_Recap, 1:3, sum)
# 
# # Cohort lookup
# cohort_info <- as.data.frame(two_rg$data$tag_release_indicator)
# colnames(cohort_info) <- c("release_region", "release_yr")
# cohort_info$cohort <- seq_len(nrow(cohort_info))
# 
# # summarize array
# tag_df <- bind_rows(
#   melt_recap(obs_collapsed,  "Observed"),
#   melt_recap(pred_collapsed, "Predicted")
# ) %>% 
#   left_join(cohort_info, by = "cohort") %>% 
#   mutate(recap_yr = release_yr + liberty - 1)  # adjust if needed
# 
# # Summarise by recap year, release region, and obs/pred
# plot_df <- tag_df %>% 
#   group_by(recap_yr, release_region, type) %>% 
#   summarise(total = sum(value), .groups = "drop") %>% 
#   mutate(release_region = factor(release_region),
#          recap_yr = recap_yr + 1959)
# 
# ggplot(plot_df, aes(x = recap_yr)) +
#   geom_col(data = filter(plot_df, type == "Observed"),
#            aes(y = total, fill = release_region), alpha = 0.3) +
#   geom_line(data = filter(plot_df, type == "Predicted"),
#             aes(y = total, color = release_region), linewidth = 0.8) +
#   scale_fill_brewer(palette = "Set1")  +
#   scale_color_brewer(palette = "Set1") +
#   facet_wrap(~release_region, scales = "free_y",
#              labeller = labeller(release_region = c("1" = "BS+AI", "2" = "GOA")) ) +
#   labs(x = "Recapture year", y = "Tag recaptures",
#        fill = "Release region", color = "Release region") +
#   theme_bw(base_size = 20) +
#   theme(legend.position = 'none')

# Three region model
obs_collapsed  <- apply(three_rg$data$Obs_Tag_Recap,  1:3, sum)
pred_collapsed <- apply(three_rg$rep$Pred_Tag_Recap, 1:3, sum)

# Cohort lookup
cohort_info <- as.data.frame(three_rg$data$tag_release_indicator)
colnames(cohort_info) <- c("release_region", "release_yr")
cohort_info$cohort <- seq_len(nrow(cohort_info))

# summarize array
tag_df <- bind_rows(
  melt_recap(obs_collapsed,  "Observed"),
  melt_recap(pred_collapsed, "Predicted")
) %>% 
  left_join(cohort_info, by = "cohort") %>% 
  mutate(recap_yr = release_yr + liberty - 1)  # adjust if needed

# Summarise by recap year, release region, and obs/pred
plot_df <- tag_df %>% 
  group_by(recap_yr, release_region, type) %>% 
  summarise(total = sum(value), .groups = "drop") %>% 
  mutate(release_region = factor(release_region),
         recap_yr = recap_yr + 1959)

pdf(here("Mods", "Spat_Cont", "tagfits_3.pdf"), width = 15, height = 10)
ggplot(plot_df %>% filter(recap_yr < 2025), aes(x = recap_yr)) +
  geom_point(data = filter(plot_df, type == "Observed", recap_yr < 2025), aes(y = total), alpha = 0.3) +
  geom_line(data = filter(plot_df, type == "Predicted", recap_yr < 2025),
            aes(y = total, color = release_region), linewidth = 0.8) +
  scale_fill_brewer(palette = "Set1")  +
  scale_color_brewer(palette = "Set1") +
  facet_wrap(~release_region, scales = "free_y",
             labeller = labeller(release_region = c("1" = "BS", "2" = "AI", "3" = "GOA")) ) +
  labs(x = "Recapture year", y = "Tag recaptures",
       fill = "Release region", color = "Release region") +
  theme_bw(base_size = 20) +
  theme(legend.position = 'none')
dev.off()

# five region model
obs_collapsed  <- apply(five_rg$data$Obs_Tag_Recap,  1:3, sum)
pred_collapsed <- apply(five_rg$rep$Pred_Tag_Recap, 1:3, sum)

# Cohort lookup
cohort_info <- as.data.frame(five_rg$data$tag_release_indicator)
colnames(cohort_info) <- c("release_region", "release_yr")
cohort_info$cohort <- seq_len(nrow(cohort_info))

# summarize array
tag_df <- bind_rows(
  melt_recap(obs_collapsed,  "Observed"),
  melt_recap(pred_collapsed, "Predicted")
) %>% 
  left_join(cohort_info, by = "cohort") %>% 
  mutate(recap_yr = release_yr + liberty - 1)  # adjust if needed

# Summarise by recap year, release region, and obs/pred
plot_df <- tag_df %>% 
  group_by(recap_yr, release_region, type) %>% 
  summarise(total = sum(value), .groups = "drop") %>% 
  mutate(release_region = factor(release_region),
         recap_yr = recap_yr + 1959)

pdf(here("Mods", "Spat_Cont", "catfits_5.pdf"), width = 18, height = 8)
ggplot(plot_df %>% filter(recap_yr < 2025), aes(x = recap_yr)) +
  geom_point(data = filter(plot_df, type == "Observed", recap_yr < 2025), aes(y = total), alpha = 0.3) +
  geom_line(data = filter(plot_df, type == "Predicted", recap_yr < 2025),
            aes(y = total, color = release_region), linewidth = 0.8) +
  scale_fill_brewer(palette = "Set1")  +
  scale_color_brewer(palette = "Set1") +
  facet_wrap(~release_region, scales = "free_y", nrow = 1,
             labeller = labeller(release_region = c("1" = "BS", "2" = "AI", "3" = "WGOA", "4" = 'CGOA', '5' = 'EGOA') ) ) +
  labs(x = "Recapture year", y = "Tag recaptures",
       fill = "Release region", color = "Release region") +
  theme_bw(base_size = 20) +
  theme(legend.position = 'none')
dev.off()
