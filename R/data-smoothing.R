# description -------------------------------------------------------------

# This code is used to smooth the raw data from the monitoring experiments. 
# The raw data is published as the bcsa R package. The data is smoothed using
# Centered Moving Average of order 5. The code also counts the number of negative
# values before and after smoothing. 

# r packages --------------------------------------------------------------

library(bcsa)
library(dplyr)
library(tidyverse)

# combine the monitoring data used in this manuscript ---------------------

df_coll <- bcsa::df_collocation |> 
  mutate(exp_type = "sensor_collocation")

df_mm_aae_paper <- df_mm_road_type |> 
  filter(time_of_day == "Morning",
         settlement_id %in% c("Ndirande", "Namiwawa", "Nyambadwe"))
  
df_pm_aae_paper <- df_pm |>
  filter(settlement_id == "Ndirande")

df_main <- df_aae |>
  bind_rows(df_mm_aae_paper) |> 
  bind_rows(df_pm_aae_paper) |>
  bind_rows(df_sm) |>
  bind_rows(df_coll) |>
  drop_na(ir_bcc, blue_bcc, uv_bcc) 

# data smoothing ----------------------------------------------------------

## count number of negative values in each experiment in raw data

df_negative_raw <- df_main |>
  group_by(id, exp_type, emission_source) |>
  summarise(neg_irbcc = sum(ir_bcc < 0),
            neg_bluebcc = sum(blue_bcc < 0),
            neg_uvbcc = sum(uv_bcc < 0),
            count = n()) |> 
  mutate(neg_irbcc_per = neg_irbcc/count*100,
         neg_bluebcc_per = neg_bluebcc/count*100,
         neg_uvbcc_per = neg_uvbcc/count*100) 

### negative values - smoothing needed

## apply cma of order 3, 5 and 7 and create a temporary dataset

df_smooth_raw <- df_main

list_df_raw <- list()

for(i in (unique(df_smooth_raw$id))){

  list_df_raw[[i]] <- df_smooth_raw |>
    filter(id == i)

  list_df_raw[[i]] <- list_df_raw[[i]] |>
    mutate(ir_bcc_3 = smooth::cma(list_df_raw[[i]]$ir_bcc, order = 3, silent = TRUE)$fitted,
           blue_bcc_3 = smooth::cma(list_df_raw[[i]]$blue_bcc, order = 3, silent = TRUE)$fitted,
           uv_bcc_3 = smooth::cma(list_df_raw[[i]]$uv_bcc, order = 3, silent = TRUE)$fitted,
           ir_bcc_5 = smooth::cma(list_df_raw[[i]]$ir_bcc, order = 5, silent = TRUE)$fitted,
           blue_bcc_5 = smooth::cma(list_df_raw[[i]]$blue_bcc, order = 5, silent = TRUE)$fitted,
           uv_bcc_5 = smooth::cma(list_df_raw[[i]]$uv_bcc, order = 5, silent = TRUE)$fitted,
           ir_bcc_7 = smooth::cma(list_df_raw[[i]]$ir_bcc, order = 7, silent = TRUE)$fitted,
           blue_bcc_7 = smooth::cma(list_df_raw[[i]]$blue_bcc, order = 7, silent = TRUE)$fitted,
           uv_bcc_7 = smooth::cma(list_df_raw[[i]]$uv_bcc, order = 7, silent = TRUE)$fitted)
}

df_smooth_temp <- bind_rows(list_df_raw)

## count number of negative values after smoothing

df_negative_count <- df_smooth_temp |>
  group_by(exp_type, emission_source) |>
  summarise(neg_irbcc = sum(ir_bcc < 0)/n()*100,
            neg_irbcc_3 = sum(ir_bcc_3 < 0)/n()*100,
            neg_irbcc_5 = sum(ir_bcc_5 < 0)/n()*100,
            neg_irbcc_7 = sum(ir_bcc_7 < 0)/n()*100) |>
  mutate_if(is.numeric,
            round,
            digits = 1) 

### number of negative values decreases with increase in the order
### of smoothing. we select the order of 5, as the decrease in the
### number of negative values between order 5 and 7 is not more than 1%

## smooth df with selected cma order (now 5)

df_smooth <- df_main

list_df <- list()

for (i in (unique(df_smooth$id))) {

  list_df[[i]] <- df_smooth |>
    filter(id == i)

  list_df[[i]] <- list_df[[i]] |>
    mutate(ir_bcc = smooth::cma(list_df[[i]]$ir_bcc, order = 5, silent = TRUE)$fitted,
           blue_bcc = smooth::cma(list_df[[i]]$blue_bcc, order = 5, silent = TRUE)$fitted,
           ir_babs = smooth::cma(list_df[[i]]$ir_babs, order = 5, silent = TRUE)$fitted,
           blue_babs = smooth::cma(list_df[[i]]$blue_babs, order = 5, silent = TRUE)$fitted,
           uv_bcc = smooth::cma(list_df[[i]]$uv_bcc, order = 5, silent = TRUE)$fitted,
           uv_babs = smooth::cma(list_df[[i]]$uv_babs, order = 5, silent = TRUE)$fitted)
}

df_smooth <- bind_rows(list_df)

# form dataframes of each experiment type --------------------------------

df_collocation <- df_smooth |>
  filter(exp_type == "sensor_collocation")

df_sm <- df_smooth |>
  filter(exp_type == "stationary_monitoring")

df_aae <- df_smooth |>
  filter(exp_type %in% c("waste_burning", "cooking", "vehicles"))

df_mm <- df_smooth |>
  filter(exp_type == "mobile_monitoring")

df_pm <- df_smooth |>
  filter(exp_type == "personal_monitoring")
