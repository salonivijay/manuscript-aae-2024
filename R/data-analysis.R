# description -------------------------------------------------------------

# This code provides the code to calculate the absorption angstrom exponent values,
# statistical test to assess similarity and dissimilarity in AAE values from different 
# emission sources, and to verify the use of AAE values in source apportionment.
# The smooth data generated from data smoothing code is used to calculate the 
# absorption angstrom exponent values, and to verify it's use in source apportionment.

# Overall, the code is used to generate the results presented in the manuscript on:
# "Absorption Ångström Exponent values to identify light-absorbing aerosol sources 
# in Blantyre, Malawi."

# read data-smoothing.R file ----------------------------------------------

source(here::here("R/data-smoothing.R"))

# r packages --------------------------------------------------------------

library(dplyr)
library(tidyverse)
library(sf)
library(basemapR)
library(hms)
library(broom)

# read files --------------------------------------------------------------

table_aae_literature <- read_csv(here::here("data/aae-comparison-literature.csv"), 
                                 locale = readr::locale(encoding = "latin1"),
                                 col_types = "cccc")

table_aae_literature[is.na(table_aae_literature)] <- ""  

table_aae_methods <- read_csv(here::here("data/aae-methods.csv"), 
                              locale = readr::locale(encoding = "latin1"),
                              col_types = "cccc")

parameters <- read_csv(here::here("data/MA200-parameters.csv"))

# define a common plot theme ----------------------------------------------

theme <- theme(axis.text.y   = element_text(size=12),
               axis.text.x   = element_text(size=12),
               axis.title.y  = element_text(size=14),
               axis.title.x  = element_text(size=14),
               panel.background = element_rect(fill='transparent'), #transparent panel bg.
               plot.background = element_rect(fill='transparent', color=NA), #transparent plot bg.
               axis.line = element_line(colour = "black"),
               panel.border = element_rect(colour = "black", fill=NA, size=0.5))

# sensor collocation ------------------------------------------------------

df_sc_aae_1 <- df_collocation |>
  filter(id %in% c(124,125)) |> 
  mutate(time = case_when(
    serial_number == "MA200-0420" ~ time - 1,
    TRUE ~ time
  )) |> 
  mutate("Monitoring frequency (s)" = 30)

df_sc_aae_2 <- df_collocation |>
  filter(id %in% c(126,127)) |> 
  mutate("Monitoring frequency (s)" = 30)

df_sc <- df_sc_aae_1 |> 
  bind_rows(df_sc_aae_2) |> 
  select(date, date_start, time, ir_bcc, blue_bcc, main_exp, serial_number, `Monitoring frequency (s)`)|>
  filter(! ir_bcc < 0) 

df_sc_wide <- df_sc |> 
  pivot_wider(names_from = serial_number,
              values_from = c(ir_bcc, blue_bcc)) |> 
  drop_na() 

lmfit_ir <- df_sc_wide |> 
  group_by(main_exp, date_start) |> 
  do(glance(lm(`ir_bcc_MA200-0416` ~ `ir_bcc_MA200-0420`, data = .))) |> 
  select(date_start, main_exp, r.squared) 

## Performance in terms of R-square score

table_sc_r2 <- df_sc_wide |> 
  group_by(main_exp, `Monitoring frequency (s)`) |> 
  summarise(count = n()) |> 
  left_join(lmfit_ir, by = "main_exp") |> 
  mutate_if(is.numeric,
            round,
            digits = 2)

# calculate aae -----------------------------------------------------------

## use blue and ir wavelength

### calculate background absorption from stationary monitoring data

df_backgr <- df_sm |>
  summarise(mean_ir = mean(ir_babs),
            mean_blue = mean(blue_babs),
            mean_uv = mean(uv_babs))

### extract mean black carbon (bc) abs values 

mean_irbabs = df_backgr$mean_ir
mean_bluebabs = df_backgr$mean_blue
mean_uvbabs = df_backgr$mean_uv

### calculate bc abs by removing the background for blue-ir wavelengths

df_aae_raw_blue_ir <- df_aae |>
  filter(ir_babs  > mean_irbabs,
         blue_babs  > mean_bluebabs)|>
  mutate(ir_babs_bg_corr  = ir_babs  - mean_irbabs,
         blue_babs_bg_corr  = blue_babs - mean_bluebabs) |>
  group_by(id, emission_source, exp_type) |>
  summarise(across(c(blue_babs, ir_babs, ir_babs_bg_corr, blue_babs_bg_corr), mean))|>
  mutate(aae_wo_bg_corr = -log(blue_babs/ir_babs)/log((parameters$wavelength[2])/(parameters$wavelength[5])),
         aae_bg_corr = -log(blue_babs_bg_corr/ir_babs_bg_corr)/log((parameters$wavelength[2])/(parameters$wavelength[5])))|>
  select(-blue_babs, -ir_babs, -ir_babs_bg_corr, -blue_babs_bg_corr) |>
  arrange(by = exp_type)

### add aae with raw data for blue-ir wavelengths

aae_blue_ir <- df_aae |>
  filter(exp_type %in% c("waste_burning", "cooking", "vehicles")) |>
  group_by(id, emission_source, exp_type) |>
  summarise(across(c(blue_babs, ir_babs), mean))|>
  mutate(aae_raw_data = -log(blue_babs/ir_babs)/log((parameters$wavelength[2])/(parameters$wavelength[5])))|>
  select(-blue_babs, -ir_babs) |>
  arrange(by = exp_type) |>
  merge(df_aae_raw_blue_ir) |> 
  mutate(wavelength = "blue_ir")

### calculate bc abs by removing the background for uv-ir wavelengths

df_aae_raw_uv_ir <- df_aae |>
  filter(ir_babs  > mean_irbabs,
         uv_babs  > mean_uvbabs) |>
  mutate(ir_babs_bg_corr  = ir_babs  - mean_irbabs,
         uv_babs_bg_corr  = uv_babs - mean_uvbabs) |>
  group_by(id, emission_source, exp_type) |>
  summarise(across(c(uv_babs, ir_babs, ir_babs_bg_corr, uv_babs_bg_corr), mean))|>
  mutate(aae_wo_bg_corr = -log(uv_babs/ir_babs)/log((parameters$wavelength[1])/(parameters$wavelength[5])),
         aae_bg_corr = -log(uv_babs_bg_corr/ir_babs_bg_corr)/log((parameters$wavelength[1])/(parameters$wavelength[5])))|>
  select(-uv_babs, -ir_babs, -ir_babs_bg_corr, -uv_babs_bg_corr) |>
  arrange(by = exp_type)

### add aae with raw data for uv-ir wavelengths

aae_uv_ir <- df_aae |>
  group_by(id, emission_source, exp_type) |>
  summarise(across(c(uv_babs, ir_babs), mean))|>
  mutate(aae_raw_data = -log(uv_babs/ir_babs)/log((parameters$wavelength[1])/(parameters$wavelength[5])))|>
  select(-uv_babs, -ir_babs) |>
  arrange(by = exp_type) |>
  merge(df_aae_raw_uv_ir) |> 
  mutate(wavelength = "uv_ir")

### combine aae data for blue-ir and uv-ir wavelengths

aae_all <- aae_uv_ir |> 
  bind_rows(aae_blue_ir)   |> 
  mutate(root_source = ifelse(emission_source %in% c("Plastics", "Textiles", "Traffic road", "Diesel truck"), "ff",
                              ifelse(emission_source %in% c("Wood", "Cardboard*", "Garden waste"), "bb", "mixed"))) 

### summary statistics of aae values

aae_summary <- aae_all |> 
  group_by(emission_source, wavelength, root_source) |>
  summarise(mean_aae_raw = mean(aae_raw_data),
            sd_aae_raw = sd(aae_raw_data),
            mean_aae_wo_bg_corr = mean(aae_wo_bg_corr),
            sd_aae_wo_bg_corr = sd(aae_wo_bg_corr),
            mean_aae_bg_corr = mean(aae_bg_corr),
            sd_aae_bg_corr = sd(aae_bg_corr),
            count = n())

### summary statistics of aae values for ff and bb sources

aae_ff_bb_summary <- aae_all |> 
  group_by(wavelength, root_source) |>
  summarise(mean_aae_raw = mean(aae_raw_data),
            sd_aae_raw = sd(aae_raw_data),
            mean_aae_wo_bg_corr = mean(aae_wo_bg_corr),
            sd_aae_wo_bg_corr = sd(aae_wo_bg_corr),
            mean_aae_bg_corr = mean(aae_bg_corr),
            sd_aae_bg_corr = sd(aae_bg_corr),
            count = n()) |> 
  mutate(mean_minus_sd = mean_aae_bg_corr - sd_aae_bg_corr,
         mean_plus_sd = mean_aae_bg_corr + sd_aae_bg_corr)

## figure

wavelength.labs <- c("470/880 nm", "375/880 nm")
names(wavelength.labs) <- c("blue_ir", "uv_ir")
color_pal <- c("black", "brown", "blue")

## aae of ff

p_aae_ff <-  aae_summary |>
  pivot_longer(cols = starts_with(c("mean","sd")), names_to = c("mean_sd", "type"), names_pattern = '(mean_aae|sd_aae)_(raw|wo_bg_corr|bg_corr)', values_to = "aae")|>
  pivot_wider(names_from = mean_sd, values_from = aae)|>
  arrange(by = root_source) |>
  filter(root_source == "ff") |> 
  mutate(type = factor(type, levels = c("raw", "wo_bg_corr", "bg_corr"))) |> 
  mutate(root_source = factor(root_source, levels = c("ff", "bb", "mixed"))) |> 
  ggplot(aes(x= reorder(emission_source, mean_aae), y = mean_aae, color = type)) +
  geom_pointrange(aes(ymin = mean_aae - sd_aae, ymax = mean_aae + sd_aae), 
                  position = position_dodge(width = 0.6),
                  size = 0.2) +
  geom_hline(yintercept = 1, color = "black") +
  #geom_hline(yintercept = 2, color = "brown") + 
  #scale_color_manual(values=c("black", "brown", "blue")) +
  scale_color_discrete(labels=c('Raw data', 'Without babs below background', 'Background correction')) +
  labs(x = "Emission source",
       y = "AAE values (mean and sd)") +
  facet_grid(root_source~wavelength,
             labeller = labeller(wavelength = wavelength.labs)) +
  guides(color=guide_legend(title="AAE calculation method")) +
  theme +
  theme(axis.text.x = element_text(angle = 45, hjust=1),
        panel.grid.major.y = element_line(color = "gray",
                                          size = 0.5,
                                          linetype = 2))

p_aae_ff

## aae of bb

p_aae_bb <-  aae_summary |>
  pivot_longer(cols = starts_with(c("mean","sd")), names_to = c("mean_sd", "type"), names_pattern = '(mean_aae|sd_aae)_(raw|wo_bg_corr|bg_corr)', values_to = "aae")|>
  pivot_wider(names_from = mean_sd, values_from = aae)|>
  arrange(by = root_source) |>
  filter(root_source == "bb") |> 
  mutate(type = factor(type, levels = c("raw", "wo_bg_corr", "bg_corr"))) |> 
  mutate(root_source = factor(root_source, levels = c("ff", "bb", "mixed"))) |> 
  ggplot(aes(x= reorder(emission_source, mean_aae), y = mean_aae, color = type)) +
  geom_pointrange(aes(ymin = mean_aae - sd_aae, ymax = mean_aae + sd_aae), 
                  position = position_dodge(width = 0.6),
                  size = 0.2) +
  #geom_hline(yintercept = 1, color = "black") +
  geom_hline(yintercept = 2, color = "brown") + 
  #scale_color_manual(values=c("black", "brown", "blue")) +
  scale_color_discrete(labels=c('Raw data', 'Without babs below background', 'Background correction')) +
  labs(x = "Emission source",
       y = "AAE values (mean and sd)") +
  facet_grid(root_source~wavelength,
             labeller = labeller(wavelength = wavelength.labs)) +
  guides(color=guide_legend(title="AAE calculation method")) +
  theme +
  theme(axis.text.x = element_text(angle = 45, hjust=1),
        panel.grid.major.y = element_line(color = "gray",
                                          size = 0.5,
                                          linetype = 2),
        panel.grid.minor.y = element_line(color = "gray",
                                          size = 0.5,
                                          linetype = 2))

p_aae_bb

p_aae_mix <-  aae_summary |>
  pivot_longer(cols = starts_with(c("mean","sd")), names_to = c("mean_sd", "type"), names_pattern = '(mean_aae|sd_aae)_(raw|wo_bg_corr|bg_corr)', values_to = "aae")|>
  pivot_wider(names_from = mean_sd, values_from = aae)|>
  filter(root_source == "mixed") |> 
  mutate(type = factor(type, levels = c("raw", "wo_bg_corr", "bg_corr"))) |> 
  ggplot(aes(x= reorder(emission_source, mean_aae), y = mean_aae, color = type)) +
  geom_pointrange(aes(ymin = mean_aae - sd_aae, ymax = mean_aae + sd_aae), 
                  position = position_dodge(width = 0.6),
                  size = 0.2) +
  geom_hline(yintercept = 1, color = "black") +
  geom_hline(yintercept = 2, color = "brown") + 
  #scale_color_manual(values=c("black", "brown", "blue")) +
  scale_color_discrete(labels=c('Raw data', 'Without babs below background', 'Background correction')) +
  labs(x = "Emission source",
       y = "AAE values (mean and sd)") +
  facet_grid(~wavelength,
             labeller = labeller(wavelength = wavelength.labs)) +
  guides(color=guide_legend(title="AAE calculation method")) +
  theme +
  theme(
        panel.grid.major.y = element_line(color = "gray",
                                          size = 0.5,
                                          linetype = 2),
        panel.grid.minor.y = element_line(color = "gray",
                                          size = 0.5,
                                          linetype = 2))

p_aae_mix

# calculate coefficient of determination of AAE values of ff and bb sources -----------------------------------------

## between background corrected data and raw data

lmfit_within_wavelengths_raw_bg_corr <- aae_all |> 
  filter(!root_source %in% "mixed") |>
  group_by(wavelength) |> 
  do(glance(lm(`aae_raw_data` ~ `aae_bg_corr`, data = .))) |> 
  mutate(variable_1 = "aae_raw_data",
         variable_2 = "aae_bg_corr")

## between data with removal of babs values below background and raw data

lmfit_within_wavelengths_raw_wo_bg_corr <- aae_all |> 
  filter(!root_source %in% "mixed") |>
  group_by(wavelength) |> 
  do(glance(lm(`aae_raw_data` ~ `aae_wo_bg_corr`, data = .))) |> 
  mutate(variable_1 = "aae_raw_data",
         variable_2 = "aae_wo_bg_corr")

## between background corrected data using uv-ir and blue-ir wavelengths

### data preparation

lmfit_within_wavelengths <- lmfit_within_wavelengths_raw_bg_corr |> 
  bind_rows(lmfit_within_wavelengths_raw_wo_bg_corr) |> 
  unite(col = "variable_1", c(variable_1, wavelength), sep = "_", remove = FALSE) |> 
  unite(col = "variable_2", c(variable_2, wavelength), sep = "_") 

aae_all_long_blue_ir <- aae_all |> 
  filter(!root_source %in% "mixed",
         wavelength == "blue_ir") |>
  select(id, emission_source, aae_raw_data, aae_bg_corr, aae_wo_bg_corr) |>
  pivot_longer(cols = starts_with(c("aae")), names_to = c("aae_type"), values_to = "aae_blue_ir")

aae_all_long_uv_ir <- aae_all |> 
  filter(!root_source %in% "mixed",
         wavelength == "uv_ir") |>
  select(id, emission_source, aae_raw_data, aae_bg_corr, aae_wo_bg_corr) |>
  pivot_longer(cols = starts_with(c("aae")), names_to = c("aae_type"), values_to = "aae_uv_ir")

### calculate r2 and merge all the data

lmfit_aae_r2 <- aae_all_long_blue_ir |> 
  left_join(aae_all_long_uv_ir, by = c("id", "emission_source", "aae_type")) |> 
  group_by(aae_type) |>
  do(glance(lm(`aae_blue_ir` ~ `aae_uv_ir`, data = .))) |> 
  mutate(variable_1 = "blue_ir",
         variable_2 = "uv_ir") |> 
  unite(col = "variable_1", c(aae_type, variable_1), sep = "_", remove = FALSE) |> 
  unite(col = "variable_2", c(aae_type, variable_2), sep = "_") |> 
  rbind(lmfit_within_wavelengths) |> 
  mutate_if(is.numeric,
            round,
            digits = 2)

### rename variables

lmfit_aae_r2$variable_1 <- recode(lmfit_aae_r2$variable_1, 	
                                  aae_bg_corr_blue_ir = 'AAE- background correction (475/880 nm)',
                                  aae_raw_data_blue_ir  = 'AAE- raw data (475/880 nm)',
                                  aae_wo_bg_corr_blue_ir = 'AAE- without babs below background (475/880 nm)',
                                  aae_raw_data_uv_ir  = 'AAE- raw data (370/880 nm)')

lmfit_aae_r2$variable_2 <- recode(lmfit_aae_r2$variable_2, 	
                                  aae_bg_corr_uv_ir = 'AAE- background correction (370/880 nm)',
                                  aae_raw_data_uv_ir  = 'AAE- raw data (370/880 nm)',
                                  aae_wo_bg_corr_uv_ir = 'AAE- without babs below background (370/880 nm)',
                                  aae_raw_data_blue_ir  = 'AAE- raw data (475/880 nm)',
                                  aae_bg_corr_blue_ir = 'AAE- background correction (475/880 nm)',
                                  aae_wo_bg_corr_blue_ir = 'AAE- without babs below background (475/880 nm)')

## final table with selected variables

table_aaer2 <- lmfit_aae_r2 |> 
  rename(`Variable 1` = variable_1,
         `Variable 2` = variable_2,
         `R2` = r.squared,
         `R2-adjusted` = adj.r.squared) |>
  select(`Variable 1`, `Variable 2`, `R2`, `R2-adjusted`) 

# calculate coefficient of determination of AAE values of mixed waste -----------------------------------------

## use same code as that for ff and bb sources

lmfit_within_wavelengths_raw_bg_corr_mix <- aae_all |> 
  filter(root_source %in% "mixed") |>
  group_by(wavelength) |> 
  do(glance(lm(`aae_raw_data` ~ `aae_bg_corr`, data = .))) |> 
  mutate(variable_1 = "aae_raw_data",
         variable_2 = "aae_bg_corr")

lmfit_within_wavelengths_raw_wo_bg_corr_mix <- aae_all |> 
  filter(root_source %in% "mixed") |>
  group_by(wavelength) |> 
  do(glance(lm(`aae_raw_data` ~ `aae_wo_bg_corr`, data = .))) |> 
  mutate(variable_1 = "aae_raw_data",
         variable_2 = "aae_wo_bg_corr")

lmfit_within_wavelengths_mix <- lmfit_within_wavelengths_raw_bg_corr_mix |> 
  bind_rows(lmfit_within_wavelengths_raw_wo_bg_corr_mix) |> 
  unite(col = "variable_1", c(variable_1, wavelength), sep = "_", remove = FALSE) |> 
  unite(col = "variable_2", c(variable_2, wavelength), sep = "_") 

aae_all_long_blue_ir_mix <- aae_all |> 
  filter(root_source %in% "mixed",
         wavelength == "blue_ir") |>
  select(id, emission_source, aae_raw_data, aae_bg_corr, aae_wo_bg_corr) |>
  pivot_longer(cols = starts_with(c("aae")), names_to = c("aae_type"), values_to = "aae_blue_ir")

aae_all_long_uv_ir_mix <- aae_all |> 
  filter(root_source %in% "mixed",
         wavelength == "uv_ir") |>
  select(id, emission_source, aae_raw_data, aae_bg_corr, aae_wo_bg_corr) |>
  pivot_longer(cols = starts_with(c("aae")), names_to = c("aae_type"), values_to = "aae_uv_ir")

lmfit_aae_r2_mix <- aae_all_long_blue_ir_mix |> 
  left_join(aae_all_long_uv_ir_mix, by = c("id", "emission_source", "aae_type")) |> 
  group_by(aae_type) |>
  do(glance(lm(`aae_blue_ir` ~ `aae_uv_ir`, data = .))) |> 
  mutate(variable_1 = "blue_ir",
         variable_2 = "uv_ir") |> 
  unite(col = "variable_1", c(aae_type, variable_1), sep = "_", remove = FALSE) |> 
  unite(col = "variable_2", c(aae_type, variable_2), sep = "_") |> 
  rbind(lmfit_within_wavelengths_mix) |> 
  mutate_if(is.numeric,
            round,
            digits = 2)

lmfit_aae_r2_mix$variable_1 <- recode(lmfit_aae_r2_mix$variable_1, 	
                                      aae_bg_corr_blue_ir = 'AAE- background correction (475/880 nm)',
                                      aae_raw_data_blue_ir  = 'AAE- raw data (475/880 nm)',
                                      aae_wo_bg_corr_blue_ir = 'AAE- without babs below background (475/880 nm)',
                                      aae_raw_data_uv_ir  = 'AAE- raw data (370/880 nm)')

lmfit_aae_r2_mix$variable_2 <- recode(lmfit_aae_r2_mix$variable_2, 	
                                      aae_bg_corr_uv_ir = 'AAE- background correction (370/880 nm)',
                                      aae_raw_data_uv_ir  = 'AAE- raw data (370/880 nm)',
                                      aae_wo_bg_corr_uv_ir = 'AAE- without babs below background (370/880 nm)',
                                      aae_raw_data_blue_ir  = 'AAE- raw data (475/880 nm)',
                                      aae_bg_corr_blue_ir = 'AAE- background correction (475/880 nm)',
                                      aae_wo_bg_corr_blue_ir = 'AAE- without babs below background (475/880 nm)')

table_aaer2_mix <- lmfit_aae_r2_mix |> 
  rename("Variable 1" = variable_1,
         `Variable 2` = variable_2,
         `R2` = r.squared,
         `R2-adjusted` = adj.r.squared) |>
  select(`Variable 1`, `Variable 2`, `R2`, `R2-adjusted`) 

# t-test on AAE values ----------------------------------------------------

## calculate mean and standard deviation of AAE values for each source

aae_avg <- aae_all |> 
  group_by(root_source, wavelength) |> 
  summarise(mean = mean(aae_raw_data),
            sd = sd(aae_raw_data)) |> 
  mutate(mean_plus_sd = mean + sd,
         mean_minus_sd = mean - sd) 

## define order of sources for t-test

new_order_ttest_aae <- c("Diesel truck", "Traffic road", "Plastics", "Textiles", "Wood", "Garden waste", "Cardboard*")

aae_blue_ir <- aae_blue_ir |> 
  arrange(factor(emission_source, levels = new_order_ttest_aae))

## t-test 

combinations <- combn(unique(aae_blue_ir$emission_source), 2, simplify = FALSE)

aae_ttest <-  list()

for (i in seq_along(combinations)){
  aae_ttest[[i]] <- aae_blue_ir |>
    infer::t_test(aae_bg_corr ~ emission_source,
                  order = combinations[[i]],
                  alternative = "two-sided")
}

## prepare data for visualization

aae_ttest_results <- aae_ttest |>
  bind_rows() |>
  mutate(id = combinations) |>
  unnest_wider(col = id, names_sep = "")|>
  mutate(p_value = sprintf("%.2f", p_value)) |>
  select(id1, id2, p_value) |> 
  pivot_wider(names_from = id1, values_from = p_value) 

new_row <- tibble("id2" = "Diesel truck")

aae_ttest_result_matrix <- new_row |> 
  bind_rows(aae_ttest_results) |> 
  mutate("Mixed waste" = NA) |> 
  column_to_rownames("id2")

aae_ttest_result_matrix <- as.matrix(aae_ttest_result_matrix)

aae_ttest_result_matrix[is.na(aae_ttest_result_matrix)] <- ""

## define custom colors for values above and below 0.05

library(reshape)

melt_aae <- melt(aae_ttest_result_matrix) |> 
  mutate(similar.not_similar = ifelse(value == "", "", ifelse(value > 0.05, "Similar", "Not similar")))
  

my_colors <- ifelse(melt_aae$value == "", "black",
                    ifelse(melt_aae$value > 0.05, "green", "red"))

## plot for t-test results

p_aae_t_test <- melt_aae |> 
  ggplot(aes(x = X1, y = X2)) +
  geom_tile(aes(fill = similar.not_similar), width=0.7, height=0.7) +
  geom_text(aes(label = value), color = "white") +
  scale_fill_manual(values = c("white", "maroon", "darkgreen"), labels = c("", "<0.05", ">0.05")) +
  guides(fill=guide_legend(title="p-value")) +
  labs(x = "Emission source",
       y = "Emission source") +
  theme +
  theme(axis.text.x = element_text(angle = 45, hjust=1))

p_aae_t_test

## t-test after combining all ff sources and all bb sources

aae_blue_ir_root_source <- aae_blue_ir |> 
  mutate(root_source = ifelse(emission_source %in% c("Plastics", "Textiles", "Traffic road", "Diesel truck"), "ff",
                              ifelse(emission_source %in% c("Wood", "Cardboard*", "Garden waste"), "bb", "mixed"))) 

combinations_ff_bb <- combn(unique(aae_blue_ir_root_source$root_source), 2, simplify = FALSE)

aae_ttest_ff_bb <-  list()

for (i in seq_along(combinations_ff_bb)){
  aae_ttest_ff_bb[[i]] <- aae_blue_ir_root_source |>
    infer::t_test(aae_bg_corr ~ root_source,
                  order = combinations_ff_bb[[i]],
                  alternative = "two-sided")
}

aae_ttest_results_ff_bb <- aae_ttest_ff_bb |>
  bind_rows() |>
  mutate(id = combinations_ff_bb) |>
  unnest_wider(col = id, names_sep = "") |>
  select(id1, id2, p_value) |>
  mutate(similar.not_similar = ifelse(p_value > 0.05, "Similar", "Not similar"))

# verification of AAE values ----------------------------------------------

## mobile monitoring data

df_mm_aae_paper <- df_mm 

## personal monitoring data

df_pm_aae_paper <- df_pm |> 
  filter(id == 90) |> 
  filter(time > parse_hms(c("15:03:00"))) ### remove values in Ndirande main road 
  ### the aim is to show values within remote areas of the settlement

## calculate AAE values for mobile and personal monitoring data

df_monitoring <- df_mm_aae_paper |> 
  bind_rows(df_pm_aae_paper) |> 
  drop_na(ir_bcc) |> 
  filter(ir_bcc > 0,
         blue_bcc > 0) |> 
  mutate(aae_blue_ir = -log(blue_babs/ir_babs)/log((parameters$wavelength[2])/(parameters$wavelength[5]))) |>
  mutate(root_source = ifelse(aae_blue_ir < 1.29, "ff",
                              ifelse(aae_blue_ir > 1.63, "bb", "mixed"))) 

## count of AAE values for each settlement in different source categories

df_monitoring_sa <- df_monitoring |> 
  group_by(settlement_id,
           root_source) |>
  summarise(count = n()) |> 
  pivot_wider(names_from = root_source, values_from = count)

## data preparation for figure on verification of AAE values

df_map <- df_monitoring |> 
  filter(settlement_id %in% c("Ndirande", "Nyambadwe"))

df_settlement <- tibble(
  settlement_id = c("Nyambadwe", "Ndirande"),
  lat = c(-15.767, -15.777),
  long = c(35.02, 35.04))

## figure on verification of AAE values

map_data_sf <- st_as_sf(df_map, coords = c("long", "lat"))

p_aae_verification <- ggplot() +
  base_map(st_bbox(map_data_sf), basemap = "mapnik", increase_zoom = 3) +
  geom_path(data = df_map, 
            aes(x = long, 
                y = lat, 
                group = as.factor(id), 
                color = root_source), 
            size = 2) +
  geom_text(data = df_settlement,
            aes(x = long,
                y = lat,
                label = settlement_id),
            size = 5) +
  guides(color=guide_legend(title="AAE470/880 (Emission source)")) +
  labs(x = "Longitude",
       y = "Latitude") +
  scale_color_manual(values = c("brown", "black", "blue"), labels = c("> 1.63 (biomass-based)", "< 1.29 (fossil-fuel-based)", "1.29-1.63 (mix of two)")) +
  theme(panel.border = element_rect(color = "black", size = 1, linetype = "solid", fill = alpha("black", 0))) + 
  theme + # Add thick border
  theme(legend.background = element_rect(linetype = 1, size = 0.1, colour = 1),
        legend.position = c(0.99, 0.99),
        legend.title = element_text(size=11),
        legend.text = element_text(size=11), #change legend title font size
        legend.justification = c("right", "top"),
        legend.box.just = "right",
        legend.margin = margin(1, 1, 1, 1)) 

