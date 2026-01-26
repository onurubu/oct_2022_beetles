##### Coskun Kucukkaragoz
#### January 2026
### Veg and soil data analysis


{
  
soil_dat_02 <- data.table::fread("./soil/soil_data_2002.csv") %>% mutate(across(c(site, replicate, year), factor)) %>% rename("rock_soil" = "rock")
soil_dat_22 <-  data.table::fread("./soil/soil_data_2022_working.csv") %>% mutate(across(c(site, replicate, year), factor)) %>% rename("rock_soil" = "rock")
soil_dat <- rbind(soil_dat_02, soil_dat_22)
rm(soil_dat_02, soil_dat_22)

}

