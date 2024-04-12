
# 5 mile data cleaning

library(base)
library(tidyverse)
library(readxl)
library(benthos)
library(lubridate)


# TO RUN TO UPDATE COMPLETENESS:
# HIGHLIGHT FROM 'START' TO 'STOP' BELOW
# HIT CTRL+SHFT+ENTER
# CHECK CONSOLE FOR ERRORS - IF NONE, YOU'RE DONE

#
### START ###
#

# cover straight from CN
cover_raw <- read.csv('./Data/Percent_cover.csv')

# trait db McLean lab
nc_eb_trait_raw <- read_excel('./Data/Traits_LaCroce.xlsx') %>%
    select(-blank, -notes, -source)
eb_trait_list <- read_excel('./Data/NC EB trait list.xlsx')
color_raw <- read_excel('../Artificial Reefs/Data/Epibiota/Trait data/color_trait_raw_data.xlsx', sheet = 'color_traits')
meta_raw <- read.csv('./Data/lacroce_metadata.csv')
temp_raw <- read.csv('./Data/vmileA_bottomwater_temp_09.XX.15-12.XX.2023.csv')

mel_spp <- colnames(cover_raw)[3:length(cover_raw)] %>%
    str_replace_all("[.]", replacement = " ") %>%
    trimws() %>%
    as.data.frame() %>%
    rename("Name" = ".") %>%
    mutate(Abund = colSums(cover_raw[3:length(cover_raw)])) %>%
    filter(Name %in% nc_eb_trait_raw$Name, Abund > 0)

color <- color_raw %>%
    group_by(Organism) %>%
    mutate(Percent_complete = sum(L*0+1, na.rm = T)/15,
           across(L:b, ~ mean(.x, na.rm = T))) %>%
    ungroup() %>% 
    select(Organism, L, a, b, Percent_complete) %>%
    rename(Name = Organism) %>%
    unique() %>%
    mutate(across(L:b, ~ round(.x, 3)),
           Color_LAB = paste(L,a,b, sep = ',')) %>%
    select(-L, -a, -b) %>%
    full_join(mel_spp) %>%
    filter(Name %in% mel_spp$Name)

checked <- left_join(mel_spp, nc_eb_trait_raw) %>%
    select(-Color_LAB) %>%
    left_join(color)

checked$Complete <- complete.cases(checked)
checked$Complete_but_colorless <- complete.cases(select(checked, -Color_LAB, -Percent_complete, -Complete))

write.csv(checked, paste('Data/Data_cleaned/x', as.character(today()), 'completeness_check.csv'))
write.csv(color, paste('Data/Data_cleaned/x', as.character(today()), 'completeness_colors.csv'))
write.csv(select(checked, -Complete), paste('Data/Data_cleaned/x', as.character(today()), 'traits_5_mile_clean.csv'), row.names = F)

#
### STOP ###
#


# # temp stuff
# 
# meta <- meta_raw %>%
#     mutate(Date = ymd(Date)) %>%
#     select(Name, Date)
# 
# temp <- temp_raw %>%
#     mutate(Date.Time = mdy_hm(Date.Time),
#            Date = date(Date.Time),
#            Temp = Temperature..C.) %>%
#     select(Date, Temp) %>%
#     subset(year(Date) < 2018) %>%
#     group_by(Date) %>%
#     mutate(Temp = mean(Temp)) %>%
#     ungroup() %>%
#     unique() %>%
#     left_join(meta)
# 
# # check that the entire temp range was surveyed
# 
# ggplot(temp, aes(x = Date, y = Temp)) + 
#     geom_line() +
#     geom_point(aes(color = is.na(Name))) +
#     geom_vline(data = meta, aes(xintercept = meta$Date))
# 
#
# meta <- meta %>%
#     left_join(temp)
#
# cover <- cover_raw %>%
#     mutate(Name = Image.name) %>%
#     select(-Points, -Image.name)
#
# is.zero <- function(n){return(!!n)}
#
# cover_pa <- cover %>%
#     mutate(across(Oculina.arbuscula:Zonaria.tournefortii, is.zero))
#
# cover_temp <- cover_pa %>%
#     left_join(meta) %>%
#     mutate(across(Oculina.arbuscula:Zonaria.tournefortii, ~ .x * Temp *.x / .x)) %>%
#     select(-Name, -Temp, -Date)
#
# cover_min <- cover_temp %>%
#     mutate(across(Oculina.arbuscula:Zonaria.tournefortii, ~ min(.x, na.rm = T))) %>%
#     unique() %>%
#     t() %>%
#     as.data.frame() %>%
#     rownames_to_column() %>%
#     rename('Species' = 'rowname', 'min' = '1')
# cover_max <- cover_temp %>%
#     mutate(across(Oculina.arbuscula:Zonaria.tournefortii, ~ max(.x, na.rm = T))) %>%
#     unique() %>%
#     t() %>%
#     as.data.frame() %>%
#     rownames_to_column() %>%
#     rename('Species' = 'rowname', 'max' = '1')
# 
# temp_ranges <- full_join(cover_min, cover_max)
# temp_ranges$Occurences <- cover_pa %>%
#     left_join(meta) %>%
#     group_by(Date) %>%
#     mutate(across(Oculina.arbuscula:Zonaria.tournefortii, ~ sum(.x))) %>%
#     ungroup() %>%
#     select(-Name, -Temp) %>%
#     unique() %>%
#     mutate(across(Oculina.arbuscula:Zonaria.tournefortii, is.zero)) %>%
#     select(-Date) %>%
#     colSums()
# 
#     
# 
# write.csv(temp_ranges, file = './Data/spp_temp_ranges_lacroce.csv')
