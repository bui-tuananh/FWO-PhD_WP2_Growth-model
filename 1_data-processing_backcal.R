# info ----
# input: otl_full.rds (finalised data from WP1 - https://github.com/bui-tuananh/FWO-PhD_WP1_Growth-linear-mixed-model)
# output: otl_full_backcal.rds

# setup ----
library(tidyverse)
library(TMB)
library(boot)
library(data.table)
library(lme4)

# data ----
# read otl data
otl <- readRDS("./input/otl_full.rds")

# backcal length ----
#### data
otl_sub <- otl %>% 
    mutate(id = FishID,
           age = Age,
           aac = AgeAtCapture,
           year_sample = SamplingYear,
           year = GrowingYear,
           cohort = Cohort,
           otl_increment = AnnulusDiameterIncrement.um,
           otl_diameter = AnnulusDiameter.um,
           otl_width = OtolithWidth.um,
           body_len = Length.mm,
           month = month(SamplingDate),
           pop = IcesAreaGroup) %>%
    select(id:pop)

# remove 1 outlier ind sol_fab_0575 (see WP1) 
# remove otl with otl_width < 1500 (typo error)
otl_sub <- otl_sub %>% filter(id != "sol_fab_0575",
                              otl_width >= 1500) 

otl_sub_cpt <- otl_sub %>%
    select(id, aac, otl_width, body_len, month, pop, year_sample, cohort) %>%
    unique() %>%
    mutate(
        log.aac = log(aac),
        log.otl_width = log(otl_width),
        log.body_len = log(body_len),
        fmonth = factor(month),
        pop.month = paste0(pop, ":", month))

#### model
otl_sub_backcal <- tibble()
for(p in c("4bc", "7a", "8ab")) {
    otl_sub_cpt_pop <- otl_sub_cpt %>% filter(pop == p)
    m <- lmer(log.body_len ~ 1 + log.otl_width + (1 | fmonth), 
                  data = otl_sub_cpt_pop)
    
    #### predict
    df_par <- tibble(pop = p,
                     a = c(fixef(m)[1]),
                     b = c(fixef(m)[2]))
    otl_sub_backcal_temp <- otl_sub %>%
        filter(pop == p) %>%
        left_join(df_par) %>%
        mutate(log.otl_width = log(otl_diameter),
               log.age = log(age),
               backcal_len = exp(a + b*log.otl_width))
    
    otl_sub_backcal <- bind_rows(otl_sub_backcal, otl_sub_backcal_temp)
}

otl_sub_backcal <- otl_sub_backcal %>% select(-a, -b)
write_rds(otl_sub_backcal, "./input/otl_full_backcal.rds")
