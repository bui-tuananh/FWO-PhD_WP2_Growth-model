# info ----
# input: 
# - data/sol_gonads_2004_2023.rds (from data_processing-script)
# - 2_data-analysis_backcal_4bc/7a/8ab
# - 3_model-validation_sdreport
# output: 
# - lp50_maturity
# - lp50_maturity_boot
# - lp50_growth

# PMRN - ref code Adriaan
# function.r

# note from Grift et al., 2003
# Reaction norm midpoints (LP50, the
# length at which the probability to mature
# is 50%) were calculated by determining the
# lengths that lead to probabilities of maturing of
# 50%.

# setup ----
library(tidyverse)
library(lme4)

# data ----
## obs ----
obs <- read_rds("./input/sol_gonads_2004_2023.rds")

# select areas, remove NA, abnormal, skipped spawning stages, NA age, select female
obs_sub <- obs %>% 
  filter(HAU.IcesArea %in% c("4b", "4c", "7a", "7d", "7f", "7g", "8a", "8b"),
         is.na(SPE.MaturityStageDescription) == F,
         !SPE.MaturityStageDescription %in% c("Abnormal", "Skipped spawning"), #2 obs Abnormal, 4 obs Skipped spawning
         is.na(SPA.Age) == F,
         SPE.Sex == "F") %>%
  mutate(SPE.MaturityStageDescription = if_else(SPE.MaturityStageDescription %in% c("Maturing", "Developing"),
                                                "Developing",
                                                SPE.MaturityStageDescription))

# sampling week: majority of trips within 1 week, yet there are trips in 1-2 weeks
# -> get the middle value of week, month
# rename field and do some merging, factoring

df_area_obs <- tibble(HAU.IcesArea = c("4b", "4c", "7a", "7d", "7f", "7g", "8a", "8b"),
                      pop          = c("4bc", "4bc", "7a", "7d", "7fg", "7fg", "8ab", "8ab"))

obs_sub <- obs_sub %>%
  left_join(df_area_obs) %>%
  mutate(id = SpecimenID, 
         age = SPA.Age,
         year = TRI.Year,
         cohort = year - age,
         month = round((month(TRI.DepartureDate) + month(TRI.ReturnDate))/2),
         week = round((week(TRI.DepartureDate) + week(TRI.ReturnDate))/2),
         week_diff = week(TRI.ReturnDate) - week(TRI.DepartureDate),
         length = SPE.Length,
         body_weight = SPE.Weight,
         gonad_weight = SPE.WeightGonads,
         gsi = gonad_weight/body_weight*100,
         maturity_stage = factor(SPE.MaturityStageDescription, level = c("Immature", "Developing", "Spawning", "Spent")),
         maturity = if_else(SPE.MaturityStageDescription %in% c("Immature"), 0, 1),
         maturity_desc = if_else(SPE.MaturityStageDescription %in% c("Immature"), "Immature", "Mature"),
  ) %>%
  select(id, pop, year, cohort, month, week, age, length, body_weight, gonad_weight, gsi, maturity_stage, maturity, maturity_desc)

# filter gsi < 1000
obs_sub <- obs_sub %>% 
  filter(gsi < 100)

# obs_sub_4bc <- obs_sub %>% filter(pop == "4bc", gonad_weight != 0)
# obs_sub_7a <- obs_sub %>% filter(pop == "7a", gonad_weight != 0)

# data for analysis
# use all the data to increase sample size
# may consider filtering months later (see Mollet 2007)
data_obs <- obs_sub

data_obs %>%
  group_by(pop, maturity_desc) %>%
  summarize(n = n())

# filter 
# age >= 1 (only 1 obs age 0)
# age <= X (consider later)
# only first 6 months - peak gsi, spawning + roughly round age in age reading (fish spawn in spring, then has round age in next spring)

data_obs <- data_obs %>%
  filter(age >= 1)

data_obs_sub <- data_obs %>%
  mutate(source = "obs") %>%
  select(id, pop, year, cohort, month, age, length, maturity, source) 


## backcal otl ----
data_backcal <- read_rds("./input/otl_full_backcal.rds") 

data_backcal_sub <- data_backcal %>%
    filter(aac >= 5, 
           age <= 15,
           !id %in% c("sol_fab_1375")) %>%
  rename(length = backcal_len) %>%
  mutate(source = "otl") %>%
  select(id, pop, year, cohort, age, length, source) 

## all data ----
data_all <- bind_rows(filter(data_obs_sub, !id %in% data_backcal_sub$id), # not include obs with backcal to avoid duplication
                         data_backcal_sub)

# lp50 ----
# function to find lp50 from results glm and pmrn - maturity ogive model: maturity ~ age*length
f_lp50_age1 <- function(alpha, beta_age, beta_length, beta_age_length, age, length) {
  f_o1 = 1/(1 + exp(-(alpha + beta_age*age + beta_length*length + beta_age_length*age*length)))
  
  f_o1 - 0.5
}

f_lp50 <-  function(alpha, beta_age, beta_length, beta_age_length, age, dl, length) {
  f_o1 = 1/(1 + exp(-(alpha + beta_age*age + beta_length*length + beta_age_length*age*length)))
  f_o0 = 1/(1 + exp(-(alpha + beta_age*(age-1) + beta_length*(length - dl) + beta_age_length*(age-1)*(length-dl))))
  
  2*f_o1 - f_o0 - 1
}

# function to find lp50 from data 
func_lp50 <- function(data_mat, data_len) {
  
  # maturity ogive
  ogive.mod<-glm(maturity~age*length,family=binomial, data=data_mat)
  
  # growth model 
  growth.mod<-lm(length ~ 0 + as.factor(age), data=data_len)
  
  # delta length
  pred_dl <- tibble(age = seq(2, 6))
  pred_dl <- pred_dl %>%
    mutate(dl = predict(growth.mod, data.frame(age=pred_dl$age)) - predict(growth.mod, data.frame(age=pred_dl$age-1)))
  
  #### lp50
  pred <- tibble()
  for(a in seq(1, 6)) {
    
    if(a == 1) {
      ## age 1 - assumed no maturation at age 0 -> pred from maturity ogive model
      lp50 <- uniroot(f_lp50_age1, 
                      interval = c(0, 500), 
                      alpha = coefficients(ogive.mod)["(Intercept)"], 
                      beta_age = coefficients(ogive.mod)["age"], 
                      beta_length = coefficients(ogive.mod)["length"],
                      beta_age_length = coefficients(ogive.mod)["age:length"], 
                      age = a)$root
      
      pred_temp <- tibble(pop = p,
                               age = a,
                               lp50 = lp50)
    } else {
      ## age 2-6 
      lp50 <- uniroot(f_lp50, 
                      interval = c(0, 500),
                      alpha = coefficients(ogive.mod)["(Intercept)"], 
                      beta_age = coefficients(ogive.mod)["age"], 
                      beta_length = coefficients(ogive.mod)["length"],
                      beta_age_length = coefficients(ogive.mod)["age:length"], 
                      age = a, 
                      dl = filter(pred_dl, age == a)$dl)$root
      
      pred_temp <- tibble(pop = p,
                               age = a,
                               lp50 = lp50) #lp50 can have 2 values
    }
    pred <- bind_rows(pred, pred_temp)
  }
  
  return(pred)
}

pred_lp50 <- tibble() 
for(p in c("4bc", "7a")) {
  data_mat <- data_obs %>% filter(pop == p, age <= 6, month <= 6)
  data_len <- data_all %>% filter(pop == p, (source == "obs" & month <= 6) | source == "otl")
  #data_len <- data_all %>% filter(pop == p, age <= 6, (source == "obs" & month <= 6) | source == "otl") # the same
  
  pred_lp50_temp <- func_lp50(data_mat, data_len)
  
  pred_lp50 <- bind_rows(pred_lp50, pred_lp50_temp) 
}

write_rds(pred_lp50, "./output/lp50_maturity.rds")

# lp50 boot ----
func_lp50_boot <- function(data_mat, data_len, r = 10) {
  pred <- tibble()
  for (i in 1:r) {
    print(paste0("processing ", i))
    
    set.seed(i)
    # Sample indices with replacement
    indices_mat <- sample(1:nrow(data_mat), replace = TRUE)
    indices_len <- sample(1:nrow(data_len), replace = TRUE)
    
    # Bootstrap sample
    data_mat_sample <- data_mat[indices_mat, ]
    data_len_sample <- data_len[indices_len, ]
    
    # maturity ogive
    ogive.mod<-glm(maturity~age*length,family=binomial, data=data_mat_sample)
    
    # growth model 
    growth.mod<-lm(length ~ 0 + as.factor(age), data=data_len_sample)
    
    # delta length
    pred_dl <- tibble(age = seq(2, 6))
    pred_dl <- pred_dl %>%
      mutate(dl = predict(growth.mod, data.frame(age=pred_dl$age)) - predict(growth.mod, data.frame(age=pred_dl$age-1)))
    
    #### lp50
    for(a in seq(1, 6)) {
      
      if(a == 1) {
        if (is.null(tryCatch({
        ## age 1 - assumed no maturation at age 0 -> pred from maturity ogive model, also no data at age 0
        lp50 <-       uniroot(f_lp50_age1, 
                              interval = c(-1000, 1000), 
                              alpha = coefficients(ogive.mod)["(Intercept)"], 
                              beta_age = coefficients(ogive.mod)["age"], 
                              beta_length = coefficients(ogive.mod)["length"],
                              beta_age_length = coefficients(ogive.mod)["age:length"], 
                              age = a)$root
        
        pred_temp <- tibble(ite_id = i,
                                 pop = unique(data_mat$pop),
                                 age = a,
                                 lp50 = lp50)
        
        }, error = function(e){}))) {
          # skip if error in uniroot 
          print("error");
          next
        }
      } else {
        ## age 2-6 
        if (is.null(tryCatch({
          lp50 <- uniroot(f_lp50, 
                          interval = c(-1000, 1000),
                          alpha = coefficients(ogive.mod)["(Intercept)"], 
                          beta_age = coefficients(ogive.mod)["age"], 
                          beta_length = coefficients(ogive.mod)["length"],
                          beta_age_length = coefficients(ogive.mod)["age:length"], 
                          age = a, 
                          dl = filter(pred_dl, age == a)$dl)$root
          
        }, error = function(e){}))) {
          # skip if error in uniroot 
          print("error");
          next
        }
        
        pred_temp <- tibble(ite_id = i,
                                 pop = unique(data_mat$pop),
                                 age = a,
                                 lp50 = lp50) #lp50 can have 2 values
      }
      pred <- bind_rows(pred, pred_temp)
    }
  }
  return(pred)
}

#func_lp50_boot(data_mat, data_len, r = 2)

pred_lp50_boot <- tibble() 
for(p in c("4bc", "7a")) {
  data_mat <- data_obs %>% filter(pop == p, age <= 6, month <= 6)
  data_len <- data_all %>% filter(pop == p, (source == "obs" & month <= 6) | source == "otl")
  
  pred_lp50_boot_temp <- func_lp50_boot(data_mat, data_len, r = 1000)
  
  pred_lp50_boot <- bind_rows(pred_lp50_boot, pred_lp50_boot_temp) 
}

pred_lp50_boot <- bind_rows(pred_lp50_boot_4bc, pred_lp50_boot_7a)

write_rds(pred_lp50_boot, "./output/lp50_maturity_boot.rds")

# calculate 95ci
## remove unrealistic estimate of lp50: < 0, 1000 at age 6
pred_lp50_boot_sum <- 
  pred_lp50_boot %>%
  filter(lp50 > 0) %>%
  filter(lp50 < 1000) %>%
  group_by(pop, age) %>%
  summarize(ci95_lower = quantile(lp50, probs = c(0.025)),
            ci95_upper = quantile(lp50, probs = c(0.975)))

# lp50 growth model ----
# function to calculate lp50 growth

fn_lp50_growth <- function(data) {
  
  pred <- tibble()
  
  # maturity ogive
  ogive.mod<-glm(maturity~age*length,family=binomial,data=data)
  
  for(a in seq(1, 6)) {
    
    if (is.null(tryCatch({
      lp50 <- uniroot(f_lp50_age1,
                      interval = c(-1000, 1000),
                      alpha = coefficients(ogive.mod)["(Intercept)"],
                      beta_age = coefficients(ogive.mod)["age"],
                      beta_length = coefficients(ogive.mod)["length"],
                      beta_age_length = coefficients(ogive.mod)["age:length"],
                      age = a)$root
      
    }, error = function(e){}))) {
      # skip if error in uniroot 
      print("error");
      next
    }
    
    pred_temp <- tibble(pop = unique(data$pop),
                             age = a,
                             lp50 = lp50) #lp50 can have 2 values
    pred <- bind_rows(pred, pred_temp)
  }
  return(pred)
}

data_rep <- read_rds("output/individual_h_Tmat_a100_age0_datarep_v2.rds")

# calculate Amat
data_rep <- data_rep %>%
  mutate(Tmat_roundup = if_else(round(Tmat) > Tmat, round(Tmat), round(Tmat) + 1)) %>%
  #mutate(Amat = if_else(Tmat_roundup - Tmat >= 1/3, Tmat_roundup, Tmat_roundup + 1)) %>%
  #mutate(Amat = if_else(Tmat_roundup - Tmat >= 1/2, Tmat_roundup, Tmat_roundup + 1)) %>%
  mutate(Amat = Tmat_roundup) %>%
  mutate(maturity = if_else(age < Amat, 0, 1)) %>%
  mutate(length = vt_sim)

pred_lp50_growth <- tibble()
for(p in c("4bc", "7a")) {
  data <- data_rep %>% filter(pop == p, age <= 6, cohort >= 1998) 
  
  pred_lp50_growth_temp <- fn_lp50_growth(data)
  
  pred_lp50_growth <- bind_rows(pred_lp50_growth, pred_lp50_growth_temp)
  
}

write_rds(pred_lp50_growth, "./output/lp50_growth.rds")

# note - demonstrate o0, o1, pmrn ---
# o1: probability of being mature at age t and length l
# o0: probability of being mature already one year before at age (t-1) and length (l - dl), where dl is growth increment
# pmrn = o1 when age = 1 (assuming this is the first maturing age + no info on age 0 to calculate dl)
# pmrn = (o1-o0)/(1-o0) when age > 1

# ref: 
# Mollet et al., 2013, Journal of Sea Research 84 109-121 - https://www.sciencedirect.com/science/article/abs/pii/S1385110112001992
# Grift et al., 2003, Marine Ecology Progress Series 257 247-257 - https://www.int-res.com/abstracts/meps/v257/p247-257/

# demonstrate using 4bc data
data_mat <- data_obs %>% filter(pop == "4bc", age <= 6, month <= 6)
data_len <- data_all %>% filter(pop == "4bc", (source == "obs" & month <= 6) | source == "otl")

# maturity ogive
ogive.mod<-glm(maturity~age*length,family=binomial, data=data_mat)

# growth model 
growth.mod<-lm(length ~ 0 + as.factor(age), data=data_len)

min.ageM = 1
age_range = seq(min.ageM, 6)
length_range <- sort(unique(filter(data_len, age %in% age_range)$length))
#length_range <- seq(-1000, 1000)

pred <- tibble(pop = p,
               age = rep(age_range, each = length(length_range) ),
               length = rep(length_range, length(age_range)))

# delta length - only age min.ageM + 1 (no obs at age 0 -> do delta length for age previous 1)
pred_dl <- tibble(age = seq(min.ageM + 1, 6))
pred_dl <- pred_dl %>%
  mutate(dl = predict(growth.mod, data.frame(age=pred_dl$age)) - predict(growth.mod, data.frame(age=pred_dl$age-1)))

pred <- pred %>%
  left_join(pred_dl) 

pred <- pred %>%
  mutate(o1 = predict(ogive.mod,pred,type="response"),
         o0 = predict(ogive.mod, data.frame(length=pred$length-dl, age=pred$age-1),type="response"),
         p = if_else(age == 1, o1, (o1-o0)/(1-o0)))

# figure similar to that in Grift et al., 2013
ggplot() +
    geom_line(data = pred, aes(x = length, y = o1)) +
    geom_line(data = pred, aes(x = length, y = o0), linetype = "dashed") +
    geom_line(data = pred, aes(x = length, y = p), linewidth = 1) +
    facet_wrap(~ factor(age))

