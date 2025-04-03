# info ----
# input: 
# - 2_data-analysis
# - 3_model-validation
# output: report

# setup ----
library(tidyverse)
library(TMB)
library(boot)
library(data.table)
library(patchwork)
library(sjPlot)

#### theme
theme_set(theme_bw())

#### pop name
df_pop <- tibble(pop = c("4bc", "7a"),
                 pop_name = factor(c("North Sea", "Irish Sea"),
                                   levels = c("North Sea", "Irish Sea")))

df_pop3 <- tibble(pop = c("4bc", "7a", "8ab"),
                  pop_name = factor(c("North Sea", "Irish Sea", "Bay of Biscay"),
                                    levels = c("North Sea", "Irish Sea", "Bay of Biscay")))

#### function
# transform optimising parameters to actual values
transpar <- function(opt){
    
    itrans <- function(x)(2/(1 + exp(-2*x))-1)
    
    par_names <- names(opt$par)
    pars <- opt$par
    for(i in 1:length(pars)){
        if(grepl("log", par_names[i]) & !grepl("logit", par_names[i])){
            pars[i] <- exp(pars[i])
            par_names[i] <- unlist(strsplit(par_names[i],"_"))[2]
        }
        if(grepl("logit", par_names[i])){
            pars[i] <- boot::inv.logit(pars[i])
            par_names[i] <- unlist(strsplit(par_names[i],"_"))[2]
        }
        if(grepl("trans", par_names[i])){
            pars[i] <- itrans(pars[i])
            par_names[i] <- unlist(strsplit(par_names[i],"_"))[2]
        }
    }
    names(pars) <- par_names
    pars
}

# function to simulate data
# parameter v0 (or l0), h, alpha, chi, Tmat are based on estimates when fitting HMM-Qunice model to North Sea population data
fn_sim <- function(seed = 2018, # for replication
                   v0 = 30,
                   beta = 2/3,
                   h = 115,
                   alpha = 0.35,
                   chi = 0.75,
                   k = -log(chi),
                   
                   precision = 300,
                   sigma_h_ind = 0.2,
                   sigma_Tmat_ind = 0.2,
                   #sigma_alpha_ind = 0.2,
                   
                   rho_h_tmat = -0.7,
                   #rho_h_alpha = -0.3, # if > 0 -> warning "sigma is numerically not positive semidefinite"
                   #rho_tmat_alpha = 0.7,
                   
                   beta_linear_pars_a = 2.2, #a
                   beta_linear_pars_b = 1.8, #Tmat
                   
                   n_ind = 250,
                   max_age = 15) {
    
    set.seed(seed)
    #### true data
    id = seq(1, n_ind)
    
    #max_age <- 20
    age = seq(1,max_age)
    
    #### variables
    Tmat = beta_linear_pars_b
    
    # variance covariance matrix
    SigmaMV = matrix(0, 2, 2)
    diag(SigmaMV) = c(sigma_h_ind^2, sigma_Tmat_ind^2)
    SigmaMV[1,2]  = rho_h_tmat * sigma_h_ind * sigma_Tmat_ind
    SigmaMV[2,1]  = SigmaMV[1,2]
    # SigmaMV[1,3]  = rho_h_alpha * sigma_h_ind * sigma_alpha_ind
    # SigmaMV[3,1]  = SigmaMV[1,3]
    # SigmaMV[2,3]  = rho_tmat_alpha * sigma_Tmat_ind * sigma_alpha_ind
    # SigmaMV[3,2]  = SigmaMV[2,3]
    
    mu_random <- mvtnorm::rmvnorm(n_ind, c(0, 0), SigmaMV)
    #mu_random <- rnorm(n_ind, c(0), sigma_h_ind)
    
    random_h_ind = mu_random[,1]
    random_Tmat_ind = mu_random[,2]
    # random_alpha_ind = mu_random[,3]
    
    h_ind = exp(log(h))*exp(random_h_ind)
    Tmat_ind = Tmat*exp(random_Tmat_ind)
    # alpha_ind = inv.logit(logit(alpha) + random_alpha_ind)
    
    vT_ind = v0 + h_ind*Tmat_ind #length at last juvenile age - mm (eq 9)
    vinf_ind = vT_ind + h_ind*alpha/(1 - chi) #length at infinity - mm (eq9)
    vc_ind = vT_ind + h_ind*(alpha - chi)/(1 - chi) #theoretical "intersect size"; if vc > 0 -> hyperallometric reproduction
    
    
    #### simulated data 
    data_ind = tibble(
        id = rep(id, each = max_age),
        age = rep(age, n_ind),
        h = rep(h_ind, each = max_age),
        alpha = alpha,
        chi = chi,
        k = k,
        Tmat = rep(Tmat_ind, each = max_age),
        vT = rep(vT_ind, each = max_age),
        vinf = rep(vinf_ind, each = max_age),
        vc = rep(vc_ind, each = max_age)) %>%
        # true data
        mutate(vt_true = if_else(age <= Tmat, 
                                 v0 + h*age,
                                 vinf - (vinf - vT)*exp(-k*(age - Tmat)))) %>%
        # simulated data
        mutate(random = rnorm(n_ind*max_age, mean = 0, sd = vt_true/precision),
               vt_sim = if_else(age <= Tmat, 
                                v0 + h*age + random,
                                vinf - (vinf - vT)*exp(-k*(age - Tmat)) + random),
               g = if_else(age <= Tmat, 
                           NA,
                           h/((1-2/3)*(vinf - vc))*(1 - vc/vt_sim)))
    
    ## population data ----
    #samp.size <- c(2324, 2883, 2188, 1247, 662, 312, 172, 89, 44, 30, 16, 12, 6, 6, 2, 2, 1, 2, 1, 1) #10000 individuals
    samp.size <- c(122, 58, 32, 16, 8, 6, 3, 2, 1, 1, 1)
    names(samp.size) <- as.character(seq(5,15))
    
    df_samp.size = tibble(age = rep(seq(5, max_age), times = samp.size),
                          id = seq(1, n_ind))
    
    data_pop <- tibble()
    for(i in 1:max_age) {
        df_samp.size_sub <- df_samp.size %>% filter(age == i)
        df_temp <- data_ind %>% filter(age <= i, id %in% df_samp.size_sub$id) # population age == 1, individual age <= i
        
        data_pop <- bind_rows(data_pop, df_temp)
    }
    
    # add age at capture and only keep fish with age_at_capture >= 3
    data_pop <- data_pop %>% 
        group_by(id) %>%
        mutate(age_at_capture = max(age)) %>%
        filter(age_at_capture >= 5)
    
    p <- ggplot() +
        geom_point(data = data_pop, aes(x = age, y = vt_sim, group = id), alpha = 0.3) +
        geom_line(data = data_pop, aes(x = age, y = vt_true, group = id))
    
    return(list(data_pop, p))
}

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

## sur ----
sur_ca <- read_rds("./input/Solea_soleaCA_BTS+BTS-VIII+DYFS+SNS_in1985till2021_in3a204a4b4c7a7d7e7f7g7h7j28a8b.RDS") 
sur_ca <- sur_ca %>% 
    mutate(id = seq(1,nrow(sur_ca)))

sur_sub <- sur_ca %>% 
    filter(Area_27 %in% c("4.a", "4.b", "4.c", "7.a", "8.a", "8.b"),
           Sex == "F", 
           is.na(Age) == F)

sur_sub <- sur_sub %>% 
    mutate(id = factor(id),
           pop = if_else(Area_27 %in% c("4.a", "4.b", "4.c"), "4bc", 
                         if_else(Area_27 %in% "7.a", "7a", "8ab")),
           age = Age,
           year = Year,
           month = Month,
           cohort = year - age,
           length = LngtClass_mm) %>%
    select(id, pop, age, year, month, cohort, length)
## data rep ----
# results from model
data_rep <- read_rds("output/need to check/individual_h_Tmat_a100_age0_datarep_v2.rds") %>%
    left_join(df_pop3)

# figSI - T and reproductive investment ----
# different reproductive investment (1-p) at different Tmat 
# implicit correlation between Tmat and (1-p) in Quince model

# individual_h_Tmat vs individual_h_Tmat_alpha

# Tmat and reproductive investment correlation in Quince model vs in Mollet model
## Quince model: reproductive investment (1-p) is higher at lower Tmat - assumed in the formula of p = alpha*chi^(age - Tmat)
## Mollet model: reproductive investment c is independent of energy acquisition a (which is positive correlated to growth) 

# simulate 2 fish with different Tmat (fish 1 low, fish 2 high)
# initial par 
v0 = 30
h = 115
alpha = 0.35
chi = 0.75
k = -log(chi)

df <- tibble(fishid = rep(c(1, 2), each = 10),
             Tmat = rep(c(2, 3), each = 10),
             age = rep(seq(1,10), 2)) %>%
    group_by(fishid) %>%
    mutate(vT = v0 + h*Tmat,
           vinf = vT + h*alpha/(1 - chi),
           vc = vT + h*(alpha - chi)/(1 - chi)) %>%
    mutate(vt_sim = if_else(age <= Tmat, 
                            v0 + h*age,
                            vinf - (vinf - vT)*exp(-k*(age - Tmat)))) %>%
    mutate(p = if_else(age <= Tmat, 0, alpha*chi^(age - (Tmat + 1))),
           one_minus_p = 1-p)

# p
ggplot(data = df %>% filter(p > 0), aes(x = age, y = p, linetype = factor(Tmat))) +
    geom_line() +
    labs(x = "Age (years)",
         y = expression(italic(p[t])),
         linetype = "Age at maturity") +
    scale_linetype_manual(values = c("solid", "dashed"))
p1 <- last_plot()

# 1-p
ggplot(data = df %>% filter(p > 0), aes(x = age, y = (1-p), linetype = factor(Tmat))) +
    geom_line() +
    labs(x = "Age (years)",
         y = expression(italic((1-p[t]))),
         linetype = "Age at maturity") +
    scale_linetype_manual(values = c("solid", "dashed"))
p2 <- last_plot()

# size at age
ggplot(data = df, aes(x = age, y = vt_sim, linetype = factor(Tmat))) +
    geom_line() +
    geom_point(alpha = 0.5) +
    labs(x = "Age (years)",
         y = "Length (mm)",
         linetype = "Age at maturity") +
    scale_linetype_manual(values = c("solid", "dashed")) +
    scale_x_continuous(breaks = seq(2,10, 2))
p3 <- last_plot()

(p1 | p2 | p3) + 
    plot_layout(guides = "collect") +
    plot_annotation(tag_level = "A") &
    theme(legend.position = "bottom",
          legend.title.position = "top",
          legend.title.align = 0.5)

ggsave(filename = "./report/figSI_corTmatReproduction.png",
       plot = last_plot(),
       width = 17,
       height = 7,
       unit = "cm",
       scale = 1.8)

ggsave(last_plot(), file = "./report/figSI_corTmatReproduction.pdf",
       device = cairo_pdf,
       width =  17*1.8, height = 9*1.8,
       units = "cm") 

# tabSI - number and proportion per age 5-15 ----
#### distribution by age (age 5-15)
sur_sub_4bc <- sur_sub %>% 
    filter(pop == "4bc",
           age >= 5, age <= 15)

sur_sub_4bc_sum <- sur_sub_4bc %>%
    mutate(n_all = n()) %>%
    group_by(age, n_all) %>%
    summarize(n = n(),
              prop = n/n_all,
              n_250 = round(prop*250,0)) %>%
    unique()

sim_size <- tibble(age = seq(5,15),
                   n = c(122, 58, 32, 16, 8, 6, 3, 2, 1, 1, 1),
                   prop = n/250)

#### summary
## survey
sur_sub_4bc_sum <- sur_sub_4bc_sum %>% ungroup() %>%
    mutate(age_name = paste0("Age ", age)) %>%
    mutate(data = "survey (n = 3367)") 

# n
sur_n <- sur_sub_4bc_sum %>% 
    select(data, age_name, n) %>%
    mutate(type = "count")
sur_n <- sur_n %>% pivot_wider(names_from = age_name,
                         values_from = n)
# prop
sur_prop <- sur_sub_4bc_sum %>%
    select(data, age_name, prop) %>%
    #mutate(prop = round(prop, 3)) %>%
    mutate(type = "prop")
sur_prop <- sur_prop %>% pivot_wider(names_from = age_name,
                         values_from = prop)

#### sim data
sim_size <- sim_size %>%
    mutate(age_name = paste0("Age ", age)) %>%
    mutate(data = "simulated (n = 250)")

# n
sim_n <- sim_size %>% 
    select(data, age_name, n)  %>%
    mutate(type = "count")
sim_n <- sim_n %>% pivot_wider(names_from = age_name,
                      values_from = n)
# prop
sim_prop <- sim_size %>%
    select(data, age_name, prop) %>%
    #mutate(prop = round(prop, 3)) %>%
    mutate(type = "prop")
sim_prop <- sim_prop %>% pivot_wider(names_from = age_name,
                         values_from = prop)

# data_all
sur_sim <- rbind(sur_n, sur_prop, sim_n, sim_prop)
sur_sim
tab_df(sur_sim,
       file = "./report/tabSI_sur_sim_count_prop_age515.html",
       digits = 3)

# figSI - backcal_len vs survey data ----
## backcal ----
otl_backcal <- readRDS("./input/otl_full_backcal.rds")
otl_backcal_sub <- otl_backcal %>% filter(age <= 15,
                                          aac >= 5, 
                                          !id %in% c("sol_fab_1375")) %>%
    mutate(length = backcal_len) %>%
    select(pop, age, length) %>%
    mutate(age_name = factor(paste0("Age ", age),
                             level = paste0("Age ", seq(1,15))))

## sur ----
sur_sub_sub <- sur_sub %>%
    filter(age >= 1, age <= 15) %>%
    select(pop, age, length) %>%
    mutate(age_name = factor(paste0("Age ", age),
                             level = paste0("Age ", seq(1,15))))

## backcal vs sur plot ----
ggplot() +
    geom_density(data = otl_backcal_sub, aes(x = length, fill = "Back-calculated data"), alpha = 0.5) +
    geom_density(data = sur_sub_sub, aes(x = length, fill = "Survey data"), alpha = 0.5) +
    facet_wrap( ~ age_name, ncol = 3) +
    scale_fill_manual(values = c("white", "grey")) +
    labs(x = "Length (mm)",
         y = "Density",
         fill = NULL) +
    theme(legend.position = "bottom")

ggsave(filename = "./report/figSI_backcal_vs_sur.png",
       plot = last_plot(),
       width = 17,
       height = 17,
       unit = "cm",
       scale = 1.2)

ggsave(last_plot(), file = "./report/figSI_backcal_vs_sur.pdf",
       device = cairo_pdf,
       width =  17*1.2, height = 17*1.2,
       units = "cm") 

# plot by pop
ggplot() +
    geom_density(data = otl_backcal_sub %>% filter(pop == "4bc"), aes(x = length, fill = "Back-calculated data"), alpha = 0.5) +
    geom_density(data = sur_sub_sub  %>% filter(pop == "4bc"), aes(x = length, fill = "Survey data"), alpha = 0.5) +
    facet_wrap( ~ age_name, ncol = 3) +
    scale_fill_manual(values = c("white", "grey")) +
    labs(x = "Length (mm)",
         y = "Density",
         fill = NULL) +
    theme(legend.position = "bottom")

ggplot() +
    geom_density(data = otl_backcal_sub %>% filter(pop == "7a"), aes(x = length, fill = "Back-calculated data"), alpha = 0.5) +
    geom_density(data = sur_sub_sub  %>% filter(pop == "7a"), aes(x = length, fill = "Survey data"), alpha = 0.5) +
    facet_wrap( ~ age_name, ncol = 3) +
    scale_fill_manual(values = c("white", "grey")) +
    labs(x = "Length (mm)",
         y = "Density",
         fill = NULL) +
    theme(legend.position = "bottom")

# tabSI - count and proportion mature obs data ----
data_mat <- data_obs %>% filter(pop %in% c("4bc", "7a", "8ab"), age <= 6, month <= 6) %>%
    left_join(df_pop3)

data_mat_count_all <- data_mat %>%
    mutate(type = "Count") %>%
    group_by(pop_name, type) %>%
    summarize(total = n()) %>%
    mutate(total = sprintf("%.0f", total))

data_mat_prop_all <- data_mat %>%
    mutate(type = "Proportion mature") %>%
    group_by(pop_name, type) %>%
    mutate(n_all = n()) %>%
    filter(maturity == 1) %>%
    mutate(n = n(),
           total = round(n/n_all,2)) %>%
    select(pop_name, type, total) %>%
    unique() %>%
    mutate(total = sprintf("%.2f", total))

data_mat_count <- data_mat %>%
    mutate(type = "Count") %>%
    group_by(pop_name, type, age) %>%
    summarize(n = n()) %>%
    mutate(n = sprintf("%.0f", n)) %>%
    pivot_wider(names_from = age, values_from = n)

data_mat_prop <- data_mat %>%
    group_by(pop_name, age) %>%
    mutate(n_all = n()) %>%
    mutate(type = "Proportion mature") %>%
    group_by(pop_name, type, age, maturity) %>%
    filter(maturity == 1) %>%
    summarize(n = n(),
              prop_mat = round(n/n_all,2)) %>%
    ungroup() %>%
    select(pop_name, type, age, prop_mat) %>%
    unique() %>%
    arrange(age) %>%
    mutate(prop_mat = sprintf("%.2f", prop_mat)) %>%
    pivot_wider(names_from = age, values_from = prop_mat)

data_mat_sum <- bind_rows(data_mat_count, data_mat_prop) %>%
    left_join(bind_rows(data_mat_count_all, data_mat_prop_all)) %>%
    arrange(pop_name, type)

tab_df(data_mat_sum,
       col.header = c("Population", "Count or proportion mature",
                      "Age 1", "Age 2", "Age 3", "Age 4", "Age 5", "Age 6", "Total"),
       file = "./report/tabSI_dataMaturity.html")

# appendix - unrealistic estimates 8ab ----
df_nll <- read_rds(file.path("./output", paste0("backcal_", "8ab", "_individual_h_Tmat_prec_a100_age0_rep1000_v2_nll.rds")))
list_opt <- read_rds(file.path("./output", paste0("backcal_", "8ab", "_individual_h_Tmat_prec_a100_age0_rep1000_v2_opt.rds")))

# best model
df_nll_best <- df_nll %>% 
    #filter(message == "relative convergence (4)") %>%
    filter(nll == min(nll)) %>%
    slice(1)
m_best <- df_nll_best$m_name
print(m_best)

# estimates - best model OK
transpar(list_opt[[m_best]]) 

# figSI - pred vs obs all pops ----
ggplot() +
    geom_point(data = data_rep, aes(x = age, y = vt_sim, group = id), alpha = 0.1) +
    geom_line(data = data_rep, aes(x = age, y = mu, group = id), alpha = 0.1) +
    labs(x = "Age (years)",
         y = "Length (mm)") +
    facet_grid(~ pop_name)

ggsave(filename = "./report/figSI_predVsObs_all.png",
       plot = last_plot(),
       width = 17,
       height = 9,
       unit = "cm",
       scale = 1.5)

ggsave(last_plot(), file = "./report/figSI_predVsObs_all.pdf",
       device = cairo_pdf,
       width =  17*1.5, height = 9*1.5,
       units = "cm") 

# figSI - g_prime ----
sdrep_sum_all <- tibble()
for(p in c("4bc", "7a")) {
    sdrep_sum_temp <- read_rds(paste0("output/individual_h_Tmat_a100_age0_sdrepsum_", p, "_v2.rds")) %>%
        mutate(pop = p) 
    sdrep_sum_all <- bind_rows(sdrep_sum_all, sdrep_sum_temp) 
}
sdrep_sum_all <- sdrep_sum_all %>%
    left_join(df_pop)

# df_g
df_g <- sdrep_sum_all %>%
    filter(pars == "g_prime_mean") %>%
    mutate(age = rep(seq(1,15), 2),
           g_prime_est = Estimate,
           g_prime_se = `Std. Error`)

list_plot <- list()
for(p in c("4bc", "7a")) {
    ggplot(data = df_g %>% filter(age >= 2, pop == p), aes(x = age, g_prime_est)) +
        geom_line() +
        geom_ribbon(aes(ymin = g_prime_est - g_prime_se,
                        ymax = g_prime_est + g_prime_se),
                    alpha = 0.5) +
        facet_wrap(~ pop_name) +
        labs(x = "Age (years)",
             y = expression(italic("g'")))
    list_plot[[p]] <- last_plot()
}

(list_plot[["4bc"]] | list_plot[["7a"]]) +
    plot_annotation(tag_levels = "A")

ggsave(filename = "./report/figSI_g_prime.png",
       plot = last_plot(),
       width = 17,
       height = 7,
       unit = "cm",
       scale = 1.2)

ggsave(last_plot(), file = "./report/figSI_g_prime.pdf",
       device = cairo_pdf,
       width =  17*1.2, height = 7*1.2,
       units = "cm") 
