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
df_pop <- tibble(pop = c("4bc", "7a", "8ab"),
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
## input data ---- 
### obs ----
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


### backcal otl ----
data_backcal <- read_rds("./input/otl_full_backcal.rds") 

data_backcal_sub <- data_backcal %>%
    filter(aac >= 5, 
           age <= 15,
           !id %in% c("sol_fab_1375")) %>%
    rename(length = backcal_len) %>%
    mutate(source = "otl") %>%
    select(id, pop, year, cohort, age, length, source, datasource) 

### all data ----
data_all <- bind_rows(filter(data_obs_sub, !id %in% data_backcal_sub$id), # not include obs with backcal to avoid duplication
                      data_backcal_sub)

## output data from the model ----
data_rep <- read_rds("output/individual_h_Tmat_a100_age0_datarep_v2.rds") %>%
    left_join(df_pop)

# calculate individual p, 1-p, g_prime
data_rep <- data_rep %>%
    mutate(p = if_else(age <= Tmat, 1, alpha*chi^(age - (Tmat + 1))),
           one_minus_p = 1-p) %>%
    mutate(vc = v0 + h*Tmat + h*(alpha - chi)/(1 - chi)) %>%
    mutate(g_prime = if_else(age < Tmat, 0, (3*h/(vinf - vc))*(1 - (vc/mu))))

# summary ----
# n per archive per pop 
data_sub <- data_backcal_sub %>% 
    select(pop, datasource, id) %>%
    unique()

# add info archive
list_ifremer <- unique(filter(data_sub, str_detect(id, "RE|AL|CO") == TRUE)$id)
data_sub <- data_sub %>% mutate(archive = if_else(id %in% list_ifremer, "IFREMER", datasource))

# n per pop
table(data_sub$pop)

# n per archive per pop
data_sub %>% 
    group_by(pop, archive) %>%
    summarize(n = n())

# fig1 - concept HMM model ----
# report/fig1_dependence structure HMM.pptx

# fig2 - concept Quince model ----
fn_sim1 <- function(seed = 2018, # for replication
                   v0 = 30,
                   h = 115,
                   Tmat = 2,
                   alpha = 0.35,
                   chi = 0.75,
                   k = -log(chi),
                   max_age = 10) {
    
    set.seed(seed)

    age = seq(0,max_age)
    
    # derived parameters
    vT = v0 + h*Tmat #length at last juvenile age - mm (eq 9)
    vinf = vT + h*alpha/(1 - chi) #length at infinity - mm (eq9)
    vc = vT + h*(alpha - chi)/(1 - chi) #theoretical "intersect size"; if vc > 0 -> hyperallometric reproduction
    
    data <- tibble(
        v0 = v0,
        h = h,
        Tmat = Tmat,
        alpha = alpha,
        chi = chi,
        age = age,
        len = if_else(age <= Tmat, 
                      v0 + h*age,
                      vinf - (vinf - vT)*exp(-k*(age - Tmat))),
        len_growth = v0 + h*age
    )
    
    return(data)
}

data <- fn_sim1()

# v0, h, Tmat
ggplot() +
    geom_line(data = data %>% filter(age <= 1), aes(x = age, y = len), color = "grey") +
    geom_point(data = data %>% filter(age <= 1), aes(x = age, y = len), color = "grey") +
    geom_line(data = data %>% filter(age >= 1), aes(x = age, y = len)) +
    geom_point(data = data %>% filter(age >= 1), aes(x = age, y = len)) +
    scale_x_continuous(breaks = seq(0,10)) +
    theme(panel.grid = element_blank()) +
    geom_segment(aes(x = 2, xend = 2, y = -Inf, yend = 30 + 2*115), linetype = "dashed") +
    labs(x = expression("Age" ~ italic("t") ~ "(years)"),
         y = expression("Length" ~ italic("l")[italic("t")] ~ "(mm)")) +
    annotate("text", x = 2.2, y = 30, label = expression(italic("T")), hjust = 0, size = 3) +
    annotate("text", x = 0.2, y = 30, label = expression(italic("l")['0']), hjust = 0, size = 3) +
    annotate("text", x = 1.2, y = 202.5, label = expression(italic("h")), hjust = 0, size = 3, angle = 65)

p1 <- last_plot()

# alpha
data_alpha025 <- fn_sim1(alpha = 0.25)
data_alpha045 <- fn_sim1(alpha = 0.45)
ggplot() +
    geom_line(data = data %>% filter(age <= 1), aes(x = age, y = len), color = "grey") +
    geom_point(data = data %>% filter(age <= 1), aes(x = age, y = len), color = "grey") +
    geom_line(data = data %>% filter(age >= 1), 
              aes(x = age, y = len, linetype = "alpha = 0.35")) +
    geom_point(data = data %>% filter(age >= 1), aes(x = age, y = len)) +
    geom_line(data = data_alpha025 %>% filter(age >= 1), 
              aes(x = age, y = len, linetype = "alpha = 0.25")) +
    geom_line(data = data_alpha045 %>% filter(age >= 1), 
              aes(x = age, y = len, linetype = "alpha = 0.45")) +
    scale_x_continuous(breaks = seq(0,10)) +
    scale_linetype_manual(values = c("dashed", "solid", "dotted"),
                          labels = c(expression(italic("\u03B1")~"= 0.25"), 
                                     expression(italic("\u03B1")~"= 0.35"), 
                                     expression(italic("\u03B1")~"= 0.45"))) +
    theme(panel.grid = element_blank()) +
    labs(x = expression("Age" ~ italic("t") ~ "(years)"),
         y = expression("Length" ~ italic("l")[italic("t")] ~ "(mm)"),
         linetype = NULL) +
    theme(legend.position = c(0.8, 0.3))
p2 <- last_plot()

# chi
data_chi065 <- fn_sim1(chi = 0.65)
data_chi085 <- fn_sim1(chi = 0.85)
ggplot() +
    geom_line(data = data %>% filter(age <= 1), aes(x = age, y = len), color = "grey") +
    geom_point(data = data %>% filter(age <= 1), aes(x = age, y = len), color = "grey") +
    geom_line(data = data %>% filter(age >= 1), 
              aes(x = age, y = len, linetype = "chi = 0.75")) +
    geom_point(data = data %>% filter(age >= 1), aes(x = age, y = len)) +
    geom_line(data = data_chi065 %>% filter(age >= 1), 
              aes(x = age, y = len, linetype = "chi = 0.65")) +
    geom_line(data = data_chi085 %>% filter(age >= 1), 
              aes(x = age, y = len, linetype = "chi = 0.85")) +
    scale_x_continuous(breaks = seq(0,10)) +
    scale_linetype_manual(values = c("dashed", "solid", "dotted"),
                          labels = c(expression(italic("χ")~"= 0.65"), 
                                     expression(italic("χ")~"= 0.75"), 
                                     expression(italic("χ")~"= 0.85"))) +
    theme(panel.grid = element_blank()) +
    labs(x = expression("Age" ~ italic("t") ~ "(years)"),
         y = expression("Length" ~ italic("l")[italic("t")] ~ "(mm)"),
         linetype = NULL) +
    theme(legend.position = c(0.8, 0.3))
p3 <- last_plot()

layout <- '
AAB
AAB
'
(p1 | (p2/p3)) + 
    plot_layout(design = layout) +
    plot_annotation(tag_levels = "A")
    
ggsave(filename = "./report/fig2_concept_QuinceModel.png",
       plot = last_plot(),
       width = 17,
       height = 9,
       unit = "cm",
       scale = 1.5)

ggsave(last_plot(), file = "./report/fig2_concept_QuinceModel.pdf",
       device = cairo_pdf,
       width =  17*1.5, height = 9*1.5,
       units = "cm") 

# tab1 - in text ----
# tab2 - in text ----
# fig3 - simulation results relative errors ----
n_sample = c(50, 150, 250) 
df_nll_best <- tibble()
for(n in n_sample) {
    print(n)
    
    df_nll <- read_rds(paste0("./output/backcal_hmm_simulation_h_Tmat_v2_nll_size", n, ".rds"))
    list_opt <- read_rds(paste0("./output/backcal_hmm_simulation_h_Tmat_v2_opt_size", n, ".rds"))
    list_obj <- read_rds(paste0("./output/backcal_hmm_simulation_h_Tmat_v2_obj_size", n, ".rds"))
    
    # get best model
    df_nll_best_temp <- df_nll %>% 
        filter(message == "relative convergence (4)") %>%
        group_by(ite_name) %>%
        filter(nll == min(nll, na.rm = T)) %>%
        #slice(1) %>%
        mutate(sample_size = n)
    
    df_nll_best <- bind_rows(df_nll_best_temp, df_nll_best)
}


n_sample = c(50, 150, 250) 
df <- tibble()
for(n in n_sample) {
    print(n)
    
    # load model results
    df_nll <- read_rds(paste0("./output/backcal_hmm_simulation_h_Tmat_v2_nll_size", n, ".rds"))
    list_opt <- read_rds(paste0("./output/backcal_hmm_simulation_h_Tmat_v2_opt_size", n, ".rds"))
    list_obj <- read_rds(paste0("./output/backcal_hmm_simulation_h_Tmat_v2_obj_size", n, ".rds"))
    
    df_nll_best <- df_nll %>% 
        filter(message == "relative convergence (4)") %>%
        group_by(ite_name) %>%
        filter(nll == min(nll, na.rm = T)) %>%
        slice(1) %>%
        mutate(sample_size = n) %>%
        mutate(ite = as.numeric(substring(ite_name, 5, 10)))

    for(i in 1:nrow(df_nll_best)) {
        m_best <- df_nll_best$m_name[i]
        
        opt <- list_opt[[m_best]]
        obj <- list_opt[[m_best]]
        
        df_temp <- tibble(
            m_name = m_best,
            nll = opt$value,
            par_name = factor(c("Tmat", "v0", "h", "alpha", "chi",
                                "prec", "sd_h", "sd_tmat", "rho_h_tmat"),
                              levels = c("v0", "h", "Tmat", "alpha", "chi", "sd_h", "sd_tmat", "rho_h_tmat", "prec")),
            par_true = c(1.8, 30, 115, 0.35, 0.75,
                         300, 0.2, 0.2, -0.7),
            par_est = transpar(opt),
            bias = (par_est - par_true)/par_true*100,
            is_convergence = opt$convergence,
            sample_size = n
        )
        
        df <- bind_rows(df, df_temp)
        
    }
}

df_label <- tibble(par_name = factor(c("v0", "h", "Tmat", "alpha", "chi", "sd_h", "sd_tmat", "rho_h_tmat", "prec"),
                                     levels = c("v0", "h", "Tmat", "alpha", "chi", "sd_h", "sd_tmat", "rho_h_tmat", "prec")),
                   par_label = factor(c("italic(l)['0']", 
                                        "italic('h')", 
                                        "italic('T')",  
                                        "italic('α')", 
                                        "italic('χ')",
                                        "italic(σ['h'])",
                                        "italic(σ['T'])",
                                        "italic(ρ['(h,T)'])",
                                        "italic('p')"),
                                      levels = c("italic(l)['0']", 
                                                 "italic('h')", 
                                                 "italic('T')",
                                                 "italic('α')", 
                                                 "italic('χ')",
                                                 "italic(σ['h'])",
                                                 "italic(σ['T'])",
                                                 "italic(ρ['(h,T)'])",
                                                 "italic('p')")))

df <- df %>%
    left_join(df_label)

ggplot(data = df %>% filter(abs(bias) <= 100), aes(x = bias)) +
    geom_rect(aes(xmin = -10, xmax = 10, ymin = -Inf, ymax = Inf), 
              color = "grey50", linetype = "dashed", alpha = 0) +
    geom_histogram(bins = 150, position = "identity", alpha = 0.8, aes(fill = factor(sample_size))) +
    facet_wrap(~ par_label, labeller = label_parsed) +
    coord_cartesian(xlim = c(-20,20)) +
    labs(x = "Bias (%)",
         y = "Count",
         fill = "Sample size") +
    theme(legend.position = "bottom",
          legend.title.position = "top",
          legend.title.align = 0.5) + 
    scale_fill_brewer(palette = "Dark2")

ggsave(filename = "./report/fig3_simulation_bias.png",
       plot = last_plot(),
       width = 17,
       height = 12,
       unit = "cm",
       scale = 1.2)

ggsave(last_plot(), file = "./report/fig3_simulation_bias.pdf",
       device = cairo_pdf,
       width =  17*1.2, height = 12*1.2,
       units = "cm") 

## summary 
# 10% bias
df_sum_10 <- df %>% 
    mutate(bias_lt10 = if_else(abs(bias) < 10, 1, 0)) %>%
    group_by(sample_size, par_name, bias_lt10) %>%
    summarize(n = n())

df_sum_10 <- df_sum_10 %>% pivot_wider(names_from = bias_lt10, values_from = n) %>%
    rename(bias_lt10 = `1`,
           bias_gte10 = `0`) %>%
    select(par_name, bias_lt10, bias_gte10) %>%
    arrange(par_name)

df_sum_10 %>% 
    select(sample_size, par_name, bias_lt10) %>%
    pivot_wider(names_from = sample_size, values_from = bias_lt10)

## summary indidual relative errors ----
# get individual estimate
# note: this stage take quite some time so it is better to do the sample size one by one
## save, and the reload the results

n_sample = c(50, 150, 250) 
for(n in n_sample) {
    print(n)
    
    # load model results
    df_nll <- read_rds(paste0("./output/backcal_hmm_simulation_h_Tmat_v2_nll_size", n, ".rds"))
    list_opt <- read_rds(paste0("./output/backcal_hmm_simulation_h_Tmat_v2_opt_size", n, ".rds"))
    list_obj <- read_rds(paste0("./output/backcal_hmm_simulation_h_Tmat_v2_obj_size", n, ".rds"))
    
    df_nll_best <- df_nll %>% 
        filter(message == "relative convergence (4)") %>%
        group_by(ite_name) %>%
        filter(nll == min(nll, na.rm = T)) %>%
        slice(1) %>%
        mutate(sample_size = n) %>%
        mutate(ite = as.numeric(substring(ite_name, 5, 10)))
    
    df_ind <- tibble()
    for (i in df_nll_best$ite) {
        print(paste0("processing ", i))
        
        # simulate data again
        set.seed(i)
        #### sampled data
        data_all <- fn_sim(i)[[1]]
        data <- data_all %>% filter(id %in% sample(unique(data_all$id), n))
        data <- data %>%
            group_by(id) %>%
            mutate(year_adult = sum(!is.na(g))) #number of year with non-na g
        
        # get obj
        ite_name_i = paste0("ite_", i)
        m_name_i <- filter(df_nll_best, ite_name == ite_name_i)$m_name
        
        obj <- list_obj[[m_name_i]]
        dyn.load(TMB::dynlib("./R/quince_hmm_individual_h_Tmat_prec"))
        
        # get estimates
        # note: rep will be faster but h was not included in rep
        sdrep <- TMB::sdreport(obj)
        sdrep_sum <- summary(sdrep)
        sdrep_sum <- as.data.frame(sdrep_sum) %>% mutate(pars = row.names(sdrep_sum))
        
        # merge estimates and true values
        data_sub <- data %>% select(id, h, Tmat)
        data_est <- tibble(ite_name = ite_name_i,
                           id = unique(data$id),
                           h_est = filter(sdrep_sum, pars == "h")$Estimate,
                           Tmat_est = filter(sdrep_sum, pars == "Tmat")$Estimate) %>%
            left_join(data_sub, by = "id") %>%
            unique()
        
        df_temp <- data_est %>%
            mutate(h_error = abs(h_est - h)/h*100,
                   Tmat_error = abs(Tmat_est - Tmat)/Tmat*100) %>%
            group_by(ite_name) %>%
            summarize(h_error_mean = mean(h_error),
                      Tmat_error_mean = mean(Tmat_error))
        df_ind <- bind_rows(df_ind, df_temp)
        
        # save file
        write_rds(df_ind, paste0("./output/backcal_hmm_simulation_h_Tmat_v2_size", n, "_results_ind.rds"))
        
    }
}

# load df_ind
df_ind_all <- tibble()
for(n in n_sample) {
    print(n)
    df_ind <- read_rds(paste0("./output/backcal_hmm_simulation_h_Tmat_v2_size", n, "_results_ind.rds")) %>%
        mutate(sample_size = n)
    df_ind_all <- bind_rows(df_ind_all, df_ind)
}

# 1 and 5% bias
df_ind_all <- df_ind_all %>% 
    mutate(h_bias_lt1 = if_else(h_error_mean < 1, 1, 0),
           Tmat_bias_lt1 = if_else(Tmat_error_mean < 1, 1, 0),
           h_bias_lt5 = if_else(h_error_mean < 5, 1, 0),
           Tmat_bias_lt5 = if_else(Tmat_error_mean < 5, 1, 0))

df_ind_all %>%
    group_by(sample_size) %>%
    summarize(h_bias_lt1 = sum(h_bias_lt1),
              Tmat_bias_lt1 = sum(Tmat_bias_lt1),
              h_bias_lt5 = sum(h_bias_lt5),
              Tmat_bias_lt5 = sum(Tmat_bias_lt5),
              n = n())

# fig4 - pred vs obs ----

p_code = "4bc"
list_sample <- c() 
for(p_code in c("4bc", "7a", "8ab")) {
    set.seed(100)
    sample_temp <- sample(filter(data_rep, pop == p_code)$id, 5)
    list_sample <- c(list_sample, sample_temp)
}

list_main = list()
list_sub = list()
for(p_code in c("4bc", "7a", "8ab")) {
    ggplot(data = data_rep %>% filter(pop == p_code)) +
        geom_point(aes(x = age, y = backcal_len, group = id), alpha = 0.5) +
        geom_line(aes(x = age, y = mu, group = id), alpha = 0.5) +
        labs(x = "Age (year)",
             y = "Length (mm)") +
        # labs(x = NULL,
        #      y = "Length (mm)") +
        facet_grid(. ~ pop_name) +
        theme(legend.position = "none") +
        xlim(1,15) +
        ylim(range(data_rep$backcal_len))
    list_main[[p_code]] = last_plot()
    
    ggplot(data = data_rep %>% filter(pop == p_code, id %in% list_sample)) +
        geom_point(aes(x = age, y = backcal_len, group = id), alpha = 0.5) +
        geom_line(aes(x = age, y = mu, group = id), alpha = 1) +
        labs(x = "Age (year)",
             y = "Length (mm)") +
        facet_grid(. ~ id_int) +
        theme(legend.position = "none") +
        xlim(1,15) +
        ylim(range(data_rep$backcal_len))
    list_sub[[p_code]] = last_plot()
}

layout = "
AAABBBCCC
AAABBBCCC
AAABBBCCC
DDDEEEFFF
"

p_all <- (list_main[["4bc"]] + list_main[["7a"]] + list_main[["8ab"]] + 
        list_sub[["4bc"]] + list_sub[["7a"]] + list_sub[["8ab"]]) + 
    plot_layout(design = layout)  +
    plot_annotation(tag_levels = "A")

ggsave(filename = "./report/fig4_predVsObs.png",
       plot = p_all,
       width = 17,
       height = 9,
       unit = "cm",
       scale = 2)

ggsave(last_plot(), file = "./report/fig4_predVsObs.pdf",
       device = cairo_pdf,
       width =  17*2, height = 9*2,
       units = "cm") 

## summary - correlation pred vs obs ----
df_r2 <- tibble()
for(p in c("4bc", "7a", "8ab")) {
    lm <- lm(mu ~ vt_sim, data = data_rep %>% filter(pop == p))
    df_temp <- tibble(pop = p,
                      R2 = summary(lm)$r.squared)
    df_r2 <- bind_rows(df_r2, df_temp)
}
df_r2

# fig5 - PMRN ----
pred_lp50 <- read_rds("output/lp50_maturity.rds") %>%
    left_join(df_pop)
pred_lp50_growth <- read_rds("output/lp50_growth.rds") %>%
    left_join(df_pop)
# pred_lp50_boot <- read_rds("output/lp50_maturity_boot.rds") %>%
#     left_join(df_pop)

# summry length at age
data_len <- data_all %>% filter((source == "obs" & month <= 6) | source == "otl") %>%
    filter(pop %in% c("4bc", "7a")) %>%
    group_by(pop, age) %>%
    summarize(len_mean = mean(length),
              len_lower = quantile(length, probs = c(0.025)),
              len_upper = quantile(length, probs = c(0.975))) %>%
    left_join(df_pop)

max_age = 3
ggplot() +
    geom_line(data = data_len %>% filter(age <= max_age), aes(x = age, y = len_mean, linetype = "3"), color = "grey") +
    geom_ribbon(data = data_len %>% filter(age <= max_age),
                aes(x = age, y = len_mean, ymin = len_lower, ymax = len_upper),
                alpha = 0.1) +
    geom_line(data = pred_lp50 %>% filter(age <= max_age), aes(x = age, y = lp50, linetype = "2")) +
    geom_line(data = pred_lp50_growth %>% filter(age <= max_age, lp50 > 0), aes(x = age, y = lp50, linetype = "1")) +
    facet_grid(~ pop_name) +
    scale_x_continuous(breaks = seq(1,6)) +
    scale_linetype_manual(values = c("solid", "dashed", "solid"),  
                          labels = c(
                              expression("Growth-model-based" ~ 'L'['P50']), 
                              expression("Maturity-staging-based" ~ 'L'['P50']),
                              "Length"
                          )) +
    labs(x = "Age (year)",
         y = "Length (mm)",
         linetype = NULL) +
    theme(legend.position = "bottom",
          legend.title.position = "top",
          legend.title.align = 0.5) 

ggsave(filename = "./report/fig5_lp50.png",
       plot = last_plot(),
       width = 17,
       height = 9,
       unit = "cm",
       scale = 1.2)

ggsave(last_plot(), file = "./report/fig5_lp50.pdf",
       device = cairo_pdf,
       width =  17*1.2, height = 9*1.2,
       units = "cm") 

# tab3 - summary par estimates ----
sdrep_sum_all <- tibble()
for(p in c("4bc", "7a")) {
    sdrep_sum <-  read_rds(paste0("output/individual_h_Tmat_a100_age0_sdrepsum_", p, "_v2.rds")) %>%
        mutate(pop = p)
    
    sdrep_sum_all <- bind_rows(sdrep_sum_all, sdrep_sum)
}

unique(sdrep_sum_all$pars)

df_par <- sdrep_sum_all %>%
    filter(pars %in% c("v0", "h_mean", "Tmat_mean", "alpha", "chi", "sd_h", "sd_tmat", "rho_h_tmat", "prec", "vinf_mean")) %>%
    unique()

df_par <- df_par %>%
    mutate(pars = factor(pars, level = c("v0", "h_mean", "Tmat_mean", "alpha", "chi", "sd_h", "sd_tmat", "rho_h_tmat", "prec", "vinf_mean"))) %>%
    mutate(est_se = sprintf("%.2f (%.3f)", Estimate, `Std. Error`)) %>%
    select(pars, pop, est_se) %>%
    pivot_wider(names_from = pop, values_from = est_se)

col_header = c("Parameters", "North Sea", "Irish Sea")

tab_df(df_par,
       col.header = col_header,
       file = "./report/tab2_sdreport.html")

## difference 4bc vs 7a ----
# assign 7a as first level 
data_rep_sub <- data_rep %>% 
    filter(pop %in% c("4bc", "7a")) %>%
    select(id, pop, h, Tmat, vinf) %>%
    unique() %>%
    mutate(pop = factor(pop, levels = c("7a", "4bc")))

data_rep_sub2 <- data_rep %>% 
    filter(pop %in% c("4bc", "7a")) %>%
    filter(one_minus_p >0) %>%
    mutate(pop = factor(pop, levels = c("7a", "4bc")))

library(betareg) # for betareg

lm1 <- lm(h ~ pop, data = data_rep_sub)
lm2 <- lm(Tmat ~ pop, data = data_rep_sub)
lm3 <- lm(vinf ~ pop, data = data_rep_sub)
beta1 <- betareg(one_minus_p ~ 1 + age + pop, data =  data_rep_sub2) #link = "logit"

pred <- data_rep_sub2 %>%
    mutate(pred = predict(beta1, newdata = data_rep_sub2, type = "response"))

ggplot(data = pred, aes(x = pred, fill = pop)) +
    geom_histogram(bins = 100)
ggplot(data = pred, aes(x = one_minus_p, fill = pop)) +
    geom_histogram(bins = 100)
ggplot(data = pred, aes(x = one_minus_p, y = pred)) +
    geom_point()

ggplot() +
    geom_point(data = pred, aes(x = age, y = one_minus_p, color = "obs")) +
    geom_point(data = pred, aes(x = age, y = pred, color = "pred")) 

ggplot() +
    geom_point(data = pred, aes(x = age, y = pred, color = pop)) 


tab_model(lm1, lm2, lm3, beta1,
          show.se = TRUE, 
          collapse.se = TRUE,
          show.ci = FALSE,
          show.p = TRUE,
          p.style = "stars", 
          string.est = "Estimate (SE)",
          show.r2 = FALSE,
          show.obs = FALSE,
          dv.labels = c("Juvenile growth rate h", 
                        "Age at maturity T",
                        "Maximum size l∞",
                        "Reproductive investment (1-pt)"
                        ),
          file = "./report/tableSI_traits_NSvsIS.html")

# fig6 - distribution h, Tmat, vinf ----

ind_pars <- data_rep %>%
    filter(pop != "8ab") %>%
    filter(age == 5) %>%
    group_by(pop) %>%
    summarize(min_h = min(h),
              max_h = max(h),
              min_Tmat = min(Tmat),
              max_Tmat = max(Tmat),
              min_1_p = min(one_minus_p),
              max_1_p = max(one_minus_p),
              mean_1_p = mean(one_minus_p),
              min_vinf = min(vinf),
              max_vinf = max(vinf)) %>%
    mutate(across(where(is.numeric), round, 2)) %>%
    left_join(df_pop)

ggplot(data = data_rep %>% filter(pop != "8ab"), aes(x = h, fill = pop_name)) +
    geom_density(bins = 50, position = "identity", alpha = 0.5) + 
    scale_fill_manual(values = c("grey", "white")) +
    labs(x = expression(italic("h")~"(mm/year)"),
         y = "Density",
         fill = "Population")
p1 <- last_plot()

ggplot(data = data_rep %>% filter(pop != "8ab"), aes(x = Tmat, fill = pop_name)) +
    geom_density(bins = 50, position = "identity", alpha = 0.5) +
    scale_fill_manual(values = c("grey", "white")) +
    labs(x = expression(italic("T")~"(year)"),
         y = "Density",
         fill = "Population")
p2 <- last_plot()

ggplot(data = data_rep %>% filter(pop != "8ab"), aes(x = vinf, fill = pop_name)) +
    geom_density(bins = 50, position = "identity", alpha = 0.5) + 
    scale_fill_manual(values = c("grey", "white")) +
    labs(x = expression(italic("l"[infinity])~"(mm)"),
         y = "Density",
         fill = "Population")
p3 <- last_plot()

(p1 | p2 | p3) + 
    plot_layout(axis_title = "collect", guides = "collect") & 
    theme(legend.position = "bottom")

ggsave(filename = "./report/fig6_distribution_h_Tmat_vinf.png",
       plot = last_plot(),
       width = 17,
       height = 7,
       unit = "cm",
       scale = 1)

ggsave(last_plot(), file = "./report/fig6_distribution_h_Tmat_vinf.pdf",
       device = cairo_pdf,
       width =  17*1, height = 7*1,
       units = "cm") 

# summary Tmat = 1
data_rep_sub <- data_rep %>%
    filter(pop != "8ab")

df_T1 <- data_rep_sub %>% 
    filter(round(Tmat,2) >= 0.99, round(Tmat,2) <= 1.0) %>%
    select(id_int, pop, Tmat) %>%
    unique()
table(df_T1$pop)
# total fish: 4bc 7a (993 + 610) - 1603
round(23/993*100,2) #2.32% 
round(4/610*100,2)  #0.66%

# fig7 - 1-p and g_prime ----
sdrep_sum_all <- tibble()
for(p in c("4bc", "7a")) {
    sdrep_sum_temp <- read_rds(paste0("output/individual_h_Tmat_a100_age0_sdrepsum_", p, "_v2.rds")) %>%
        mutate(pop = p) 
    sdrep_sum_all <- bind_rows(sdrep_sum_all, sdrep_sum_temp) 
}
sdrep_sum_all <- sdrep_sum_all %>%
    left_join(df_pop)

# 1-p
df_1p <- sdrep_sum_all %>%
    filter(pars %in% c("one_minus_p_mean")) %>%
    mutate(age = rep(seq(1,15), 2),
           one_minus_p_est = Estimate,
           one_minus_p_se = `Std. Error`) %>%
    mutate(est_se = sprintf("%.2f (%.3f)", Estimate, `Std. Error`)) 

# df_g_prime
df_g <- sdrep_sum_all %>%
    filter(pars == "g_prime_mean") %>%
    mutate(age = rep(seq(1,15), 2),
           g_prime_est = Estimate,
           g_prime_se = `Std. Error`)

list_plot_1p <- list()
for(p in c("4bc", "7a")) {
    ggplot(data = df_1p %>% filter(age >= 2, pop == p), aes(x = age, one_minus_p_est)) +
        geom_line() +
        geom_ribbon(aes(ymin = one_minus_p_est - one_minus_p_se,
                        ymax = one_minus_p_est + one_minus_p_se),
                    alpha = 0.5) +
        facet_wrap(~ pop_name) +
        labs(x = "Age (year)",
             y = expression((1-italic(p[t]))))
    list_plot_1p[[p]] <- last_plot()
}

list_plot_g <- list()
for(p in c("4bc", "7a")) {
    ggplot(data = df_g %>% filter(age >= 2, pop == p), aes(x = age, g_prime_est)) +
        geom_line() +
        geom_ribbon(aes(ymin = g_prime_est - g_prime_se,
                        ymax = g_prime_est + g_prime_se),
                    alpha = 0.5) +
        facet_wrap(~ pop_name) +
        labs(x = "Age (year)",
             y = expression(italic("g'")))
    list_plot_g[[p]] <- last_plot()
}

layout <- "
AB
CD
"

(list_plot_1p[["4bc"]] | list_plot_1p[["7a"]] |
        list_plot_g[["4bc"]] | list_plot_g[["7a"]]) +
    plot_layout(axis_titles = "collect", design = layout) 

# save file
ggsave(filename = "./report/fig7_1p_g.png",
       plot = last_plot(),
       width = 17,
       height = 15,
       unit = "cm",
       scale = 1)

ggsave(last_plot(), file = "./report/fig7_1p_g.pdf",
       device = cairo_pdf,
       width =  17*1, height = 15*1,
       units = "cm") 

# fig8 - correlation h, Tmat ----
ggplot(data = data_rep %>% filter(pop != "8ab"), aes(x = h, y = Tmat)) +
    geom_point(alpha = 0.1) +
    facet_grid(~ pop_name) +
    labs(x = expression(italic("h")~"(mm/year)"),
         y = expression(italic("T")~"(year)")) +
    theme(legend.position = "bottom",
          plot.margin = margin(t = 5.5, r = 10, b = 5.5, l = 5.5), 
          panel.spacing.x = unit(1, "lines"))

ggsave(filename = "./report/fig8_correlation_h_Tmat.png",
       plot = last_plot(),
       width = 17,
       height = 9,
       unit = "cm",
       scale = 1)

ggsave(last_plot(), file = "./report/fig8_correlation_h_Tmat.pdf",
       device = cairo_pdf,
       width =  17*1, height = 9*1,
       units = "cm") 

# example ID
df_example <- data_rep %>%
    filter(pop == "4bc", 
           id_int %in% c(384, 956)) %>%
    mutate(h = sprintf("%.2f", h),
           Tmat = sprintf("%.2f", Tmat)) 

df_example %>%
    select(id_int, h, Tmat) %>%
    unique()

# fig9 - trend Tmat ----
# cohort with at least 10 ind
data_rep_ncohort <- data_rep %>%
    select(pop, cohort, id) %>%
    unique() %>%
    group_by(pop, cohort) %>%
    summarize(n = n()) %>%
    mutate(pop_cohort = paste0(pop, "_", cohort)) %>%
    filter(n >= 10)

# summarize Tmat by cohort
data_rep_cohort10 <- data_rep %>%
    filter(paste0(pop, "_", cohort) %in% data_rep_ncohort$pop_cohort) %>%
    select(pop, cohort, h, Tmat, vinf) %>%
    unique() %>%
    group_by(pop, cohort) %>%
    mutate(Tmat_mean = mean(Tmat),
           Tmat_sd = sd(Tmat)) %>%
    left_join(df_pop)

# add trend label: 
# after 1970 4bc and 1980 7a - used to assess trend
# before that not

ggplot(data = data_rep_cohort10 %>% filter(pop %in% c("4bc", "7a")), aes(x = cohort, y = Tmat_mean)) +
    geom_smooth(data = data_rep_cohort10 %>% 
                    filter((pop == "4bc" & cohort > 1970)), 
                method = "lm",
                color = "black",
                linetype = "solid",
                linewidth = 0.5,
                se = F) +
    # geom_smooth(data = data_rep_age %>% 
    #                 filter((pop == "4bc" & cohort > 1970) | (pop == "7a" & cohort > 1980)), 
    #             method = "loess",
    #             color = "black",
    #             linetype = "dashed",
    #             linewidth = 0.5,
    #             se = F) +
    geom_point() +
    geom_linerange(aes(ymin = Tmat_mean - Tmat_sd, 
                       ymax = Tmat_mean + Tmat_sd)) +
    facet_grid(~ pop_name) +
    labs(x = "Cohort",
         y = expression("Cohort-mean" ~ italic("T") ~ "(year)")) 

ggsave(filename = "./report/fig9_trend_Tmat.png",
       plot = last_plot(),
       width = 17,
       height = 7,
       unit = "cm",
       scale = 1.5)

ggsave(last_plot(), file = "./report/fig9_trend_Tmat.pdf",
       device = cairo_pdf,
       width =  17*1.5, height = 9*1.5,
       units = "cm") 

## summary - lm trend Tmat ----
df_age <- tibble()
for(p in c("4bc", "7a")) {
    if(p == "4bc") {
        data_sub <- data_rep_cohort10 %>% filter(pop == p, cohort > 1970) 
    } else {
        data_sub <- data_rep_cohort10 %>% filter(pop == p, cohort > 1980) 
    }
    
    lm1 <- lm(Tmat ~ cohort, data = data_sub) 
    lm2 <- lm(Tmat_mean ~ cohort, data = data_sub) 
    
    
    df_temp <- tibble(pop = p,
                      cohort_effect_all = round(lm1$coefficients["cohort"],3),
                      p_value_all = round(summary(lm1)$coefficients["cohort",4],2),
                      cohort_effect_mean = round(lm2$coefficients["cohort"],3),
                      p_value_mean = round(summary(lm2)$coefficients["cohort",4],2))
    
    df_age <- bind_rows(df_age, df_temp)
}
df_age
