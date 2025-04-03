# info ----
# input: 
# - 2_data-analysis_backcal_4bc/7a/8ab
# output: output/individual_h_Tmat_a100_age0_datarep_v2.rds

# setup ----
library(tidyverse)
library(TMB)
library(boot)
library(data.table)

# data ----
otl_backcal <- readRDS("./input/otl_full_backcal.rds")
otl_backcal_sub <- otl_backcal %>% filter(age <= 15,
                                      aac >= 5, 
                                      !id %in% c("sol_fab_1375")) 
    
# sdreport ----
# get best fit for each population
df_nll_best <- tibble()
for(p in c("4bc", "7a", "8ab")) {
    df_nll <- read_rds(paste0("./output/", "backcal_", p, "_individual_h_Tmat_prec_a100_age0_rep1000_v2_nll.rds"))
    df_nll_best_temp <- df_nll %>% 
        filter(message == "relative convergence (4)") %>%
        filter(nll == min(nll)) %>%
        mutate(pop = p)
    
    df_nll_best <- bind_rows(df_nll_best, df_nll_best_temp)
}

# loop to get data_rep
data_rep <- tibble()
for(i in 1:3) {
    
    # note
    print(paste0("processing ", df_nll_best$pop[i]))
    
    ## data ----
    data <- otl_backcal_sub %>% 
        filter(pop == df_nll_best$pop[i]) %>%
        mutate(id_int = as.numeric(factor(id))-1)

    # some data prep
    setDT(data)               # convert to datatable class
    setorder(data, id, age)   # rearrange by id and age
    dat.ls <- split(data, f = data$id)
    dat.ls <- lapply(dat.ls, function(x){
        x0 <- x[1,]
        x0$age = 0
        x0$backcal_len = NA
        rbind(x0,x)
    })
    data <- do.call("rbind", dat.ls)
    
    ids <- c(1,1+cumsum(table(data$id))) # index corresponding to the first observation of an individual fish growth trajectory
    ids <- ids[-length(ids)]             # figure out where the data starts
    
    ### tmb ----
    #dyn.unload(TMB::dynlib("./R/quince_hmm_individual_h_Tmat_prec"))
    #TMB::compile("./R/quince_hmm_individual_h_Tmat_prec.cpp")
    dyn.load(TMB::dynlib("./R/quince_hmm_individual_h_Tmat_prec"))
    
    # note TMB:
    # keyVar has no +1 (0,0)
    # id has -1 -> starts from 0
    n_states = 2
    
    ## report ----
    df_nll <- read_rds(file.path("./output", paste0("backcal_", df_nll_best$pop[i], "_individual_h_Tmat_prec_a100_age0_rep1000_v2_nll.rds")))
    list_opt <- read_rds(file.path("./output", paste0("backcal_", df_nll_best$pop[i], "_individual_h_Tmat_prec_a100_age0_rep1000_v2_opt.rds")))
    list_obj <- read_rds(file.path("./output", paste0("backcal_", df_nll_best$pop[i], "_individual_h_Tmat_prec_a100_age0_rep1000_v2_obj.rds")))
    # saved rep and sdrep below to avoid load list_obj here
    
    #### report
    obj <- list_obj[[df_nll_best$m_name[i]]]
    
    rep <- obj$report()
    write_rds(rep, paste0("output/individual_h_Tmat_a100_age0_rep_", df_nll_best$pop[i], "_v2.rds"))
    rep  <-  read_rds(paste0("output/individual_h_Tmat_a100_age0_rep_", df_nll_best$pop[i], "_v2.rds"))
    # note: save rep to avoid load obj
    
    #### sdreport
    sdrep <- TMB::sdreport(obj)
    sdrep_sum <- summary(sdrep)
    sdrep_sum <- as.data.frame(sdrep_sum) %>% mutate(pars = row.names(sdrep_sum))
    any(is.na(sdrep_sum))
    
    # note: save sdrep_sum to avoid running sdreport
    write_rds(sdrep_sum, paste0("output/individual_h_Tmat_a100_age0_sdrepsum_", df_nll_best$pop[i], "_v2.rds"))
    sdrep_sum <-  read_rds(paste0("output/individual_h_Tmat_a100_age0_sdrepsum_", df_nll_best$pop[i], "_v2.rds"))
    
    
    df_ind <- tibble(id = unique(data$id),
                     v0 = filter(sdrep_sum, pars == "v0")$Estimate[1],
                     h = filter(sdrep_sum, pars == "h")$Estimate,
                     Tmat = filter(sdrep_sum, pars == "Tmat")$Estimate,
                     alpha = filter(sdrep_sum, pars == "alpha")$Estimate,
                     chi = filter(sdrep_sum, pars == "chi")$Estimate,
                     vinf = filter(sdrep_sum, pars == "vinf")$Estimate)
    
    data_rep_temp <- data %>% 
        left_join(df_ind) %>%
        mutate(mu_juv = rep$mu[,1],
               mu_adu = rep$mu[,2],
               mu = if_else(age <= Tmat, mu_juv, mu_adu)) %>%
        filter(age >= 1)
    
    data_rep <- bind_rows(data_rep, data_rep_temp)
    
}

# save data_rep
write_rds(data_rep, "output/individual_h_Tmat_a100_age0_datarep_v2.rds")
