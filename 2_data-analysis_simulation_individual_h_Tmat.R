# info ----
# input: NA
# output: output/backcal_simulation

### note
# do 100 simulation round
## each simulation round fit model 10 times, with random initial pars
# get the model with the highest nll

# setup ----
library(tidyverse)
library(TMB)
library(boot)
library(data.table)

#### function
# constrain value (-1,1)
itrans <- function(x)(2/(1 + exp(-2*x))-1)

# reverse itrans
itrans_rev <- function(x)(log(2/(x+1)-1)/(-2))

# data ----
# function to simulate data
# parameter v0 (or l0), h, alpha, chi, Tmat are based on estimates when fitting HMM-Qunice model to North Sea population data
fn_sim <- function(seed = 2018, # seed for random generator, allowing reproducibility 
                   v0 = 30, # size at age 0 (l0)
                   h = 115, # potential annual growth rate
                   Tmat = 1.8, # age at maturity (last age at juvenile state)
                   alpha = 0.35, # proportion of the growing season in the first adult year devoted to growth
                   chi = 0.75, # indicator of the annual decrease rate of α
                   k = -log(chi), # rate of deceleration of growth 
                   
                   precision = 300, # precision (mu/sd)
                   sigma_h_ind = 0.2, # standard deviation of individual random effect in h
                   sigma_Tmat_ind = 0.2, # standard deviation of individual random effect in Tmat
                   rho_h_tmat = -0.7, # correlation between individual random effects in h and T 

                   n_ind = 250, # number of individual
                   max_age = 15 # maximum age
                   ) {
  
  set.seed(seed) # set seed to allow reproducibility
    
  #### true data
  id = seq(1, n_ind)        # fish id
  age = seq(1,max_age)      # age
  
  # variance covariance matrix
  SigmaMV = matrix(0, 2, 2)
  diag(SigmaMV) = c(sigma_h_ind^2, sigma_Tmat_ind^2)
  SigmaMV[1,2]  = rho_h_tmat * sigma_h_ind * sigma_Tmat_ind
  SigmaMV[2,1]  = SigmaMV[1,2]
 
  # random effects of h and Tmat
  mu_random <- mvtnorm::rmvnorm(n_ind, c(0, 0), SigmaMV) 
  random_h_ind = mu_random[,1]
  random_Tmat_ind = mu_random[,2]

  h_ind = exp(log(h))*exp(random_h_ind) # individual h
  Tmat_ind = Tmat*exp(random_Tmat_ind)  # individual Tmat

  # derived parameters
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
           g_prime = if_else(age <= Tmat, 
                       NA,
                       ((3*h)/(vinf - vc))*(1 - (vc/vt_sim)))
           )
  
  ## population data ----
  samp.size <- c(122, 58, 32, 16, 8, 6, 3, 2, 1, 1, 1) # sample size by age at capture (5-15)
  names(samp.size) <- as.character(seq(5,15))
  # assign age at capture for each fish
  df_samp.size = tibble(age = rep(seq(5, max_age), times = samp.size),
                        id = seq(1, n_ind))
  
  # population data, where fish have age at capture 5-15 and sample size as in samp.size
  data_pop <- tibble()
  for(i in 5:max_age) {
    df_samp.size_sub <- df_samp.size %>% filter(age == i)
    df_temp <- data_ind %>% filter(age <= i, id %in% df_samp.size_sub$id) 
    
    data_pop <- bind_rows(data_pop, df_temp)
  }
  
  # add age at capture 
  data_pop <- data_pop %>% 
    group_by(id) %>%
    mutate(age_at_capture = max(age))
  
  # plot
  p <- ggplot() +
    geom_point(data = data_pop, aes(x = age, y = vt_sim, group = id), alpha = 0.3) +
    geom_line(data = data_pop, aes(x = age, y = vt_true, group = id)) +
    labs(x = "Age (year)",
         y = "Length (mm)") +
    theme_
  
  # return data and plot
  return(list(data_pop, p))
}


# model ----
df_nll <- tibble()
list_opt <- list()
list_obj <- list()

for(i in 1:1000) {
  print(paste0("processing ", i))
  set.seed(i)
  #### sampled data
  #id_sample <- sample(data_pop$id, size = sample_size[s], replace = F)
  
  data <- fn_sim(i)[[1]]
  data <- data %>% 
    group_by(id) %>%
    mutate(year_adult = sum(!is.na(g))) #number of year with non-na g
  
  if(min(data$year_adult) < 1) {
    print(paste("minimum year adult <", 1, sep = " "))
    next
  } 
  
  ggplot() +
    geom_point(data = data, aes(x = age, y = vt_sim, group = id), alpha = 0.3) +
    geom_line(data = data, aes(x = age, y = vt_true, group = id))
  
  # add age 0 observation with NA
  # calculate growth increments
  setDT(data)               # convert to datatable class
  setorder(data, id, age)   # rearrange by id and age 
  
  data[ , len.incr:= vt_sim - shift(vt_sim, 1,type = "lag", fill = NA) , by = id]
  setDF(data)
  # fill in length at age 1
  data$len.incr[is.na(data$len.incr)] <- data$vt_sim[is.na(data$len.incr)]
  
  dat.ls <- split(data, f = data$id)
  dat.ls <- lapply(dat.ls, function(x){
    x0 <- x[1,]
    x0$age = 0
    x0$vt_sim = NA
    x0$len.incr = NA
    rbind(x0,x)
  })
  
  data <- do.call("rbind", dat.ls)
  
  # some data prep
  ids <- c(1,1+cumsum(table(data$id))) # index corresponding to the first observation of an individual fish growth trajectory 
  ids <- ids[-length(ids)]             # figure out where the data starts
  
  #### random initial pars - set 1000
  set.seed(i)
  v0_med <- median(data$vt_sim[data$age==1] - (data$vt_sim[data$age==2]-data$vt_sim[data$age==1]))
  h_med <- median(data$vt_sim[data$age==2]-data$vt_sim[data$age==1])
  
  ini_loga <- round(log(100),1)
  list_tmat <- runif(1000, min = 1, max = 3) #tmat 1-3
  list_v0 <- runif(1000, min = v0_med - 0.5*v0_med, max = v0_med + 0.5*v0_med)
  list_h <- runif(1000, min = h_med - 0.5*h_med, max = h_med + 0.5*h_med)
  list_alpha <- runif(1000, min = 0.3, max = 0.7) #assume alpha within 0.3-0.7
  list_chi <- runif(1000, min = 0.3, max = 0.7) #assume chi within 0.3-0.7
  # prec <- 2 - high value for sigma (numerically more stable)
  list_sd <- runif(1000, min = 0.1, max = 1) # assume sd within 0-1
  list_rho <- runif(1000, min = -1, max = 1) # rho -1 to 1
  
  #### tmb 
  #dyn.unload(TMB::dynlib("./R/quince_hmm_individual_h_Tmat_prec"))
  #TMB::compile("./R/quince_hmm_individual_h_Tmat_prec.cpp")
  dyn.load(TMB::dynlib("./R/quince_hmm_individual_h_Tmat_prec"))
  
  # note TMB:
  # id has -1 -> starts from 0
  n_states = 2
  dat <- list(
    "n_states" = n_states,
    "growth_first_idx" = ids,                                     # integer vector
    "beta_linear_covs" = cbind(rep(1,nrow(data)), data$age),      # matrix with covariates for the transition state (intercept = 1 and age)
    "age" = c(data$age),                                          # vector of age covariates (for the observations)
    "len" = matrix(data$vt_sim, ncol = 1),                        # vector of observations
    #"keyVar" = (numeric(n_states)),                              # initial state probabilities (can be estimated in the model, but here we assume fish are immature at the start of their life)
    "id" = as.numeric(as.factor(data$id)) - 1
  )
  
  #### each iteration run model 10 times
  df_nll_ite <- tibble()
  for(j in 1:100) {
    print(paste("processing random par", j, sep = " "))
    
    par <- list(
      "log_beta_linear_pars" = matrix(c(ini_loga, log(sample(list_tmat[j],1))),nrow = 2),         # parameter matrix correpsonding to the transition probability from juvenile to adult (slope and Tmat)
      "v0"= sample(list_v0[j],1),
      "log_h" = log(sample(list_h[j],1)),                                                
      "logit_alpha" = logit(sample(list_alpha[j],1)),    
      "logit_chi" = logit(sample(list_chi[j],1)),      
      "log_prec" = log(2),  # high value for sigma (numerically more stable)
      "log_sd_h" = log(sample(list_sd[j],1)),
      "log_sd_tmat" = log(sample(list_sd[j],1)),
      "trans_rho" = c(itrans_rev(sample(list_rho[j],1))), #rho_h_tmat
      "mu_random" = matrix(0, nrow = 2, ncol = length(unique(dat$id)))  # row 1 = log_h, row 2 = log_Tmat
    )
    
    # create fitting object with data, initial parameters, indication of random effects, and model
    obj <- TMB::MakeADFun(dat, 
                          par,  
                          random = "mu_random",
                          map = list("log_beta_linear_pars" = factor(c(NA, 0))), #fix loga = log(100)
                          DLL="quince_hmm_individual_h_Tmat_prec")
    
    # next/break when there is error in opt
    if ( is.null(tryCatch({
      opt <- do.call("optim",obj)
      
      ite_name <- paste("ite", i, sep = "_")
      m_name <- paste("ite", i, "rep", j, sep = "_")
      
      list_opt[[m_name]] <- opt
      list_obj[[m_name]] <- obj
      
      df_nll_temp <- tibble(#ite_id = i,
        ite_name = ite_name, 
        m_name = m_name,
        nll = opt$value) 
      df_nll_ite <- bind_rows(df_nll_ite, df_nll_temp)
      df_nll <- bind_rows(df_nll, df_nll_temp)
      
    }, error = function(e){}))) {
      print("error");
      next
    }
    
    write_rds(df_nll, "./output/backcal_hmm_simulation_h_Tmat_v2_nll.rds")
    write_rds(list_opt, "./output/backcal_hmm_simulation_h_Tmat_v2_opt.rds")
    write_rds(list_obj, "./output/backcal_hmm_simulation_h_Tmat_v2_obj.rds")
    
    
    if (nrow(df_nll_ite) == 10) {
      print("run 10 times, next loop")
      break
    }
    
  }
  
  if (nrow(df_nll) == 1000) {
    print("run 1000 times, end")
    break
  }
  
}
  