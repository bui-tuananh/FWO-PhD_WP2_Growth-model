# info ----
# input: 1_data-procesing_backcal.R
# output: ./output/backcal_7a_XXX

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
otl_backcal <- readRDS("./input/otl_full_backcal.rds")

# clean data (each cohort has at least age 1-5)
# sol_fab_1375 - strange growth curve with very low age 2 growth
# ggplot(data = otl_backcal %>% filter(id == "sol_fab_1375"),
#        aes(x = age, y = otl_diameter)) +
#     geom_point()

data <- otl_backcal %>% filter(age <= 15,
                        aac >= 5, 
                        pop == "7a",
                        !id %in% c("sol_fab_1375")) 

#### process data pre-fitting 

n_distinct(data$id)
# add id_int (0-(n_id-1)); cohort_int (0-(n_cohort-1))
data <- data %>%
    mutate(id_int = as.numeric(factor(id))-1,
           cohort_int = as.numeric(factor(cohort))-1)

# add age 0 observation with NA
# calculate growth increments
setDT(data)               # convert to datatable class
setorder(data, id, age)   # rearrange by id and age

data[ , len.incr:= backcal_len - shift(backcal_len, 1,type = "lag", fill = NA) , by = id]
setDF(data)
# fill in length at age 1
data$len.incr[is.na(data$len.incr)] <- data$backcal_len[is.na(data$len.incr)]

dat.ls <- split(data, f = data$id)
dat.ls <- lapply(dat.ls, function(x){
    x0 <- x[1,]
    x0$age = 0
    x0$backcal_len = NA
    x0$len.incr = NA
    rbind(x0,x)
})

data <- do.call("rbind", dat.ls)

# some data prep
ids <- c(1,1+cumsum(table(data$id))) # index corresponding to the first observation of an individual fish growth trajectory
ids <- ids[-length(ids)]             # figure out where the data starts

# tmb ----
dyn.load(TMB::dynlib("./R/quince_hmm_individual_h_Tmat_prec"))

n_states = 2
dat <- list(
    "n_states" = n_states,
    "growth_first_idx" = ids,                                     # integer vector
    "beta_linear_covs" = cbind(rep(1,nrow(data)), data$age),      # matrix with covariates for the transition state (intercept = 1 and age)
    "age" = c(data$age),                                          # vector of age covariates (for the observations)
    "len" = matrix(data$backcal_len, ncol = 1),                   # vector of observations
    #"keyVar" = (numeric(n_states)),                              # initial state probabilities (can be estimated in the model, but here we assume fish are immature at the start of their life)
    "id" = data$id_int
)

## initial pars ----
n_states = 2

# 1000 random initial pars 
set.seed(77777)

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

## optimisation ----
df_nll <- tibble()
list_opt <- list()
list_obj <- list()
# df_nll <- read_rds("./output/backcal_7a_individual_h_Tmat_prec_a100_age0_rep1000_v2_nll.rds")
# list_opt <- read_rds("./output/backcal_7a_individual_h_Tmat_prec_a100_age0_rep1000_v2_opt.rds")
# list_obj <- read_rds("./output/backcal_7a_individual_h_Tmat_prec_a100_age0_rep1000_v2_obj.rds")

for(i in 1:10000) {
    
    print(paste("processing", i, sep = " "))
    
    set.seed(i)
    par <- list(
        "log_beta_linear_pars" = matrix(c(ini_loga, log(sample(list_tmat,1))),nrow = 2),         # parameter matrix correpsonding to the transition probability from juvenile to adult (slope and Tmat)
        "v0"= sample(list_v0,1),
        "log_h" = log(sample(list_h,1)),                                                
        "logit_alpha" = logit(sample(list_alpha,1)),    
        "logit_chi" = logit(sample(list_chi,1)),      
        "log_prec" = log(2),  # high value for sigma (numerically more stable)
        "log_sd_h" = log(sample(list_sd,1)),
        "log_sd_tmat" = log(sample(list_sd,1)),
        "trans_rho" = c(itrans_rev(sample(list_rho,1))), #rho_h_tmat
        "mu_random" = matrix(0, nrow = 2, ncol = length(unique(dat$id)))  # row 1 = log_h, row 2 = log_Tmat
    )
    
    obj <- TMB::MakeADFun(dat, 
                          par,  
                          random = "mu_random",
                          map = list("log_beta_linear_pars" = factor(c(NA, 0))), #fix loga = log(100)
                          DLL="quince_hmm_individual_h_Tmat_prec")
    
    # next/break when there is error in opt
    if ( is.null(tryCatch({
        
        opt <- nlminb(obj$par, obj$fn, obj$gr,
                      control = list(iter.max = 10000, eval.max = 10000))
        
        #opt <- do.call("optim",obj)
        
        m_name <- paste("rep", i, sep = "_")
        
        list_opt[[m_name]] <- opt
        list_obj[[m_name]] <- obj
        
        df_nll_temp <- tibble(#ite_id = i,
            m_name = m_name,
            #nll = opt$value,
            nll = opt$objective, #nlminb - objective instead of value
            convergence = opt$convergence,
            message = opt$message) 
        
        df_nll <- bind_rows(df_nll, df_nll_temp)
        
        
    }, error = function(e){}))) {
        print("error");
        next
    }
    
    write_rds(df_nll, "./output/backcal_7a_individual_h_Tmat_prec_a100_age0_rep1000_v2_nll.rds")
    write_rds(list_opt, "./output/backcal_7a_individual_h_Tmat_prec_a100_age0_rep1000_v2_opt.rds")
    write_rds(list_obj, "./output/backcal_7a_individual_h_Tmat_prec_a100_age0_rep1000_v2_obj.rds")
    
    if (nrow(df_nll) == 1000) {
        print("run 1000 times, end")
        break
    }
    
}
