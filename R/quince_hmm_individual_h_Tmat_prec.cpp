#include <TMB.hpp>

// helper functions for NA handling
template<class Type>
bool isNA(Type x){
  return R_IsNA(asDouble(x));
}

template<class Type>
bool isNAINT(int x){
  return NA_INTEGER==x;
}

template <class Type>
Type itrans(Type x){ // scaled ilogit
  return Type(2)/(Type(1) + exp(-Type(2) * x)) - Type(1);
}


// function to transform initial state vector
// n-1 parameters required for n states (sum == 1)
template <class Type>
vector<Type> delta_transform(vector<Type> delta_in)
{
  vector<Type> res_unscaled(delta_in.size()+1);
  res_unscaled(0) = 1;
  for(int i = 1; i<res_unscaled.size();i++){
    res_unscaled(i) = exp(delta_in(i-1));
  }
  return res_unscaled/res_unscaled.sum();
}

template<class Type>
Type objective_function<Type>::operator() ()
{
  /// I. DATA & PARAMETER SECTION ///
  
  // 1. global variables
  DATA_INTEGER(n_states);                       // number of states
  DATA_IVECTOR(growth_first_idx);               // first index of each growth trajectory in the data
  
  // 2. state transition model
  // linear covariates
  DATA_MATRIX(beta_linear_covs);                // Intercept and linear covariate for the transition probability matrix
  
  // 3. covariates/observations & indexes
  // linear effects
  DATA_VECTOR(age);                               // covariates in the growth models
  // the observations
  DATA_MATRIX(len);                               // vector with length-at-age observations      
  //DATA_IVECTOR(keyVar);                         // vector that allows to define the variance parameter by state
  DATA_IVECTOR(id);                               // vector fish id
  
  // 4. initial state parameters
  //PARAMETER_VECTOR(log_delta);                  // initial state probability parameters
  
  // 5. transition matrix parameters
  PARAMETER_MATRIX(log_beta_linear_pars);          // parameters of linear covariates that govern the state transition probablities
  
  // 6. parameters of the Quince growth model
  PARAMETER(v0);             // size at age 0 
  PARAMETER(log_h);          // potential annual growth rate, h > 0
  PARAMETER(logit_alpha);    // proportion of the growing season in the first adult year devoted to growth (0 < alpha < 1)
  PARAMETER(logit_chi);      // indicator of the annual decrease rate of alpha (0 < chi < 1)
  // Tmat is defined below

  
  // 7. observation variance
  PARAMETER(log_prec);        // precision
  PARAMETER(log_sd_h);        // standard deviation h
  PARAMETER(log_sd_tmat);     // standard deviation Tmat
  
  PARAMETER_VECTOR(trans_rho); // correlation h and Tmat
  PARAMETER_MATRIX(mu_random); // matrix of random effects: log_h, log_Tmat
  
  /// II. MODEL SECTION ///
  
  // 1. parameter transformation: working scale to natural scale
  
  // slope maturing probability and Tmat
  matrix<Type> beta_linear_pars = (log_beta_linear_pars.array());   
  beta_linear_pars(0,0) = exp(beta_linear_pars(0,0)); // slope maturing prob > 0
  beta_linear_pars(1,0) = exp(beta_linear_pars(1,0)); // Tmat > 0
  
  //variance covariance matrix
  // sd_h^2                    rho_h_Tmat*sd_h*sd_Tmat         
  // rho_h_Tmat*sd_h*sd_Tmat   sd_Tmat^2                       
  matrix<Type> SigmaMV(2, 2);
  SigmaMV.setZero();
  // diagonal
  Type rho_h_tmat = itrans(trans_rho(0)); //transform trans_rho to rho [-1,1]
  SigmaMV(0,0) = exp(log_sd_h)*exp(log_sd_h);
  SigmaMV(1,1) = exp(log_sd_tmat)*exp(log_sd_tmat);
  // off-diagonal
  SigmaMV(0,1) = rho_h_tmat*exp(log_sd_h)*exp(log_sd_tmat);
  SigmaMV(1,0) = SigmaMV(0,1);
  
  // individual h and Tmat
  vector<Type> h(mu_random.cols()); 
  vector<Type> Tmat(mu_random.cols()); 
  for(int i = 0; i < mu_random.cols(); i++) {
    h(i) = exp(log_h)*exp(mu_random(0,i));
    Tmat(i) = beta_linear_pars(1,0)*exp(mu_random(1,i));
  }
  
  // precision, alpha, chi, k
  Type prec = exp(log_prec);
  Type alpha = invlogit(logit_alpha);
  Type chi = invlogit(logit_chi);
  Type k = -log(chi);
  
  // 2. derived parameters & variables
  vector<Type> vT = v0 + h*Tmat;                    //length at last juvenile age 
  vector<Type> vinf = vT+((h*alpha)/(1-chi));       //asymptotic length
  vector<Type> vc = vT + h*(alpha - chi)/(1 - chi); //theoretical 'intersect size' 
  
  int n_obs = len.size(); // number of observations

  // initial state distribution (Pr(juvenile) = 1, Pr(adult) = 0)
  vector<Type> delta(n_states);
  delta  = delta *0;
  delta(0)=1;
  
  // 3. Compute likelihood of observations for each state

  // matrix to store the expected length for each observation for each stage
  matrix<Type> mu(n_obs,n_states);
  mu.setZero();
  
  for(int i = 0; i < n_obs; i++){
    for(int j = 0; j < n_states; j++){
      if(j == 0){
        // immature growth;
        mu(i,j) = v0 + age(i) * h(id(i));                   
        
      }
      if(j == 1){
        // adult growth
        mu(i,j) = vinf(id(i)) - (vinf(id(i)) - vT(id(i)))*exp(-k*(age(i)-Tmat(id(i))));  
      }
    }
  }
  
  // calculate the probabilities of the data across all data streams. 
  // Here we only have 1 data stream: length-at-age but the model can be extended for multiple streams
  matrix<Type> allProbs(n_obs,n_states);
  allProbs.setOnes();
  
  int nr_data_streams = len.cols();  // we only have a single observation at each time step
  
  array<Type> dataProb(nr_data_streams,n_obs,n_states);
  for(int d = 0; d < nr_data_streams; d++){  
    vector<Type> observed = len.col(d);
    for(int state = 0; state < n_states; state++){
      for(int i=0;i<n_obs;i++) {
        
        if (isNA(observed(i))){
          dataProb(d,i,state) = 1;
        } else {
          dataProb(d,i,state) = dnorm(observed(i), mu(i,state), mu(i,state)/prec, false);  // gamma/lognormal distribution may be more appropriate (strictly positive)
        }
        allProbs(i,state) = allProbs(i,state) * dataProb(d,i,state);
      }
    }
  }
  
  // 4. Computation of transition probability matrix trMat 
  
  //-------------------------------------------------------------
  // calculate the log probabilities
  //-------------------------------------------------------------
  
  int nr_transitions = n_states * (n_states-1); 
  matrix<Type> g_mat(n_obs,nr_transitions);
  
  // initialize a transition probability matrix for each state (3d array)
  array<Type> trMat(n_states,n_states,n_obs);   
  trMat.setZero();                              
  Type pMat = 0;
  
  for(int o = 0; o < n_obs; o++) {
    pMat = 1/(1+exp(beta_linear_pars(0,0)*beta_linear_covs(o,0)*(Tmat(id(o)) - beta_linear_covs(o,1))));
    trMat(0,0,o) = 1 - pMat;
    trMat(0,1,o) = pMat;
    trMat(1,0,o) = 0;
    trMat(1,1,o) = 1;
  }
  
  // nll
  // initialization 
  Type nll = 0;

  //random effects
  // mu_random
  int n_random = mu_random.cols(); //ncol(mu_random) 
  
  using namespace density;
  MVNORM_t<Type> neg_log_dmvnorm(SigmaMV);
  
  for(int i=0; i < n_random; i++)
  {
      nll += neg_log_dmvnorm(mu_random.col(i));
  }  
 
  // 6. Forward algorithm  ***bring it all together*** 
  
  matrix<Type> Gamma(n_states,n_states);           // store values of the transition matrix for observation i      
  matrix<Type> lalpha(n_obs,n_states);             // matrix that stores output for residuals
  
  // some objects to store data while looping
  int z=0;                               
  vector<Type> alfa(n_states);          // alpha was already declared, used alfa instead
  vector<Type> alpha_vec(n_obs);
  matrix<Type> alpha_mat(1,n_states);
  matrix<Type> alpha_gamma(1,n_states);
  vector<Type> alpha_gamma_vec(n_states);
  vector<Type> probvec(n_states); 
  
  for(int i = 0; i < n_obs; i++) {
    
    // take the right slice from the transition state probability matrix (probably, a short-cut code exists to do so)
    for(int nrows = 0; nrows < n_states; nrows++){
      for(int ncols = 0; ncols < n_states; ncols++){
        Gamma(nrows,ncols) = trMat(nrows,ncols,i);
      }
    }
    
    if(z < growth_first_idx.size() && i == (growth_first_idx(z)-1)) {
      // if 'i' is the 'z'-th element of 'growth_first_idx', switch to the next time series
      z++;
      probvec = allProbs.row(i);
      alfa = delta * probvec;  //elementwise product: probabilities first state * state probabilities
    } 
    else {
      alpha_mat.row(0) = alfa ;
      alpha_gamma = alpha_mat * Gamma;
      alpha_gamma_vec = alpha_gamma.row(0);  // matrix to vector
      probvec = allProbs.row(i);
      alfa = alpha_gamma_vec * probvec;
    }
    alpha_vec(i) = alfa.sum();
    nll -= log(alfa.sum());
    alfa = alfa/alfa.sum();
    lalpha.row(i) = log(alfa) - nll;
  }
  
  /// III. REPORT SECTION ///
  
  // REPORT
  REPORT(v0);
  REPORT(h);
  REPORT(Tmat); 
  REPORT(alpha);
  REPORT(chi);
  REPORT(vT) ;
  REPORT(vinf);
  REPORT(vc);
  REPORT(mu);
  REPORT(mu_random);
  REPORT(rho_h_tmat);
  
  REPORT(delta);
  REPORT(allProbs);
  REPORT(dataProb);
  REPORT(trMat);
  REPORT(g_mat);
  
  // ADREPORT
  // ADREPORT(delta);
  ADREPORT(v0);
  Type h_mean = exp(log_h);
  ADREPORT(h_mean);
  ADREPORT(h);
  Type Tmat_mean = beta_linear_pars(1,0);
  ADREPORT(Tmat_mean);
  ADREPORT(Tmat);
  ADREPORT(alpha);
  ADREPORT(chi);
  ADREPORT(vT);
  ADREPORT(vinf);
  ADREPORT(vc);
  Type sd_h = exp(log_sd_h);
  Type sd_tmat = exp(log_sd_tmat);
  ADREPORT(sd_h);
  ADREPORT(sd_tmat);
  ADREPORT(rho_h_tmat);
  ADREPORT(prec);
  
  // Prediction population mean  
  Type vT_mean = v0 + h_mean*Tmat_mean;                    // length at last juvenile age 
  Type vinf_mean = vT_mean + ((h_mean*alpha)/(1-chi));     // asymptotic length
  Type vc_mean = vT_mean + h_mean*(alpha - chi)/(1 - chi); // theoretical 'intersect size' 
  
  int max_age = 15;  // maximum age for prediction
  vector<Type> pred_mean_juv(max_age);
  vector<Type> pred_mean_adu(max_age);
  vector<Type> g_prime_mean(max_age);
  for(int i = 0; i < max_age; i++){
        Type a = i+1; // age (+1 as index starts from 0)
    
        // immature growth;
        pred_mean_juv(i) = v0 + a*h_mean;

        // adult growth
        pred_mean_adu(i) = vinf_mean - (vinf_mean - vT_mean)*exp(-k*(a-Tmat_mean));

        // g_prime
        g_prime_mean(i) = (3*h_mean/(vinf_mean - vc_mean))*(1 - (vc_mean/pred_mean_adu(i)));
    }

  ADREPORT(vT_mean);
  ADREPORT(vinf_mean);
  ADREPORT(vc_mean);
  ADREPORT(pred_mean_juv);
  ADREPORT(pred_mean_adu);
  ADREPORT(g_prime_mean);
  
  // nll
  return nll;
}
