

data 
{
int<lower=0> k_area; 
int<lower=0> t_year; 
int<lower=0> n_covariate; 
int Y_kt[k_area, t_year];
real E_kt[k_area, t_year];
real area_AP_one[k_area, t_year];
// real area_AP_two[k_area, t_year];

// matrix[k_area, n_covariate] covariates;

real covariate_one[k_area, t_year];
real covariate_two[k_area, t_year];
real covariate_three[k_area, t_year];
real covariate_four[k_area, t_year];

int<lower=0> N_edges;
int<lower=1, upper=k_area> node1[N_edges];  // node1[i] adjacent to node2[i]
int<lower=1, upper=k_area> node2[N_edges];  // and node1[i] < node2[i]

}


transformed data
{
   int Y_kt_new[k_area*t_year];
  real E_kt_new[k_area*t_year];
  real area_AP_one_new[k_area*t_year];
  // real area_AP_two_new[k_area*t_year];
  real covariate_one_new[k_area*t_year];
  real covariate_two_new[k_area*t_year];
  real covariate_three_new[k_area*t_year];
  real covariate_four_new[k_area*t_year];
  
  for(n in 1:t_year)
  {
    Y_kt_new[(1+(n-1)*k_area):n*k_area]=Y_kt[,n];
    E_kt_new[(1+(n-1)*k_area):n*k_area]=E_kt[,n];
     area_AP_one_new[(1+(n-1)*k_area):n*k_area]=area_AP_one[,n];
     // area_AP_two_new[(1+(n-1)*k_area):n*k_area]=area_AP_two[,n];
    covariate_one_new[(1+(n-1)*k_area):n*k_area]=covariate_one[,n];
    covariate_two_new[(1+(n-1)*k_area):n*k_area]=covariate_two[,n];
    covariate_three_new[(1+(n-1)*k_area):n*k_area]=covariate_three[,n];
    covariate_four_new[(1+(n-1)*k_area):n*k_area]=covariate_four[,n];
  }
}


parameters 
{
vector[n_covariate] alpha;
real lambda1;
// real lambda2;
vector[k_area] phi_kt;
vector[k_area*t_year] theta;
real<lower=0.0001> nu2;
real<lower=0> tau_theta; // precision of heterogeneous effects

}

transformed parameters {
vector[k_area*t_year] phi_kt_new;
real<lower=0> sigma_theta = inv(sqrt(tau_theta)); // convert precision to sigma

for(n in 1:t_year)
{
  phi_kt_new[(1+(n-1)*k_area):n*k_area]=phi_kt;
}



}

model 
{
  alpha ~ normal(0,10);
  lambda1 ~ normal(0,10);
  // lambda2 ~ normal(0,10);
     nu2 ~ exponential(0.5);
  
// Y_kt ~ poisson_log(to_vector(log(E_kt)) + phi_kt+lambda*to_vector(area_AP)+lambda2*to_vector(area_AP2)+theta * sigma_theta + covariates*alpha);
// vector operator*(matrix x, vector y)

Y_kt_new ~ poisson_log(to_vector(log(E_kt_new)) 
                    + phi_kt_new
                    +lambda1*to_vector(area_AP_one_new)
                    // +lambda2*to_vector(area_AP_two_new)
                    +theta * sigma_theta
                    +to_vector(covariate_one_new)*alpha[1]
                    +to_vector(covariate_two_new)*alpha[2]
                    +to_vector(covariate_three_new)*alpha[3]
                    +to_vector(covariate_four_new)*alpha[4]);
 
 
target += -0.5 * (1.0/nu2^2) * dot_self(phi_kt[node1] - phi_kt[node2])-0.5*k_area*log(nu2^2);
  sum(phi_kt) ~ normal(0, 0.001 * k_area);  // equivalent to mean(phi) ~ normal(0,0.001)
  
  theta ~ normal(0, 1);
  tau_theta ~ gamma(3.2761, 1.81); // Carlin WinBUGS priors
  
}
