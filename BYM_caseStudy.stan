

data 
{
int<lower=0> k_area; 
int<lower=0> n_covariate; 
int Y_kt[k_area];
real E_kt[k_area];
real area_AP[k_area];

matrix[k_area, n_covariate] covariates;

int<lower=0> N_edges;
int<lower=1, upper=k_area> node1[N_edges];  // node1[i] adjacent to node2[i]
int<lower=1, upper=k_area> node2[N_edges];  // and node1[i] < node2[i]

}


parameters 
{
vector[n_covariate] alpha;
real lambda;
vector[k_area] phi_kt;
vector[k_area] theta;
real<lower=0.0001> nu2;
real<lower=0> tau_theta; // precision of heterogeneous effects

}

transformed parameters {
real<lower=0> sigma_theta = inv(sqrt(tau_theta)); // convert precision to sigma
}

model 
{
  alpha ~ normal(0,10);
  lambda ~ normal(0,10);
     nu2 ~ exponential(0.5);
  
Y_kt ~ poisson_log(to_vector(log(E_kt)) + phi_kt+lambda*to_vector(area_AP)+theta * sigma_theta + covariates*alpha);
// vector operator*(matrix x, vector y)
 
target += -0.5 * (1.0/nu2^2) * dot_self(phi_kt[node1] - phi_kt[node2])-0.5*k_area*log(nu2^2);
  sum(phi_kt) ~ normal(0, 0.001 * k_area);  // equivalent to mean(phi) ~ normal(0,0.001)
  
  theta ~ normal(0, 1);
  tau_theta ~ gamma(3.2761, 1.81); // Carlin WinBUGS priors
  
}
