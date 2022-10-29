functions {
  vector phi_inv(vector mu, int K) {
    vector[K] e_mu;

    e_mu = exp(mu);
    return append_row(e_mu / (1 + sum(e_mu)), 1 / (1 + sum(e_mu)));
  }
}

data {
  int<lower=1> n_species; // number of species - 1 
  int<lower=1> n_timepoints;// number of time points
  array[n_timepoints, n_species+1] int<lower=0> y; //microbiom abundance //change to array
  vector[n_timepoints] alpha;
}

parameters {
  vector<lower=0>[n_species]sigma;
  matrix[2, n_species] beta; //coefficient for in/out clusters
}


model {
  vector[n_species + 1] p;
  for (l in 1:2) {
    for (k in 1:n_species) {
      // beta[l, k] ~ normal(0, 1);
      beta[l, k] ~ normal(0, sigma[k]);// mean and variance can be changed based on parameters
    }
  }
  vector[n_species] mu_;
  for (i in 1:n_timepoints) {
    mu_ = (beta' * to_vector({alpha[i], 1-alpha[i]}));
    p = phi_inv(mu_, n_species);
    y[i] ~ multinomial(to_vector(p));
  }
}

