functions {
  vector phi_inv(vector mu, int K) {
    vector[K] e_mu;

    e_mu = exp(mu);
    return append_row(e_mu / (1 + sum(e_mu)), 1 / (1 + sum(e_mu)));
  }
}

data {
  int<lower=1> n_species; // number of species - 1 //Q1: why I need to use species - 1 sometimes
  int<lower=1> n_steps;// number of time points
  int<lower=1> n_clust;//number of clusters
  int<lower=1> n_person;//number of persons
 // int<lower=1> n;// number of samples(without replicates)
  array[(n_steps*n_person), n_species+1] int<lower=0> y; //microbiom abundance //change to array
  //int<lower=1> depth;//sequence depth
}

parameters {
 //vector[n_species] mu;
 vector<lower=0>[n_species]sigma;
 //vector[n_species]mu;//May give too much variant for the spiecies
 real<lower=0.1,upper=10> dirich;
 real<lower = 0, upper = 1> scale;
  //real<lower=0> sigma;
  matrix[n_clust, n_species] beta; //coefficient for in/out clusters
  //matrix[n_clust, n_clust]; //transition matrix
  array[(n_steps*n_person)] simplex[n_clust] theta; // in/out clusters
  //vector[n_person] person_effect;
}


model {
  vector[n_species + 1] p;
  for (l in 1:n_clust) {
    for (k in 1:n_species) {
      // beta[l, k] ~ normal(0, 1);
      beta[l, k] ~ normal(0, sigma[k]);// mean and variance can be changed based on parameters
    }
  }

  // for (i in 1:n_species){
  //   beta[i] ~ multi_normal(mu, Sigma);
  // }
  dirich ~ gamma(5,scale);//can be changed to 10 // can make it a user choice
  for (i in 1:(n_steps*n_person)) {
    theta[i] ~ dirichlet(rep_vector(dirich, n_clust));//set to 0.1 will be easier to threshould
  }

  vector[n_species] mu_;
  // Can add person's effect
  // for(j in 1:n_person){
  //   person_effect[j] ~ normal(0,1);
  // }
  for (i in 1:(n_steps*n_person)) {
    mu_ = (beta' * to_vector(theta[i]));//Q2: how the dimension actually works here
    p = phi_inv(mu_, n_species);
    y[i] ~ multinomial(to_vector(p));
  }
}
//Q3: how we simulate the transition among states





//can be used to evaluate the model.
// generated quantities {
//   // simulate samples based on current treatment
//   vector[n_species + 1] p_sim;
//   matrix[n_clust, n_species] mu_;
//   array[n_clust, n_species + 1] int<lower=0> y_sim;
//
//   for (j in 1:n) {
//     mu_[j] =  (beta' * to_vector(theta[j]))';
//   }
  //
  // for (i in 1:n) {
  //   p_sim = phi_inv(mu_[j], n_species);
  //   y_sim[i] = multinomial_rng(to_vector(p_sim), depth);
  // }

//}
