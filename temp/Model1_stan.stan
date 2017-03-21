data {
    int<lower=0> n_max;
    int<lower=0> t_max;
    matrix[t_max,n_max] z;
}

parameters {
  vector<lower=0,upper=1>[n_max] lambda;
  vector[t_max] w;	
}

transformed parameters {
  vector[n_max] sigma_tilde;
  for (i in 1:n_max)
	sigma_tilde[i] <- sqrt(1-lambda[i]^2);
}

model {
  for (i in 1:n_max) 
	lambda[i] ~ uniform(0.0001, 0.9999);

  for (t in 1:t_max) 
	w[t] ~ normal(0,1);

  for (t in 1:t_max)
	for (j in 1:n_max)
		z[t,j] ~ normal(w[t]*lambda[j], sigma_tilde[j]);

}

