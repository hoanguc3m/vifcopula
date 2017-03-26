functions {

    // double bicop_normal_log(vector u, vector v, double rho) {
    //
    //   return 1;
    // }



}


data {
    int<lower=0> n_max;
    int<lower=0> t_max;
    vector[2] u[t_max];
    vector[2] v[t_max];
    int copula_type;
}

transformed data {
    vector[2] mu;
    mu[1] = 0;
    mu[2] = 0;
}

parameters {
    real<lower = -1, upper = 1> theta;
    real<lower = 2, upper = 30> theta2;
}

transformed parameters {
    matrix[2, 2] Sigma;
    Sigma[1, 1] = 1;
    Sigma[1, 2] = theta2;
    Sigma[2, 1] = theta2;
    Sigma[2, 2] = 1;
}

model {
    theta ~ uniform(-1,1);
    theta2 ~ uniform(2,30);

    for (j in 1:n_max){
        target += multi_normal_lpdf(u[j] | mu, Sigma);
    }

}


