functions {

    // double bicop_normal_log(vector u, vector v, double rho) {
    //
    //   return 1;
    // }



}


data {
    int<lower=0> n_max;
    int<lower=0> t_max;
    matrix[t_max,n_max] u;
    int<lower=0> k;
    int gid[n_max];
    int copula_type[n_max];
}

parameters {
  vector<lower=0,upper=1>[t_max] v;
  vector<lower=-1,upper=1>[n_max] theta;
}



model {
    v ~ uniform(-1,1);
    theta ~ uniform(-1,1);

    for (j in 1:n_max){
        target += normal_lpdf( col(u,j) | 0, theta[j]  ) ;
    }

}

