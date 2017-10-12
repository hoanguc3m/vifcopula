library(devtools)
install_github("hoanguc3m/vifcopula")
#setwd("/home/hoanguc3m/Dropbox/WP2/")
library(vifcopula)
t_max = 1000
n_max = 100
gauss_init <- matrix(1, nrow = n_max, ncol = 1)

load("/home/hoanguc3m/Dropbox/WP2/code/vifcopula/vignettes/sim/data_fin.RData")

n_max = dim(fin_u_mat)[2]
t_max = dim(fin_u_mat)[1]
gid = fin_SIC_group
k_max = max(gid)+1
structfactor = 3

copfamily = matrix(sample(c(1,2,3,4,5,6),size = t_max, replace = T),ncol=1)
latent_copfamily = matrix(sample(c(1,2,3,4,5,6),size = k_max-1, replace = T),ncol=1)
gauss_init <- rep(1, t_max)
gauss_latent_init <- rep(1, k_max-1)


data <- list(u = fin_u_mat,
             n_max = n_max,
             t_max = t_max,
             k_max = k_max,
             gid = gid,
             structfactor = structfactor)
init <- list(copula_type = gauss_init,
             latent_copula_type = gauss_latent_init,
             v = matrix(0.5,nrow = t_max, ncol = 1),
             par = rep(0.5,t_max),
             par2 = rep(0.5,t_max))
other <- list(seed = 126, core = 8, iter = 1000,
              n_monte_carlo_grad = 1, n_monte_carlo_elbo = 10,
              eval_elbo = 100, adapt_bool = F, adapt_val = 1,
              adapt_iterations = 50, tol_rel_obj = 0.1, copselect = T)
vi_mix <- vifcopula::vifcop(data,init,other)
