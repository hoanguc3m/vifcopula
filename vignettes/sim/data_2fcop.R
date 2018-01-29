library(devtools)
install_github("hoanguc3m/vifcopula")
#setwd("/home/hoanguc3m/Dropbox/WP2/")
library(vifcopula)
load("/home/hoanguc3m/Dropbox/WP2/code/vifcopula/vignettes/sim/data_clean.RData")

t_max = 1000
n_max = 100

copfamily = matrix(sample(c(1,2,3,4,5,6),size = t_max, replace = T),ncol=1)

gid <- rep(1,t_max)
k_max = max(gid)+1
structfactor = 2

gauss_init <- rep(1, n_max)
gauss_latent_init <- rep(1, n_max)


data <- list(u = new_u_mat,
             n_max = n_max,
             t_max = t_max,
             k_max = k_max,
             gid = gid,
             structfactor = structfactor)
init <- list(copula_type = gauss_init,
             latent_copula_type = gauss_latent_init)

other <- list(seed = 126, core = 8, iter = 1000,
              n_monte_carlo_grad = 1, n_monte_carlo_elbo = 10,
              eval_elbo = 100, adapt_bool = F, adapt_val = 1,
              adapt_iterations = 50, tol_rel_obj = 0.1, copselect = T)
vi_mix_2f <- vifcopula::vifcop(data,init,other)
plot.vifcop(vi_mix_2f)
sum(vi_mix_2f$cop_type == 2)
sum(vi_mix_2f$latent_copula_type == 0)
save.image("/media/hoanguc3m/Data/wp2/data_2f.Rdata")
plot.vifcop(vi_mix_2f)

plot.vifcop(vi_mix_2f)
round(vi_mix_2f$time)
round(vi_mix_2f$criteria)
vi_mix_2f$iteration
sum(vi_mix_2f$cop_type == 1) + sum(vi_mix_2f$cop_type == 21)
sum(vi_mix_2f$cop_type == 2) + sum(vi_mix_2f$cop_type == 22)
sum(vi_mix_2f$cop_type == 3) + sum(vi_mix_2f$cop_type == 13) + sum(vi_mix_2f$cop_type == 23) + sum(vi_mix_2f$cop_type == 33)
sum(vi_mix_2f$cop_type == 4) + sum(vi_mix_2f$cop_type == 14) + sum(vi_mix_2f$cop_type == 24) + sum(vi_mix_2f$cop_type == 34)
sum(vi_mix_2f$cop_type == 5) + sum(vi_mix_2f$cop_type == 15) + sum(vi_mix_2f$cop_type == 25) + sum(vi_mix_2f$cop_type == 35)
sum(vi_mix_2f$cop_type == 6)

sum(vi_mix_2f$latent_copula_type == 1) + sum(vi_mix_2f$latent_copula_type == 21)
sum(vi_mix_2f$latent_copula_type == 2) + sum(vi_mix_2f$latent_copula_type == 22)
sum(vi_mix_2f$latent_copula_type == 3) + sum(vi_mix_2f$latent_copula_type == 13) + sum(vi_mix_2f$latent_copula_type == 23) + sum(vi_mix_2f$latent_copula_type == 33)
sum(vi_mix_2f$latent_copula_type == 4) + sum(vi_mix_2f$latent_copula_type == 14) + sum(vi_mix_2f$latent_copula_type == 24) + sum(vi_mix_2f$latent_copula_type == 34)
sum(vi_mix_2f$latent_copula_type == 5) + sum(vi_mix_2f$latent_copula_type == 15) + sum(vi_mix_2f$latent_copula_type == 25) + sum(vi_mix_2f$latent_copula_type == 35)
sum(vi_mix_2f$latent_copula_type == 6)
sum(vi_mix_2f$latent_copula_type == 0)


