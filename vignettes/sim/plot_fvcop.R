library(vifcopula)
set.seed(0)
t_max = 1000
n_max = 10
k_max = 1
g_max = 1
gid = sample(1:g_max,n_max,replace = T)

datagen_gauss <- fcopsim(t_max = t_max, n_max = n_max, k_max = k_max, gid = gid,
                         family = 1, family_vine = 1, seed_num = 100,
                         structfactor = 4)
datagen <- datagen_gauss

data <- list(u = datagen$u,
             n_max = datagen$n_max,
             t_max = datagen$t_max,
             k_max = datagen$k_max,
             gid = datagen$gid,
             structfactor = datagen$structfactor)
init <- list(copula_type = datagen$family,
             vine_copula_type = datagen$family_vine,
             vine_edges = datagen$egdes)
other <- list(seed = 100909, core = 8, iter = 1000,
              n_monte_carlo_grad = 1, n_monte_carlo_elbo = 100,
              eval_elbo = 100, adapt_bool = F, adapt_val = 1,
              adapt_iterations = 50, tol_rel_obj = 0.01, copselect = F, modelselect = T)
vi_gauss <- vifcopula::vifcop(data,init,other)
plot(datagen$v, head(vi_gauss$mean_vi,t_max))
abline(a = 0, b = 1, col = "red")
plot(datagen$theta, vi_gauss$mean_vi[(t_max+1):(t_max+n_max)])
abline(a = 0, b = 1, col = "red")
plot(datagen$theta_vine, tail(vi_gauss$mean_vi,n_max - max(gid)))
abline(a = 0, b = 1, col = "red")

plot(datagen$theta, vi_gauss$mean_vi[seq(t_max+1,t_max+2*n_max, by = 2)])
abline(a = 0, b = 1, col = "red")
plot(datagen$theta2, vi_gauss$mean_vi[seq(t_max+2,t_max+2*n_max, by = 2)])
abline(a = 0, b = 1, col = "red")
#############################################################################
library(vifcopula)
set.seed(0)
t_max = 1000
n_max = 100
k_max = 1
g_max = 10
gid = sample(1:g_max,n_max,replace = T)

datagen_gauss <- fcopsim(t_max = t_max, n_max = n_max, k_max = k_max, gid = gid,
                         family = 1, family_vine = 1, seed_num = 100,
                         structfactor = 4)
datagen <- datagen_gauss

data <- list(u = datagen$u,
             n_max = datagen$n_max,
             t_max = datagen$t_max,
             k_max = datagen$k_max,
             gid = datagen$gid,
             structfactor = datagen$structfactor)
init <- list(copula_type = datagen$family,
             vine_copula_type = datagen$family_vine,
             vine_edges = datagen$egdes)
other <- list(seed = 100909, core = 8, iter = 1000,
              n_monte_carlo_grad = 1, n_monte_carlo_elbo = 100,
              eval_elbo = 100, adapt_bool = F, adapt_val = 1,
              adapt_iterations = 50, tol_rel_obj = 0.01, copselect = F, modelselect = T)
vi_gauss <- vifcopula::vifcop(data,init,other)


################################################################################
