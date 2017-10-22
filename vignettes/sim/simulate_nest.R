library(devtools)
install_github("hoanguc3m/vifcopula")
#setwd("/home/hoanguc3m/Dropbox/WP2/")
library(vifcopula)
set.seed(0)
t_max = 1000
n_max = 100
k_max = 6
gid = sample(1:(k_max-1),n_max,replace = T)

gauss_init <- matrix(1, nrow = n_max, ncol = 1)
gauss_latent_init <- matrix(1, nrow = k_max-1, ncol = 1)

copfamily_init <- sample(c(1,2,3,4,5,6),size = 100, replace = T)
copfamily_latent_init <- sample(c(1,2,3,4,5,6),size = k_max-1, replace = T)

datagen_gauss <- fcopsim(t_max = t_max, n_max = n_max, k_max = k_max,
                         family = 1, family_latent = 1,
                         gid = gid, structfactor = 3, seed_num = 0)

datagen <- datagen_gauss
data <- list(u = datagen$u,
             n_max = datagen$n_max,
             t_max = datagen$t_max,
             k_max = datagen$k_max,
             gid = datagen$gid,
             structfactor = datagen$structfactor)
init <- list(copula_type = datagen$family,
             latent_copula_type = datagen$family_latent)
other <- list(seed = 126, core = 8, iter = 1000,
              n_monte_carlo_grad = 1, n_monte_carlo_elbo = 10,
              eval_elbo = 100, adapt_bool = F, adapt_val = 1,
              adapt_iterations = 50, tol_rel_obj = 0.1, copselect = F)

vi_gauss <- vifcopula::vifcop(data,init,other)
comparefcop(datagen, vi_gauss)

init <- list(copula_type = copfamily_init,
             latent_copula_type = copfamily_latent_init)
other <- list(seed = 126, core = 8, iter = 1000,
              n_monte_carlo_grad = 1, n_monte_carlo_elbo = 10,
              eval_elbo = 50, adapt_bool = F, adapt_val = 1,
              adapt_iterations = 50, tol_rel_obj = 0.1, copselect = T)
vi_gauss_rng <- vifcopula::vifcop(data,init,other)

sum(vi_gauss_rng$cop_type == datagen_gauss$family)
sum(vi_gauss_rng$latent_copula_type == datagen_gauss$family_latent)
comparefcop(datagen, vi_gauss_rng)

#save.image("/media/hoanguc3m/Data/wp2/simnf_gauss.Rdata")

#################################################################################
datagen_Student <- fcopsim(t_max = 1000, n_max = 100, k_max = k_max, gid = gid,
                           family = 2, family_latent = 2,  seed_num = 100,
                           structfactor = 3)
datagen <- datagen_Student
data <- list(u = datagen$u,
             n_max = datagen$n_max,
             t_max = datagen$t_max,
             k_max = datagen$k_max,
             gid = datagen$gid,
             structfactor = datagen$structfactor)
init <- list(copula_type = datagen$family,
             latent_copula_type = datagen$family_latent)
other <- list(seed = 126, core = 8, iter = 1000,
              n_monte_carlo_grad = 1, n_monte_carlo_elbo = 10,
              eval_elbo = 100, adapt_bool = F, adapt_val = 1,
              adapt_iterations = 50, tol_rel_obj = 0.1, copselect = F)
vi_student <- vifcopula::vifcop(data,init,other)
comparefcop(datagen, vi_student)

plot(datagen$theta, abs(tail(vi_student$mean_iv,2*n_max)[seq(1,200, by = 2)]), xlab = expression(theta[t]), ylab = expression(theta[approximated]))
plot(datagen$theta2, tail(vi_student$mean_iv,2*n_max)[seq(2,200, by = 2)], xlab = expression(theta[t]), ylab = expression(theta[approximated]))
plot(datagen$theta_latent, vi_student$mean_iv[seq(t_max*k_max+1,t_max*k_max+k_max*2-2,by = 2)])
plot(datagen$theta2_latent, vi_student$mean_iv[seq(t_max*k_max+2,t_max*k_max+k_max*2-1,by = 2)])


abline(a= 0, b=1, col="red")

init <- list(copula_type = copfamily_init,
             latent_copula_type = copfamily_latent_init)
other <- list(seed = 126, core = 8, iter = 1000,
              n_monte_carlo_grad = 1, n_monte_carlo_elbo = 10,
              eval_elbo = 100, adapt_bool = F, adapt_val = 1,
              adapt_iterations = 50, tol_rel_obj = 0.1, copselect = T)
vi_student_rng <- vifcopula::vifcop(data,init,other)
comparefcop(datagen, vi_student_rng)
sum(vi_student_rng$cop_type == datagen_Student$family)
sum(vi_student_rng$latent_copula_type == datagen_Student$family_latent)
#save.image("/media/hoanguc3m/Data/wp2/simnf_student.Rdata")
#################################################################################

#################################################################################
datagen_clayton <- fcopsim(t_max = 1000, n_max = 100, k_max = k_max, gid = gid,
                           family = 3, family_latent = 3, seed_num = 10,
                           structfactor = 3)
datagen <- datagen_clayton
data <- list(u = datagen$u,
             n_max = datagen$n_max,
             t_max = datagen$t_max,
             k_max = datagen$k_max,
             gid = datagen$gid,
             structfactor = datagen$structfactor)
init <- list(copula_type = datagen$family,
             latent_copula_type = datagen$family_latent)
other <- list(seed = 126, core = 8, iter = 1000,
              n_monte_carlo_grad = 1, n_monte_carlo_elbo = 10,
              eval_elbo = 100, adapt_bool = F, adapt_val = 1,
              adapt_iterations = 50, tol_rel_obj = 0.1, copselect = F)
vi_clayton <- vifcopula::vifcop(data,init,other)
comparefcop(datagen, vi_clayton)

init <- list(copula_type = copfamily_init,
             latent_copula_type = copfamily_latent_init)
other <- list(seed = 126, core = 8, iter = 1000,
              n_monte_carlo_grad = 1, n_monte_carlo_elbo = 10,
              eval_elbo = 100, adapt_bool = F, adapt_val = 1,
              adapt_iterations = 50, tol_rel_obj = 0.1, copselect = T)
vi_clayton_rng <- vifcopula::vifcop(data,init,other)
comparefcop(datagen, vi_clayton_rng)

sum(vi_clayton_rng$cop_type == datagen_clayton$family)
sum(vi_clayton_rng$latent_copula_type == datagen_clayton$family_latent)
#save.image("/media/hoanguc3m/Data/wp2/simnf_clayton.Rdata")

#################################################################################
datagen_gumbel <- fcopsim(t_max = 1000, n_max = 100, k_max = k_max, gid = gid,
                          family = 4, family_latent = 4, seed_num = 100,
                          structfactor = 3)
datagen <- datagen_gumbel
data <- list(u = datagen$u,
             n_max = datagen$n_max,
             t_max = datagen$t_max,
             k_max = datagen$k_max,
             gid = datagen$gid,
             structfactor = datagen$structfactor)
init <- list(copula_type = datagen$family,
             latent_copula_type = datagen$family_latent)
other <- list(seed = 126, core = 8, iter = 1000,
              n_monte_carlo_grad = 1, n_monte_carlo_elbo = 10,
              eval_elbo = 100, adapt_bool = F, adapt_val = 1,
              adapt_iterations = 50, tol_rel_obj = 0.1, copselect = F)
vi_gumbel <- vifcopula::vifcop(data,init,other)
comparefcop(datagen, vi_gumbel)

init <- list(copula_type = copfamily_init,
             latent_copula_type = copfamily_latent_init)
other <- list(seed = 126, core = 8, iter = 1000,
              n_monte_carlo_grad = 1, n_monte_carlo_elbo = 10,
              eval_elbo = 100, adapt_bool = F, adapt_val = 1,
              adapt_iterations = 50, tol_rel_obj = 0.1, copselect = T)
vi_gumbel_rng <- vifcopula::vifcop(data,init,other)
comparefcop(datagen, vi_gumbel_rng)

sum(vi_gumbel_rng$cop_type == datagen$family)
sum(vi_gumbel_rng$latent_copula_type == datagen_gumbel$family_latent)
#save.image("/media/hoanguc3m/Data/wp2/simnf_gumbel.Rdata")

#################################################################################

datagen_frank <- fcopsim(t_max = 1000, n_max = 100, k_max = k_max, gid = gid,
                         family = 5, family_latent = 5, seed_num = 100,
                         structfactor = 3)
datagen <- datagen_frank
data <- list(u = datagen$u,
             n_max = datagen$n_max,
             t_max = datagen$t_max,
             k_max = datagen$k_max,
             gid = datagen$gid,
             structfactor = datagen$structfactor)
init <- list(copula_type = datagen$family,
             latent_copula_type = datagen$family_latent)
other <- list(seed = 126, core = 8, iter = 1000,
              n_monte_carlo_grad = 1, n_monte_carlo_elbo = 10,
              eval_elbo = 100, adapt_bool = F, adapt_val = 1,
              adapt_iterations = 50, tol_rel_obj = 0.1, copselect = F)
vi_frank <- vifcopula::vifcop(data,init,other)
comparefcop(datagen, vi_frank)

init <- list(copula_type = copfamily_init,
             latent_copula_type = copfamily_latent_init)
other <- list(seed = 126, core = 8, iter = 1000,
              n_monte_carlo_grad = 1, n_monte_carlo_elbo = 10,
              eval_elbo = 100, adapt_bool = F, adapt_val = 1,
              adapt_iterations = 50, tol_rel_obj = 0.1, copselect = T)
vi_frank_rng <- vifcopula::vifcop(data,init,other)

sum(vi_frank_rng$cop_type == datagen$family)
sum(vi_frank_rng$latent_copula_type == datagen$family_latent)
comparefcop(datagen, vi_frank_rng)
#save.image("/media/hoanguc3m/Data/wp2/simnf_frank.Rdata")

#################################################################################
datagen_joe <- fcopsim(t_max = 1000, n_max = 100, k_max = k_max, gid = gid,
                       family = 6, family_latent = 6,  seed_num = 100,
                       structfactor = 3)
datagen <- datagen_joe
data <- list(u = datagen$u,
             n_max = datagen$n_max,
             t_max = datagen$t_max,
             k_max = datagen$k_max,
             gid = datagen$gid,
             structfactor = datagen$structfactor)
init <- list(copula_type = datagen$family,
             latent_copula_type = datagen$family_latent)
other <- list(seed = 126, core = 8, iter = 1000,
              n_monte_carlo_grad = 1, n_monte_carlo_elbo = 10,
              eval_elbo = 100, adapt_bool = F, adapt_val = 1,
              adapt_iterations = 50, tol_rel_obj = 0.1, copselect = F)
vi_joe <- vifcopula::vifcop(data,init,other)
comparefcop(datagen, vi_joe)

init <- list(copula_type = copfamily_init,
             latent_copula_type = copfamily_latent_init)
other <- list(seed = 126, core = 8, iter = 1000,
              n_monte_carlo_grad = 1, n_monte_carlo_elbo = 10,
              eval_elbo = 100, adapt_bool = F, adapt_val = 1,
              adapt_iterations = 50, tol_rel_obj = 0.1, copselect = T)
vi_joe_rng <- vifcopula::vifcop(data,init,other)
comparefcop(datagen, vi_joe_rng)
sum(vi_joe_rng$cop_type == datagen$family)
sum(vi_joe_rng$latent_copula_type == datagen$family_latent)
#save.image("/media/hoanguc3m/Data/wp2/simnf_joe.Rdata")
#################################################################################
copfamily = sample(c(1,2,3,4,5,6),size = n_max, replace = T)
latentcopfamily = sample(c(1,2,3,4,5,6),size = k_max-1, replace = T)

datagen_mix <- fcopsim(t_max = 1000, n_max = 100, k_max = k_max, gid = gid,
                       family = copfamily, family_latent = latentcopfamily,
                       seed_num = 100, structfactor = 3)
datagen <- datagen_mix
data <- list(u = datagen$u,
             n_max = datagen$n_max,
             t_max = datagen$t_max,
             k_max = datagen$k_max,
             gid = datagen$gid,
             structfactor = datagen$structfactor)
init <- list(copula_type = datagen$family,
             latent_copula_type = datagen$family_latent)
other <- list(seed = 126, core = 8, iter = 1000,
              n_monte_carlo_grad = 1, n_monte_carlo_elbo = 10,
              eval_elbo = 100, adapt_bool = F, adapt_val = 1,
              adapt_iterations = 50, tol_rel_obj = 0.1, copselect = F)
vi_mix <- vifcopula::vifcop(data,init,other)
comparefcop(datagen, vi_mix)


init <- list(copula_type = copfamily_init,
             latent_copula_type = copfamily_latent_init)
other <- list(seed = 126, core = 8, iter = 1000,
              n_monte_carlo_grad = 1, n_monte_carlo_elbo = 10,
              eval_elbo = 100, adapt_bool = F, adapt_val = 1,
              adapt_iterations = 50, tol_rel_obj = 0.1, copselect = T)
vi_mix_rng <- vifcopula::vifcop(data,init,other)
comparefcop(datagen, vi_mix_rng)

sum(vi_mix_rng$cop_type == datagen$family)
sum(vi_mix_rng$latent_copula_type == datagen$family_latent)
#save.image("/media/hoanguc3m/Data/wp2/simnf_mix.Rdata")

#################################################################################
#################################################################################
