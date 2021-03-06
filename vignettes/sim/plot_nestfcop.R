# library(devtools)
# install_github("hoanguc3m/vifcopula")
setwd("/home/hoanguc3m/Dropbox/WP2")
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
                         gid = gid, structfactor = 3, seed_num = 10)

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
              n_monte_carlo_grad = 1, n_monte_carlo_elbo = 100,
              eval_elbo = 100, adapt_bool = F, adapt_val = 1,
              adapt_iterations = 50, tol_rel_obj = 0.01, copselect = F, modelselect = T)

vi_gauss <- vifcopula::vifcop(data,init,other)
compare_sim_vi(datagen, vi_gauss)

# init <- list(copula_type = copfamily_init,
#              latent_copula_type = copfamily_latent_init)
# other <- list(seed = 126, core = 8, iter = 1000,
#               n_monte_carlo_grad = 1, n_monte_carlo_elbo = 100,
#               eval_elbo = 50, adapt_bool = F, adapt_val = 1,
#               adapt_iterations = 50, tol_rel_obj = 0.01, copselect = T)
# vi_gauss_rng <- vifcopula::vifcop(data,init,other)
#
# sum(vi_gauss_rng$cop_type == datagen_gauss$family)
# sum(vi_gauss_rng$latent_copula_type == datagen_gauss$family_latent)
# compare_sim_vi(datagen, vi_gauss_rng)

#################################################################################
datagen_Student <- fcopsim(t_max = 1000, n_max = 100, k_max = k_max, gid = gid,
                           family = 2, seed_num = 0,
                           structfactor = 3, tau_latent_range = c(0.6,0.8),
                           family_latent = sample(c(1,3,4,5,6),size = k_max-1, replace = T) )
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
              n_monte_carlo_grad = 1, n_monte_carlo_elbo = 100,
              eval_elbo = 100, adapt_bool = F, adapt_val = 1,
              adapt_iterations = 50, tol_rel_obj = 0.01, copselect = F)
vi_student <- vifcopula::vifcop(data,init,other)
compare_sim_vi(datagen, vi_student)

# init <- list(copula_type = copfamily_init,
#              latent_copula_type = copfamily_latent_init)
# other <- list(seed = 126, core = 8, iter = 1000,
#               n_monte_carlo_grad = 1, n_monte_carlo_elbo = 100,
#               eval_elbo = 100, adapt_bool = F, adapt_val = 1,
#               adapt_iterations = 50, tol_rel_obj = 0.01, copselect = T)
# vi_student_rng <- vifcopula::vifcop(data,init,other)
# plot.vifcop(vi_student_rng)
# compare_sim_vi(datagen, vi_student_rng)
# sum(vi_student_rng$cop_type == datagen_Student$family)
# sum(vi_student_rng$latent_copula_type == datagen_Student$family_latent)
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
              n_monte_carlo_grad = 1, n_monte_carlo_elbo = 100,
              eval_elbo = 100, adapt_bool = F, adapt_val = 1,
              adapt_iterations = 50, tol_rel_obj = 0.01, copselect = F)
vi_clayton <- vifcopula::vifcop(data,init,other)
# compare_sim_vi(datagen, vi_clayton)
#
# init <- list(copula_type = copfamily_init,
#              latent_copula_type = copfamily_latent_init)
# other <- list(seed = 126, core = 8, iter = 1000,
#               n_monte_carlo_grad = 1, n_monte_carlo_elbo = 100,
#               eval_elbo = 100, adapt_bool = F, adapt_val = 1,
#               adapt_iterations = 50, tol_rel_obj = 0.01, copselect = T)
# vi_clayton_rng <- vifcopula::vifcop(data,init,other)
# compare_sim_vi(datagen, vi_clayton_rng)
#
# sum(vi_clayton_rng$cop_type == datagen_clayton$family)
# sum(vi_clayton_rng$latent_copula_type == datagen_clayton$family_latent)

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
              n_monte_carlo_grad = 1, n_monte_carlo_elbo = 100,
              eval_elbo = 100, adapt_bool = F, adapt_val = 1,
              adapt_iterations = 50, tol_rel_obj = 0.01, copselect = F)
vi_gumbel <- vifcopula::vifcop(data,init,other)
# compare_sim_vi(datagen, vi_gumbel)
#
# init <- list(copula_type = copfamily_init,
#              latent_copula_type = copfamily_latent_init)
# other <- list(seed = 126, core = 8, iter = 1000,
#               n_monte_carlo_grad = 1, n_monte_carlo_elbo = 100,
#               eval_elbo = 100, adapt_bool = F, adapt_val = 1,
#               adapt_iterations = 50, tol_rel_obj = 0.01, copselect = T)
# vi_gumbel_rng <- vifcopula::vifcop(data,init,other)
# compare_sim_vi(datagen, vi_gumbel_rng)
#
# sum(vi_gumbel_rng$cop_type == datagen_gumbel$family)
# sum(vi_gumbel_rng$latent_copula_type == datagen_gumbel$family_latent)

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
              n_monte_carlo_grad = 1, n_monte_carlo_elbo = 100,
              eval_elbo = 100, adapt_bool = F, adapt_val = 1,
              adapt_iterations = 50, tol_rel_obj = 0.01, copselect = F, modelselect = T)
vi_frank <- vifcopula::vifcop(data,init,other)
# compare_sim_vi(datagen, vi_frank)
#
# init <- list(copula_type = copfamily_init,
#              latent_copula_type = copfamily_latent_init)
# other <- list(seed = 126, core = 8, iter = 1000,
#               n_monte_carlo_grad = 1, n_monte_carlo_elbo = 100,
#               eval_elbo = 100, adapt_bool = F, adapt_val = 1,
#               adapt_iterations = 50, tol_rel_obj = 0.01, copselect = T)
# vi_frank_rng <- vifcopula::vifcop(data,init,other)
#
# sum(vi_frank_rng$cop_type == datagen_frank$family)
# sum(vi_frank_rng$latent_copula_type == datagen_frank$family_latent)
# compare_sim_vi(datagen, vi_frank_rng)

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
              n_monte_carlo_grad = 1, n_monte_carlo_elbo = 100,
              eval_elbo = 100, adapt_bool = F, adapt_val = 1,
              adapt_iterations = 50, tol_rel_obj = 0.01, copselect = F)
vi_joe <- vifcopula::vifcop(data,init,other)
# compare_sim_vi(datagen_joe, vi_joe)
#
# init <- list(copula_type = copfamily_init,
#              latent_copula_type = copfamily_latent_init)
# other <- list(seed = 126, core = 8, iter = 1000,
#               n_monte_carlo_grad = 1, n_monte_carlo_elbo = 100,
#               eval_elbo = 100, adapt_bool = F, adapt_val = 1,
#               adapt_iterations = 50, tol_rel_obj = 0.01, copselect = T)
# vi_joe_rng <- vifcopula::vifcop(data,init,other)
# compare_sim_vi(datagen_joe, vi_joe_rng)
# sum(vi_joe_rng$cop_type == datagen_joe$family)
# sum(vi_joe_rng$latent_copula_type == datagen_joe$family_latent)

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
              n_monte_carlo_grad = 1, n_monte_carlo_elbo = 100,
              eval_elbo = 100, adapt_bool = F, adapt_val = 1,
              adapt_iterations = 50, tol_rel_obj = 0.01, copselect = F)
vi_mix <- vifcopula::vifcop(data,init,other)
# compare_sim_vi(datagen, vi_mix)
#
#
# init <- list(copula_type = copfamily_init,
#              latent_copula_type = copfamily_latent_init)
# other <- list(seed = 126, core = 8, iter = 1000,
#               n_monte_carlo_grad = 1, n_monte_carlo_elbo = 100,
#               eval_elbo = 100, adapt_bool = F, adapt_val = 1,
#               adapt_iterations = 50, tol_rel_obj = 0.01, copselect = T)
# vi_mix_rng <- vifcopula::vifcop(data,init,other)
# compare_sim_vi(datagen, vi_mix_rng)
#
# sum(vi_mix_rng$cop_type == datagen_mix$family)
# sum(vi_mix_rng$latent_copula_type == datagen_mix$family_latent)

#################################################################################
#################################################################################
#
# time <- c(vi_gauss$time, vi_student$time, vi_clayton$time, vi_gumbel$time, vi_frank$time, vi_joe$time, vi_mix$time)
# print(time, digits = 0)
#
# ELBO_init <- c(vi_gauss$criteria[1], vi_student$criteria[1], vi_clayton$criteria[1], vi_gumbel$criteria[1], vi_frank$criteria[1], vi_joe$criteria[1], vi_mix$criteria[1])
# print(ELBO_init, digits = 1)
#
# iter_num <- c(vi_gauss_rng$iteration, vi_student_rng$iteration, vi_clayton_rng$iteration, vi_gumbel_rng$iteration, vi_frank_rng$iteration, vi_joe_rng$iteration, vi_mix_rng$iteration)
# print(iter_num, digits = 0)
#
# correct_percent <- c(sum(vi_gauss_rng$cop_type == 1),
#     sum(vi_student_rng$cop_type == 2),
#     sum(vi_clayton_rng$cop_type == 3),
#     sum(vi_gumbel_rng$cop_type == 4),
#     sum(vi_frank_rng$cop_type == 5),
#     sum(vi_joe_rng$cop_type == 6),
#     sum(vi_mix_rng$cop_type == datagen_mix$family))
# print(correct_percent, digits = 0)
#
# correct_latent_percent <- c(sum(vi_gauss_rng$latent_copula_type == 1),
#     sum(vi_student_rng$latent_copula_type == 2),
#     sum(vi_clayton_rng$latent_copula_type == 3),
#     sum(vi_gumbel_rng$latent_copula_type == 4),
#     sum(vi_frank_rng$latent_copula_type == 5),
#     sum(vi_joe_rng$latent_copula_type == 6),
#     sum(vi_mix_rng$latent_copula_type == datagen_mix$family_latent))
# print(correct_latent_percent, digits = 0)
#
# time_rng <- c(vi_gauss_rng$time, vi_student_rng$time, vi_clayton_rng$time, vi_gumbel_rng$time, vi_frank_rng$time, vi_joe_rng$time, vi_mix_rng$time)
# print(time_rng, digits = 0)
#
# ELBO_rng <- c(vi_gauss_rng$criteria[1], vi_student_rng$criteria[1], vi_clayton_rng$criteria[1], vi_gumbel_rng$criteria[1], vi_frank_rng$criteria[1], vi_joe_rng$criteria[1], vi_mix_rng$criteria[1])
# AIC_rng <- c(vi_gauss_rng$criteria[2], vi_student_rng$criteria[2], vi_clayton_rng$criteria[2], vi_gumbel_rng$criteria[2], vi_frank_rng$criteria[2], vi_joe_rng$criteria[2], vi_mix_rng$criteria[2])
# BIC_rng <- c(vi_gauss_rng$criteria[3], vi_student_rng$criteria[3], vi_clayton_rng$criteria[3], vi_gumbel_rng$criteria[3], vi_frank_rng$criteria[3], vi_joe_rng$criteria[3], vi_mix_rng$criteria[3])
# logP_rng <- c(vi_gauss_rng$criteria[4], vi_student_rng$criteria[4], vi_clayton_rng$criteria[4], vi_gumbel_rng$criteria[4], vi_frank_rng$criteria[4], vi_joe_rng$criteria[4], vi_mix_rng$criteria[4])
#
#
# rng_tab <- rbind(iter_num, correct_percent, ELBO_rng, AIC_rng, BIC_rng, logP_rng)
#
# #print(ELBO_rng, digits = 1)
# print(xtable(rng_tab, digits = 0))


#############################################################################

pdf(file='img/Nestedfactor.pdf', width = 19, height = 12)
par(mfrow =c(4,7))
par(mar=c(5,5,3,1))
cex_main = 1.25
cex_lab = 2
cex_axis = 0.75

plot(datagen_gauss$v[,1], get_v0(vi_gauss),
    xlab = expression(v0[true]), ylab = expression(v0[vi]),
    main = " Gaussian nested factor copula", cex.main=cex_main, cex.lab=cex_lab, cex.axis=cex_axis)
abline(a= 0, b=1, col="red")

plot(datagen_Student$v[,1], get_v0(vi_student),
    xlab = expression(v0[true]), ylab = expression(v0[vi]),
    main = " Student nested factor copula", cex.main=cex_main, cex.lab=cex_lab, cex.axis=cex_axis)
abline(a= 0, b=1, col="red")

plot(datagen_clayton$v[,1], get_v0(vi_clayton),
    xlab = expression(v0[true]), ylab = expression(v0[vi]),
    main = " Clayton nested factor copula", cex.main=cex_main, cex.lab=cex_lab, cex.axis=cex_axis)
abline(a= 0, b=1, col="red")

plot(datagen_gumbel$v[,1], get_v0(vi_gumbel),
    xlab = expression(v0[true]), ylab = expression(v0[vi]),
    main = " Gumbel nested factor copula", cex.main=cex_main, cex.lab=cex_lab, cex.axis=cex_axis)
abline(a= 0, b=1, col="red")

plot(datagen_frank$v[,1], get_v0(vi_frank),
    xlab = expression(v0[true]), ylab = expression(v0[vi]),
    main = " Frank nested factor copula", cex.main=cex_main, cex.lab=cex_lab, cex.axis=cex_axis)
abline(a= 0, b=1, col="red")

plot(datagen_joe$v[,1], get_v0(vi_joe),
    xlab = expression(v0[true]), ylab = expression(v0[vi]),
    main = " Joe nested factor copula", cex.main=cex_main, cex.lab=cex_lab, cex.axis=cex_axis)
abline(a= 0, b=1, col="red")

plot(datagen_BB1$v[,1], get_v0(vi_BB1),
     xlab = expression(v0[true]), ylab = expression(v0[vi]),
     main = " BB1 nested factor copula", cex.main=cex_main, cex.lab=cex_lab, cex.axis=cex_axis)
abline(a= 0, b=1, col="red")


plot(datagen_gauss$v[,2:k_max], get_v(vi_gauss),
    xlab = expression(vg[true]), ylab = expression(vg[vi]),
    main = "", cex.main=cex_main, cex.lab=cex_lab, cex.axis=cex_axis)
#main = " Gaussian nested factor copula")
abline(a= 0, b=1, col="red")

plot(datagen_Student$v[,2:k_max], get_v(vi_student),
    xlab = expression(vg[true]), ylab = expression(vg[vi]),
    main = "", cex.main=cex_main, cex.lab=cex_lab, cex.axis=cex_axis)
#main = " Student nested factor copula")
abline(a= 0, b=1, col="red")

plot(datagen_clayton$v[,2:k_max], get_v(vi_clayton),
    xlab = expression(vg[true]), ylab = expression(vg[vi]),
    main = "", cex.main=cex_main, cex.lab=cex_lab, cex.axis=cex_axis)
#main = " Clayton nested factor copula")
abline(a= 0, b=1, col="red")

plot(datagen_gumbel$v[,2:k_max], get_v(vi_gumbel),
    xlab = expression(vg[true]), ylab = expression(vg[vi]),
    main = "", cex.main=cex_main, cex.lab=cex_lab, cex.axis=cex_axis)
#main = " Gumbel nested factor copula")
abline(a= 0, b=1, col="red")

plot(datagen_frank$v[,2:k_max], get_v(vi_frank),
    xlab = expression(vg[true]), ylab = expression(vg[vi]),
    main = "", cex.main=cex_main, cex.lab=cex_lab, cex.axis=cex_axis)
#main = " Frank nested factor copula")
abline(a= 0, b=1, col="red")

plot(datagen_joe$v[,2:k_max], get_v(vi_joe),
    xlab = expression(vg[true]), ylab = expression(vg[vi]),
    main = "", cex.main=cex_main, cex.lab=cex_lab, cex.axis=cex_axis)
#main = " Joe nested factor copula")
abline(a= 0, b=1, col="red")

plot(datagen_BB1$v[,2:k_max], get_v(vi_BB1),
     xlab = expression(vg[true]), ylab = expression(vg[vi]),
     main = "", cex.main=cex_main, cex.lab=cex_lab, cex.axis=cex_axis)
#main = " BB1 nested factor copula")
abline(a= 0, b=1, col="red")



plot(datagen_gauss$theta, get_theta(vi_gauss) ,
    xlab = expression(theta[true]), ylab = expression(theta[vi]),
    main = "", cex.main=cex_main, cex.lab=cex_lab, cex.axis=cex_axis)
#main = " Gaussian nested factor copula")
abline(a= 0, b=1, col="red")

plot(datagen_Student$theta, get_theta(vi_student),
    xlab = expression(theta[true]), ylab = expression(theta[vi]),
    main = "", cex.main=cex_main, cex.lab=cex_lab, cex.axis=cex_axis)
#main = " Student nested factor copula")
abline(a= 0, b=1, col="red")

plot(datagen_clayton$theta, get_theta(vi_clayton),
    xlab = expression(theta[true]), ylab = expression(theta[vi]),
    main = "", cex.main=cex_main, cex.lab=cex_lab, cex.axis=cex_axis)
#main = " Clayton nested factor copula")
abline(a= 0, b=1, col="red")

plot(datagen_gumbel$theta, get_theta(vi_gumbel),
    xlab = expression(theta[true]), ylab = expression(theta[vi]),
    main = "", cex.main=cex_main, cex.lab=cex_lab, cex.axis=cex_axis)
#main = " Gumbel nested factor copula")
abline(a= 0, b=1, col="red")

plot(datagen_frank$theta, get_theta(vi_frank),
    xlab = expression(theta[true]), ylab = expression(theta[vi]),
    main = "", cex.main=cex_main, cex.lab=cex_lab, cex.axis=cex_axis)
#main = " Frank nested factor copula")
abline(a= 0, b=1, col="red")

plot(datagen_joe$theta, get_theta(vi_joe),
    xlab = expression(theta[true]), ylab = expression(theta[vi]),
    main = "", cex.main=cex_main, cex.lab=cex_lab, cex.axis=cex_axis)
#main = " Joe nested factor copula")
abline(a= 0, b=1, col="red")

plot(datagen_BB1$theta, get_theta(vi_BB1),
     xlab = expression(theta[true]), ylab = expression(theta[vi]),
     main = "", cex.main=cex_main, cex.lab=cex_lab, cex.axis=cex_axis)
#main = " BB1 nested factor copula")
abline(a= 0, b=1, col="red")

plot.new()

plot(datagen_Student$theta2, get_theta2(vi_student), xlim = c(2,20),ylim = c(2,20),
    xlab = expression(nu[true]), ylab = expression(nu[vi]),
    main = "", cex.main=cex_main, cex.lab=cex_lab, cex.axis=cex_axis)
#main = " Student nested factor copula")
abline(a= 0, b=1, col="red")

plot.new()
plot.new()
plot.new()
plot.new()

plot(datagen_BB1$theta2, get_theta2(vi_BB1),
     xlab = expression(delta[true]), ylab = expression(delta[vi]),
     main = "", cex.main=cex_main, cex.lab=cex_lab, cex.axis=cex_axis)
#main = " BB1 nested factor copula")
abline(a= 0, b=1, col="red")

dev.off()
#save.image("vi_nestfcop_plot.RData")

