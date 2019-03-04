library(vifcopula)
set.seed(0)
t_max = 1000
n_max = 100
k_max = 1
g_max = 10
gid = sample(1:g_max,n_max,replace = T)

################################################################################

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
compare_sim_vi(datagen_gauss, vi_gauss)

plot(datagen$v, head(vi_gauss$mean_vi,t_max))
abline(a = 0, b = 1, col = "red")
plot(datagen$theta, vi_gauss$mean_vi[(t_max+1):(t_max+n_max)])
abline(a = 0, b = 1, col = "red")
plot(datagen$theta_vine, tail(vi_gauss$mean_vi,n_max - max(gid)))
abline(a = 0, b = 1, col = "red")

#############################################################################

datagen_student <- fcopsim(t_max = t_max, n_max = n_max, k_max = k_max, gid = gid,
                           family = 2, seed_num = 100,
                           structfactor = 4, tau_range = c(0.2,0.8),
                           tau_latent_range = c(0.2,0.8), family_vine = sample(c(1,3,4,5,6),size = n_max, replace = T))
datagen <- datagen_student
data <- list(u = datagen$u,
             n_max = datagen$n_max,
             t_max = datagen$t_max,
             k_max = datagen$k_max,
             gid = datagen$gid,
             structfactor = datagen$structfactor)
init <- list(copula_type = datagen$family,
             vine_copula_type = datagen$family_vine,
             vine_edges = datagen$egdes)
other <- list(seed = 126, core = 8, iter = 1000,
              n_monte_carlo_grad = 1, n_monte_carlo_elbo = 100,
              eval_elbo = 100, adapt_bool = F, adapt_val = 1,
              adapt_iterations = 50, tol_rel_obj = 0.01, copselect = F)
vi_student <- vifcopula::vifcop(data,init,other)

plot(datagen$theta, vi_gauss$mean_vi[seq(t_max+1,t_max+2*n_max, by = 2)])
abline(a = 0, b = 1, col = "red")
plot(datagen$theta2, vi_gauss$mean_vi[seq(t_max+2,t_max+2*n_max, by = 2)])
abline(a = 0, b = 1, col = "red")

# compare_sim_vi(datagen_student, vi_student)
# plot.vifcop(vi_student)

# init <- list(copula_type = copfamily_rng,
#              latent_copula_type = latentcopfamily_rng)
# other <- list(seed = 126, core = 8, iter = 1000,
#               n_monte_carlo_grad = 1, n_monte_carlo_elbo = 100,
#               eval_elbo = 100, adapt_bool = F, adapt_val = 1,
#               adapt_iterations = 50, tol_rel_obj = 0.01, copselect = T)
# vi_student_rng <- vifcopula::vifcop(data,init,other)
# compare_sim_vi(datagen_student, vi_student_rng)
#
# sum(vi_student_rng$cop_type == datagen_student$family)
# sum(vi_student_rng$latent_copula_type == datagen_student$family_vine)
#

#################################################################################

datagen_clayton <- fcopsim(t_max = 1000, n_max = 100, k_max = k_max, gid = gid,
                           family = 3, family_vine = 3, seed_num = 0,
                           structfactor = 4)
datagen <- datagen_clayton
data <- list(u = datagen$u,
             n_max = datagen$n_max,
             t_max = datagen$t_max,
             k_max = datagen$k_max,
             gid = datagen$gid,
             structfactor = datagen$structfactor)
init <- list(copula_type = datagen$family,
             vine_copula_type = datagen$family_vine,
             vine_edges = datagen$egdes)
other <- list(seed = 126, core = 8, iter = 1000,
              n_monte_carlo_grad = 1, n_monte_carlo_elbo = 100,
              eval_elbo = 100, adapt_bool = F, adapt_val = 1,
              adapt_iterations = 50, tol_rel_obj = 0.01, copselect = F)
vi_clayton <- vifcopula::vifcop(data,init,other)
# compare_sim_vi(datagen_clayton, vi_clayton)

# init <- list(copula_type = copfamily_rng,
#             latent_copula_type = latentcopfamily_rng)
# other <- list(seed = 126, core = 8, iter = 1000,
#               n_monte_carlo_grad = 1, n_monte_carlo_elbo = 100,
#               eval_elbo = 100, adapt_bool = F, adapt_val = 1,
#               adapt_iterations = 50, tol_rel_obj = 0.01, copselect = T)
# vi_clayton_rng <- vifcopula::vifcop(data,init,other)
# compare_sim_vi(datagen_clayton, vi_clayton_rng)
# sum(vi_clayton_rng$cop_type == datagen_clayton$family)
# sum(vi_clayton_rng$latent_copula_type == datagen_clayton$family_vine)

#################################################################################
datagen_gumbel <- fcopsim(t_max = 1000, n_max = 100, k_max = k_max, gid = gid,
                          family = 4, family_vine = 4, seed_num = 100,
                          structfactor = 4)
datagen <- datagen_gumbel
data <- list(u = datagen$u,
             n_max = datagen$n_max,
             t_max = datagen$t_max,
             k_max = datagen$k_max,
             gid = datagen$gid,
             structfactor = datagen$structfactor)
init <- list(copula_type = datagen$family,
             vine_copula_type = datagen$family_vine,
             vine_edges = datagen$egdes)
other <- list(seed = 126, core = 8, iter = 1000,
              n_monte_carlo_grad = 1, n_monte_carlo_elbo = 100,
              eval_elbo = 100, adapt_bool = F, adapt_val = 1,
              adapt_iterations = 50, tol_rel_obj = 0.01, copselect = F)
vi_gumbel <- vifcopula::vifcop(data,init,other)
# compare_sim_vi(datagen_gumbel, vi_gumbel)
#
# init <- list(copula_type = copfamily_rng,
#              latent_copula_type = latentcopfamily_rng)
# other <- list(seed = 126, core = 8, iter = 1000,
#               n_monte_carlo_grad = 1, n_monte_carlo_elbo = 100,
#               eval_elbo = 100, adapt_bool = F, adapt_val = 1,
#               adapt_iterations = 50, tol_rel_obj = 0.01, copselect = T)
# vi_gumbel_rng <- vifcopula::vifcop(data,init,other)
# compare_sim_vi(datagen_gumbel, vi_gumbel_rng)
#
# sum(vi_gumbel_rng$cop_type == datagen_gumbel$family)
# sum(vi_gumbel_rng$latent_copula_type == datagen_gumbel$family_vine)


#################################################################################

datagen_frank <- fcopsim(t_max = 1000, n_max = 100, k_max = k_max, gid = gid,
                         family = 5, family_vine = 5, seed_num = 100,
                         structfactor = 4)
datagen <- datagen_frank
data <- list(u = datagen$u,
             n_max = datagen$n_max,
             t_max = datagen$t_max,
             k_max = datagen$k_max,
             gid = datagen$gid,
             structfactor = datagen$structfactor)
init <- list(copula_type = datagen$family,
             vine_copula_type = datagen$family_vine,
             vine_edges = datagen$egdes)
other <- list(seed = 126, core = 8, iter = 1000,
              n_monte_carlo_grad = 1, n_monte_carlo_elbo = 100,
              eval_elbo = 100, adapt_bool = F, adapt_val = 1,
              adapt_iterations = 50, tol_rel_obj = 0.01, copselect = F, modelselect = T)
vi_frank <- vifcopula::vifcop(data,init,other)
# compare_sim_vi(datagen_frank, vi_frank)
#
# init <- list(copula_type = copfamily_rng,
#              latent_copula_type = latentcopfamily_rng)
# other <- list(seed = 112321, core = 8, iter = 1000,
#               n_monte_carlo_grad = 1, n_monte_carlo_elbo = 100,
#               eval_elbo = 100, adapt_bool = F, adapt_val = 1,
#               adapt_iterations = 50, tol_rel_obj = 0.01, copselect = T)
# vi_frank_rng <- vifcopula::vifcop(data,init,other)
# compare_sim_vi(datagen_frank, vi_frank_rng)
#
# sum(vi_frank_rng$cop_type == datagen_frank$family)
# sum(vi_frank_rng$latent_copula_type == datagen_frank$family_vine)

#################################################################################
datagen_joe <- fcopsim(t_max = 1000, n_max = 100, k_max = k_max, gid = gid,
                       family = 6, family_vine = 6, seed_num = 100,
                       structfactor = 4)
datagen <- datagen_joe
data <- list(u = datagen$u,
             n_max = datagen$n_max,
             t_max = datagen$t_max,
             k_max = datagen$k_max,
             gid = datagen$gid,
             structfactor = datagen$structfactor)
init <- list(copula_type = datagen$family,
             vine_copula_type = datagen$family_vine,
             vine_edges = datagen$egdes)
other <- list(seed = 126, core = 8, iter = 1000,
              n_monte_carlo_grad = 1, n_monte_carlo_elbo = 100,
              eval_elbo = 100, adapt_bool = F, adapt_val = 1,
              adapt_iterations = 50, tol_rel_obj = 0.01, copselect = F)
vi_joe <- vifcopula::vifcop(data,init,other)
# compare_sim_vi(datagen_joe, vi_joe)

# init <- list(copula_type = copfamily_rng,
#              latent_copula_type = latentcopfamily_rng)
# other <- list(seed = 126, core = 8, iter = 1000,
#               n_monte_carlo_grad = 1, n_monte_carlo_elbo = 100,
#               eval_elbo = 100, adapt_bool = F, adapt_val = 1,
#               adapt_iterations = 50, tol_rel_obj = 0.01, copselect = T)
# vi_joe_rng <- vifcopula::vifcop(data,init,other)
# compare_sim_vi(datagen_joe, vi_joe_rng)
#
# sum(vi_joe_rng$cop_type == datagen_joe$family)
# sum(vi_joe_rng$latent_copula_type == datagen_joe$family_vine)

#################################################################################
copfamily = sample(c(1,2,3,4,5,6), size = n_max, replace = T)
vinecopfamily = sample(c(1,3,4,5,6),size = n_max - g_max, replace = T)

datagen_mix <- fcopsim(t_max = 1000, n_max = 100, k_max = k_max, gid = gid,
                       family = copfamily, family_vine = vinecopfamily,
                       seed_num = 100, structfactor = 4)

datagen <- datagen_mix
data <- list(u = datagen$u,
             n_max = datagen$n_max,
             t_max = datagen$t_max,
             k_max = datagen$k_max,
             gid = datagen$gid,
             structfactor = datagen$structfactor)
init <- list(copula_type = datagen$family,
             vine_copula_type = datagen$family_vine,
             vine_edges = datagen$egdes)
other <- list(seed = 126, core = 8, iter = 1000,
              n_monte_carlo_grad = 1, n_monte_carlo_elbo = 100,
              eval_elbo = 100, adapt_bool = F, adapt_val = 1,
              adapt_iterations = 50, tol_rel_obj = 0.01, copselect = F)
vi_mix <- vifcopula::vifcop(data,init,other)

plot(datagen$v, head(vi_mix$mean_vi,t_max))
abline(a = 0, b = 1, col = "red")
plot(datagen$theta, vi_mix$mean_vi[(t_max+1):(t_max+n_max)])
abline(a = 0, b = 1, col = "red")
plot(datagen$theta_vine, tail(vi_mix$mean_vi,n_max - max(gid)))
abline(a = 0, b = 1, col = "red")


# compare_sim_vi(datagen_mix, vi_mix)
#
#
# init <- list(copula_type = copfamily_rng,
#              latent_copula_type = latentcopfamily_rng)
#
# other <- list(seed = 126, core = 8, iter = 1000,
#               n_monte_carlo_grad = 1, n_monte_carlo_elbo = 100,
#               eval_elbo = 100, adapt_bool = F, adapt_val = 1,
#               adapt_iterations = 50, tol_rel_obj = 0.01, copselect = T)
# vi_mix_rng <- vifcopula::vifcop(data,init,other)
#
# sum(vi_mix_rng$cop_type == datagen$family)
# sum(vi_mix_rng$latent_copula_type == datagen$family_vine)
# compare_sim_vi(datagen_mix, vi_mix_rng)
#############################################################################

pdf(file='img/truncatedfv.pdf', width = 18, height = 12)
par(mfrow =c(4,6))
par(mar=c(5,5,3,1))

plot(datagen_gauss$v[,1], get_v0(vi_gauss),
     xlab = expression(v[gen]), ylab = expression(v[vi]),
     main = " Gaussian truncated factor vine copula")
abline(a= 0, b=1, col="red")

plot(datagen_student$v[,1], get_v0(vi_student),
     xlab = expression(v[gen]), ylab = expression(v[vi]),
     main = " Student truncated factor vine copula")
abline(a= 0, b=1, col="red")

plot(datagen_clayton$v[,1], get_v0(vi_clayton),
     xlab = expression(v[gen]), ylab = expression(v[vi]),
     main = " Clayton truncated factor vine copula")
abline(a= 0, b=1, col="red")

plot(datagen_gumbel$v[,1], get_v0(vi_gumbel),
     xlab = expression(v[gen]), ylab = expression(v[vi]),
     main = " Gumbel truncated factor vine copula")
abline(a= 0, b=1, col="red")

plot(datagen_frank$v[,1], get_v0(vi_frank),
     xlab = expression(v[gen]), ylab = expression(v[vi]),
     main = " Frank truncated factor vine copula")
abline(a= 0, b=1, col="red")

plot(datagen_joe$v[,1], get_v0(vi_joe),
     xlab = expression(v[gen]), ylab = expression(v[vi]),
     main = " Joe truncated factor vine copula")
abline(a= 0, b=1, col="red")




plot(datagen_gauss$v[,2:k_max], get_v(vi_gauss),
     xlab = expression(v[gen]), ylab = expression(v[vi]),
     main = " Gaussian truncated factor vine copula")
abline(a= 0, b=1, col="red")

plot(datagen_student$v[,2:k_max], get_v(vi_student),
     xlab = expression(v[gen]), ylab = expression(v[vi]),
     main = " Student truncated factor vine copula")
abline(a= 0, b=1, col="red")

plot(datagen_clayton$v[,2:k_max], get_v(vi_clayton),
     xlab = expression(v[gen]), ylab = expression(v[vi]),
     main = " Clayton truncated factor vine copula")
abline(a= 0, b=1, col="red")

plot(datagen_gumbel$v[,2:k_max], get_v(vi_gumbel),
     xlab = expression(v[gen]), ylab = expression(v[vi]),
     main = " Gumbel truncated factor vine copula")
abline(a= 0, b=1, col="red")

plot(datagen_frank$v[,2:k_max], get_v(vi_frank),
     xlab = expression(v[gen]), ylab = expression(v[vi]),
     main = " Frank truncated factor vine copula")
abline(a= 0, b=1, col="red")

plot(datagen_joe$v[,2:k_max], get_v(vi_joe),
     xlab = expression(v[gen]), ylab = expression(v[vi]),
     main = " Joe truncated factor vine copula")
abline(a= 0, b=1, col="red")



plot(datagen_gauss$theta, get_theta(vi_gauss) ,
     xlab = expression(theta[gen]), ylab = expression(theta[vi]),
     main = " Gaussian truncated factor vine copula")
abline(a= 0, b=1, col="red")

plot(datagen_student$theta, get_theta(vi_student),
     xlab = expression(theta[gen]), ylab = expression(theta[vi]),
     main = " Student truncated factor vine copula")
abline(a= 0, b=1, col="red")

plot(datagen_clayton$theta, get_theta(vi_clayton),
     xlab = expression(theta[gen]), ylab = expression(theta[vi]),
     main = " Clayton truncated factor vine copula")
abline(a= 0, b=1, col="red")

plot(datagen_gumbel$theta, get_theta(vi_gumbel),
     xlab = expression(theta[gen]), ylab = expression(theta[vi]),
     main = " Gumbel truncated factor vine copula")
abline(a= 0, b=1, col="red")

plot(datagen_frank$theta, get_theta(vi_frank),
     xlab = expression(theta[gen]), ylab = expression(theta[vi]),
     main = " Frank truncated factor vine copula")
abline(a= 0, b=1, col="red")

plot(datagen_joe$theta, get_theta(vi_joe),
     xlab = expression(theta[gen]), ylab = expression(theta[vi]),
     main = " Joe truncated factor vine copula")
abline(a= 0, b=1, col="red")



plot(datagen_gauss$theta_latent, get_latent_theta(vi_gauss) ,
     xlab = expression(theta[gen]), ylab = expression(theta[vi]),
     main = " Gaussian truncated factor vine copula")
abline(a= 0, b=1, col="red")

plot(datagen_student$theta_latent, get_latent_theta(vi_student),
     xlab = expression(theta[gen]), ylab = expression(theta[vi]),
     main = " Student truncated factor vine copula")
abline(a= 0, b=1, col="red")

plot(datagen_clayton$theta_latent, get_latent_theta(vi_clayton),
     xlab = expression(theta[gen]), ylab = expression(theta[vi]),
     main = " Clayton truncated factor vine copula")
abline(a= 0, b=1, col="red")

plot(datagen_gumbel$theta_latent, get_latent_theta(vi_gumbel),
     xlab = expression(theta[gen]), ylab = expression(theta[vi]),
     main = " Gumbel truncated factor vine copula")
abline(a= 0, b=1, col="red")

plot(datagen_frank$theta_latent, get_latent_theta(vi_frank),
     xlab = expression(theta[gen]), ylab = expression(theta[vi]),
     main = " Frank truncated factor vine copula")
abline(a= 0, b=1, col="red")

plot(datagen_joe$theta_latent, get_latent_theta(vi_joe),
     xlab = expression(theta[gen]), ylab = expression(theta[vi]),
     main = " Joe truncated factor vine copula")
abline(a= 0, b=1, col="red")


# plot.new()
#
# plot(datagen_student$theta2, get_theta2(vi_student), xlim = c(2,20),ylim = c(2,20),
#     xlab = expression(nu[gen]), ylab = expression(nu[vi]),
#     main = " Student truncated factor vine copula")
# abline(a= 0, b=1, col="red")
#
# plot.new()
# plot.new()
# plot.new()
# plot.new()
#
dev.off()

save.image("vi_truncated_fv_plot.RData")

