# library(devtools)
# install_github("hoanguc3m/vifcopula")
setwd("/home/hoanguc3m/Dropbox/WP2/")
library(vifcopula)
set.seed(0)
t_max = 1000
n_max = 100
gauss_init <- rep(1, n_max)
copfamily_init <- sample(c(1,2,3,4,5,6),size = 100, replace = T)

datagen_gauss <- fcopsim(t_max = t_max, n_max = n_max, family = 1, seed_num = 0)
datagen <- datagen_gauss
data <- list(u = datagen$u,
             n_max = datagen$n_max,
             t_max = datagen$t_max,
             k_max = datagen$k_max,
             gid = datagen$gid,
             structfactor = datagen$structfactor)
init <- list(copula_type = datagen$family)
other <- list(seed = 126, core = 8, iter = 1000,
              n_monte_carlo_grad = 1, n_monte_carlo_elbo = 100,
              eval_elbo = 100, adapt_bool = F, adapt_val = 1,
              adapt_iterations = 50, tol_rel_obj = 0.01, copselect = F)
vi_gauss <- vifcopula::vifcop(data,init,other)

pdf(file='img/Gaussian1.pdf', width = 9, height = 4.5)
par(mfrow =c(1,2))
par(mar=c(5,5,3,1))
plot(datagen$v, get_v0(vi_gauss),
     xlab = expression(v[gen]), ylab = expression(v[vi]),
     main = " Gaussian one factor copula")
abline(a= 0, b=1, col="red")
plot(datagen$theta, get_theta(vi_gauss) ,
     xlab = expression(theta[gen]), ylab = expression(theta[vi]),
     main = " Gaussian one factor copula")
abline(a= 0, b=1, col="red")
dev.off()

# init <- list(copula_type = copfamily_init)
# other <- list(seed = 126, core = 8, iter = 1000,
#                 n_monte_carlo_grad = 1, n_monte_carlo_elbo = 100,
#                 eval_elbo = 100, adapt_bool = F, adapt_val = 1,
#                 adapt_iterations = 50, tol_rel_obj = 0.01, copselect = T)
# vi_gauss_rng <- vifcopula::vifcop(data,init,other)
# sum(vi_gauss_rng$cop_type == datagen_gauss$family)


################################################################################

datagen_student <- fcopsim(t_max = 1000, n_max = 100, family = 2)
datagen <- datagen_student
data <- list(u = datagen$u,
             n_max = datagen$n_max,
             t_max = datagen$t_max,
             k_max = datagen$k_max,
             gid = datagen$gid,
             structfactor = datagen$structfactor)
init <- list(copula_type = datagen$family)
other <- list(seed = 126, core = 8, iter = 1000,
              n_monte_carlo_grad = 1, n_monte_carlo_elbo = 100,
              eval_elbo = 100, adapt_bool = F, adapt_val = 1,
              adapt_iterations = 50, tol_rel_obj = 0.01, copselect = F)
vi_student <- vifcopula::vifcop(data,init,other)


# init <- list(copula_type = copfamily_init)
# other <- list(seed = 126, core = 8, iter = 1000,
#               n_monte_carlo_grad = 1, n_monte_carlo_elbo = 100,
#               eval_elbo = 100, adapt_bool = F, adapt_val = 1,
#               adapt_iterations = 50, tol_rel_obj = 0.01, copselect = T)
# vi_student_rng <- vifcopula::vifcop(data,init,other)
# sum(vi_student_rng$cop_type == datagen_student$family)
# plot.vifcop(vi_student_rng)
pdf(file='img/Student1.pdf', width = 15, height = 5)
par(mfrow =c(1,3))
par(mar=c(5,5,3,1))

plot(datagen$v, get_v0(vi_student),
     xlab = expression(v[gen]), ylab = expression(v[vi]),
     main = " Student one factor copula")
abline(a= 0, b=1, col="red")

plot(datagen$theta, get_theta(vi_student), xlim = c(0.2,1),ylim = c(0.2,1),
     xlab = expression(theta[gen]), ylab = expression(theta[vi]),
     main = " Student one factor copula")
abline(a= 0, b=1, col="red")

plot(datagen$theta2, get_theta2(vi_student), xlim = c(2,20),ylim = c(2,20),
     xlab = expression(nu[gen]), ylab = expression(nu[vi]),
     main = " Student one factor copula")
abline(a= 0, b=1, col="red")
dev.off()
#save.image("/media/hoanguc3m/Data/wp2/sim1f_student.Rdata")

################################################################################

datagen_clayton <- fcopsim(t_max = 1000, n_max = 100, family = 3)
datagen <- datagen_clayton
data <- list(u = datagen$u,
             n_max = datagen$n_max,
             t_max = datagen$t_max,
             k_max = datagen$k_max,
             gid = datagen$gid,
             structfactor = datagen$structfactor)
init <- list(copula_type = datagen$family)
other <- list(seed = 126, core = 8, iter = 1000,
              n_monte_carlo_grad = 1, n_monte_carlo_elbo = 100,
              eval_elbo = 100, adapt_bool = F, adapt_val = 1,
              adapt_iterations = 50, tol_rel_obj = 0.01, copselect = F)
vi_clayton <- vifcopula::vifcop(data,init,other)

pdf(file='img/Clayton1.pdf', width = 9, height = 4.5)
par(mfrow =c(1,2))
par(mar=c(5,5,3,1))
plot(datagen$v, get_v0(vi_clayton),
     xlab = expression(v[gen]), ylab = expression(v[vi]),
     main = " Clayton one factor copula")

abline(a= 0, b=1, col="red")
plot(datagen$theta, get_theta(vi_clayton),
     xlab = expression(theta[gen]), ylab = expression(theta[vi]),
     main = " Clayton one factor copula")
abline(a= 0, b=1, col="red")
dev.off()

# init <- list(copula_type = copfamily_init)
# other <- list(seed = 126, core = 8, iter = 1000,
#               n_monte_carlo_grad = 1, n_monte_carlo_elbo = 100,
#               eval_elbo = 100, adapt_bool = F, adapt_val = 1,
#               adapt_iterations = 50, tol_rel_obj = 0.01, copselect = T)
# vi_clayton_rng <- vifcopula::vifcop(data,init,other)
#
# sum(vi_clayton_rng$cop_type == datagen_clayton$family)

################################################################################
datagen_gumbel <- fcopsim(t_max = 1000, n_max = 100, family = 4)
datagen <- datagen_gumbel
data <- list(u = datagen$u,
             n_max = datagen$n_max,
             t_max = datagen$t_max,
             k_max = datagen$k_max,
             gid = datagen$gid,
             structfactor = datagen$structfactor)
init <- list(copula_type = datagen$family)
other <- list(seed = 126, core = 8, iter = 1000,
              n_monte_carlo_grad = 1, n_monte_carlo_elbo = 100,
              eval_elbo = 100, adapt_bool = F, adapt_val = 1,
              adapt_iterations = 50, tol_rel_obj = 0.01, copselect = F)
vi_gumbel <- vifcopula::vifcop(data,init,other)

pdf(file='img/Gumbel1.pdf', width = 9, height = 4.5)
par(mfrow =c(1,2))
par(mar=c(5,5,3,1))
plot(datagen$v, get_v0(vi_gumbel),
     xlab = expression(v[gen]), ylab = expression(v[vi]),
     main = " Gumbel one factor copula")

abline(a= 0, b=1, col="red")
plot(datagen$theta, get_theta(vi_gumbel),
     xlab = expression(theta[gen]), ylab = expression(theta[vi]),
     main = " Gumbel one factor copula")
abline(a= 0, b=1, col="red")
dev.off()

# init <- list(copula_type = copfamily_init)
# other <- list(seed = 126, core = 8, iter = 1000,
#               n_monte_carlo_grad = 1, n_monte_carlo_elbo = 100,
#               eval_elbo = 100, adapt_bool = F, adapt_val = 1,
#               adapt_iterations = 50, tol_rel_obj = 0.01, copselect = T)
# vi_gumbel_rng <- vifcopula::vifcop(data,init,other)
# sum(vi_gumbel_rng$cop_type == datagen_gumbel$family)


################################################################################

datagen_frank <- fcopsim(t_max = 1000, n_max = 100, family = 5)
datagen <- datagen_frank
data <- list(u = datagen$u,
             n_max = datagen$n_max,
             t_max = datagen$t_max,
             k_max = datagen$k_max,
             gid = datagen$gid,
             structfactor = datagen$structfactor)
init <- list(copula_type = datagen$family)
other <- list(seed = 126, core = 8, iter = 1000,
              n_monte_carlo_grad = 1, n_monte_carlo_elbo = 100,
              eval_elbo = 100, adapt_bool = F, adapt_val = 1,
              adapt_iterations = 50, tol_rel_obj = 0.01, copselect = F)
vi_frank <- vifcopula::vifcop(data,init,other)

pdf(file='img/Frank1.pdf', width = 9, height = 4.5)
par(mfrow =c(1,2))
par(mar=c(5,5,3,1))
plot(datagen$v, get_v0(vi_frank),
     xlab = expression(v[gen]), ylab = expression(v[vi]),
     main = " Frank one factor copula")

abline(a= 0, b=1, col="red")
plot(datagen$theta, get_theta(vi_frank),
     xlab = expression(theta[gen]), ylab = expression(theta[vi]),
     main = " Frank one factor copula")
abline(a= 0, b=1, col="red")
dev.off()

# init <- list(copula_type = copfamily_init)
# other <- list(seed = 126, core = 8, iter = 1000,
#               n_monte_carlo_grad = 1, n_monte_carlo_elbo = 100,
#               eval_elbo = 100, adapt_bool = F, adapt_val = 1,
#               adapt_iterations = 50, tol_rel_obj = 0.01, copselect = T)
# vi_frank_rng <- vifcopula::vifcop(data,init,other)
# sum(vi_frank_rng$cop_type == datagen_frank$family)

###############################################################################

datagen_joe <- fcopsim(t_max = 1000, n_max = 100, family = 6)
datagen <- datagen_joe
data <- list(u = datagen$u,
             n_max = datagen$n_max,
             t_max = datagen$t_max,
             k_max = datagen$k_max,
             gid = datagen$gid,
             structfactor = datagen$structfactor)
init <- list(copula_type = datagen$family)
other <- list(seed = 126, core = 8, iter = 1000,
              n_monte_carlo_grad = 1, n_monte_carlo_elbo = 100,
              eval_elbo = 100, adapt_bool = F, adapt_val = 1,
              adapt_iterations = 50, tol_rel_obj = 0.01, copselect = F)
vi_joe <- vifcopula::vifcop(data,init,other)

pdf(file='img/Joe1.pdf', width = 9, height = 4.5)
par(mfrow =c(1,2))
par(mar=c(5,5,3,1))
plot(datagen$v, get_v0(vi_joe),
     xlab = expression(v[gen]), ylab = expression(v[vi]),
     main = " Joe one factor copula")

abline(a= 0, b=1, col="red")
plot(datagen$theta, get_theta(vi_joe),
     xlab = expression(theta[gen]), ylab = expression(theta[vi]),
     main = " Joe one factor copula")
abline(a= 0, b=1, col="red")
dev.off()

# init <- list(copula_type = copfamily_init)
# other <- list(seed = 126, core = 8, iter = 1000,
#               n_monte_carlo_grad = 1, n_monte_carlo_elbo = 100,
#               eval_elbo = 100, adapt_bool = F, adapt_val = 1,
#               adapt_iterations = 50, tol_rel_obj = 0.01, copselect = T)
# vi_joe_rng <- vifcopula::vifcop(data,init,other)
# sum(vi_joe_rng$cop_type == datagen_joe$family)


###############################################################################

copfamily = matrix(sample(c(1,2,3,4,5,6),size = 100, replace = T),ncol=1)
datagen_mix <- fcopsim(t_max = 1000, n_max = 100, family = copfamily, family_latent = 0 )
datagen <- datagen_mix
data <- list(u = datagen$u,
             n_max = datagen$n_max,
             t_max = datagen$t_max,
             k_max = datagen$k_max,
             gid = datagen$gid,
             structfactor = datagen$structfactor)
init <- list(copula_type = datagen_mix$family)
other <- list(seed = 126, core = 8, iter = 1000,
              n_monte_carlo_grad = 1, n_monte_carlo_elbo = 100,
              eval_elbo = 100, adapt_bool = F, adapt_val = 1,
              adapt_iterations = 50, tol_rel_obj = 0.01, copselect = F)
vi_mix <- vifcopula::vifcop(data,init,other)

sum(vi_mix$cop_type == datagen_mix$family)

pdf(file='img/Mix1.pdf', width = 9, height = 4.5)
par(mfrow =c(1,2))
par(mar=c(5,5,3,1))
plot(datagen$v, get_v0(vi_mix),
     xlab = expression(v[gen]), ylab = expression(v[vi]),
     main = " Mixed one factor copula")

abline(a= 0, b=1, col="red")
plot(datagen$theta, get_theta(vi_mix),
     xlab = expression(theta[gen]), ylab = expression(theta[vi]),
     main = " Mixed one factor copula")
abline(a= 0, b=1, col="red")
dev.off()

# init <- list(copula_type = copfamily_init)
# other <- list(seed = 126, core = 8, iter = 1000,
#               n_monte_carlo_grad = 1, n_monte_carlo_elbo = 100,
#               eval_elbo = 100, adapt_bool = F, adapt_val = 1,
#               adapt_iterations = 50, tol_rel_obj = 0.01, copselect = T)
# vi_mix_rng <- vifcopula::vifcop(data,init,other)
#
# sum(vi_mix_rng$cop_type == datagen_mix$family)

###############################################################################

# datagen_rotate <- fcopsim(t_max = 1000, n_max = 100, family = 36)
# datagen <- datagen_rotate
# data <- list(u = datagen$u,
#     n_max = datagen$n_max,
#     t_max = datagen$t_max,
#     k_max = datagen$k_max,
#     gid = datagen$gid,
#     structfactor = datagen$structfactor)
# init <- list(copula_type = datagen$family)
# other <- list(seed = 126, core = 8, iter = 1000,
#     n_monte_carlo_grad = 1, n_monte_carlo_elbo = 100,
#     eval_elbo = 100, adapt_bool = F, adapt_val = 1,
#     adapt_iterations = 50, tol_rel_obj = 0.01, copselect = F)
# datagen_rotate <- vifcopula::vifcop(data,init,other)
#
# par(mfrow =c(1,2))
# par(mar=c(5,5,3,1))
# plot(datagen$v, get_v0(datagen_rotate),
#     xlab = expression(v[gen]), ylab = expression(v[vi]),
#     main = " Rotated one factor copula")
#
# abline(a= 0, b=1, col="red")
# plot(datagen$theta, get_theta(datagen_rotate),
#     xlab = expression(theta[gen]), ylab = expression(theta[vi]),
#     main = " Rotated one factor copula")
# abline(a= 0, b=1, col="red")
# dev.off()

#############################################################################
#
# time <- c(vi_gauss$time, vi_student$time, vi_clayton$time, vi_gumbel$time, vi_frank$time, vi_joe$time, vi_mix$time)
# print(time, digits = 0)
#
# ELBO_init <- c(vi_gauss$ELBO, vi_student$ELBO, vi_clayton$ELBO, vi_gumbel$ELBO, vi_frank$ELBO, vi_joe$ELBO, vi_mix$ELBO)
# print(ELBO_init, digits = 1)
#
# iter_num <- c(vi_gauss_rng$iteration, vi_student_rng$iteration, vi_clayton_rng$iteration, vi_gumbel_rng$iteration, vi_frank_rng$iteration, vi_joe_rng$iteration, vi_mix_rng$iteration)
# print(iter_num, digits = 0)
#
# correct_percent <- c(sum(vi_gauss_rng$cop_type == 1),
#                      sum(vi_student_rng$cop_type == 2),
#                      sum(vi_clayton_rng$cop_type == 3),
#                      sum(vi_gumbel_rng$cop_type == 4),
#                      sum(vi_frank_rng$cop_type == 5),
#                      sum(vi_joe_rng$cop_type == 6),
#                      sum(vi_mix_rng$cop_type == datagen_mix$family))
# print(correct_percent, digits = 0)
#
# time_rng <- c(vi_gauss_rng$time, vi_student_rng$time, vi_clayton_rng$time, vi_gumbel_rng$time, vi_frank_rng$time, vi_joe_rng$time, vi_mix_rng$time)
# print(time_rng, digits = 0)
#
# ELBO_rng <- c(vi_gauss_rng$ELBO, vi_student_rng$ELBO, vi_clayton_rng$ELBO, vi_gumbel_rng$ELBO, vi_frank_rng$ELBO, vi_joe_rng$ELBO, vi_mix_rng$ELBO)
# print(ELBO_rng, digits = 1)
#
# pdf(file='img/OnefactorCompare1.pdf', width = 7, height = 2.5)
# par(mfrow =c(3,6))
# par(mar=c(5,5,3,1))
# compare_sim_vi(datagen_mix,vi_mix)
# dev.off()
#
# pdf(file='img/OnefactorCompare2.pdf', width = 7, height = 2.5)
# par(mfrow =c(3,6))
# par(mar=c(5,5,3,1))
# compare_sim_vi(datagen_mix,vi_mix_rng)
# dev.off()


#############################################################################

pdf(file='img/SimOnefactor.pdf', width = 18, height = 9)
par(mfrow =c(3,6))
par(mar=c(5,5,3,1))

plot(datagen_gauss$v, get_v0(vi_gauss),
    xlab = expression(v[gen]), ylab = expression(v[vi]),
    main = " Gaussian one factor copula")
abline(a= 0, b=1, col="red")

plot(datagen_student$v, get_v0(vi_student),
    xlab = expression(v[gen]), ylab = expression(v[vi]),
    main = " Student one factor copula")
abline(a= 0, b=1, col="red")

plot(datagen_clayton$v, get_v0(vi_clayton),
    xlab = expression(v[gen]), ylab = expression(v[vi]),
    main = " Clayton one factor copula")
abline(a= 0, b=1, col="red")

plot(datagen_gumbel$v, get_v0(vi_gumbel),
    xlab = expression(v[gen]), ylab = expression(v[vi]),
    main = " Gumbel one factor copula")
abline(a= 0, b=1, col="red")

plot(datagen_frank$v, get_v0(vi_frank),
    xlab = expression(v[gen]), ylab = expression(v[vi]),
    main = " Frank one factor copula")
abline(a= 0, b=1, col="red")

plot(datagen_joe$v, get_v0(vi_joe),
    xlab = expression(v[gen]), ylab = expression(v[vi]),
    main = " Joe one factor copula")
abline(a= 0, b=1, col="red")





plot(datagen_gauss$theta, get_theta(vi_gauss) ,
    xlab = expression(theta[gen]), ylab = expression(theta[vi]),
    main = " Gaussian one factor copula")
abline(a= 0, b=1, col="red")

plot(datagen_student$theta, get_theta(vi_student),
    xlab = expression(theta[gen]), ylab = expression(theta[vi]),
    main = " Student one factor copula")
abline(a= 0, b=1, col="red")

plot(datagen_clayton$theta, get_theta(vi_clayton),
    xlab = expression(theta[gen]), ylab = expression(theta[vi]),
    main = " Clayton one factor copula")
abline(a= 0, b=1, col="red")

plot(datagen_gumbel$theta, get_theta(vi_gumbel),
    xlab = expression(theta[gen]), ylab = expression(theta[vi]),
    main = " Gumbel one factor copula")
abline(a= 0, b=1, col="red")

plot(datagen_frank$theta, get_theta(vi_frank),
    xlab = expression(theta[gen]), ylab = expression(theta[vi]),
    main = " Frank one factor copula")
abline(a= 0, b=1, col="red")

plot(datagen_joe$theta, get_theta(vi_joe),
    xlab = expression(theta[gen]), ylab = expression(theta[vi]),
    main = " Joe one factor copula")
abline(a= 0, b=1, col="red")

plot.new()

plot(datagen_student$theta2, get_theta2(vi_student), xlim = c(2,20),ylim = c(2,20),
    xlab = expression(nu[gen]), ylab = expression(nu[vi]),
    main = " Student one factor copula")
abline(a= 0, b=1, col="red")

plot.new()
plot.new()
plot.new()
plot.new()
dev.off()
save.image("vi_onefcop_plot.RData")
