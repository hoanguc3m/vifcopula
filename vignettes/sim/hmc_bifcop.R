setwd("/home/hoanguc3m/Dropbox/WP2/")
library(vifcopula)
seed_num = 5623
set.seed(seed_num)
t_max = 1000
n_max = 100
k_max = 6
gid = sample(1:(k_max-1),n_max,replace = T)

copfamily_rng = sample(c(1,3,4,5,6), size = n_max, replace = T)
latentcopfamily_rng = sample(c(1,3,4,5,6),size = n_max, replace = T)

gauss_init <- matrix(1, nrow = n_max, ncol = 1)
gauss_latent_init <- matrix(1, nrow = n_max, ncol = 1)

datagen_gauss <- fcopsim(t_max = t_max, n_max = n_max, k_max = k_max, gid = gid,
                         family = 1, family_latent = 1, seed_num = seed_num,
                         structfactor = 2)
datagen <- datagen_gauss
data <- list(u = datagen$u,
             n_max = datagen$n_max,
             t_max = datagen$t_max,
             k_max = datagen$k_max,
             gid = datagen$gid,
             structfactor = datagen$structfactor)
init <- list(copula_type = datagen$family,
             latent_copula_type = datagen$family_latent)
other <- list(seed = 126, core = 8, num_warmup = 500, num_samples = 1000, copselect = F)
hmc_gauss <- vifcopula::hmcfcop(data,init,other)
vi_gauss <- vifcopula::vifcop(data,init,other)
#
# plot(hmc_gauss)
# compare_vi_hmc(vi_gauss, hmc_gauss)
################################################################################

# datagen_student <- fcopsim(t_max = t_max, n_max = n_max, k_max = k_max, gid = gid,
#                            family = 2, seed_num = seed_num,
#                            structfactor = 2, tau_range = c(0.2,0.8),
#                            tau_latent_range = c(0.2,0.8), family_latent = sample(c(1,3,4,5,6),size = n_max, replace = T))
# datagen <- datagen_student
# data <- list(u = datagen$u,
#              n_max = datagen$n_max,
#              t_max = datagen$t_max,
#              k_max = datagen$k_max,
#              gid = datagen$gid,
#              structfactor = datagen$structfactor)
# init <- list(copula_type = datagen$family,
#              latent_copula_type = datagen$family_latent)
# other <- list(seed = 126, core = 8, num_warmup = 500, num_samples = 1000, copselect = F)
# hmc_student <- vifcopula::hmcfcop(data,init,other)
# vi_student <- vifcopula::vifcop(data,init,other)

################################################################################

datagen_clayton <- fcopsim(t_max = 1000, n_max = 100, k_max = k_max, gid = gid,
                           family = 3, family_latent = 3, seed_num = seed_num,
                           tau_range = c(0.2,0.7), tau_latent_range = c(0.2,0.7),
                           structfactor = 2)
datagen <- datagen_clayton
data <- list(u = datagen$u,
             n_max = datagen$n_max,
             t_max = datagen$t_max,
             k_max = datagen$k_max,
             gid = datagen$gid,
             structfactor = datagen$structfactor)
init <- list(copula_type = datagen$family,
             latent_copula_type = datagen$family_latent)
other <- list(seed = 126, core = 8, num_warmup = 500, num_samples = 1000, copselect = F)
hmc_clayton <- vifcopula::hmcfcop(data,init,other)
vi_clayton <- vifcopula::vifcop(data,init,other)

################################################################################
datagen_gumbel <- fcopsim(t_max = 1000, n_max = 100, k_max = k_max, gid = gid,
                          family = 4, family_latent = 4, seed_num = seed_num,
                          structfactor = 2)
datagen <- datagen_gumbel
data <- list(u = datagen$u,
             n_max = datagen$n_max,
             t_max = datagen$t_max,
             k_max = datagen$k_max,
             gid = datagen$gid,
             structfactor = datagen$structfactor)
init <- list(copula_type = datagen$family,
             latent_copula_type = datagen$family_latent)
other <- list(seed = 126, core = 8, num_warmup = 500, num_samples = 1000, copselect = F)
hmc_gumbel <- vifcopula::hmcfcop(data,init,other)
vi_gumbel <- vifcopula::vifcop(data,init,other)

################################################################################

datagen_frank <- fcopsim(t_max = 1000, n_max = 100, k_max = k_max, gid = gid,
                         family = 5, family_latent = 5, seed_num = seed_num,
                         structfactor = 2)
datagen <- datagen_frank
data <- list(u = datagen$u,
             n_max = datagen$n_max,
             t_max = datagen$t_max,
             k_max = datagen$k_max,
             gid = datagen$gid,
             structfactor = datagen$structfactor)
init <- list(copula_type = datagen$family,
             latent_copula_type = datagen$family_latent)
other <- list(seed = 126, core = 8, num_warmup = 500, num_samples = 1000, copselect = F)
hmc_frank <- vifcopula::hmcfcop(data,init,other)
vi_frank <- vifcopula::vifcop(data,init,other)

###############################################################################

datagen_joe <- fcopsim(t_max = 1000, n_max = 100, k_max = k_max, gid = gid,
                       family = 6, family_latent = 6, seed_num = seed_num,
                       structfactor = 2)
datagen <- datagen_joe
data <- list(u = datagen$u,
             n_max = datagen$n_max,
             t_max = datagen$t_max,
             k_max = datagen$k_max,
             gid = datagen$gid,
             structfactor = datagen$structfactor)
init <- list(copula_type = datagen$family,
             latent_copula_type = datagen$family_latent)
other <- list(seed = 126, core = 8, num_warmup = 500, num_samples = 1000, copselect = F)
hmc_joe <- vifcopula::hmcfcop(data,init,other)
vi_joe <- vifcopula::vifcop(data,init,other)

###############################################################################
latentcopfamily = sample(c(1,3,4,5,6),size = n_max, replace = T)

datagen_BB1 <- fcopsim(t_max = 1000, n_max = 100, k_max = k_max, gid = gid,
                       family = 7, family_latent = latentcopfamily, seed_num = seed_num,
                       structfactor = 2)
datagen <- datagen_BB1
data <- list(u = datagen$u,
             n_max = datagen$n_max,
             t_max = datagen$t_max,
             k_max = datagen$k_max,
             gid = datagen$gid,
             structfactor = datagen$structfactor)
init <- list(copula_type = datagen$family,
             latent_copula_type = datagen$family_latent)
other <- list(seed = 126, core = 8, num_warmup = 500, num_samples = 1000, copselect = F)
hmc_BB1 <- vifcopula::hmcfcop(data,init,other)
vi_BB1 <- vifcopula::vifcop(data,init,other)

###############################################################################

copfamily = sample(c(1,2,3,4,5,6), size = n_max, replace = T)
latentcopfamily = sample(c(1,3,4,5,6),size = n_max, replace = T)

datagen_mix <- fcopsim(t_max = 1000, n_max = 100, k_max = k_max, gid = gid,
                       family = copfamily, family_latent = latentcopfamily,
                       seed_num = seed_num, structfactor = 2)

datagen <- datagen_mix
data <- list(u = datagen$u,
             n_max = datagen$n_max,
             t_max = datagen$t_max,
             k_max = datagen$k_max,
             gid = datagen$gid,
             structfactor = datagen$structfactor)
init <- list(copula_type = datagen$family,
             latent_copula_type = datagen$family_latent)
other <- list(seed = 126, core = 8, num_warmup = 500, num_samples = 1000, copselect = F)
hmc_mix <- vifcopula::hmcfcop(data,init,other)
vi_mix <- vifcopula::vifcop(data,init,other)


###############################################################################

# Mix selected
###############################################################################

copfamily = c(rep(1, 40), rep(3, 10), rep(4, 10), rep(5, 10), rep(6, 10), rep(2, 10), rep(7, 10))
latentcopfamily = sample(c(1,3,4,5),size = n_max, replace = T)

datagen_mix <- fcopsim(t_max = 1000, n_max = 100, k_max = k_max, gid = gid,
                       family = copfamily, family_latent = latentcopfamily,
                       seed_num = 100, structfactor = 2)
tail(BiCopPar2Tau(family = datagen_mix$family, par = datagen_mix$theta, par2 = datagen_mix$theta2),10)

datagen <- datagen_mix
data <- list(u = datagen$u,
             n_max = datagen$n_max,
             t_max = datagen$t_max,
             k_max = datagen$k_max,
             gid = datagen$gid,
             structfactor = datagen$structfactor)
init <- list(copula_type = datagen$family,
             latent_copula_type = datagen$family_latent)
other <- list(seed = 126, core = 8, num_warmup = 2000, num_samples = 1000, copselect = F,
              tol_rel_obj = 0.005)
hmc_mix <- vifcopula::hmcfcop(data,init,other)
vi_mix <- vifcopula::vifcop(data,init,other)


###############################################################################

time_vi <- c(vi_gauss$time, vi_student$time, vi_clayton$time, vi_gumbel$time, vi_frank$time, vi_joe$time, vi_BB1$time, vi_mix$time)
time_hmc <- c(hmc_gauss$time, hmc_student$time, hmc_clayton$time, hmc_gumbel$time, hmc_frank$time, hmc_joe$time, hmc_BB1$time, hmc_mix$time)

time_tab <- rbind(time_vi, time_hmc)
print(xtable(time_tab, digits = 0))

compare_vi_hmc(vi_gauss, hmc_gauss)
compare_vi_hmc(vi_student, hmc_student)
compare_vi_hmc(vi_clayton, hmc_clayton)
compare_vi_hmc(vi_gumbel, hmc_gumbel)
compare_vi_hmc(vi_frank, hmc_frank)
compare_vi_hmc(vi_joe, hmc_joe)
compare_vi_hmc(vi_mix, hmc_mix)

#############################################################################

pdf(file='img/VIHMCbifcop.pdf', width = 19, height = 12)
par(mfrow =c(4,7))
par(mar=c(5,5,3,1))
cex_main = 1.25
cex_lab = 2
cex_axis = 0.75

plot(get_v0_sd(hmc_gauss), get_v0_sd(vi_gauss),
     xlab = expression(sd(v0[hmc]) ), ylab = expression(sd (v0[vi])),
     main = " Gaussian bi-factor copula", cex.main=cex_main, cex.lab=cex_lab, cex.axis=cex_axis)
abline(a= 0, b=1, col="red")

plot(get_v0_sd(hmc_student), get_v0_sd(vi_student),
     xlab = expression(sd(v0[hmc]) ), ylab = expression(sd (v0[vi])),
     main = " Student bi-factor copula", cex.main=cex_main, cex.lab=cex_lab, cex.axis=cex_axis)
abline(a= 0, b=1, col="red")

plot(get_v0_sd(hmc_clayton), get_v0_sd(vi_clayton),
     xlab = expression(sd(v0[hmc]) ), ylab = expression(sd (v0[vi])),
     main = " Clayton bi-factor copula", cex.main=cex_main, cex.lab=cex_lab, cex.axis=cex_axis)
abline(a= 0, b=1, col="red")

plot(get_v0_sd(hmc_gumbel), get_v0_sd(vi_gumbel),
     xlab = expression(sd(v0[hmc]) ), ylab = expression(sd (v0[vi])),
     main = " Gumbel bi-factor copula", cex.main=cex_main, cex.lab=cex_lab, cex.axis=cex_axis)
abline(a= 0, b=1, col="red")

plot(get_v0_sd(hmc_frank), get_v0_sd(vi_frank),
     xlab = expression(sd(v0[hmc]) ), ylab = expression(sd (v0[vi])),
     main = " Frank bi-factor copula", cex.main=cex_main, cex.lab=cex_lab, cex.axis=cex_axis)
abline(a= 0, b=1, col="red")

plot(get_v0_sd(hmc_joe), get_v0_sd(vi_joe),
     xlab = expression(sd(v0[hmc]) ), ylab = expression(sd (v0[vi])),
     main = " Joe bi-factor copula", cex.main=cex_main, cex.lab=cex_lab, cex.axis=cex_axis)
abline(a= 0, b=1, col="red")

plot(get_v0_sd(hmc_BB1), get_v0_sd(vi_BB1),
     xlab = expression(sd(v0[hmc]) ), ylab = expression(sd (v0[vi])),
     main = " BB1 bi-factor copula", cex.main=cex_main, cex.lab=cex_lab, cex.axis=cex_axis)
abline(a= 0, b=1, col="red")

plot(get_v_sd(hmc_gauss), get_v_sd(vi_gauss),
     xlab = expression(sd(vg[hmc]) ), ylab = expression(sd (vg[vi])),
     main = "", cex.main=cex_main, cex.lab=cex_lab, cex.axis=cex_axis)
#main = " Gaussian bi-factor copula")
abline(a= 0, b=1, col="red")

plot(get_v_sd(hmc_student), get_v_sd(vi_student),
     xlab = expression(sd(vg[hmc]) ), ylab = expression(sd (vg[vi])),
     main = "", cex.main=cex_main, cex.lab=cex_lab, cex.axis=cex_axis)
#main = " Student bi-factor copula")
abline(a= 0, b=1, col="red")

plot(get_v_sd(hmc_clayton), get_v_sd(vi_clayton),
     xlab = expression(sd(vg[hmc]) ), ylab = expression(sd (vg[vi])),
     main = "", cex.main=cex_main, cex.lab=cex_lab, cex.axis=cex_axis)
#main = " Clayton bi-factor copula")
abline(a= 0, b=1, col="red")

plot(get_v_sd(hmc_gumbel), get_v_sd(vi_gumbel),
     xlab = expression(sd(vg[hmc]) ), ylab = expression(sd (vg[vi])),
     main = "", cex.main=cex_main, cex.lab=cex_lab, cex.axis=cex_axis)
#main = " Gumbel bi-factor copula")
abline(a= 0, b=1, col="red")

plot(get_v_sd(hmc_frank), get_v_sd(vi_frank),
     xlab = expression(sd(vg[hmc]) ), ylab = expression(sd (vg[vi])),
     main = "", cex.main=cex_main, cex.lab=cex_lab, cex.axis=cex_axis)
#main = " Frank bi-factor copula")
abline(a= 0, b=1, col="red")

plot(get_v_sd(hmc_joe), get_v_sd(vi_joe),
     xlab = expression(sd(vg[hmc]) ), ylab = expression(sd (vg[vi])),
     main = "", cex.main=cex_main, cex.lab=cex_lab, cex.axis=cex_axis)
#main = " Joe bi-factor copula")
abline(a= 0, b=1, col="red")

plot(get_v_sd(hmc_BB1), get_v_sd(vi_BB1),
     xlab = expression(sd(vg[hmc]) ), ylab = expression(sd (vg[vi])),
     main = "", cex.main=cex_main, cex.lab=cex_lab, cex.axis=cex_axis)
#main = " BB1 bi-factor copula")
abline(a= 0, b=1, col="red")

plot(get_theta_sd(hmc_gauss), get_theta_sd(vi_gauss),
     xlab = expression(sd(theta[hmc]) ), ylab = expression(sd (theta[vi])),
     main = "", cex.main=cex_main, cex.lab=cex_lab, cex.axis=cex_axis)
#main = " Gaussian bi-factor copula")
abline(a= 0, b=1, col="red")

plot(get_theta_sd(hmc_student), get_theta_sd(vi_student),
     xlab = expression(sd(theta[hmc]) ), ylab = expression(sd (theta[vi])),
     main = "", cex.main=cex_main, cex.lab=cex_lab, cex.axis=cex_axis)
#main = " Student bi-factor copula")
abline(a= 0, b=1, col="red")

plot(get_theta_sd(hmc_clayton), get_theta_sd(vi_clayton),
     xlab = expression(sd(theta[hmc]) ), ylab = expression(sd (theta[vi])),
     main = "", cex.main=cex_main, cex.lab=cex_lab, cex.axis=cex_axis)
#main = " Clayton bi-factor copula")
abline(a= 0, b=1, col="red")

plot(get_theta_sd(hmc_gumbel), get_theta_sd(vi_gumbel),
     xlab = expression(sd(theta[hmc]) ), ylab = expression(sd (theta[vi])),
     main = "", cex.main=cex_main, cex.lab=cex_lab, cex.axis=cex_axis)
#main = " Gumbel bi-factor copula")
abline(a= 0, b=1, col="red")

plot(get_theta_sd(hmc_frank), get_theta_sd(vi_frank),
     xlab = expression(sd(theta[hmc]) ), ylab = expression(sd (theta[vi])),
     main = "", cex.main=cex_main, cex.lab=cex_lab, cex.axis=cex_axis)
#main = " Frank bi-factor copula")
abline(a= 0, b=1, col="red")

plot(get_theta_sd(hmc_joe), get_theta_sd(vi_joe),
     xlab = expression(sd(theta[hmc]) ), ylab = expression(sd (theta[vi])),
     main = "", cex.main=cex_main, cex.lab=cex_lab, cex.axis=cex_axis)
#main = " Joe bi-factor copula")
abline(a= 0, b=1, col="red")


plot(get_theta_sd(hmc_BB1), get_theta_sd(vi_BB1),
     xlab = expression(sd(theta[hmc]) ), ylab = expression(sd (theta[vi])),
     main = "", cex.main=cex_main, cex.lab=cex_lab, cex.axis=cex_axis)
#main = " BB1 bi-factor copula")
abline(a= 0, b=1, col="red")

plot(get_latent_theta_sd(hmc_gauss), get_latent_theta_sd(vi_gauss) ,
     xlab = expression(sd(theta[gtrue])), ylab = expression(sd (theta[gvi])),
     main = "", cex.main=cex_main, cex.lab=cex_lab, cex.axis=cex_axis)
#main = " Gaussian bi-factor copula")
abline(a= 0, b=1, col="red")


plot(get_latent_theta_sd(hmc_student), get_latent_theta_sd(vi_student) ,
     xlab = expression(sd(theta[gtrue])), ylab = expression(sd (theta[gvi])),
     main = "", cex.main=cex_main, cex.lab=cex_lab, cex.axis=cex_axis)
#main = " Student bi-factor copula")
abline(a= 0, b=1, col="red")

plot(get_latent_theta_sd(hmc_clayton), get_latent_theta_sd(vi_clayton) ,
     xlab = expression(sd(theta[gtrue])), ylab = expression(sd (theta[gvi])),
     main = "", cex.main=cex_main, cex.lab=cex_lab, cex.axis=cex_axis)
#main = " Clayton bi-factor copula")
abline(a= 0, b=1, col="red")

plot(get_latent_theta_sd(hmc_gumbel), get_latent_theta_sd(vi_gumbel) ,
     xlab = expression(sd(theta[gtrue])), ylab = expression(sd (theta[gvi])),
     main = "", cex.main=cex_main, cex.lab=cex_lab, cex.axis=cex_axis)
#main = " Gumbel bi-factor copula")
abline(a= 0, b=1, col="red")

plot(get_latent_theta_sd(hmc_frank), get_latent_theta_sd(vi_frank) ,
     xlab = expression(sd(theta[gtrue])), ylab = expression(sd (theta[gvi])),
     main = "", cex.main=cex_main, cex.lab=cex_lab, cex.axis=cex_axis)
#main = " Frank bi-factor copula")
abline(a= 0, b=1, col="red")

plot(get_latent_theta_sd(hmc_joe), get_latent_theta_sd(vi_joe) ,
     xlab = expression(sd(theta[gtrue])), ylab = expression(sd (theta[gvi])),
     main = "", cex.main=cex_main, cex.lab=cex_lab, cex.axis=cex_axis)
#main = " Joe bi-factor copula")
abline(a= 0, b=1, col="red")

plot(get_latent_theta_sd(hmc_BB1), get_latent_theta_sd(vi_BB1) ,
     xlab = expression(sd(theta[gtrue])), ylab = expression(sd (theta[gvi])),
     main = "", cex.main=cex_main, cex.lab=cex_lab, cex.axis=cex_axis)
#main = " BB1 bi-factor copula")
abline(a= 0, b=1, col="red")
dev.off()

###############################################################################
hmc_col <- ncol(hmc_mix$sample_hmc)
hmc_sample <- hmc_mix$sample_hmc[,8:hmc_col]
v_idx <- c(100,700, 1000)
theta_idx <- c(10,101, 102)

paranum <- cumsum(rep(1,100) + (datagen_mix$family == 2) + (datagen_mix$family == 7))
a1 <- which( paranum == theta_idx[1]);  datagen_mix$family[a1] # check here
a2 <- which( paranum == theta_idx[2]);  datagen_mix$family[a2] # check here
a3 <- which( paranum == theta_idx[3]);  datagen_mix$family[a3] # check here
theta_sample_idx <- k_max * t_max + paranum[c(a1,a2,a3)] # check here

theta_sample_idx <- k_max * t_max + theta_idx # check here

hmc_sample_select <- hmc_sample[,c(v_idx,theta_sample_idx)]
vi_sample_select <- vi_mix$sample_vi[,c(v_idx,theta_sample_idx)]

labels_ii = c(expression('v'[100]),expression('v'[500]), expression('v'[900]),
              expression(theta[10]^G[p]),expression(theta[90]^BB1), expression(delta[90]^BB1))


true_val <- c(datagen_mix$v[v_idx], datagen_mix$theta[a1], datagen_mix$theta[a2], datagen_mix$theta[a3])
true_val <- c(datagen_mix$v[v_idx], datagen_mix$theta[a1], datagen_mix$theta[a3], datagen_mix$theta2[a3])
true_val
apply(hmc_sample_select, MARGIN = 2, FUN = mean)
apply(vi_sample_select, MARGIN = 2, FUN = mean)

# pdf(file='img/scatterHMCVI.pdf', width = 7, height = 7)
# pairs(rbind(hmc_sample_select, vi_sample_select),
#       pch=rep(c(1,2), c(nrow(hmc_sample_select), nrow(vi_sample_select))),
#       col=rep(c(rgb(red = 0, green = 0, blue = 0, alpha = 0.4),
#                 rgb(red = 1, green = 0, blue = 0, alpha = 0.2)), c(nrow(hmc_sample_select), nrow(vi_sample_select))),
#       labels = c(expression('v'[10]),expression('v'[100]), expression('v'[1000]),
#                  expression(theta[1]),expression(theta[10]), expression(theta[100])),
#       upper.panel = NULL)
#
# dev.off()

pdf(file='img/scatterHMCVI.pdf', width = 8, height = 7)


library("KernSmooth")
i = 1; j = 2
maxplot = 6
par(mfrow = c(maxplot,maxplot))
par(mar=c(2,2,0,0))
library(RColorBrewer)
k <- 11
colfunc_hmc <- colorRampPalette(c(rgb(0,0.2,1,0.5), rgb(0,0.2,1,1)), alpha = TRUE)
colfunc_vi <- colorRampPalette(c(rgb(1,0,0,0.25), rgb(1,0,0,0.75)), alpha = TRUE)
my.cols <- colfunc_hmc(k)
my.cols.vi <- colfunc_vi(k)


for (i in c(1:maxplot)){
        for (j in c(1:maxplot)){
                if (i == j) {
                        plot(c(0, 1), c(0, 1), ann = F, bty = 'n', type = 'n', xaxt = 'n', yaxt = 'n')
                        text(x = 0.4, y = 0.25, labels_ii[i],
                             cex = 2.5, col = "black")
                }
                if (i < j) {
                        plot.new()
                }
                if (i > j) {
                        xaxt_n = 'n'
                        yaxt_n = 'n'
                        if (i == maxplot) {xaxt_n = 's'}
                        if (j == 1) {yaxt_n = 's'}
                        x <- data.frame(X = hmc_sample_select[,j], Y = hmc_sample_select[,i])
                        est <- bkde2D(x, bandwidth=sapply(x,dpik)*1.4,
                                      range.x=list(c(min(x[,1]),max(x[,1])),
                                                   c(min(x[,2]),max(x[,2]))))

                        contour(est$x1, est$x2, est$fhat,
                                drawlabels = F, nlevels = 7, lty = 2, lwd = 1.5,
                                yaxt=yaxt_n, xaxt=xaxt_n,
                                col = my.cols)
                        points(true_val[j], true_val[i], pch = 8, col = "black", cex = 2)

                        x <- data.frame(X = vi_sample_select[,j], Y = vi_sample_select[,i])
                        est <- bkde2D(x, bandwidth=sapply(x,dpik)*1.4,
                                      range.x=list(c(min(x[,1]),max(x[,1])),
                                                   c(min(x[,2]),max(x[,2]))))
                        contour(est$x1, est$x2, est$fhat,
                                drawlabels = F, nlevels = 7, col=my.cols.vi,
                                add=TRUE, lwd = 1.5)



                }
        }
}
dev.off()

par(mfrow=c(1,1))

####################################################################
# Gumbel and Student
####################################################################
hmc_col <- ncol(hmc_mix$sample_hmc)
hmc_sample <- hmc_mix$sample_hmc[,8:hmc_col]
v_idx <- c(100,700, 1000)
theta_idx <- c(1,2, 3)

paranum <- cumsum(rep(1,100) + (datagen_mix$family == 2) + (datagen_mix$family == 7))
a1 <- which( paranum == theta_idx[1]);  datagen_mix$family[a1] # check here
a2 <- which( paranum == theta_idx[2]);  datagen_mix$family[a2] # check here
a3 <- which( paranum == theta_idx[3]);  datagen_mix$family[a3] # check here
theta_sample_idx <- k_max * t_max + paranum[c(a1,a2,a3)] # check here

theta_sample_idx <- k_max * t_max + theta_idx # check here

hmc_sample_select <- hmc_sample[,c(v_idx,theta_sample_idx)]
vi_sample_select <- vi_mix$sample_vi[,c(v_idx,theta_sample_idx)]

labels_ii = c(expression('v'[100]),expression('v'[500]), expression('v'[900]),
              expression(theta[50]),expression(theta[80]^1), expression(theta[80]^2))


true_val <- c(datagen_mix$v[v_idx], datagen_mix$theta[a1], datagen_mix$theta[a2], datagen_mix$theta[a3])
true_val <- c(datagen_mix$v[v_idx], datagen_mix$theta[a1], datagen_mix$theta[a3], datagen_mix$theta2[a3])
true_val
apply(hmc_sample_select, MARGIN = 2, FUN = mean)
apply(vi_sample_select, MARGIN = 2, FUN = mean)

# pdf(file='img/scatterHMCVIv2.pdf', width = 8, height = 7)


library("KernSmooth")
i = 1; j = 2
maxplot = 6
par(mfrow = c(maxplot,maxplot))
par(mar=c(2,2,0,0))
library(RColorBrewer)
k <- 11
colfunc_hmc <- colorRampPalette(c(rgb(0,0.2,1,0.5), rgb(0,0.2,1,1)), alpha = TRUE)
colfunc_vi <- colorRampPalette(c(rgb(1,0,0,0.25), rgb(1,0,0,0.75)), alpha = TRUE)
my.cols <- colfunc_hmc(k)
my.cols.vi <- colfunc_vi(k)


for (i in c(1:maxplot)){
        for (j in c(1:maxplot)){
                if (i == j) {
                        plot(c(0, 1), c(0, 1), ann = F, bty = 'n', type = 'n', xaxt = 'n', yaxt = 'n')
                        text(x = 0.4, y = 0.25, labels_ii[i],
                             cex = 2.5, col = "black")
                }
                if (i < j) {
                        plot.new()
                }
                if (i > j) {
                        xaxt_n = 'n'
                        yaxt_n = 'n'
                        if (i == maxplot) {xaxt_n = 's'}
                        if (j == 1) {yaxt_n = 's'}
                        x <- data.frame(X = hmc_sample_select[,j], Y = hmc_sample_select[,i])
                        est <- bkde2D(x, bandwidth=sapply(x,dpik)*1.4,
                                      range.x=list(c(min(x[,1]),max(x[,1])),
                                                   c(min(x[,2]),max(x[,2]))))

                        contour(est$x1, est$x2, est$fhat,
                                drawlabels = F, nlevels = 7, lty = 2, lwd = 1.5,
                                yaxt=yaxt_n, xaxt=xaxt_n,
                                col = my.cols)
                        points(true_val[j], true_val[i], pch = 8, col = "black", cex = 2)

                        x <- data.frame(X = vi_sample_select[,j], Y = vi_sample_select[,i])
                        est <- bkde2D(x, bandwidth=sapply(x,dpik)*1.4,
                                      range.x=list(c(min(x[,1]),max(x[,1])),
                                                   c(min(x[,2]),max(x[,2]))))
                        contour(est$x1, est$x2, est$fhat,
                                drawlabels = F, nlevels = 7, col=my.cols.vi,
                                add=TRUE, lwd = 1.5)



                }
        }
}
dev.off()


###############################################################################
# pdf(file='img/VIvsHMCbifcop.pdf', width = 18, height = 6)
# par(mfrow =c(2,6))
# par(mar=c(5,5,3,1))
#
# plot(get_theta(hmc_gauss), get_theta(vi_gauss) ,
#      xlab = expression(theta[gen]), ylab = expression(theta[vi]),
#      main = " Gaussian bi-factor copula")
# abline(a= 0, b=1, col="red")
#
# plot(get_theta(hmc_student), get_theta(vi_student),
#      xlab = expression(theta[gen]), ylab = expression(theta[vi]),
#      main = " Student bi-factor copula")
# abline(a= 0, b=1, col="red")
#
# plot(get_theta(hmc_clayton), get_theta(vi_clayton),
#      xlab = expression(theta[gen]), ylab = expression(theta[vi]),
#      main = " Clayton bi-factor copula")
# abline(a= 0, b=1, col="red")
#
# plot(get_theta(hmc_gumbel), get_theta(vi_gumbel),
#      xlab = expression(theta[gen]), ylab = expression(theta[vi]),
#      main = " Gumbel bi-factor copula")
# abline(a= 0, b=1, col="red")
#
# plot(get_theta(hmc_frank), get_theta(vi_frank),
#      xlab = expression(theta[gen]), ylab = expression(theta[vi]),
#      main = " Frank bi-factor copula")
# abline(a= 0, b=1, col="red")
#
# plot(get_theta(hmc_joe), get_theta(vi_joe),
#      xlab = expression(theta[gen]), ylab = expression(theta[vi]),
#      main = " Joe bi-factor copula")
# abline(a= 0, b=1, col="red")
#
#
# plot(get_theta_sd(hmc_gauss), get_theta_sd(vi_gauss),
#      xlab = expression(sd(theta[hmc]) ), ylab = expression(sd (theta[vi])),
#      main = " Gaussian bi-factor copula")
# abline(a= 0, b=1, col="red")
#
# plot(get_theta_sd(hmc_student), get_theta_sd(vi_student),
#      xlab = expression(sd(theta[hmc]) ), ylab = expression(sd (theta[vi])),
#      main = " Student bi-factor copula")
# abline(a= 0, b=1, col="red")
#
# plot(get_theta_sd(hmc_clayton), get_theta_sd(vi_clayton),
#      xlab = expression(sd(theta[hmc]) ), ylab = expression(sd (theta[vi])),
#      main = " Clayton bi-factor copula")
# abline(a= 0, b=1, col="red")
#
# plot(get_theta_sd(hmc_gumbel), get_theta_sd(vi_gumbel),
#      xlab = expression(sd(theta[hmc]) ), ylab = expression(sd (theta[vi])),
#      main = " Gumbel bi-factor copula")
# abline(a= 0, b=1, col="red")
#
# plot(get_theta_sd(hmc_frank), get_theta_sd(vi_frank),
#      xlab = expression(sd(theta[hmc]) ), ylab = expression(sd (theta[vi])),
#      main = " Frank bi-factor copula")
# abline(a= 0, b=1, col="red")
#
# plot(get_theta_sd(hmc_joe), get_theta_sd(vi_joe),
#      xlab = expression(sd(theta[hmc]) ), ylab = expression(sd (theta[vi])),
#      main = " Joe bi-factor copula")
# abline(a= 0, b=1, col="red")
#
# dev.off()
#
