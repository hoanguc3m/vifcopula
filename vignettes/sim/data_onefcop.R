library(devtools)
install_github("hoanguc3m/vifcopula")
#setwd("/home/hoanguc3m/Dropbox/WP2/")
library(vifcopula)
t_max = 1000
n_max = 100
gauss_init <- matrix(1, nrow = n_max, ncol = 1)

load("/home/hoanguc3m/Dropbox/WP2/code/vifcopula/vignettes/sim/data_fin.RData")

copfamily = matrix(sample(c(1,2,3,4,5,6),size = t_max, replace = T),ncol=1)
n_max = dim(fin_u_mat)[2]
t_max = dim(fin_u_mat)[1]
gid = rep(1, n_max)
k_max = 1
structfactor = 1
data <- list(u = fin_u_mat,
             n_max = n_max,
             t_max = t_max,
             k_max = k_max,
             gid = gid,
             structfactor = structfactor)
init <- list(copula_type = copfamily,
             v = matrix(0.5,nrow = t_max, ncol = 1),
             par = rep(0.5,t_max),
             par2 = rep(0.5,t_max))
other <- list(seed = 126, core = 8, iter = 1000,
              n_monte_carlo_grad = 1, n_monte_carlo_elbo = 10,
              eval_elbo = 100, adapt_bool = F, adapt_val = 1,
              adapt_iterations = 50, tol_rel_obj = 0.1, copselect = T)
vi_mix <- vifcopula::vifcop(data,init,other)

sum(vi_mix$cop_type == 2)

pdf(file='img/Mix1.pdf', width = 9, height = 4.5)
par(mfrow =c(1,2))
par(mar=c(5,5,3,1))
plot(datagen$v, vi_mix$mean_iv[1:t_max], xlab = expression(v[t]), ylab = expression(v[approximated]))
abline(a= 0, b=1, col="red")
plot(datagen$theta, vi_mix$mean_iv[(t_max+1):(t_max+n_max)], xlab = expression(theta[t]), ylab = expression(theta[approximated]))
abline(a= 0, b=1, col="red")
dev.off()
