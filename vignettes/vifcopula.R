## ---- eval=FALSE---------------------------------------------------------
#  library(vifcopula)
#  t_max = 1000
#  n_max = 100
#  datagen_gauss <- fcopsim(t_max = t_max, n_max = n_max, family = 1)
#  datagen <- datagen_gauss
#  data <- list(u = datagen$u,
#               n_max = datagen$n_max,
#               t_max = datagen$t_max,
#               k_max = datagen$k_max,
#               gid = datagen$gid,
#               structfactor = datagen$structfactor)
#  init <- list(copula_type = datagen$family,
#               v = datagen$v,
#               par = datagen$theta,
#               par2 = datagen$theta2)
#  other <- list(seed = 126, core = 8, iter = 1000,
#                n_monte_carlo_grad = 1, n_monte_carlo_elbo = 10,
#                eval_elbo = 50, adapt_bool = T, adapt_val = 1,
#                adapt_iterations = 50, tol_rel_obj = 0.1, copselect = T)
#  vi_gauss <- vifcopula::vifcop(data,init,other)
#  
#  plot(datagen$v, vi_gauss$mean_vi[1:t_max])
#  abline(a= 0, b=1, col="red")
#  plot(datagen$theta, vi_gauss$mean_vi[(t_max+1):(t_max+n_max)])
#  abline(a= 0, b=1, col="red")
#  vi_gauss$cop_vec_new
#  rho_vi = vi_gauss$mean_vi[(t_max+1):(t_max+n_max)]
#  selec <- as.logical(vi_gauss$cop_vec_new)
#  pos = seq(1,n_max)
#  rho_vi[pos[!selec]]
#  
#  init <- list(copula_type = matrix(vi_gauss$cop_vec_new, ncol=1),
#               v = datagen$v,
#               par = datagen$theta,
#               par2 = datagen$theta2)
#  vi_gauss <- vifcopula::vifcop(data,init,other)

## ---- eval=FALSE---------------------------------------------------------
#  datagen_student <- fcopsim(t_max = 1000, n_max = 100, family = 2)
#  datagen <- datagen_student
#  data <- list(u = datagen$u,
#               n_max = datagen$n_max,
#               t_max = datagen$t_max,
#               k_max = datagen$k_max,
#               gid = datagen$gid,
#               structfactor = datagen$structfactor)
#  init <- list(copula_type = datagen$family,
#               v = datagen$v,
#               par = datagen$theta,
#               par2 = datagen$theta2)
#  other <- list(seed = 126, core = 8, iter = 1000,
#                n_monte_carlo_grad = 1, n_monte_carlo_elbo = 10,
#                eval_elbo = 100, adapt_bool = T, adapt_val = 1,
#                adapt_iterations = 50, tol_rel_obj = 0.1, copselect = F)
#  vi_student <- vifcopula::vifcop(data,init,other)
#  
#  plot(datagen$v, vi_student$mean_vi[1:t_max])
#  abline(a= 0, b=1, col="red")
#  plot(datagen$theta, vi_student$mean_vi[seq(from = t_max+1, to = t_max+2*n_max, by = 2 )])
#  plot(datagen$theta2, vi_student$mean_vi[seq(from = t_max+2, to = t_max+n_max, by = 2 )])
#  abline(a= 0, b=1, col="red")

## ---- eval=FALSE---------------------------------------------------------
#  datagen_clayton <- fcopsim(t_max = 1000, n_max = 100, family = 3)
#  datagen <- datagen_clayton
#  data <- list(u = datagen$u,
#               n_max = datagen$n_max,
#               t_max = datagen$t_max,
#               k_max = datagen$k_max,
#               gid = datagen$gid,
#               structfactor = datagen$structfactor)
#  init <- list(copula_type = datagen$family,
#               v = datagen$v,
#               par = datagen$theta,
#               par2 = datagen$theta2)
#  other <- list(seed = 126, core = 8, iter = 1000,
#                n_monte_carlo_grad = 1, n_monte_carlo_elbo = 10,
#                eval_elbo = 100, adapt_bool = T, adapt_val = 1,
#                adapt_iterations = 50, tol_rel_obj = 0.1, copselect = T)
#  vi_clayton <- vifcopula::vifcop(data,init,other)
#  
#  plot(datagen$v, vi_clayton$mean_vi[1:t_max])
#  abline(a= 0, b=1, col="red")
#  plot(datagen$theta, vi_clayton$mean_vi[(t_max+1):(t_max+n_max)])
#  abline(a= 0, b=1, col="red")
#  vi_clayton$cop_vec_new

## ---- eval=FALSE---------------------------------------------------------
#  datagen_gumbel <- fcopsim(t_max = 1000, n_max = 100, family = 4)
#  datagen <- datagen_gumbel
#  data <- list(u = datagen$u,
#               n_max = datagen$n_max,
#               t_max = datagen$t_max,
#               k_max = datagen$k_max,
#               gid = datagen$gid,
#               structfactor = datagen$structfactor)
#  init <- list(copula_type = datagen$family,
#               v = datagen$v,
#               par = datagen$theta,
#               par2 = datagen$theta2)
#  other <- list(seed = 126, core = 8, iter = 1000,
#                n_monte_carlo_grad = 1, n_monte_carlo_elbo = 10,
#                eval_elbo = 100, adapt_bool = T, adapt_val = 1,
#                adapt_iterations = 50, tol_rel_obj = 0.1, copselect = T)
#  vi_gumbel <- vifcopula::vifcop(data,init,other)
#  
#  plot(datagen$v, vi_gumbel$mean_vi[1:t_max])
#  abline(a= 0, b=1, col="red")
#  plot(datagen$theta, vi_gumbel$mean_vi[(t_max+1):(t_max+n_max)])
#  abline(a= 0, b=1, col="red")
#  vi_gumbel$cop_vec_new
#  

## ---- eval=FALSE---------------------------------------------------------
#  datagen_frank <- fcopsim(t_max = 1000, n_max = 100, family = 5)
#  datagen <- datagen_frank
#  data <- list(u = datagen$u,
#               n_max = datagen$n_max,
#               t_max = datagen$t_max,
#               k_max = datagen$k_max,
#               gid = datagen$gid,
#               structfactor = datagen$structfactor)
#  init <- list(copula_type = datagen$family,
#               v = datagen$v,
#               par = datagen$theta,
#               par2 = datagen$theta2)
#  other <- list(seed = 126, core = 8, iter = 1000,
#                n_monte_carlo_grad = 1, n_monte_carlo_elbo = 10,
#                eval_elbo = 100, adapt_bool = T, adapt_val = 1,
#                adapt_iterations = 50, tol_rel_obj = 0.1, copselect = T)
#  vi_frank <- vifcopula::vifcop(data,init,other)
#  
#  plot(datagen$v, vi_frank$mean_vi[1:t_max])
#  abline(a= 0, b=1, col="red")
#  plot(datagen$theta, vi_frank$mean_vi[(t_max+1):(t_max+n_max)])
#  abline(a= 0, b=1, col="red")
#  vi_frank$cop_vec_new
#  

## ---- eval=FALSE---------------------------------------------------------
#  datagen_joe <- fcopsim(t_max = 1000, n_max = 100, family = 6)
#  datagen <- datagen_joe
#  data <- list(u = datagen$u,
#               n_max = datagen$n_max,
#               t_max = datagen$t_max,
#               k_max = datagen$k_max,
#               gid = datagen$gid,
#               structfactor = datagen$structfactor)
#  init <- list(copula_type = datagen$family,
#               v = datagen$v,
#               par = datagen$theta,
#               par2 = datagen$theta2)
#  other <- list(seed = 126, core = 8, iter = 1000,
#                n_monte_carlo_grad = 1, n_monte_carlo_elbo = 10,
#                eval_elbo = 100, adapt_bool = T, adapt_val = 1,
#                adapt_iterations = 50, tol_rel_obj = 0.1, copselect = T)
#  vi_joe <- vifcopula::vifcop(data,init,other)
#  
#  plot(datagen$v, vi_joe$mean_vi[1:t_max])
#  abline(a= 0, b=1, col="red")
#  plot(datagen$theta, vi_joe$mean_vi[(t_max+1):(t_max+n_max)])
#  abline(a= 0, b=1, col="red")
#  vi_joe$cop_vec_new
#  sum(vi_joe$cop_vec_new == datagen_joe$family)
#  

## ---- eval=FALSE---------------------------------------------------------
#  copfamily = matrix(sample(c(1,3,4,5,6),size = 100, replace = T),ncol=1)
#  datagen_mix <- fcopsim(t_max = 1000, n_max = 100, family = copfamily)
#  datagen <- datagen_mix
#  data <- list(u = datagen$u,
#               n_max = datagen$n_max,
#               t_max = datagen$t_max,
#               k_max = datagen$k_max,
#               gid = datagen$gid,
#               structfactor = datagen$structfactor)
#  init <- list(copula_type = datagen$family,
#               v = datagen$v,
#               par = datagen$theta,
#               par2 = datagen$theta2)
#  other <- list(seed = 126, core = 8, iter = 1000,
#                n_monte_carlo_grad = 1, n_monte_carlo_elbo = 10,
#                eval_elbo = 100, adapt_bool = T, adapt_val = 1,
#                adapt_iterations = 50, tol_rel_obj = 0.1, copselect = T)
#  vi_mix <- vifcopula::vifcop(data,init,other)
#  
#  plot(datagen$v, vi_mix$mean_vi[1:t_max])
#  abline(a= 0, b=1, col="red")
#  plot(datagen$theta, vi_mix$mean_vi[(t_max+1):(t_max+n_max)])
#  abline(a= 0, b=1, col="red")
#  sum(vi_mix$cop_vec_new == datagen$family)
#  
#  init <- list(copula_type = matrix(1,nrow = n_max,ncol=1),
#               v = datagen$v,
#               par = datagen$theta,
#               par2 = datagen$theta2)
#  vi_mix <- vifcopula::vifcop(data,init,other)
#  sum(vi_mix$cop_vec_new == datagen_mix$family)
#  

