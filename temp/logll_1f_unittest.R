# simnestfact(nn,grsize,cop,param)
# simbifact(nn,grsize,cop,param)
# f90str1nllk(param,dstruct,iprfn=F) # nested-factor
# f90str2nllk(param,dstruct,iprfn=F) # bi-factor

gl=gausslegendre(25)
grsize=c(6,6,6,6,6)
d=sum(grsize)
n=500
# nested-factor copula
mgrp=length(grsize)
set.seed(123)
parne=c(rep(4,6),rep(6,6),rep(6.5,6),rep(6.5,6),rep(6.5,6),rep(7,5))
udatne=simnestfact(n,grsize,cop=5,parne)
dstrfrk=list(data=udatne,copname="frank",quad=gl,repar=0,grsize=grsize)
npar=mgrp+d
nllk=f90str1nllk(parne,dstrfrk)
outn= pdhessminb(rep(3,npar),f90str1nllk, ifixed=rep(FALSE,npar), dstrfrk,
                 LB=rep(0,npar), UB=rep(30,npar), mxiter=30, eps=5.e-5,iprint=TRUE)





data <- list(u = udatne,
             n_max = d,
             t_max = n,
             k_max = length(grsize)+1,
             gid = c(rep(1, grsize[1]),rep(2, grsize[2]),rep(3, grsize[3]),
                     rep(4, grsize[4]),rep(5, grsize[5]) ),
             structfactor = 3)
init <- list(copula_type = rep(3,d),
             latent_copula_type = rep(3,length(grsize)))
other <- list(seed = 126, core = 8, iter = 1000,
              n_monte_carlo_grad = 1, n_monte_carlo_elbo = 100,
              eval_elbo = 100, adapt_bool = F, adapt_val = 1,
              adapt_iterations = 50, tol_rel_obj = 0.01, copselect = F, modelselect = T)

vi_frank <- vifcopula::vifcop(data,init,other)




##############################################################
gl=gausslegendre(25)
grsize=c(4,4,3)
d=sum(grsize)
n=500
# nested-factor copula
mgrp=length(grsize)
set.seed(123)
parne=c(rep(1.5,3),rep(2.6,4),rep(3.7,4),rep(2.7,3))
udatne=simnestfact(n,grsize,cop=3,parne)
dstrfrk=list(data=udatne,copname="gumbel",quad=gl,repar=0,grsize=grsize)
npar=mgrp+d
outn= pdhessminb(rep(3,npar),f90str1nllk, ifixed=rep(FALSE,npar), dstrfrk,
                 LB=rep(1,npar), UB=rep(30,npar), mxiter=30, eps=5.e-5,iprint=TRUE)




data <- list(u = udatne,
             n_max = d,
             t_max = n,
             k_max = length(grsize)+1,
             gid = c(rep(1, grsize[1]),rep(2, grsize[2]),rep(3, grsize[3]) ),
             structfactor = 3)
init <- list(copula_type = rep(4,d),
             latent_copula_type = rep(4,length(grsize)))
other <- list(seed = 126, core = 8, iter = 1000,
              n_monte_carlo_grad = 1, n_monte_carlo_elbo = 100,
              eval_elbo = 100, adapt_bool = F, adapt_val = 1,
              adapt_iterations = 50, tol_rel_obj = 0.01, copselect = F, modelselect = F)

vi_frank <- vifcopula::vifcop(data,init,other)


##############################################################
gl=gausslegendre(25)
grsize=c(4,4,3)
d=sum(grsize)
n=1000

set.seed(123)
parbi=c(rep(4,11),rep(6,4),rep(6.5,4),rep(7,3))
udatbi=simbifact(n,grsize,cop=5,parbi)
npar=2*d
dstrfrk=list(data=udatbi,copname="frank",quad=gl,repar=0,grsize=grsize,pdf=0)
nllk=f90str2nllk(parbi,dstrfrk)
outb=pdhessminb(c(rep(2,d),rep(3,d)),f90str2nllk,ifixed=rep(FALSE,npar),dstrfrk,
                LB=rep(0,npar), UB=rep(20,npar), mxiter=30, eps=5.e-5,iprint=TRUE)

ml1fact()

data <- list(u = udatbi,
             n_max = d,
             t_max = n,
             k_max = length(grsize)+1,
             gid = c(rep(1, grsize[1]),rep(2, grsize[2]),rep(3, grsize[3]) ),
             structfactor = 2)
init <- list(copula_type = rep(5,d),
             latent_copula_type = rep(5,d))
other <- list(seed = 126, core = 8, iter = 1000,
              n_monte_carlo_grad = 1, n_monte_carlo_elbo = 100,
              eval_elbo = 100, adapt_bool = F, adapt_val = 1,
              adapt_iterations = 50, tol_rel_obj = 0.01, copselect = F, modelselect = T)

vi_frank <- vifcopula::vifcop(data,init,other)



###################################################

gl=gausslegendre(25)
grsize=c(4,4,3)
d=sum(grsize)
n=1000

set.seed(123)
parbi=c(rep(4,11),rep(6,4),rep(6.5,4),rep(7,3))
udatbi=simbifact(n,grsize,cop=3,parbi)
npar=2*d
dstrfrk=list(data=udatbi,copname="gumbel",quad=gl,repar=0,grsize=grsize,pdf=0)
nllk=f90str2nllk(parbi,dstrfrk)


data <- list(u = udatbi,
             n_max = d,
             t_max = n,
             k_max = length(grsize)+1,
             gid = c(rep(1, grsize[1]),rep(2, grsize[2]),rep(3, grsize[3]) ),
             structfactor = 2)
init <- list(copula_type = rep(4,d),
             latent_copula_type = rep(4,d))
other <- list(seed = 126, core = 8, iter = 1000,
              n_monte_carlo_grad = 1, n_monte_carlo_elbo = 100,
              eval_elbo = 100, adapt_bool = F, adapt_val = 1,
              adapt_iterations = 50, tol_rel_obj = 0.01, copselect = F, modelselect = T)

vi_gumbel <- vifcopula::vifcop(data,init,other)
