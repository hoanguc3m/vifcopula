#' @export
get_v0 <- function(vi) {
    head(vi$mean_iv,vi$t_max)
}

#' @export
get_v <- function(vi) {
    t_max <- vi$t_max
    n_max <- vi$n_max
    k_max <- vi$k_max
    v <- matrix(head(vi$mean_iv,t_max*k_max), nrow = t_max)
    v[,2:k_max]
}

#' @export
get_theta <- function(vi) {
    t_max <- vi$t_max
    n_max <- vi$n_max
    k_max <- vi$k_max
    structfactor <- vi$structfactor

    all_theta <- vi$mean_iv[(t_max * k_max +1):length(vi$mean_iv)]
    theta <- rep(0,n_max)
    count <- 0
    if (structfactor == 1){
        for (i in 1:n_max){
            if (vi$cop_type[i] > 0){
                count = count + 1
                theta[i] <- all_theta[count]
                if (vi$cop_type[i] == 2){
                    count = count + 1
                }
            }
        }
    }

    if (structfactor == 2){
        for (i in 1:n_max){
            if (vi$cop_type[i] > 0){
                count = count + 1
                theta[i] <- all_theta[count]
                if (vi$cop_type[i] == 2){
                    count = count + 1
                }
            }
        }
    }

    if (structfactor == 3){
        count = k_max - 1 + sum(vi$latent_copula_type == 2) - sum(vi$latent_copula_type == 0)
        for (i in 1:n_max){
            if (vi$cop_type[i] > 0){
                count = count + 1
                theta[i] <- all_theta[count]
                if (vi$cop_type[i] == 2){
                    count = count + 1
                }
            }
        }
    }

    theta
}
#' @export
get_theta2 <- function(vi) {
    t_max <- vi$t_max
    n_max <- vi$n_max
    k_max <- vi$k_max
    structfactor <- vi$structfactor

    all_theta <- vi$mean_iv[(t_max * k_max +1):length(vi$mean_iv)]
    theta2 <- rep(0,n_max)
    count <- 0
    if (structfactor == 1){
        for (i in 1:n_max){
            if (vi$cop_type[i] > 0){
                count = count + 1
                if (vi$cop_type[i] == 2){
                    count = count + 1
                    theta2[i] <- all_theta[count]
                }
            }
        }
    }

    if (structfactor == 2){
        for (i in 1:n_max){
            if (vi$cop_type[i] > 0){
                count = count + 1
                if (vi$cop_type[i] == 2){
                    count = count + 1
	            theta2[i] <- all_theta[count]
                }
            }
        }
    }

    if (structfactor == 3){
        count = k_max - 1 + sum(vi$latent_copula_type == 2) - sum(vi$latent_copula_type == 0)
        for (i in 1:n_max){
            if (vi$cop_type[i] > 0){
                count = count + 1

                if (vi$cop_type[i] == 2){
                    count = count + 1
                    theta2[i] <- all_theta[count]
                }
            }
        }
    }
    theta2
}
#' @export
get_latent_theta <- function(vi) {
    t_max <- vi$t_max
    n_max <- vi$n_max
    k_max <- vi$k_max
    structfactor <- vi$structfactor


    all_theta <- vi$mean_iv[(t_max * k_max +1):length(vi$mean_iv)]
    latent_theta <- NULL
    count <- 0
    if (structfactor == 1){
	latent_theta = NULL
	}

    if (structfactor == 2){
    	latent_theta <- rep(0,n_max)
        count = n_max + sum(vi$cop_type == 2) - sum(vi$cop_type == 0)
        for (i in 1:n_max){
            if (vi$latent_copula_type[i] > 0){
                count = count + 1
	        latent_theta[i] <- all_theta[count]

                if (vi$latent_copula_type[i] == 2){
                    count = count + 1
                }
            }
        }
    }

    if (structfactor == 3){
    	latent_theta <- rep(0,k_max-1)
        for (i in 1:(k_max-1)){
            if (vi$latent_copula_type[i] > 0){
                count = count + 1
                latent_theta[i] <- all_theta[count]
                if (vi$latent_copula_type[i] == 2){
                    count = count + 1
                }
            }
        }
    }
    latent_theta
}
#' @export
get_latent_theta2 <- function(vi) {
    t_max <- vi$t_max
    n_max <- vi$n_max
    k_max <- vi$k_max
    structfactor <- vi$structfactor


    all_theta <- vi$mean_iv[(t_max * k_max +1):length(vi$mean_iv)]
    latent_theta2 <- NULL

    count <- 0

    if (structfactor == 2){	
    	latent_theta2 <- rep(0,n_max)
        count = n_max + sum(vi$cop_type == 2) - sum(vi$cop_type == 0)
        for (i in 1:n_max){
            if (vi$latent_copula_type[i] > 0){
                count = count + 1
                if (vi$latent_copula_type[i] == 2){
                    count = count + 1
	            latent_theta2[i] <- all_theta[count]
                }
            }
        }
    }

    if (structfactor == 3){
    	latent_theta2 <- rep(0,k_max-1)
        count = n_max + sum(vi$latent_copula_type == 2) - sum(vi$latent_copula_type == 0)
        for (i in 1:(k_max-1)){
            if (vi$latent_copula_type[i] > 0){
                count = count + 1

                if (vi$latent_copula_type[i] == 2){
                    count = count + 1
                    latent_theta2[i] <- all_theta[count]
                }
            }
        }
    }
    latent_theta2
}
#' @export
num_param <- function(vi) {
    t_max <- vi$t_max
    n_max <- vi$n_max
    k_max <- vi$k_max
    structfactor <- vi$structfactor
    count <- 0

    if (structfactor == 1){
	count = t_max + n_max + sum(vi$cop_type == 2) - sum(vi$cop_type == 0)
	}

    if (structfactor == 2){	
	count = t_max*k_max + n_max + sum(vi$latent_copula_type == 2) - sum(vi$latent_copula_type == 0) +
		 		n_max + sum(vi$cop_type == 2) - sum(vi$cop_type == 0)
    }

    if (structfactor == 3){	
	count = t_max*k_max + (k_max - 1) + sum(vi$latent_copula_type == 2) - sum(vi$latent_copula_type == 0) +
		 		n_max + sum(vi$cop_type == 2) - sum(vi$cop_type == 0)
    }
    count
}

# cbind(vi$cop_type, get_theta(vi), get_theta2(vi), BiCopPar2Tau(vi$cop_type, get_theta(vi), get_theta2(vi)) )
# cbind(vi$latent_copula_type, get_latent_theta(vi), get_latent_theta2(vi), BiCopPar2Tau(vi$latent_copula_type, get_latent_theta(vi), get_latent_theta2(vi)))

# cbind(vi$latent_copula_type, get_latent_theta(vi), get_latent_theta2(vi), BiCopPar2Tau(vi$latent_copula_type, get_latent_theta(vi), get_latent_theta2(vi)) ,
#     vi_gauss$latent_copula_type, get_latent_theta(vi_gauss), get_latent_theta2(vi_gauss), BiCopPar2Tau(vi_gauss$latent_copula_type, get_latent_theta(vi_gauss), get_latent_theta2(vi_gauss)))
#
# cbind(vi$cop_type, get_theta(vi), get_theta2(vi), BiCopPar2Tau(vi$cop_type, get_theta(vi), get_theta2(vi)) ,
#     vi_gauss$cop_type, get_theta(vi_gauss), get_theta2(vi_gauss),
#     BiCopPar2Tau(vi_gauss$cop_type, get_theta(vi_gauss), get_theta2(vi_gauss)))
