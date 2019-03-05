#' @export
get_v0 <- function(x) UseMethod("get_v0",x)
#' @export
get_v  <- function(x) UseMethod("get_v",x)
#' @export
get_theta <- function(x) UseMethod("get_theta",x)
#' @export
get_theta2  <- function(x) UseMethod("get_theta2",x)
#' @export
get_latent_theta  <- function(x) UseMethod("get_latent_theta",x)
#' @export
get_latent_theta2  <- function(x) UseMethod("get_latent_theta2",x)
#' @export
get_vine_theta  <- function(x) UseMethod("get_vine_theta",x)
#' @export
get_vine_theta2  <- function(x) UseMethod("get_vine_theta2",x)
#' @export
get_v0_sd <- function(x) UseMethod("get_v0_sd",x)
#' @export
get_v_sd  <- function(x) UseMethod("get_v_sd",x)
#' @export
get_theta_sd <- function(x) UseMethod("get_theta_sd",x)
#' @export
get_theta2_sd  <- function(x) UseMethod("get_theta2_sd",x)
#' @export
get_latent_theta_sd  <- function(x) UseMethod("get_latent_theta_sd",x)
#' @export
get_latent_theta2_sd  <- function(x) UseMethod("get_latent_theta2_sd",x)
#' @export
get_vine_theta_sd  <- function(x) UseMethod("get_vine_theta_sd",x)
#' @export
get_vine_theta2_sd  <- function(x) UseMethod("get_vine_theta2_sd",x)
#' @export
num_param  <- function(x) UseMethod("num_param",x)


###############################################################################

#' @export
get_v0.vifcop <- function(vi) {
    head(vi$mean_vi,vi$t_max)
}

#' @export
get_v.vifcop <- function(vi) {
    v_out <- NULL
    if (vi$structfactor == 1 | vi$structfactor == 11 | vi$structfactor == 4 | vi$structfactor == 14){
        v_out <- head(vi$mean_vi,vi$t_max)
    } else {
        t_max <- vi$t_max
        n_max <- vi$n_max
        k_max <- vi$k_max
        v <- matrix(head(vi$mean_vi,t_max*k_max), nrow = t_max)
        v_out <- v[,2:k_max]
    }
    v_out
}

#' @export
get_theta.vifcop <- function(vi) {
    t_max <- vi$t_max
    n_max <- vi$n_max
    k_max <- vi$k_max
    structfactor <- vi$structfactor

    all_theta <- vi$mean_vi[(t_max * k_max +1):length(vi$mean_vi)]
    theta <- rep(0,n_max)
    count <- 0
    if (structfactor == 1 | structfactor == 4 | structfactor == 2){
        count <- 0
    } else {
        if (structfactor == 3){
            count = k_max - 1 + sum(vi$latent_copula_type == 2) - sum(vi$latent_copula_type == 0)
        }
    }

    for (i in 1:n_max){
        if (vi$cop_type[i] > 0){
            count = count + 1
            theta[i] <- all_theta[count]
            if (vi$cop_type[i] == 2){
                count = count + 1
            }
        }
    }
    theta
}
#' @export
get_theta2.vifcop <- function(vi) {
    t_max <- vi$t_max
    n_max <- vi$n_max
    k_max <- vi$k_max
    structfactor <- vi$structfactor

    all_theta <- vi$mean_vi[(t_max * k_max +1):length(vi$mean_vi)]
    theta2 <- rep(0,n_max)
    count <- 0
    if (structfactor == 1 | structfactor == 4 | structfactor == 2){
        count <- 0
    } else {
        if (structfactor == 3){
            count = k_max - 1 + sum(vi$latent_copula_type == 2) - sum(vi$latent_copula_type == 0)
        }
    }

    for (i in 1:n_max){
        if (vi$cop_type[i] > 0){
            count = count + 1
            if (vi$cop_type[i] == 2){
                count = count + 1
                theta2[i] <- all_theta[count]
            }
        }
    }
    theta2
}
#' @export
get_latent_theta.vifcop <- function(vi) {
    t_max <- vi$t_max
    n_max <- vi$n_max
    k_max <- vi$k_max
    structfactor <- vi$structfactor


    all_theta <- vi$mean_vi[(t_max * k_max +1):length(vi$mean_vi)]
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
get_latent_theta2.vifcop <- function(vi) {
    t_max <- vi$t_max
    n_max <- vi$n_max
    k_max <- vi$k_max
    structfactor <- vi$structfactor


    all_theta <- vi$mean_vi[(t_max * k_max +1):length(vi$mean_vi)]
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
        count = 0
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
get_vine_theta.vifcop <- function(vi) {
    t_max <- vi$t_max
    n_max <- vi$n_max
    k_max <- vi$k_max
    g_max <- max(vi$gid)
    structfactor <- vi$structfactor


    all_theta <- vi$mean_vi[(t_max * k_max +1):length(vi$mean_vi)]
    vine_theta <- NULL
    count <- 0
    if (structfactor != 4){
        vine_theta = NULL
    }

    vine_theta <- rep(0, n_max - g_max)
    count = n_max + sum(vi$cop_type == 2) - sum(vi$cop_type == 0)
    for (i in 1:length(vi$vine_copula_type)){
        if (vi$vine_copula_type[i] > 0){
            count = count + 1
            vine_theta[i] <- all_theta[count]

            if (vi$vine_copula_type[i] == 2){
                count = count + 1
            }
        }
    }

    vine_theta
}
#' @export
get_vine_theta2.vifcop <- function(vi) {
    t_max <- vi$t_max
    n_max <- vi$n_max
    k_max <- vi$k_max
    g_max <- max(vi$gid)
    structfactor <- vi$structfactor


    all_theta <- vi$mean_vi[(t_max * k_max +1):length(vi$mean_vi)]
    vine_theta2 <- NULL

    count <- 0
    if (structfactor != 4){
        vine_theta2 <- NULL
    }

    if (structfactor == 4){
        vine_theta2 <- rep(0,n_max - g_max)
        count = n_max + sum(vi$cop_type == 2) - sum(vi$cop_type == 0)
        for (i in 1:length(vi$vine_copula_type)){
            if (vi$vine_copula_type[i] > 0){
                count = count + 1
                if (vi$vine_copula_type[i] == 2){
                    count = count + 1
                    vine_theta2[i] <- all_theta[count]
                }
            }
        }
    }
    vine_theta2
}
#' @export
num_param.vifcop <- function(vi) {
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
    if (structfactor == 4){
        g_max <- max(vi$gid)
        count = t_max*k_max + sum(vi$vine_copula_type > 0) + sum(vi$vine_copula_type == 2) - sum(vi$vine_copula_type == 0) +
            n_max + sum(vi$cop_type == 2) - sum(vi$cop_type == 0)
    }
    count
}

###############################################################################

#' @export
get_v0.hmcfcop <- function(hmc) {
    # No need
    # "lp__"          "accept_stat__" "stepsize__"    "treedepth__"   "n_leapfrog__"  "divergent__"   "energy__"
    head(hmc$mean_hmc,hmc$t_max + 7)[-c(1:7)]
}

#' @export
get_v.hmcfcop <- function(hmc) {
    v_out <- NULL
    if (hmc$structfactor == 1 | hmc$structfactor == 11 | hmc$structfactor == 4 | hmc$structfactor == 14){
        v_out <- get_v0(hmc)
    } else {
        t_max <- hmc$t_max
        n_max <- hmc$n_max
        k_max <- hmc$k_max
        v <- matrix(head(hmc$mean_hmc,t_max*k_max+7)[-c(1:7)], nrow = t_max)
        v_out <- v[,2:k_max]
    }
    v_out
}

#' @export
get_theta.hmcfcop <- function(hmc) {
    t_max <- hmc$t_max
    n_max <- hmc$n_max
    k_max <- hmc$k_max
    structfactor <- hmc$structfactor

    all_theta <- hmc$mean_hmc[(t_max * k_max +8):length(hmc$mean_hmc)] # Eliminate 7 at beginning
    theta <- rep(0,n_max)
    count <- 0
    if (structfactor == 1 | structfactor == 4 | structfactor == 2){
        count <- 0
    } else {
        if (structfactor == 3){
            count = k_max - 1 + sum(hmc$latent_copula_type == 2) - sum(hmc$latent_copula_type == 0)
        }
    }

    for (i in 1:n_max){
        if (hmc$cop_type[i] > 0){
            count = count + 1
            theta[i] <- all_theta[count]
            if (hmc$cop_type[i] == 2){
                count = count + 1
            }
        }
    }

    theta
}


#' @export
get_theta2.hmcfcop <- function(hmc) {
    t_max <- hmc$t_max
    n_max <- hmc$n_max
    k_max <- hmc$k_max
    structfactor <- hmc$structfactor

    all_theta <- hmc$mean_hmc[(t_max * k_max +8):length(hmc$mean_hmc)] # Eliminate 7 at beginning
    theta2 <- rep(0,n_max)
    count <- 0
    if (structfactor == 1 | structfactor == 4 | structfactor == 2){
        count <- 0
    } else {
        if (structfactor == 3){
            count = k_max - 1 + sum(hmc$latent_copula_type == 2) - sum(hmc$latent_copula_type == 0)
        }
    }

    for (i in 1:n_max){
        if (hmc$cop_type[i] > 0){
            count = count + 1
            if (hmc$cop_type[i] == 2){
                count = count + 1
                theta2[i] <- all_theta[count]
            }
        }
    }

    theta2
}
#' @export
get_latent_theta.hmcfcop <- function(hmc) {
    t_max <- hmc$t_max
    n_max <- hmc$n_max
    k_max <- hmc$k_max
    structfactor <- hmc$structfactor


    all_theta <- hmc$mean_hmc[(t_max * k_max +8):length(hmc$mean_hmc)] # Eliminate 7 at beginning
    latent_theta <- NULL
    count <- 0
    if (structfactor == 1){
        latent_theta = NULL
    }

    if (structfactor == 2){
        latent_theta <- rep(0,n_max)
        count = n_max + sum(hmc$cop_type == 2) - sum(hmc$cop_type == 0)
        for (i in 1:n_max){
            if (hmc$latent_copula_type[i] > 0){
                count = count + 1
                latent_theta[i] <- all_theta[count]

                if (hmc$latent_copula_type[i] == 2){
                    count = count + 1
                }
            }
        }
    }

    if (structfactor == 3){
        latent_theta <- rep(0,k_max-1)
        for (i in 1:(k_max-1)){
            if (hmc$latent_copula_type[i] > 0){
                count = count + 1
                latent_theta[i] <- all_theta[count]
                if (hmc$latent_copula_type[i] == 2){
                    count = count + 1
                }
            }
        }
    }
    latent_theta
}
#' @export
get_latent_theta2.hmcfcop <- function(hmc) {
    t_max <- hmc$t_max
    n_max <- hmc$n_max
    k_max <- hmc$k_max
    structfactor <- hmc$structfactor


    all_theta <- hmc$mean_hmc[(t_max * k_max +8):length(hmc$mean_hmc)] # Eliminate 7 at beginning
    latent_theta2 <- NULL

    count <- 0

    if (structfactor == 2){
        latent_theta2 <- rep(0,n_max)
        count = n_max + sum(hmc$cop_type == 2) - sum(hmc$cop_type == 0)
        for (i in 1:n_max){
            if (hmc$latent_copula_type[i] > 0){
                count = count + 1
                if (hmc$latent_copula_type[i] == 2){
                    count = count + 1
                    latent_theta2[i] <- all_theta[count]
                }
            }
        }
    }

    if (structfactor == 3){
        latent_theta2 <- rep(0,k_max-1)
        count = 0
        for (i in 1:(k_max-1)){
            if (hmc$latent_copula_type[i] > 0){
                count = count + 1

                if (hmc$latent_copula_type[i] == 2){
                    count = count + 1
                    latent_theta2[i] <- all_theta[count]
                }
            }
        }
    }
    latent_theta2
}

#' @export
get_vine_theta.hmcfcop <- function(hmc) {
    t_max <- hmc$t_max
    n_max <- hmc$n_max
    k_max <- hmc$k_max
    g_max <- max(hmc$gid)
    structfactor <- hmc$structfactor


    all_theta <- hmc$mean_hmc[(t_max * k_max +8):length(hmc$mean_hmc)] # Eliminate 7 at beginning
    vine_theta <- NULL
    count <- 0
    if (structfactor != 4){
        vine_theta = NULL
    }
    vine_theta <- rep(0,n_max - g_max)
    count = n_max + sum(hmc$cop_type == 2) - sum(hmc$cop_type == 0)

    for (i in 1:length(hmc$vine_copula_type)){
        if (hmc$vine_copula_type[i] > 0){
            count = count + 1
            vine_theta[i] <- all_theta[count]

            if (hmc$vine_copula_type[i] == 2){
                count = count + 1
            }
        }
    }
    vine_theta
}
#' @export
get_vine_theta2.hmcfcop <- function(hmc) {
    t_max <- hmc$t_max
    n_max <- hmc$n_max
    k_max <- hmc$k_max
    g_max <- max(hmc$gid)
    structfactor <- hmc$structfactor


    all_theta <- hmc$mean_hmc[(t_max * k_max +8):length(hmc$mean_hmc)] # Eliminate 7 at beginning
    vine_theta2 <- NULL

    count <- 0
    if (structfactor != 4){
        vine_theta2 <- NULL
    }

    if (structfactor == 4){
        vine_theta2 <- rep(0,n_max - g_max)
        count = n_max + sum(hmc$cop_type == 2) - sum(hmc$cop_type == 0)
        for (i in 1:length(hmc$vine_copula_type)){
            if (hmc$vine_copula_type[i] > 0){
                count = count + 1
                if (hmc$vine_copula_type[i] == 2){
                    count = count + 1
                    vine_theta2[i] <- all_theta[count]
                }
            }
        }
    }
    vine_theta2
}

#' @export
num_param.hmcfcop <- function(hmc) {
    t_max <- hmc$t_max
    n_max <- hmc$n_max
    k_max <- hmc$k_max
    structfactor <- hmc$structfactor
    count <- 0

    if (structfactor == 1){
        count = t_max + n_max + sum(hmc$cop_type == 2) - sum(hmc$cop_type == 0)
    }

    if (structfactor == 2){
        count = t_max*k_max + n_max + sum(hmc$latent_copula_type == 2) - sum(hmc$latent_copula_type == 0) +
            n_max + sum(hmc$cop_type == 2) - sum(hmc$cop_type == 0)
    }

    if (structfactor == 3){
        count = t_max*k_max + (k_max - 1) + sum(hmc$latent_copula_type == 2) - sum(hmc$latent_copula_type == 0) +
            n_max + sum(hmc$cop_type == 2) - sum(hmc$cop_type == 0)
    }
    if (structfactor == 4){
        g_max <- max(hmc$gid)
        count = t_max*k_max + sum(hmc$vine_copula_type > 0) + sum(hmc$vine_copula_type == 2) - sum(hmc$vine_copula_type == 0) +
            n_max + sum(hmc$cop_type == 2) - sum(hmc$cop_type == 0)
    }
    count
}

###############################################################################

#' @export
get_v0_sd.vifcop <- function(vi) {
    all_param_sd = apply(vi$sample_vi, MARGIN = 2, FUN = sd)
    head(all_param_sd,vi$t_max)
}

#' @export
get_v_sd.vifcop <- function(vi) {
    sd_out <- NULL
    all_param_sd = apply(vi$sample_vi, MARGIN = 2, FUN = sd)

    if (vi$structfactor > 1){
        t_max <- vi$t_max
        n_max <- vi$n_max
        k_max <- vi$k_max
        sd_mat <- matrix(head(all_param_sd,t_max*k_max), nrow = t_max)
        sd_out <- sd_mat[,2:k_max]
    }
    sd_out
}

#' @export
get_theta_sd.vifcop <- function(vi) {
    t_max <- vi$t_max
    n_max <- vi$n_max
    k_max <- vi$k_max
    structfactor <- vi$structfactor
    all_param_sd = apply(vi$sample_vi, MARGIN = 2, FUN = sd)

    all_theta <- all_param_sd[(t_max * k_max +1):length(vi$mean_vi)]
    theta <- rep(0,n_max)
    count <- 0

    if (structfactor == 1 | structfactor == 4 | structfactor == 2){
        count <- 0
    } else {
        if (structfactor == 3){
            count = k_max - 1 + sum(vi$latent_copula_type == 2) - sum(vi$latent_copula_type == 0)
        }
    }


    for (i in 1:n_max){
        if (vi$cop_type[i] > 0){
            count = count + 1
            theta[i] <- all_theta[count]
            if (vi$cop_type[i] == 2){
                count = count + 1
            }
        }
    }

    theta
}

#' @export
get_theta2_sd.vifcop <- function(vi) {
    t_max <- vi$t_max
    n_max <- vi$n_max
    k_max <- vi$k_max
    structfactor <- vi$structfactor
    all_param_sd = apply(vi$sample_vi, MARGIN = 2, FUN = sd)

    all_theta <- all_param_sd[(t_max * k_max +1):length(vi$mean_vi)]
    theta2 <- rep(0,n_max)
    count <- 0

    if (structfactor == 1 | structfactor == 4 | structfactor == 2){
        count <- 0
    } else {
        if (structfactor == 3){
            count = k_max - 1 + sum(vi$latent_copula_type == 2) - sum(vi$latent_copula_type == 0)
        }
    }

    for (i in 1:n_max){
        if (vi$cop_type[i] > 0){
            count = count + 1
            if (vi$cop_type[i] == 2){
                count = count + 1
                theta2[i] <- all_theta[count]
            }
        }
    }
    theta2
}

#' @export
get_latent_theta_sd.vifcop <- function(vi) {
    t_max <- vi$t_max
    n_max <- vi$n_max
    k_max <- vi$k_max
    structfactor <- vi$structfactor
    all_param_sd = apply(vi$sample_vi, MARGIN = 2, FUN = sd)

    all_theta <- all_param_sd[(t_max * k_max +1):length(vi$mean_vi)]
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
get_latent_theta2_sd.vifcop <- function(vi) {
    t_max <- vi$t_max
    n_max <- vi$n_max
    k_max <- vi$k_max
    structfactor <- vi$structfactor
    all_param_sd = apply(vi$sample_vi, MARGIN = 2, FUN = sd)

    all_theta <- all_param_sd[(t_max * k_max +1):length(vi$mean_vi)]
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
        count = 0
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
get_vine_theta_sd.vifcop <- function(vi) {
    t_max <- vi$t_max
    n_max <- vi$n_max
    k_max <- vi$k_max
    g_max <- max(vi$gid)
    structfactor <- vi$structfactor
    all_param_sd = apply(vi$sample_vi, MARGIN = 2, FUN = sd)

    all_theta <- all_param_sd[(t_max * k_max +1):length(vi$mean_vi)]
    vine_theta <- NULL
    count <- 0
    if (structfactor != 4){
        vine_theta = NULL
    }

    if (structfactor == 4){
        vine_theta <- rep(0,n_max - g_max)
        count = n_max + sum(vi$cop_type == 2) - sum(vi$cop_type == 0)
        for (i in 1:length(vi$vine_copula_type)){
            if (vi$vine_copula_type[i] > 0){
                count = count + 1
                vine_theta[i] <- all_theta[count]

                if (vi$vine_copula_type[i] == 2){
                    count = count + 1
                }
            }
        }
    }
    vine_theta
}
#' @export
get_vine_theta2_sd.vifcop <- function(vi) {
    t_max <- vi$t_max
    n_max <- vi$n_max
    k_max <- vi$k_max
    g_max <- max(vi$gid)
    structfactor <- vi$structfactor
    all_param_sd = apply(vi$sample_vi, MARGIN = 2, FUN = sd)

    all_theta <- all_param_sd[(t_max * k_max +1):length(vi$mean_vi)]
    vine_theta2 <- NULL

    count <- 0
    if (structfactor != 4){
        vine_theta2 = NULL
    }

    if (structfactor == 4){
        vine_theta2 <- rep(0,n_max - g_max)
        count = n_max + sum(vi$cop_type == 2) - sum(vi$cop_type == 0)
        for (i in 1:length(vi$vine_copula_type)){
            if (vi$vine_copula_type[i] > 0){
                count = count + 1
                if (vi$vine_copula_type[i] == 2){
                    count = count + 1
                    vine_theta2[i] <- all_theta[count]
                }
            }
        }
    }
    vine_theta2
}
###############################################################################


#' @export
get_v0_sd.hmcfcop <- function(hmc) {
    all_param_sd = apply(hmc$sample_hmc, MARGIN = 2, FUN = sd)[-c(1:7)]
    head(all_param_sd,hmc$t_max)
}

#' @export
get_v_sd.hmcfcop <- function(hmc) {
    sd_out <- NULL
    if (hmc$structfactor == 1 | hmc$structfactor == 11 | hmc$structfactor == 4 | hmc$structfactor == 14){
        sd_out <- get_v0_sd(hmc)
    } else {
        all_param_sd = apply(hmc$sample_hmc, MARGIN = 2, FUN = sd)[-c(1:7)]
        t_max <- hmc$t_max
        n_max <- hmc$n_max
        k_max <- hmc$k_max
        sd_mat <- matrix(head(all_param_sd,t_max*k_max), nrow = t_max)
        sd_out <- sd_mat[,2:k_max]
    }
    sd_out
}

#' @export
get_theta_sd.hmcfcop <- function(hmc) {
    t_max <- hmc$t_max
    n_max <- hmc$n_max
    k_max <- hmc$k_max
    structfactor <- hmc$structfactor
    all_param_sd = apply(hmc$sample_hmc, MARGIN = 2, FUN = sd)[-c(1:7)]

    all_theta <- all_param_sd[(t_max * k_max +1):length(hmc$mean_hmc)]
    theta <- rep(0,n_max)
    count <- 0
    if (structfactor == 1 | structfactor == 4 | structfactor == 2){
        count <- 0
    } else {
        if (structfactor == 3){
            count = k_max - 1 + sum(hmc$latent_copula_type == 2) - sum(hmc$latent_copula_type == 0)
        }
    }


    for (i in 1:n_max){
        if (hmc$cop_type[i] > 0){
            count = count + 1
            theta[i] <- all_theta[count]
            if (hmc$cop_type[i] == 2){
                count = count + 1
            }
        }
    }
    theta
}

#' @export
get_theta2_sd.hmcfcop <- function(hmc) {
    t_max <- hmc$t_max
    n_max <- hmc$n_max
    k_max <- hmc$k_max
    structfactor <- hmc$structfactor
    all_param_sd = apply(hmc$sample_hmc, MARGIN = 2, FUN = sd)[-c(1:7)]

    all_theta <- all_param_sd[(t_max * k_max +1):length(hmc$mean_hmc)]
    theta2 <- rep(0,n_max)
    count <- 0
    if (structfactor == 1 | structfactor == 4 | structfactor == 2){
        count <- 0
    } else {
        if (structfactor == 3){
            count = k_max - 1 + sum(hmc$latent_copula_type == 2) - sum(hmc$latent_copula_type == 0)
        }
    }


    for (i in 1:n_max){
        if (hmc$cop_type[i] > 0){
            count = count + 1
            if (hmc$cop_type[i] == 2){
                count = count + 1
                theta2[i] <- all_theta[count]
            }
        }
    }
    theta2
}

#' @export
get_latent_theta_sd.hmcfcop <- function(hmc) {
    t_max <- hmc$t_max
    n_max <- hmc$n_max
    k_max <- hmc$k_max
    structfactor <- hmc$structfactor
    all_param_sd = apply(hmc$sample_hmc, MARGIN = 2, FUN = sd)[-c(1:7)]

    all_theta <- all_param_sd[(t_max * k_max +1):length(hmc$mean_hmc)]
    latent_theta <- NULL
    count <- 0
    if (structfactor == 1){
        latent_theta = NULL
    }

    if (structfactor == 2){
        latent_theta <- rep(0,n_max)
        count = n_max + sum(hmc$cop_type == 2) - sum(hmc$cop_type == 0)
        for (i in 1:n_max){
            if (hmc$latent_copula_type[i] > 0){
                count = count + 1
                latent_theta[i] <- all_theta[count]

                if (hmc$latent_copula_type[i] == 2){
                    count = count + 1
                }
            }
        }
    }

    if (structfactor == 3){
        latent_theta <- rep(0,k_max-1)
        for (i in 1:(k_max-1)){
            if (hmc$latent_copula_type[i] > 0){
                count = count + 1
                latent_theta[i] <- all_theta[count]
                if (hmc$latent_copula_type[i] == 2){
                    count = count + 1
                }
            }
        }
    }
    latent_theta
}
#' @export
get_latent_theta2_sd.hmcfcop <- function(hmc) {
    t_max <- hmc$t_max
    n_max <- hmc$n_max
    k_max <- hmc$k_max
    structfactor <- hmc$structfactor
    all_param_sd = apply(hmc$sample_hmc, MARGIN = 2, FUN = sd)[-c(1:7)]

    all_theta <- all_param_sd[(t_max * k_max +1):length(hmc$mean_hmc)]
    latent_theta2 <- NULL

    count <- 0

    if (structfactor == 2){
        latent_theta2 <- rep(0,n_max)
        count = n_max + sum(hmc$cop_type == 2) - sum(hmc$cop_type == 0)
        for (i in 1:n_max){
            if (hmc$latent_copula_type[i] > 0){
                count = count + 1
                if (hmc$latent_copula_type[i] == 2){
                    count = count + 1
                    latent_theta2[i] <- all_theta[count]
                }
            }
        }
    }

    if (structfactor == 3){
        latent_theta2 <- rep(0,k_max-1)
        count = 0
        for (i in 1:(k_max-1)){
            if (hmc$latent_copula_type[i] > 0){
                count = count + 1

                if (hmc$latent_copula_type[i] == 2){
                    count = count + 1
                    latent_theta2[i] <- all_theta[count]
                }
            }
        }
    }
    latent_theta2
}
#' @export
get_vine_theta_sd.hmcfcop <- function(hmc) {
    t_max <- hmc$t_max
    n_max <- hmc$n_max
    k_max <- hmc$k_max
    structfactor <- hmc$structfactor
    all_param_sd = apply(hmc$sample_hmc, MARGIN = 2, FUN = sd)[-c(1:7)]

    all_theta <- all_param_sd[(t_max * k_max +1):length(hmc$mean_hmc)]
    vine_theta <- NULL
    count <- 0
    if (structfactor != 4){
        vine_theta = NULL
    }

    if (structfactor == 4){
        g_max = max(hmc$gid)
        vine_theta <- rep(0,n_max - g_max)
        count = n_max + sum(hmc$cop_type == 2) - sum(hmc$cop_type == 0)
        for (i in 1:length(hmc$vine_copula_type)){
            if (hmc$vine_copula_type[i] > 0){
                count = count + 1
                vine_theta[i] <- all_theta[count]

                if (hmc$vine_copula_type[i] == 2){
                    count = count + 1
                }
            }
        }
    }

    vine_theta
}
#' @export
get_vine_theta2_sd.hmcfcop <- function(hmc) {
    t_max <- hmc$t_max
    n_max <- hmc$n_max
    k_max <- hmc$k_max
    structfactor <- hmc$structfactor
    all_param_sd = apply(hmc$sample_hmc, MARGIN = 2, FUN = sd)[-c(1:7)]

    all_theta <- all_param_sd[(t_max * k_max +1):length(hmc$mean_hmc)]
    vine_theta2 <- NULL

    count <- 0
    if (structfactor != 4){
        vine_theta2 = NULL
    }

    if (structfactor == 4){
        g_max = max(hmc$gid)
        vine_theta2 <- rep(0,n_max - g_max)
        count = n_max + sum(hmc$cop_type == 2) - sum(hmc$cop_type == 0)
        for (i in 1:length(hmc$vine_copula_type)){
            if (hmc$vine_copula_type[i] > 0){
                count = count + 1
                if (hmc$vine_copula_type[i] == 2){
                    count = count + 1
                    vine_theta2[i] <- all_theta[count]
                }
            }
        }
    }
    vine_theta2
}

###############################################################################
