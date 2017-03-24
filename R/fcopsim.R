#' Sum of vector elements.
#'
#' \code{sum} returns the sum of all the values present in its arguments.
#'
#' This is a generic function: methods can be defined for it directly
#' or via the \code{\link{Summary}} group generic. For this to work properly,
#' the arguments \code{...} should be unnamed, and dispatch is on the
#' first argument.
#'
#' @param ... Numeric, complex, or logical vectors.
#' @param na.rm A logical scalar. Should missing values (including NaN)
#'   be removed?
#' @return If all inputs are integer and logical, then the output
#'   will be an integer. If integer overflow
#'   \url{http://en.wikipedia.org/wiki/Integer_overflow} occurs, the output
#'   will be NA with a warning. Otherwise it will be a length-one numeric or
#'   complex vector.
#'
#'   Zero-length vectors have sum 0 by definition. See
#'   \url{http://en.wikipedia.org/wiki/Empty_sum} for more details.
#' @examples
#' sum(1:10)
#' sum(1:5, 6:10)
#' sum(F, F, F, T, T)
#'
#' sum(.Machine$integer.max, 1L)
#' sum(.Machine$integer.max, 1)
#'
#' \dontrun{
#' sum("a")
#' }
rtheta <-  function(family, theta = TRUE) {
    if (theta){
        if (family == 1) theta_gen = runif(1, min = 0, max = 1)
        if (family == 2) theta_gen = runif(1, min = 0, max = 1)
        if (family == 3) theta_gen = runif(1, min = 0.1, max = 2)
        if (family == 4) theta_gen = runif(1, min = 1.1, max = 3)
        if (family == 5) theta_gen = runif(1, min = 1, max = 4)
        if (family == 6) theta_gen = runif(1, min = 1.3, max = 2.5)

    } else {
        if (family == 1) theta_gen = 0
        if (family == 2) theta_gen = runif(1, min = 2, max = 10)
        if (family == 3) theta_gen = 0
        if (family == 4) theta_gen = 0
        if (family == 5) theta_gen = 0
        if (family == 6) theta_gen = 0
    }
    theta_gen
}

fcopsim <- function(t_max, n_max, k_max = 1, family, gid = 0, structfactor = 1, seed_num = 0) {
    set.seed(seed_num)

    if (! all.equal(t_max, as.integer(t_max)))
        stop("'t_max' has to be a integer")
    if (! all.equal(n_max, as.integer(n_max)))
        stop("'n_max' has to be a integer")
    if (! all.equal(k_max, as.integer(k_max)))
        stop("'n_max' has to be a integer")
    if (! (t_max > 0) & (n_max > 0) & (k_max > 0) )
        stop("'t_max', 'n_max', 'k_max' has to be greater than 0")

    if (!(length(family) %in% c(1, n_max)))
        stop("'family' has to be a single number or a size n_max vector")

    if (gid == 0) {
            gid = rep(1,n_max)
        }
    if (length(family) == 1){
        family = matrix(family, nrow = n_max, ncol = k_max)
        }

    u <- matrix(0,nrow=t_max, ncol = n_max)
    theta <- matrix(0,nrow=n_max, ncol = k_max)
    theta2 <- matrix(0,nrow=n_max, ncol = k_max)

    if (structfactor == 1)
    {
        v <- matrix(runif(t_max*k_max), nrow = t_max, ncol = k_max)
        if (k_max == 1)
        {
            for (i in 1:n_max){
                theta[i,1] <- rtheta(family[i,1], TRUE)
                theta2[i,1] <- rtheta(family[i,1], FALSE)

                obj <- BiCop(family = family[i,1], par = theta[i], par2 = theta2[i])
                u[,i] <- BiCopCondSim(t_max, cond.val = v[,1], cond.var = 1, obj)
            }
        }

    }
    list( u = u,
          v = v,
          t_max = t_max,
          n_max = n_max,
          k_max = k_max,
          family = family,
          theta = theta,
          theta2 = theta2,
          gid = gid,
          structfactor = structfactor)
}
