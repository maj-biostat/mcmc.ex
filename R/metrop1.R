library(data.table)
library(ggplot2)

# Metropolis Hastings MCMC
# For the MCMC example, there is a data generation function,
# an assumed model (here aligned with the data generation process)
# a likelihood and a prior specification.

get_data <- function(
  N = 31, b0 = 0, b1 = 5, se = 10
  ){

  # mean centering gives better sampling
  x <- (-(N-1)/2):((N-1)/2)
  y <- b0 + b1 * x + rnorm(N, 0, se)

  data.table(x, y)

}


lik <- function(d, lpar){
  mu = lpar[["b0"]] + lpar[["b1"]] * d$x
  # individual contributions to the likelihood
  # (on log scale to avoid overflow)
  individ_lik <- dnorm(d$y, mu, lpar[["se"]], log = T)
  sumll <- sum(individ_lik)
  sumll
}

pri <- function(lpar){
  b0pri <- dnorm(lpar[["b0"]], 0, 10, log = T)
  b1pri <- dnorm(lpar[["b1"]], 0, 10, log = T)
  sepri <- dexp(lpar[["se"]], 1, log = T)
  b0pri + b1pri + sepri
}

postr <- function(d, lpar){
  lik(d, lpar) + pri(lpar)
}

# The actual workhorse that constructs samples from the
# posterior. It can be used to sample from any integrable
# function. The goal is to jump around the parameter space
# in a way that the probability to be at a point is proportional
# to the function we sample from (the target). In our case, the
# target is the posterior defined above.
#
# Start from a random place (well chosen)
# Choose a new param close to the old one based on a proposal dist
# Jump to this new point with prob p(new)/p(old) where the p
# is the target function and p>1 means jumping as well
#
# lpar0 initial values
# B number of smpls
mcmc_mh <- function(d, lpar0 = c(0, 1, 1), B = 1000){
  # proposal function
  propl <-function(lpar){
    # b0, b1, se
    pp <- rnorm(3, mean = unlist(lpar), sd = c(0.5, 0.5, 0.3))
    l <- list(b0 = pp[1], b1 = pp[2], se = pp[3])
    l
  }
  smpl <- array(dim= c(B+1, 3))
  colnames(smpl) <- c("b0", "b1", "se")
  smpl[1, ] <- unlist(lpar0)
  i <- 1
  for(i in 1:B){
    proposal <- propl(smpl[i, ])

    # because we are working in logarithms, the probab is obtained
    # using p1/p2 = exp(log(p1) - log(p2))

    probab = exp(postr(d, proposal) - postr(d, smpl[i, ]))
    if(runif(1) < probab){
      # accept
      smpl[i+1, ] <- unlist(proposal)
    } else {
      # reject
      smpl[i+1, ] <- smpl[i, ]
    }
  }
  smpl
}


main <- function(){

  d <- get_data()
  smpls <- data.table(mcmc_mh(d, B = 100000))
  smpls[, id := 1:.N]

  burnIn <- 5000
  acceptance <- 1 - mean(duplicated(smpls[id>burnIn, .(b0, b1, se)]))
  message("Target accept rates around 20-30%, actual was: ", acceptance)
  message("Acceptance rates influenced by proposal function.")
  message("Very high rates imply that the alg is getting stuck at the same point.")

  setcolorder(smpls, "id")
  dfig <- melt(smpls, measure.vars = names(smpls)[-1])
  ggplot(dfig, aes(x = id, y = value)) +
    geom_line() +
    facet_wrap(~variable, scales = "free")
#

  # summary stats (discarding burn-in)
  dfig[id > burnIn , .(
    val = sprintf("%.2f (%.2f,%.2f)",
                  mean(value),
                  quantile(value, 0.025),
                  quantile(value, 0.975))
    ), keyby = variable]

  # compare to linear model results
  f1 <- lm(y~x, data = d)
  cbind(coef(f1), confint(f1))
  sigma(f1)
}
