library(data.table)
library(splines)
library(ggplot2)

# Metropolis Hastings MCMC
#
# Logistic regression model.
#
# For the MCMC example, there is a data generation function,
# an assumed model (here aligned with the data generation process)
# a likelihood and a prior specification.

get_data <- function(
    N = 31,
    p = c(0.15, 0.25, 0.25, 0.25, 0.45, 0.7, 0.7, 0.7),
    xp = c(seq(0, 400, len = 5), 490, 750, 1000)

){
  spl1 <- splinefun(xp, p, method = "hyman")
  xpred <- c(0, sort(sample(1:999, N-2)), 1000)
  ppred <- spl1(xpred)
  # plot(xpred, ppred)
  trials <- sample(50:100, N, replace = T)
  y <- rbinom(N, trials, ppred)
  d <- data.table(
    x = xpred,
    y = y,
    n = trials,
    p_tru = ppred
  )
  d
}

lik <- function(y, n, X, lpar){
  b0 <- lpar[[1]]
  b <- unlist(lapply(2:length(lpar), function(z) lpar[[z]]))
  p = plogis(X %*% c(b0, b))
  # individual contributions to the likelihood
  # (on log scale to avoid overflow)
  individ_lik <- sapply(seq_along(y), function(ii){
    dbinom(y[ii], n[ii], p[ii], log = T)
  })
  sumll <- sum(individ_lik)
  sumll
}

pri <- function(lpar, hyp){
  b0 <- lpar[[1]]
  b <- unlist(lapply(2:length(lpar), function(z) lpar[[z]]))
  if(is.null(hyp)){
    # default priors
    b0pri <- dnorm(b0, 0, 10, log = T)
    bpri <- sapply(b, function(z) dnorm(z, 0, 10, log = T))
  } else {
    b0pri <- dnorm(b0, hyp$mu0, hyp$sd0, log = T)
    bpri <-  sapply(b, function(z) dnorm(z, hyp$mu, hyp$sd, log = T))
  }
  b0pri + sum(bpri)
}

postr <- function(y, n, X, lpar, hyp){
  lik(y, n, X, lpar) + pri(lpar, hyp)
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
mcmc_mh <- function(
    y, n, X,
    lpar0 = NULL,
    hyp = NULL,
    iter = 100,
    jumps = NULL
){
  # proposal function
  propl <-function(lpar, jumps){
    pp <- rnorm(length(lpar), mean = unlist(lpar), sd = jumps)
    pp
  }
  smpl <- array(dim= c(iter+1, ncol(X)))
  colnames(smpl) <- paste0("b", 0:(ncol(X)-1))
  smpl[1, ] <- unlist(lpar0)
  i <- 1
  for(i in 1:iter){
    proposal <- propl(smpl[i, ], jumps)
    # because we are working in logarithms, the probab is obtained
    # using p1/p2 = exp(log(p1) - log(p2))
    probab = exp(postr(y, n, X, proposal, hyp) - postr(y, n, X, smpl[i, ], hyp))
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

  p = c(0.15, 0.25, 0.25, 0.25, 0.45, 0.7, 0.7, 0.7)
  xp = c(seq(0, 400, len = 5), 490, 750, 1000)

  d <- get_data(N = 31, p, xp)
  # set up design matrix
  knots_inner = c(150, 300, 500, 700 ,850)
  B <- splines::bs(d$x, knots=knots_inner, degree=3, intercept = F)
  X <- cbind(1, B)
  y = d$y
  n = d$n
  iter <- 10000

  smpls <- mcmc_mh(
      y, n, X,
      lpar0 = c(-1, rep(0, ncol(B))),
      hyp = list(mu0=-1, sd0=3, mu=0, sd=1),
      iter = iter,
      jumps = rep(0.1, ncol(X))
    )

  burnIn <- 5000
  acceptance <- 1 - mean(duplicated(smpls[burnIn:nrow(smpls), ]))
  message("Target accept rates around 20-30%, actual was: ", acceptance)
  message("Acceptance rates influenced by proposal function.")
  message("Very high rates imply that the alg is getting stuck at the same point.")

  dall <- data.table(id = 1:nrow(smpls), smpls)
  dall <- melt(dall, measure.vars = colnames(smpls))
  ggplot(dall, aes(x = id, y = value)) +
    geom_line() +
    facet_wrap(~variable, scales = "free")
  #

  # summary stats (discarding burn-in and thin)
  dsmpl <- copy(dall[id > burnIn & id %% 2 == 0])
  dsmpl[ , .(
    val = sprintf("%.2f (%.2f,%.2f)",
                  mean(value),
                  quantile(value, 0.025),
                  quantile(value, 0.975))
  ), keyby = variable]

  # prediction
  Xpred <- cbind(1, splines::bs(0:1000, knots=knots_inner, degree=3, intercept = F))
  eta_hat <- smpls[burnIn:nrow(smpls), ] %*% t(Xpred)
  p_hat <- apply(eta_hat, 2, plogis)

  # compare to glm results
  f1 <- glm(cbind(y,n-y) ~ -1 + X, family = binomial)
  d[, p_freq := predict(f1,type = "response")]

  # figure
  dfig1 <- data.table(
    x = 0:1000,
    mu = colMeans(p_hat),
    mu_lb = apply(p_hat, 2, quantile, 0.025),
    mu_ub = apply(p_hat, 2, quantile, 0.975)
  )

  ggplot(dfig1, aes(x = x, y= mu)) +
    geom_ribbon(aes(ymin = mu_lb, ymax = mu_ub), alpha = 0.2) +
    geom_line() +
    geom_point(
      data = d,
      aes(x = x, y = y/n), inherit.aes = F, col = 1
    ) +
    geom_point(
      data = d,
      aes(x = x, y = p_tru), inherit.aes = F, col = 2, pch = 8
    ) +
    geom_line(
      data = d,
      aes(x = x, y = p_freq), inherit.aes = F, lty = 2
    )
#



}
