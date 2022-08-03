library(data.table)
library(ggplot2)

# Metropolis Hastings MCMC
#
# Logistic regression model.
#
# For the MCMC example, there is a data generation function,
# an assumed model (here aligned with the data generation process)
# a likelihood and a prior specification.

get_data <- function(
    N = 31, b0 = -2.2, b1 = 1.35
){

  # in real life you will need to
  # scale the covariates that have a large range otherwise sampling
  # will be dreadful.
  # for example, if the covariate is 0 to 1000, then maybe conver it to
  # -1 to 1.
  # note that you will need to adjust the parameter as well
  # but that should just be a linear transform on the log odds scale.
  x <- sort(runif(N, -1, 1))
  eta <- b0 + b1 * x
  # plot(x, plogis(eta))
  trials <- sample(50:100, N, replace = T)
  y <- rbinom(N, trials, plogis(eta))

  d <- data.table(
    x = x,
    y = y,
    n = trials
  )
  d
}

lik <- function(d, lpar){
  b0 <- lpar[[1]]
  b1 <- lpar[[2]]
  p = plogis(b0 + b1 * d$x)
  # individual contributions to the likelihood
  # (on log scale to avoid overflow)
  individ_lik <- sapply(seq_along(d$y), function(ii){
    dbinom(d$y[ii], d$n[ii], p[ii], log = T)
  })
  sumll <- sum(individ_lik)
  sumll
}

pri <- function(lpar, hyp){
  b0 <- lpar[[1]]
  b1 <- lpar[[2]]
  if(is.null(hyp)){
    # default priors
    b0pri <- dnorm(b0, 0, 10, log = T)
    b1pri <- dnorm(b1, 0, 10, log = T)
  } else {
    b0pri <- dnorm(b0, hyp$mu0, hyp$sd0, log = T)
    b1pri <- dnorm(b1, hyp$mu1, hyp$sd1, log = T)
  }
  b0pri + b1pri
}

postr <- function(d, lpar, hyp){
  lik(d, lpar) + pri(lpar, hyp)
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
    d,
    lpar0 = c(-1, 0),
    hyp = NULL,
    B = 100,
    jumps = c(0.5, 0.5)
    ){
  # proposal function
  propl <-function(lpar, jumps){
    pp <- rnorm(2, mean = unlist(lpar), sd = jumps)
    l <- list(b0 = pp[1], b1 = pp[2])
    l
  }
  smpl <- array(dim= c(B+1, 2))
  colnames(smpl) <- c("b0", "b1")
  smpl[1, ] <- unlist(lpar0)
  i <- 1
  for(i in 1:B){
    proposal <- propl(smpl[i, ], jumps)
    # because we are working in logarithms, the probab is obtained
    # using p1/p2 = exp(log(p1) - log(p2))
    probab = exp(postr(d, proposal, hyp) - postr(d, smpl[i, ], hyp))
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
  smpls <- data.table(
    mcmc_mh(
      d,
      lpar0 = c(-1, 0),
      hyp = list(mu0=-1, sd0=3, mu1=0, sd1=1),
      B = 10000,
      jumps = c(0.2, 0.2)
      )
  )
  smpls[, id := 1:.N]

  burnIn <- 5000
  acceptance <- 1 - mean(duplicated(smpls[id>burnIn, .(b0, b1)]))
  message("Target accept rates around 20-30%, actual was: ", acceptance)
  message("Acceptance rates influenced by proposal function.")
  message("Very high rates imply that the alg is getting stuck at the same point.")

  setcolorder(smpls, "id")
  dall <- melt(smpls, measure.vars = names(smpls)[-1])
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

  # compare to glm results
  f1 <- glm(cbind(y,n-y)~x, data = d, family = binomial)
  cbind(coef(f1), confint(f1))
}


