# This script aims to simulate data using ecomix, and then hopefully 
# recover the parameters using a greta model


# Simulate data following ecomix example ----------------------------------

library(ecomix)

set.seed(42)
sam_form <-
    stats::as.formula(paste0('cbind(', paste(paste0('spp', 1:20),
                                             collapse = ','), ")~x1+x2"))
sp_form <- ~ 1
beta_true <- matrix(c(-2.9, -3.6, -0.9, 1, .9, 1.9), 3, 2, byrow = TRUE)
dat <- data.frame(
    y = rep(1, 100),
    x1 = stats::runif(100, 0, 2.5),
    x2 = stats::rnorm(100, 0, 2.5)
)
dat[, -1] <- scale(dat[, -1])
simulated_data <-
    species_mix.simulate(
        archetype_formula = sam_form,
        species_formula = sp_form,
        data = dat,
        beta = beta_true,
        family = "bernoulli"
    )




# greta -------------------------------------------------------------------

library(greta)

# data
y <- as.matrix(simulated_data[, 1:20])
x <- dat[, -1]

# indices
n <- nrow(y)  # number of sites
p <- ncol(y)  # number of species
m <- ncol(x)  # number of covariates
k <- 3        # number of archetypes

# slopes and intercepts
beta  <- normal(0, 1, dim = c(k, m))  # need to be ordered variable?
beta0 <- normal(-1, 0.5, dim = p)   # using the known input in ecomix, for now

# probability of archetype assignment
# not sure about the prior, using a weak prior like the ones in 
# stable isotope models
z <- dirichlet(t(rep(1/k, k)), n_realisations = p)

# linear predictor
eta <- x %*% t(beta) %*% t(z)
eta <- sweep(eta, 2, beta0, "+")
prob <- ilogit(eta)

# likelihood and model
distribution(y) <- bernoulli(prob)

mod <- model(beta, beta0, z)

# optimisation / draws
draws <- mcmc(mod, 
              sampler = hmc(15, 20),
              # initial_values = initials(beta = beta_true),
              chain = 1)




# Diagnostics -------------------------------------------------------------

# look for label switching across chains
library(bayesplot)
mcmc_trace(draws) 

mcmc_intervals(draws, regex_pars = "beta\\[")
beta_sim <- calculate(beta, values = draws, nsim = 1000)
apply(beta_sim$beta, 2:3, mean); beta_true  # label switched but in reality we don't really care?
