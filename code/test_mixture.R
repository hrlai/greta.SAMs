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
source("code/util.R")

# data
y <- as.matrix(simulated_data[, 1:20])
x <- dat[, -1]

# indices
n <- nrow(y)  # number of sites
p <- ncol(y)  # number of species
m <- ncol(x)  # number of covariates
k <- 3        # number of archetypes

# archetype slopes
# beta  <- normal(0, 1, dim = c(k, m))  # need to be ordered variable?
beta <- triangular_ordered_matrix(k, m)

# species intercepts
beta0 <- normal(-1, 0.5, dim = p)   # using the known input in ecomix, for now

# site intercepts
sd_alpha <- exponential(1)
z_alpha  <- normal(0, 1, dim = n)
alpha    <- z_alpha * sd_alpha

# probability of archetype assignment
# not sure about the prior, using a weak prior like the ones in 
# stable isotope models
z <- dirichlet(t(rep(1/k, k)), n_realisations = p)

# From hereon we turn the response matrix to long format
# to fit a SAM using the mixture function in greta
y_long <- as.vector(y)
x_long <- x[as.vector(row(y)), ]

# indices for the long-format data
sites <- as.vector(row(y))
sp    <- as.vector(col(y))

# linear predictor
intercepts <- beta0[sp] + alpha[sites]
eta <- x_long %*% t(beta)
eta <- sweep(eta, 1, intercepts, "+")
prob <- ilogit(eta)

# likelihood and model
distribution(y_long) <- mixture(bernoulli(prob[, 1]),
                                bernoulli(prob[, 2]),
                                bernoulli(prob[, 3]),
                                weights = t(z[sp, ]))

mod <- model(beta, beta0, alpha, sd_alpha, z)

draws <- mcmc(mod, 
              sampler = hmc(15, 20),
              # initial_values = initials(beta = beta_true),  # potentially using kmeans here
              chain = 1)  # I presume we need one chain to avoid label switching?



# Diagnostics -------------------------------------------------------------

# look for label switching across chains
library(bayesplot)
mcmc_trace(draws, regex_pars = "beta\\[") 
mcmc_trace(draws, regex_pars = "z") 

mcmc_intervals(draws, regex_pars = "beta\\[")
beta_sim <- calculate(beta, values = draws, nsim = 1000)   # doesn't work, prob related to recent issues that Nick Tierney faced
apply(beta_sim$beta, 2:3, median); beta_true  # label switched but in reality we don't really care?
