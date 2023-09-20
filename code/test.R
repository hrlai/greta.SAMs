library(ecomix)

set.seed(42)
sam_form <-
    stats::as.formula(paste0('cbind(', paste(paste0('spp', 1:20),
                                             collapse = ','), ")~x1+x2"))
sp_form <- ~ 1
beta <- matrix(c(-2.9, -3.6, -0.9, 1, .9, 1.9), 3, 2, byrow = TRUE)
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
        beta = beta,
        family = "bernoulli"
    )




# greta -------------------------------------------------------------------

library(greta)

y <- as.matrix(simulated_data[, 1:20])
x <- dat[, -1]

n <- nrow(y)
p <- ncol(y)
m <- ncol(x)
k <- 3

beta  <- normal(0, 1, dim = c(k, m))
beta0 <- normal(0, 1, dim = p) 

z <- dirichlet(t(rep(0.1, k)), n_realisations = p)

eta <- x %*% t(beta) %*% t(z)
eta <- sweep(eta, 2, beta0, "+")
prob <- ilogit(eta)

distribution(y) <- bernoulli(prob)

mod <- model(beta, beta0, z)

o <- opt(mod)
# draws <- mcmc(mod, chains = 1)
# bayesplot::mcmc_trace(draws)
