df <- expression(S * M)
df_fmm <- expression(S * (K-1) + K * M)

S <- 10    # number of species
M <- 1:10  # number of covariates
K <- 3     # number of archetypes

plot(M, eval(df), type = "b")
lines(M, eval(df_fmm), type = "b", col = "steelblue")
