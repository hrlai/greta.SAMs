df <- expression(S * M)
df_fmm <- expression(S * (K-1) + K * M)

S <- 10
M <- 1:10
K <- 3

plot(M, eval(df), type = "b")
lines(M, eval(df_fmm), type = "b", col = "steelblue")
