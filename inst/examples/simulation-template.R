library(psmqte)

set.seed(2026)

taus <- c(0.25, 0.5, 0.75)
M_grid <- m_grid_diverging(512)

one_rep <- function(n, seed) {
  dat <- simulate_psmqte_data(n, seed = seed)
  fit <- psmqte(
    Y ~ X1 + X2,
    treatment = "D",
    data = dat,
    taus = taus,
    M = M_grid,
    se = FALSE
  )
  summary(fit)
}

mc <- do.call(rbind, lapply(seq_len(10), function(r) one_rep(512, 1000 + r)))
aggregate(cbind(qte, qtt) ~ M + tau, data = mc, FUN = mean)
