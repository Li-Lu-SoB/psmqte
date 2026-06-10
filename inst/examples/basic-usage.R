library(psmqte)

dat <- simulate_psmqte_data(1000, seed = 1)

qfit <- psmqte(
  Y ~ X1 + X2,
  treatment = "D",
  data = dat,
  taus = seq(0.1, 0.9, by = 0.1),
  M = c(1, 4, 8),
  se = TRUE,
  seed = 1
)

print(qfit)
head(summary(qfit))
qfit$diagnostics$reuse

afit <- psmate(
  Y ~ X1 + X2,
  treatment = "D",
  data = dat,
  M = c(1, 4, 8)
)

print(afit)

if (requireNamespace("causaldata", quietly = TRUE)) {
  nhefs <- psmqte_example_data("nhefs")
  nhefs_fit <- psmqte(
    Y ~ age + sex + race + school + smokeintensity + smokeyrs +
      exercise + active + wt71,
    treatment = "D",
    data = nhefs,
    taus = c(0.25, 0.5, 0.75),
    M = c(1, 2),
    se = FALSE
  )
  summary(nhefs_fit)
} else {
  message("Install causaldata to run the NHEFS example: install.packages('causaldata')")
}
