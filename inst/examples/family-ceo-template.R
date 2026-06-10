library(haven)
library(psmqte)

ceo <- haven::read_dta("psmqte.dta")

# Expected variables:
#   chxs_opass: outcome
#   heir:       treatment, coded 1 for family-related incoming CEO
#   uxs_opass:  prior OROA
#   bown:       board ownership
#   highfamd:   family directors indicator

qfit <- psmqte(
  chxs_opass ~ uxs_opass + bown + highfamd,
  treatment = "heir",
  data = ceo,
  taus = seq(0.1, 0.9, by = 0.1),
  M = c(1, 2, 4, 6, 8, 13),
  multiplier = 499,
  seed = 20260531
)

qte_qtt_table <- summary(qfit)
reuse_table <- qfit$diagnostics$reuse
overlap_table <- qfit$diagnostics$overlap$by_arm

afit <- psmate(
  chxs_opass ~ uxs_opass + bown + highfamd,
  treatment = "heir",
  data = ceo,
  M = c(1, 2, 4, 6, 8, 13)
)

ate_att_table <- summary(afit)
