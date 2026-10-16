## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(collapse = TRUE, comment = "#>", fig.width = 7, fig.height = 4)
library(BreedReliabR)

## ----import-------------------------------------------------------------------
data(maize_demo)
d <- maize_demo
head(d$training_environments)
head(d$training_phenotypes)
d$candidates
pred <- do.call(br_prediction, d$prediction_inputs)
pred

## ----tpe----------------------------------------------------------------------
d$tpe_states
tpe <- br_tpe(d$tpe_basis, weights = d$tpe_states$weight,
              state_id = d$tpe_states$state_id)
tpe

## ----support------------------------------------------------------------------
support <- br_support(tpe, d$training_basis)
support
support$states
support$ranges
plot(support)

## ----omega--------------------------------------------------------------------
om <- br_omega(pred, tpe, H = c(1, 3, 5, 10, 20))
om
head(om$candidate)

## ----comparison---------------------------------------------------------------
comparison <- br_compare(pred, tpe, H = 5, k = 3, nsim = 2000, seed = 410)
summary(comparison)
plot(comparison)

## ----advancement--------------------------------------------------------------
joint <- comparison$analyses$M3
joint$candidates
joint$expected_merit_set
joint$probability_set
stopifnot(abs(sum(joint$candidates$p_advance) - 3) < 1e-10)

## ----curve--------------------------------------------------------------------
curve <- br_curve(pred, tpe, H = c(1, 3, 5, 10, 20), k = 3,
                  nsim = 2000, seed = 410)
summary(curve)
plot(curve)
plot(curve, type = "omega")

## ----sensitivity--------------------------------------------------------------
# A transparent alternative scenario, not an automatic correction of the TPE.
w <- d$tpe_states$weight
w[d$tpe_states$state_id == "extreme"] <- 0
tpe_alternative <- br_tpe(d$tpe_basis, w, d$tpe_states$state_id)
br_support(tpe_alternative, d$training_basis)
summary(br_compare(pred, tpe_alternative, H = 5, k = 3,
                   nsim = 2000, seed = 410))

