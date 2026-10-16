## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(collapse = TRUE, comment = "#>", fig.width = 6,
                      fig.height = 4)

## ----objects------------------------------------------------------------------
library(BreedReliabR)
data(maize_demo)
prediction <- do.call(br_prediction, maize_demo$prediction_inputs)
tpe <- br_tpe(maize_demo$tpe_basis,
              weights = maize_demo$tpe_states$weight,
              state_id = maize_demo$tpe_states$state_id)
prediction
tpe

## ----support-omega------------------------------------------------------------
support <- br_support(tpe, maize_demo$training_basis)
support$summary
omega <- br_omega(prediction, tpe, H = c(1, 3, 5, 10, 21))
omega$summary

## ----compare------------------------------------------------------------------
comparison <- br_compare(prediction, tpe, H = 5, k = 3,
                         nsim = 2000, seed = 410)
comparison$summary
head(comparison$candidates)

## ----curve--------------------------------------------------------------------
curve <- br_curve(prediction, tpe, H = c(1, 3, 5, 10), k = 3,
                 nsim = 2000, seed = 410)
curve$summary
plot(curve, type = "reliability")
plot(curve, type = "omega")

