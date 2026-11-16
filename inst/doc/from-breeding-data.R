## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(collapse = TRUE, comment = "#>", fig.width = 7, fig.height = 4)
library(BreedReliabR)

## ----data---------------------------------------------------------------------
set.seed(2026)

n_gen <- 8
n_env <- 6
n_rep <- 2

genotypes <- paste0("G", sprintf("%02d", seq_len(n_gen)))
environments <- paste0("E", sprintf("%02d", seq_len(n_env)))

env <- data.frame(
  Environment = environments,
  Heat = c(24, 27, 29, 31, 33, 35),
  Water = c(72, 88, 60, 80, 55, 75)
)

heat_mean <- mean(env$Heat)
heat_sd <- sd(env$Heat)
water_mean <- mean(env$Water)
water_sd <- sd(env$Water)

env$heat_z <- (env$Heat - heat_mean) / heat_sd
env$water_z <- (env$Water - water_mean) / water_sd
env$intercept <- 1

basis_names <- c("intercept", "heat_z", "water_z")
training_basis <- as.matrix(env[, basis_names])
rownames(training_basis) <- env$Environment

## ----phenotype----------------------------------------------------------------
true_coef <- cbind(
  intercept = rnorm(n_gen, 5.0, 0.30),
  heat_z = rnorm(n_gen, -0.25, 0.12),
  water_z = rnorm(n_gen, 0.25, 0.12)
)
rownames(true_coef) <- genotypes

dat <- expand.grid(
  Genotype = genotypes,
  Environment = environments,
  Rep = seq_len(n_rep),
  KEEP.OUT.ATTRS = FALSE,
  stringsAsFactors = FALSE
)

dat <- merge(dat, env, by = "Environment", sort = FALSE)
dat <- dat[
  order(match(dat$Genotype, genotypes),
        match(dat$Environment, environments),
        dat$Rep),
]
rownames(dat) <- NULL

g_index <- match(dat$Genotype, genotypes)
dat$Yield <-
  true_coef[g_index, "intercept"] +
  true_coef[g_index, "heat_z"] * dat$heat_z +
  true_coef[g_index, "water_z"] * dat$water_z +
  rnorm(nrow(dat), 0, 0.25)

head(dat[, c("Genotype", "Environment", "Rep", "Yield")])

## ----fit----------------------------------------------------------------------
p <- length(basis_names)

X <- matrix(0, nrow = nrow(dat), ncol = n_gen * p)
for (i in seq_along(genotypes)) {
  rows <- dat$Genotype == genotypes[i]
  cols <- ((i - 1) * p + 1):(i * p)
  X[rows, cols] <- as.matrix(dat[rows, basis_names])
}

coefficient_order <- unlist(
  lapply(genotypes, function(g) paste(g, basis_names, sep = "__")),
  use.names = FALSE
)
colnames(X) <- coefficient_order

fit <- lm.fit(X, dat$Yield)
stopifnot(fit$rank == ncol(X))

beta_hat <- fit$coefficients
coef_mean <- matrix(
  beta_hat,
  nrow = n_gen,
  ncol = p,
  byrow = TRUE,
  dimnames = list(genotypes, basis_names)
)

sigma2 <- sum(fit$residuals^2) / fit$df.residual
coef_vcov <- sigma2 * solve(crossprod(X))
rownames(coef_vcov) <- coefficient_order
colnames(coef_vcov) <- coefficient_order

dim(coef_mean)
dim(coef_vcov)

## ----checks-------------------------------------------------------------------
expected_order <- unlist(
  lapply(genotypes, function(g) paste(g, basis_names, sep = "__")),
  use.names = FALSE
)

stopifnot(
  identical(rownames(coef_vcov), expected_order),
  identical(colnames(coef_vcov), expected_order),
  all(is.finite(coef_mean)),
  all(is.finite(coef_vcov)),
  isSymmetric(coef_vcov, tol = 1e-8)
)

pev_eigen <- eigen(coef_vcov, symmetric = TRUE, only.values = TRUE)$values
pev_scale <- max(1, max(abs(pev_eigen)))
stopifnot(min(pev_eigen) >= -1e-8 * pev_scale)

## ----target-------------------------------------------------------------------
tpe_raw <- data.frame(
  State = paste0("TPE", 1:4),
  Heat = c(26, 30, 33, 35),
  Water = c(82, 72, 63, 56),
  Weight = c(0.20, 0.30, 0.30, 0.20)
)

tpe_raw$heat_z <- (tpe_raw$Heat - heat_mean) / heat_sd
tpe_raw$water_z <- (tpe_raw$Water - water_mean) / water_sd
tpe_raw$intercept <- 1

tpe_basis <- as.matrix(tpe_raw[, basis_names])
rownames(tpe_basis) <- tpe_raw$State
tpe_weights <- tpe_raw$Weight

## ----preflight----------------------------------------------------------------
input_check <- br_check_inputs(
  candidates = genotypes,
  coef_mean = coef_mean,
  coef_vcov = coef_vcov,
  basis_names = basis_names,
  tpe_basis = tpe_basis,
  tpe_weights = tpe_weights,
  training_basis = training_basis,
  k = 3,
  H = 5
)
input_check

## ----objects------------------------------------------------------------------
prediction <- br_prediction(
  candidates = genotypes,
  coef_mean = coef_mean,
  coef_vcov = coef_vcov,
  basis_names = basis_names,
  source = "Pedagogical genotype-specific reaction-norm model"
)

tpe <- br_tpe(
  basis = tpe_basis,
  weights = tpe_weights,
  state_id = tpe_raw$State
)

prediction
tpe

## ----support------------------------------------------------------------------
support <- br_support(tpe, training_basis)
support
plot(support)

## ----omega--------------------------------------------------------------------
omega <- br_omega(prediction, tpe, H = c(1, 3, 5, 10))
omega$summary
head(omega$candidate)

## ----advancement--------------------------------------------------------------
comparison <- br_compare(
  prediction, tpe, H = 5, k = 3,
  nsim = 1000, seed = 2026
)
summary(comparison)

advancement <- br_advance(
  prediction, tpe, H = 5, k = 3,
  nsim = 2000, method = "M3", seed = 2026
)
advancement$candidates
sum(advancement$candidates$p_advance)

## ----curve--------------------------------------------------------------------
curve_results <- br_curve(
  prediction, tpe,
  H = c(1, 3, 5, 10),
  k = 3,
  nsim = 1000,
  seed = 2026
)
summary(curve_results)
plot(curve_results, type = "reliability")
plot(curve_results, type = "omega")

## ----export-------------------------------------------------------------------
template_dir <- file.path(tempdir(), "BreedReliabR_beginner_example")
dir.create(template_dir, showWarnings = FALSE, recursive = TRUE)

candidate_coefficients <- data.frame(
  candidate = rownames(coef_mean),
  coef_mean,
  check.names = FALSE
)

write.csv(
  candidate_coefficients,
  file.path(template_dir, "candidate_coefficients.csv"),
  row.names = FALSE
)

saveRDS(
  coef_vcov,
  file.path(template_dir, "candidate_joint_PEV.rds")
)

write.csv(
  data.frame(
    state = tpe_raw$State,
    tpe_basis,
    weight = tpe_weights,
    check.names = FALSE
  ),
  file.path(template_dir, "target_environments.csv"),
  row.names = FALSE
)

write.csv(
  data.frame(
    environment = rownames(training_basis),
    training_basis,
    check.names = FALSE
  ),
  file.path(template_dir, "training_environments.csv"),
  row.names = FALSE
)

