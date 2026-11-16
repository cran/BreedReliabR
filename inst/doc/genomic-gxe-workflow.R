## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(collapse = TRUE, comment = "#>", fig.width = 7, fig.height = 4)
library(BreedReliabR)

## ----dimensions---------------------------------------------------------------
set.seed(2026)

n_train <- 24
n_candidate <- 8
n_total <- n_train + n_candidate
n_env <- 8
n_rep <- 2
n_markers <- 150

training_genotypes <- paste0("T", sprintf("%02d", seq_len(n_train)))
candidate_genotypes <- paste0("C", sprintf("%02d", seq_len(n_candidate)))
all_genotypes <- c(training_genotypes, candidate_genotypes)

## ----environment--------------------------------------------------------------
env <- data.frame(
  Environment = paste0("E", sprintf("%02d", seq_len(n_env))),
  Heat = c(24, 27, 29, 31, 33, 35, 28, 32),
  Water = c(72, 88, 60, 80, 55, 75, 65, 90)
)

cor(env$Heat, env$Water)

heat_mean <- mean(env$Heat)
heat_sd <- sd(env$Heat)
water_mean <- mean(env$Water)
water_sd <- sd(env$Water)

env$heat_z <- (env$Heat - heat_mean) / heat_sd
env$water_z <- (env$Water - water_mean) / water_sd
env$intercept <- 1

basis_names <- c("intercept", "heat_z", "water_z")
p <- length(basis_names)
training_basis <- as.matrix(env[, basis_names])
rownames(training_basis) <- env$Environment

## ----markers------------------------------------------------------------------
marker_frequency <- runif(n_markers, min = 0.10, max = 0.50)
M <- sapply(
  marker_frequency,
  function(freq) rbinom(n_total, size = 2, prob = freq)
)
M <- as.matrix(M)
rownames(M) <- all_genotypes
colnames(M) <- paste0("SNP", seq_len(n_markers))
stopifnot(all(dim(M) == c(n_total, n_markers)))

p_train <- colMeans(M[training_genotypes, , drop = FALSE]) / 2
W <- sweep(M, MARGIN = 2, STATS = 2 * p_train, FUN = "-")
relationship_denominator <- 2 * sum(p_train * (1 - p_train))
stopifnot(is.finite(relationship_denominator), relationship_denominator > 0)

K_raw <- tcrossprod(W) / relationship_denominator
K_raw <- (K_raw + t(K_raw)) / 2
rownames(K_raw) <- all_genotypes
colnames(K_raw) <- all_genotypes

ridge <- 1e-6 * max(1, mean(diag(K_raw)))
K_model <- K_raw + diag(ridge, n_total)
rownames(K_model) <- all_genotypes
colnames(K_model) <- all_genotypes

relationship_rcond <- rcond(K_model)
stopifnot(is.finite(relationship_rcond), relationship_rcond > 0)

## ----genomic-covariance-------------------------------------------------------
G_beta <- matrix(
  c(
     0.20,   0.02,  -0.01,
     0.02,   0.06,  -0.015,
    -0.01,  -0.015,  0.05
  ),
  nrow = p,
  byrow = TRUE,
  dimnames = list(basis_names, basis_names)
)

G_beta_eigen <- eigen(G_beta, symmetric = TRUE, only.values = TRUE)$values
stopifnot(min(G_beta_eigen) > 0)

G_u <- kronecker(K_model, G_beta)
stopifnot(all(dim(G_u) == c(n_total * p, n_total * p)))

## ----simulate-----------------------------------------------------------------
chol_G <- chol(G_u)
u_true <- as.vector(t(chol_G) %*% rnorm(n_total * p))

beta_true <- c(
  intercept = 5.00,
  heat_z = -0.20,
  water_z = 0.25
)

dat <- expand.grid(
  Genotype = training_genotypes,
  Environment = env$Environment,
  Rep = seq_len(n_rep),
  KEEP.OUT.ATTRS = FALSE,
  stringsAsFactors = FALSE
)
dat <- merge(dat, env, by = "Environment", sort = FALSE)
dat <- dat[
  order(match(dat$Genotype, training_genotypes),
        match(dat$Environment, env$Environment),
        dat$Rep),
]
rownames(dat) <- NULL

stopifnot(
  nrow(dat) == n_train * n_env * n_rep,
  !any(dat$Genotype %in% candidate_genotypes)
)

X <- as.matrix(dat[, basis_names])
Z <- matrix(0, nrow = nrow(dat), ncol = n_total * p)
for (i in seq_len(n_train)) {
  rows <- dat$Genotype == training_genotypes[i]
  cols <- ((i - 1) * p + 1):(i * p)
  Z[rows, cols] <- as.matrix(dat[rows, basis_names])
}

sigma_e2 <- 0.25^2
dat$Yield <- as.vector(
  X %*% beta_true + Z %*% u_true +
    rnorm(nrow(dat), mean = 0, sd = sqrt(sigma_e2))
)

## ----henderson----------------------------------------------------------------
G_u_inv <- solve(G_u)
R_inv_scalar <- 1 / sigma_e2

C11 <- crossprod(X) * R_inv_scalar
C12 <- crossprod(X, Z) * R_inv_scalar
C21 <- t(C12)
C22 <- crossprod(Z) * R_inv_scalar + G_u_inv
C <- rbind(cbind(C11, C12), cbind(C21, C22))

stopifnot(all(dim(C) == c(p + n_total * p, p + n_total * p)))
C_rcond <- rcond(C)
stopifnot(is.finite(C_rcond), C_rcond > 0)

rhs <- c(
  as.vector(crossprod(X, dat$Yield) * R_inv_scalar),
  as.vector(crossprod(Z, dat$Yield) * R_inv_scalar)
)
theta_hat <- solve(C, rhs)

## ----targets------------------------------------------------------------------
n_targets <- n_candidate * p
D <- matrix(0, nrow = n_targets, ncol = length(theta_hat))
target_names <- character(n_targets)
candidate_index <- (n_train + 1):n_total

row_id <- 1L
for (g in candidate_index) {
  for (j in seq_len(p)) {
    D[row_id, j] <- 1
    random_position <- p + (g - 1) * p + j
    D[row_id, random_position] <- 1
    target_names[row_id] <- paste(all_genotypes[g], basis_names[j], sep = "__")
    row_id <- row_id + 1L
  }
}
rownames(D) <- target_names

candidate_coef_vector <- as.vector(D %*% theta_hat)
coef_mean <- matrix(
  candidate_coef_vector,
  nrow = n_candidate,
  ncol = p,
  byrow = TRUE,
  dimnames = list(candidate_genotypes, basis_names)
)
stopifnot(all(dim(coef_mean) == c(n_candidate, p)))
round(coef_mean, 4)

## ----pev----------------------------------------------------------------------
C_inv_Dt <- solve(C, t(D))
coef_vcov <- D %*% C_inv_Dt
rownames(coef_vcov) <- target_names
colnames(coef_vcov) <- target_names

stopifnot(all(dim(coef_vcov) == c(n_targets, n_targets)))

pev_scale <- max(1, max(abs(coef_vcov)))
pev_asymmetry <- max(abs(coef_vcov - t(coef_vcov)))
stopifnot(pev_asymmetry <= 1e-8 * pev_scale)

pev_eigen <- eigen(coef_vcov, symmetric = TRUE, only.values = TRUE)$values
eigen_scale <- max(1, max(abs(pev_eigen)))
stopifnot(min(pev_eigen) >= -1e-8 * eigen_scale)

expected_order <- unlist(
  lapply(candidate_genotypes, function(g) paste(g, basis_names, sep = "__")),
  use.names = FALSE
)
stopifnot(
  identical(rownames(coef_vcov), expected_order),
  identical(colnames(coef_vcov), expected_order)
)

## ----tpe-data-----------------------------------------------------------------
tpe_raw <- data.frame(
  State = paste0("TPE", seq_len(4)),
  Heat = c(27, 30, 33, 35),
  Water = c(82, 72, 63, 56),
  Weight = c(0.20, 0.30, 0.30, 0.20)
)

tpe_raw$heat_z <- (tpe_raw$Heat - heat_mean) / heat_sd
tpe_raw$water_z <- (tpe_raw$Water - water_mean) / water_sd
tpe_raw$intercept <- 1

tpe_basis <- as.matrix(tpe_raw[, basis_names])
rownames(tpe_basis) <- tpe_raw$State
tpe_weights <- tpe_raw$Weight

stopifnot(identical(colnames(coef_mean), colnames(tpe_basis)))

## ----preflight----------------------------------------------------------------
input_check <- br_check_inputs(
  candidates = candidate_genotypes,
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
  candidates = candidate_genotypes,
  coef_mean = coef_mean,
  coef_vcov = coef_vcov,
  basis_names = basis_names,
  source = "Transparent genomic reaction-norm reference model",
  metadata = list(
    training_genotypes = n_train,
    candidate_genotypes = n_candidate,
    markers = n_markers,
    environments = n_env,
    candidate_phenotypes = 0
  )
)

tpe <- br_tpe(
  basis = tpe_basis,
  weights = tpe_weights,
  state_id = tpe_raw$State
)

## ----analysis-----------------------------------------------------------------
support <- br_support(tpe, training_basis)
omega <- br_omega(prediction, tpe, H = c(1, 3, 5, 10))

comparison <- br_compare(
  prediction, tpe, H = 5, k = 3,
  nsim = 1000, seed = 2026
)

advancement <- br_advance(
  prediction, tpe, H = 5, k = 3,
  nsim = 2000, method = "M3", seed = 2026
)

support$summary
omega$summary
summary(comparison)
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
template_dir <- file.path(tempdir(), "BreedReliabR_genomic_example")
dir.create(template_dir, showWarnings = FALSE, recursive = TRUE)

write.csv(
  data.frame(candidate = candidate_genotypes, coef_mean, check.names = FALSE),
  file.path(template_dir, "candidate_coefficients.csv"),
  row.names = FALSE
)
saveRDS(coef_vcov, file.path(template_dir, "candidate_joint_PEV.rds"))
write.csv(
  data.frame(state = tpe_raw$State, tpe_basis, weight = tpe_weights,
             check.names = FALSE),
  file.path(template_dir, "target_environments.csv"),
  row.names = FALSE
)
write.csv(
  data.frame(environment = rownames(training_basis), training_basis,
             check.names = FALSE),
  file.path(template_dir, "training_environments.csv"),
  row.names = FALSE
)

