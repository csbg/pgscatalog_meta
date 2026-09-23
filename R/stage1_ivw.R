# Stage-1 inverse-variance weighting on the logit AUC scale.

clamp01 <- function(a, eps = 1e-6) pmin(pmax(a, eps), 1 - eps)

auc_ci_to_logit <- function(auc, lo, hi, eps_se = 1e-8) {
  a   <- clamp01(as.numeric(auc))
  lo1 <- clamp01(as.numeric(lo))
  hi1 <- clamp01(as.numeric(hi))
  eta <- qlogis(a)
  se  <- (qlogis(hi1) - qlogis(lo1)) / (2 * 1.96)
  se  <- ifelse(is.finite(se) & se > 0, se, eps_se)
  tibble::tibble(eta = eta, se = se)
}

ivw_pool_logit <- function(eta, se) {
  ok <- is.finite(eta) & is.finite(se) & (se > 0)
  eta <- eta[ok]
  se <- se[ok]
  k <- length(eta)
  if (k == 0) return(tibble::tibble(eta = NA_real_, se = NA_real_, Q = NA_real_, I2 = NA_real_, k_eval = 0L))
  if (k == 1) return(tibble::tibble(eta = eta, se = se, Q = NA_real_, I2 = NA_real_, k_eval = 1L))

  w <- 1 / (se^2)
  w <- ifelse(is.finite(w), w, 0)
  sw <- sum(w)
  if (!is.finite(sw) || sw <= 0) {
    return(tibble::tibble(eta = NA_real_, se = NA_real_, Q = NA_real_, I2 = NA_real_, k_eval = k))
  }

  eta_hat <- sum(w * eta) / sw
  se_hat  <- sqrt(1 / sw)

  Q  <- sum(w * (eta - eta_hat)^2)
  df <- k - 1
  I2 <- if (isTRUE(Q > 0)) max(0, (Q - df) / Q) * 100 else 0
  assert_i2_percent(I2)

  tibble::tibble(eta = eta_hat, se = se_hat, Q = Q, I2 = I2, k_eval = k)
}
