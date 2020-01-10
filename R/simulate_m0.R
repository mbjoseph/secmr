#' Simulate encounter-capture-recapture data.
#'
#' @param n The true population size.
#' @param pr_encounter The probability of encountering an individual on one survey.
#' @param pr_cap The probability of capture conditional on encounter.
#' @param n_timesteps The number of repeat surveys.
#' @param n_aug The number of capture histories that augment observed histories.
#' @return A list of simulated values
#' @examples
#' simulate_m0(n = 20, pr_encounter = 0.5, pr_cap = 0.5, n_timesteps = 2, n_aug = 100)
#' @export
simulate_m0 <- function(n = 100, pr_encounter = 0.8, pr_cap = 0.8, n_timesteps = 3,
                        n_aug = 100) {
  e_full <- matrix(nrow = n, ncol = n_timesteps)
  for (t in 1:n_timesteps) {
    e_full[, t] <- stats::rbinom(n, size = 1, prob = pr_encounter)
  }

  y_full <- matrix(nrow = n, ncol = n_timesteps)
  for (t in 1:n_timesteps) {
    y_full[, t] <- stats::rbinom(n = n, size = 1, prob = e_full[, t] * pr_cap)
  }
  y_encountered_not_captured <- e_full & !y_full
  ever_captured <- apply(y_full, 1, max)
  C <- sum(ever_captured)
  y_obs <- y_full[ever_captured == 1, ]
  y_aug <- rbind(y_obs, matrix(0, nrow = n_aug, ncol = n_timesteps))

  n_encountered <- apply(e_full, 2, sum)
  n_encountered_not_captured <- apply(y_encountered_not_captured, 2, sum)
  stopifnot(all(n_encountered >= n_encountered_not_captured))
  stopifnot(all(n_encountered_not_captured >= 0))

  list(n = n,
       pr_encounter = pr_encounter,
       pr_cap = pr_cap,
       C = C,
       n_timesteps = n_timesteps,
       y_full = y_full,
       y_aug = y_aug,
       y_cmr = apply(y_aug, 1, sum),
       n_encountered = n_encountered,
       n_encountered_not_captured = n_encountered_not_captured)
}
