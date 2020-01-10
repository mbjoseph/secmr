#' Run model experiments.
#'
#' This runs the experiments where data sets are simulated, and a basic
#' capture-recapture model (M0) is compared against a model that includes
#' encounter data leading to failed captures.
#' @importFrom magrittr %>%
#' @importFrom rlang .data
#' @export
lower_bound_analysis <- function() {
  outfile <- here::here('analysis', 'figures', 'sim-verification.png')
  if (file.exists(outfile)) {
    return(NA)
  }

  N <- 100
  eta <- .8
  kappa <- .2
  K <- 5

  # Maxima from counts ------------------------------------------------------
  vals <- 0:N
  Fmax <- rep(NA, length(vals))
  dmax <- rep(NA, length(vals))
  Fmax[1] <- stats::pbinom(vals[1], size = N, prob = eta)^K

  for (i in 2:length(vals)) {
    Fmax[i] <- stats::pbinom(vals[i], size = N, prob = eta)^K
    dmax[i] <- Fmax[i] - Fmax[i - 1]
  }

  param_grid <- data.frame(
    K = c(3, 6, 9),
    eta = seq(0.01, 0.99, by = .01),
    kappa = seq(0.01, 0.99, by = .01),
    N = c(10, 50, 100)) %>%
    tidyr::complete(.data$K, .data$eta, .data$kappa, .data$N) %>%
    tibble::as_tibble() %>%
    dplyr::mutate(E_cmr = .data$N * (1 - (1 - .data$eta * .data$kappa)^.data$K)) %>%
    dplyr::rowwise() %>%
    dplyr::mutate(E_ves = expected_ves(.data$N, .data$eta, .data$K)) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(E_diff = .data$E_cmr - .data$E_ves,
                  N_label = paste("N =", .data$N))

  p <- param_grid %>%
    tidyr::pivot_longer(cols = dplyr::starts_with("E_")) %>%
    dplyr::filter(.data$name == "E_diff") %>%
    dplyr::mutate(Klab = paste0("K = ", .data$K, " surveys")) %>%
    ggplot2::ggplot(ggplot2::aes(.data$eta, .data$kappa,
                                 fill = .data$value)) +
    ggplot2::facet_grid(reorder(.data$N_label,
                                .data$N) ~ reorder(.data$Klab,
                                                   .data$K)) +
    ggplot2::geom_tile() +
    ggplot2::scale_fill_gradient2(expression(paste(E(c[min] - n[min])))) +
    ggplot2::stat_contour(ggplot2::aes(z = .data$value),
                 color = 'black', breaks = c(0)) +
    ggplot2::stat_contour(ggplot2::aes(z = .data$value),
                 color = 'black', size = .1, alpha = .5,
                 breaks = seq(-100, 100, by = 2)) +
    ggplot2::xlab(expression(paste("Encounter probability: ", eta))) +
    ggplot2::ylab(expression(paste("Capture probability: ", kappa))) +
    ggplot2::theme_minimal() +
    ggplot2::coord_equal() +
    ggplot2::theme(panel.grid = ggplot2::element_blank()) +
    ggplot2::scale_x_continuous(labels = c(0, .25, .5, .75, 1)) +
    ggplot2::scale_y_continuous(labels = c(0, .25, .5, .75, 1))
  p

  ggplot2::ggsave(here::here('analysis', 'figures', 'm0-maxima.png'),
                  plot = p, width = 6, height = 5)



  # Verify with simulation --------------------------------
  sim_grid <- tibble::as_tibble(
    expand.grid(K = c(3, 6, 9),
                eta = seq(0.01, 0.99, by = .01),
                kappa = seq(0.01, 0.99, by = .01),
                N = c(10, 50, 100),
                iter = 1:5)
    ) %>%
    dplyr::mutate(E_cmr = .data$N * (1 - (1 - .data$eta * .data$kappa)^.data$K)) %>%
    dplyr::rowwise() %>%
    dplyr::mutate(E_ves = expected_ves(.data$N, .data$eta, .data$K)) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(E_diff = .data$E_cmr - .data$E_ves)

  sim_grid$n_min <- NA
  sim_grid$c_min <- NA
  pb <- utils::txtProgressBar(max = nrow(sim_grid))
  for (i in 1:nrow(sim_grid)) {
    d <- simulate_m0(n = sim_grid$N[i], pr_encounter = sim_grid$eta[i],
                     pr_cap = sim_grid$kappa[i], n_timesteps = sim_grid$K[i],
                     n_aug = 1)
    sim_grid$n_min[i] <- max(d$n_encountered)
    sim_grid$c_min[i] <- d$C
    utils::setTxtProgressBar(pb, i)
  }
  close(pb)

  sim_grid <- sim_grid %>%
    dplyr::mutate(id = 1:dplyr::n())

  c_data <- sim_grid %>%
    dplyr::select(.data$id, .data$E_cmr, .data$c_min) %>%
    dplyr::rename(expected = .data$E_cmr, observed = .data$c_min) %>%
    dplyr::mutate(src = "Capture lower bounds")
  e_data <- sim_grid %>%
    dplyr::select(.data$id, .data$E_ves, .data$n_min) %>%
    dplyr::rename(expected = .data$E_ves,
                  observed = .data$n_min) %>%
    dplyr::mutate(src = "Encounter lower bounds")

  sim_df <- dplyr::full_join(c_data, e_data) %>%
    dplyr::left_join(dplyr::select(sim_grid,
                                   .data$iter, .data$id,
                                   .data$K, .data$eta,
                                   .data$kappa, .data$N))

  sim_verify <- sim_df %>%
    dplyr::group_by(.data$K, .data$eta, .data$kappa, .data$N, .data$src) %>%
    dplyr::summarize(expected = mean(.data$expected),
                     observed = mean(.data$observed)) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(facet = paste0(.data$src, ": ", paste("N =", .data$N)),
                  facet = forcats::fct_reorder(.data$facet, .data$N)) %>%
    ggplot2::ggplot(ggplot2::aes(.data$expected, .data$observed)) +
    ggplot2::geom_point(alpha = .01, size = .1) +
    ggplot2::facet_wrap(dplyr::vars(.data$facet),
                        scales = "free", ncol = 2) +
    ggplot2::xlab("Theoretical expectation") +
    ggplot2::ylab("Empirical mean") +
    ggplot2::theme_minimal() +
    ggplot2::theme(panel.grid.minor = ggplot2::element_blank())
  sim_verify
  ggplot2::ggsave(outfile,
                  plot = sim_verify, width = 5, height = 4.5)
}


# Compute expected maximum encounters over grid of values
expected_ves <- function(N, eta, K) {
  vals <- 0:N
  Fmax <- rep(NA, length(vals))
  dmax <- rep(NA, length(vals))
  Fmax[1] <- stats::pbinom(vals[1], size = N, prob = eta)^K

  for (i in 2:length(vals)) {
    Fmax[i] <- stats::pbinom(vals[i], size = N, prob = eta)^K
    dmax[i] <- Fmax[i] - Fmax[i - 1]
  }

  expected_max_ves <- sum(vals * dmax, na.rm = TRUE)
  expected_max_ves
}
