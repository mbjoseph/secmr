#' Run M0 model experiments.
#'
#' This runs the experiments where data sets are simulated, and a basic
#' capture-recapture model (M0) is compared against a model that includes
#' encounter data leading to failed captures.
#' @importFrom magrittr %>%
#' @importFrom rlang .data
#' @export
m0_experiments <- function() {
  set.seed(123)

  param_grid <- data.frame(pr_encounter = .5,
                           pr_cap = .2,
                           n_timesteps = c(3, 6, 9)) %>%
    tidyr::complete(.data$pr_encounter, .data$pr_cap, .data$n_timesteps) %>%
    dplyr::mutate(outfile = here::here(
      "analysis", "output",
      paste0('det_', .data$pr_encounter,
             '_cap_', .data$pr_cap,
             '_timesteps_', .data$n_timesteps, "_.rds")),
      exists = file.exists(.data$outfile),
      n_aug = dplyr::case_when(
        .data$n_timesteps == 3 ~ 1200,
        .data$n_timesteps == 6 ~ 300,
        .data$n_timesteps == 9 ~ 100
        )
      )

  for (i in 1:nrow(param_grid)) {
    fit_model(param_grid[i, ])
    print(i)
  }

  # For each model, read the data and generate a data frame for plotting
  post_dfs <- list()
  max_pts <- list()
  for (i in 1:nrow(param_grid)) {
    out <- readr::read_rds(param_grid$outfile[i])
    m0 <- out$m0
    me <- out$me
    post_dfs[[i]] <- tibble::tibble(
      N_m0 = m0$sims.list$N,
      p_m0 = m0$sims.list$p,
      N_m = me$sims.list$N,
      p_m = me$sims.list$p_encounter * me$sims.list$p_capture
      ) %>%
      dplyr::mutate(iter = 1:dplyr::n()) %>%
      tidyr::pivot_longer(tidyselect::contains("_")) %>%
      tidyr::separate(.data$name, into = c("par", "mod")) %>%
      tidyr::pivot_wider(names_from = .data$par, values_from = .data$value) %>%
      dplyr::mutate(Model = ifelse(.data$mod == "m", "Encounter-capture-recapture", "Capture-recapture")) %>%
      dplyr::arrange(.data$Model) %>%
      dplyr::mutate(n_timesteps = paste("K =", out$d$n_timesteps, "surveys"))
    max_pts[[i]] <- tibble::tibble(C = c(out$d$C, max(out$d$n_encountered)),
                           label = c("c[min]", "n[min]"),
                           n_timesteps = paste("K =", out$d$n_timesteps, "surveys"),
                           ypos = c(0.03, -0.02),
                           xadj = c(-5, 4))
  }

  post_dfs <- dplyr::bind_rows(post_dfs)
  max_pts <- dplyr::bind_rows(max_pts)

  joint_post_plot <- post_dfs %>%
    ggplot2::ggplot(ggplot2::aes(.data$N, .data$p, color = .data$Model)) +
    ggplot2::facet_wrap(~.data$n_timesteps) +
    ggplot2::geom_point(size = .1, alpha = .4) +
    ggplot2::theme_minimal() +
    ggplot2::geom_vline(ggplot2::aes(xintercept = .data$C),
                        linetype = "dashed", size = .2, data = max_pts) +
    ggplot2::xlab("Abundance") +
    ggplot2::ylab("Capture probability") +
    ggplot2::theme(legend.position = "bottom",
          panel.grid.minor = ggplot2::element_blank(),
          plot.margin = grid::unit(c(.1, 1, .1, .5), "cm")) +
    ggplot2::scale_color_manual(values = c("black", "dodgerblue"), "Model:") +
    ggplot2::guides(colour = ggplot2::guide_legend(override.aes = list(size = 2, alpha = 1))) +
    ggplot2::scale_x_log10(limits = c(5, max(post_dfs$N))) +
    ggplot2::geom_text(ggplot2::aes(x = .data$C + .data$xadj, y = ypos, label = .data$label),
                       parse = TRUE,
              data = max_pts, inherit.aes = FALSE, size = 3, nudge_x = .1) +
    ggplot2::theme(strip.background = ggplot2::element_blank(),
                   strip.text.x = ggplot2::element_blank())
  joint_post_plot

  n_dens <- post_dfs %>%
    ggplot2::ggplot(ggplot2::aes(.data$N, color = .data$Model)) +
    ggplot2::facet_wrap(~.data$n_timesteps) +
    ggplot2::stat_density(geom = "line", position = "identity", trim = TRUE) +
    ggplot2::theme_minimal() +
    ggplot2::geom_vline(ggplot2::aes(xintercept = .data$C),
                        linetype = "dashed", size = .2, data = max_pts) +
    ggplot2::xlab("") +
    ggplot2::ylab("Posterior density") +
    ggplot2::theme(legend.position = "none",
          axis.text.x = ggplot2::element_blank(),
          panel.grid.minor = ggplot2::element_blank(),
          plot.margin = grid::unit(c(.1, 1, .1, .5), "cm")) +
    ggplot2::scale_color_manual(values = c("black", "dodgerblue"), "Model:") +
    ggplot2::scale_fill_manual(values = c("black", "dodgerblue"), "Model:") +
    ggplot2::scale_x_log10(limits = c(5, max(post_dfs$N))) +
    ggplot2::geom_text(ggplot2::aes(x = .data$C + .data$xadj, y = ypos * 5, label = .data$label),
                       parse = TRUE, data = max_pts,
                       inherit.aes = FALSE, size = 3, nudge_x = .1)

  p <- patchwork::wrap_plots(n_dens, joint_post_plot, ncol = 1)
  p
  ggplot2::ggsave(here::here("analysis", "figures", "joint-post.pdf"),
                  plot = p, width = 6, height = 4)
}



fit_model <- function(df) {
  if (file.exists(df$outfile)) return(df$outfile)
  pr_encounter <- df$pr_encounter
  pr_cap <- df$pr_cap
  outfile <- df$outfile

  d <- simulate_m0(n = 50,
                   n_aug = df$n_aug,
                   pr_encounter = pr_encounter,
                   pr_cap = pr_cap,
                   n_timesteps = df$n_timesteps)

  m0_d <- list(M = nrow(d$y_aug),
               n_timesteps = d$n_timesteps,
               z = ifelse(apply(d$y_aug, 1, max), 1, NA),
               y_cmr = d$y_cmr)

  encounter_mat <- d$y_aug
  encounter_mat[d$y_aug == 0] <- NA
  encounter_df <- reshape2::melt(encounter_mat,
                                 varnames = c("indiv", "t")) %>%
    tibble::as_tibble() %>%
    dplyr::filter(.data$value == 1)
  z_obs <- ifelse(apply(d$y_aug, 1, max) == 1, 1, NA)
  n_min <- max(c(d$C, d$n_encountered))

  n_iter <- 300000

  #CMR
  system.time({
    m_par <- jagsUI::jags(data = m0_d,
                          model.file = here::here("analysis", "models", "m0.txt"),
                          n.chains = 6,
                          n.adapt = n_iter / 2, n.iter = n_iter,
                          n.burnin = n_iter / 10,
                          n.thin = ifelse(n_iter > 1000, n_iter / 1000, 1),
                          parameters.to.save = c('z', 'p', 'N', 'omega'),
                          parallel = TRUE,
                          verbose = TRUE)
  })
  # traceplot(m_par, parameters = "N", ask = FALSE)

  # integrated model
  y_cat <- array(dim = c(length(d$y_cmr), d$n_timesteps, 3))
  y_star_init <- array(dim = c(length(d$y_cmr), d$n_timesteps))
  for (i in 1:length(d$y_cmr)) {
    for (k in 1:d$n_timesteps) {
      y_cat[i, k, 3] <- d$y_aug[i, k]
      if (d$y_aug[i, k] == 1) {
        y_cat[i, k, 1:2] <- 0
        y_star_init[i, k] <- 3
      }
    }
  }
  ves_d <- list(n_timesteps = d$n_timesteps,
                M = length(d$y_cmr),
                y = d$y_aug,
                failures = d$n_encountered_not_captured,
                z = z_obs,
                y_star = y_star_init)

  z_init <- ifelse(is.na(z_obs), stats::rbinom(ves_d$M, 1, .5), NA)

  init_list <- function() list(z = z_init)

  system.time({
    m_ves <- jagsUI::jags(data = ves_d,
                          inits = init_list,
                          parameters.to.save = c('z', 'p_encounter',
                                                 'p_capture', 'N',
                                                 'n_missed', 'omega'),
                          model.file = here::here("analysis", "models", "m0_multinom.txt"),
                          n.chains = 6, n.iter = n_iter,
                          n.burnin = n_iter / 10, n.adapt = n_iter / 2,
                          n.thin = ifelse(n_iter > 1000, n_iter / 1000, 1),
                          parallel = TRUE)
  })
  # traceplot(m_ves, "N")
  # traceplot(m_ves, "n_missed")

  list(d = d,
       m0 = m_par,
       me = m_ves) %>%
    readr::write_rds(outfile)
}
