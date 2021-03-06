% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/simulate_m0.R
\name{simulate_m0}
\alias{simulate_m0}
\title{Simulate encounter-capture-recapture data.}
\usage{
simulate_m0(
  n = 100,
  pr_encounter = 0.8,
  pr_cap = 0.8,
  n_timesteps = 3,
  n_aug = 100
)
}
\arguments{
\item{n}{The true population size.}

\item{pr_encounter}{The probability of encountering an individual on one survey.}

\item{pr_cap}{The probability of capture conditional on encounter.}

\item{n_timesteps}{The number of repeat surveys.}

\item{n_aug}{The number of capture histories that augment observed histories.}
}
\value{
A list of simulated values
}
\description{
Simulate encounter-capture-recapture data.
}
\examples{
simulate_m0(n = 20, pr_encounter = 0.5, pr_cap = 0.5, n_timesteps = 2, n_aug = 100)
}
