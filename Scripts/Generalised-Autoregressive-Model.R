library(nlme)
library(mnormt)

get.temporal.correlation.given.time.difference <- function(temporal.correlation, time.difference) {
  result <- temporal.correlation ^ abs(time.difference)
  return(result)
}

get.sd.shock.given.time.difference <- function(temporal.correlation, sd.shock, time.difference) {
  temporal.correlation.given.time.difference <- get.temporal.correlation.given.time.difference(temporal.correlation, time.difference)
  result <- sd.shock * sqrt((1 - temporal.correlation.given.time.difference ^ 2) / (1 - temporal.correlation ^ 2))
  return(result)
}

generate.next.value <- function(time.difference, previous.value, temporal.correlation, sd.shock) {
  temporal.correlation.given.time.difference <- get.temporal.correlation.given.time.difference(temporal.correlation, time.difference)
  sd.shock.given.time.difference <- get.sd.shock.given.time.difference(temporal.correlation, sd.shock, time.difference)
  result <- temporal.correlation.given.time.difference * previous.value + sd.shock.given.time.difference * rnorm(1)
  return(result)
}

get.sd.correlated.error <- function(temporal.correlation, sd.shock) {
  result <- get.sd.shock.given.time.difference(temporal.correlation, sd.shock, Inf)
  return(result)
}

generate.population.with.generalised.autoregressive.values <- function(n.participants, n.time.steps, time.step.size, sd.random.intercept, temporal.correlation, sd.shock, sd.uncorrelated.error) {
  result <- expand.grid(
    participant.id = seq(n.participants),
    wave = seq(n.time.steps)
  )
  result$time <- time.step.size * (result$wave - 1)
  result$random.intercept <- sd.random.intercept * rep(rnorm(n.participants), n.time.steps)
  sd.correlated.error <- get.sd.correlated.error(temporal.correlation, sd.shock)
  result$autoregressive <- sd.correlated.error * rnorm(n.participants * n.time.steps)
  for (i in seq(2, n.time.steps)) {
    previous.value <- result[result$wave == i - 1, "autoregressive"]
    result[result$wave == i, "autoregressive"] <- sapply(previous.value, generate.next.value, time.difference = time.step.size, temporal.correlation = temporal.correlation, sd.shock = sd.shock)
  }
  result$uncorrelated.error <- sd.uncorrelated.error * rnorm(n.participants * n.time.steps)
  result$value <- result$random.intercept + result$autoregressive + result$uncorrelated.error
  return(result)
}

get.time.constant <- function(temporal.correlation) {
  result <- -1 / log(temporal.correlation)
  return(result)
}

get.total.variance <- function(sd.random.intercept, temporal.correlation, sd.shock, sd.uncorrelated.error) {
  sd.correlated.error <- get.sd.correlated.error(temporal.correlation, sd.shock)
  result <- sd.random.intercept ^ 2 + sd.correlated.error ^ 2 + sd.uncorrelated.error ^ 2
  return(result)
}

get.covariance.given.time.difference <- function(time.difference, sd.random.intercept, temporal.correlation, sd.shock, sd.uncorrelated.error) {
  sd.correlated.error <- get.sd.correlated.error(temporal.correlation, sd.shock)
  temporal.correlation.given.time.difference <- get.temporal.correlation.given.time.difference(temporal.correlation, abs(time.difference))
  result <- sd.random.intercept ^ 2 + sd.correlated.error ^ 2 * temporal.correlation.given.time.difference + sd.uncorrelated.error ^ 2 * (time.difference == 0)
  return(result)
}

get.covariance.matrix <- function(times, sd.random.intercept, temporal.correlation, sd.shock, sd.uncorrelated.error) {
  time.differences <- outer(times, times, `-`)
  result <- get.covariance.given.time.difference(
    time.difference = time.differences,
    sd.random.intercept = sd.random.intercept,
    temporal.correlation = temporal.correlation,
    sd.shock = sd.shock,
    sd.uncorrelated.error = sd.uncorrelated.error
  )
  return(result)
}

generate.generalised.autoregressive.values <- function(n, times, sd.random.intercept, temporal.correlation, sd.shock, sd.uncorrelated.error) {
  varcov <- get.covariance.matrix(
    times = times,
    sd.random.intercept = sd.random.intercept,
    temporal.correlation = temporal.correlation,
    sd.shock = sd.shock,
    sd.uncorrelated.error = sd.uncorrelated.error
  )
  result <- rmnorm(
    n = n,
    varcov = varcov
  )
  return(result)
}

generate.population.with.values.simultaneously <- function(n.participants, n.time.steps, time.step.size, sd.random.intercept, temporal.correlation, sd.shock, sd.uncorrelated.error) {
  result <- expand.grid(
    participant.id = seq(n.participants),
    wave = seq(n.time.steps)
  )
  result$time <- time.step.size * (result$wave - 1)
  result$value <- c(
    generate.generalised.autoregressive.values(
      n = n.participants,
      times = seq(0, (n.time.steps - 1) * time.step.size, time.step.size),
      sd.random.intercept = sd.random.intercept,
      temporal.correlation = temporal.correlation,
      sd.shock = sd.shock,
      sd.uncorrelated.error = sd.uncorrelated.error
    )
  )
  return(result)
}

get.transition.probability <- function(start.lower, start.upper, end.lower, end.upper, time.difference, sd.random.intercept, temporal.correlation, sd.shock, sd.uncorrelated.error) {
  variance <- get.total.variance(
    sd.random.intercept = sd.random.intercept,
    temporal.correlation = temporal.correlation,
    sd.shock = sd.shock,
    sd.uncorrelated.error = sd.uncorrelated.error
  )
  covariance <- get.covariance.given.time.difference(
    time.difference = time.difference,
    sd.random.intercept = sd.random.intercept,
    temporal.correlation = temporal.correlation,
    sd.shock = sd.shock,
    sd.uncorrelated.error = sd.uncorrelated.error
  )
  result <- biv.nt.prob(
    df = Inf,
    lower = c(start.lower, end.lower),
    upper = c(start.upper, end.upper),
    mean = c(0, 0),
    S = matrix(c(variance, covariance, covariance, variance), 2, 2)
  ) / (pnorm(start.upper, 0, sqrt(variance)) - pnorm(start.lower, 0, sqrt(variance)))
  return(result)
}

fit.generalised.autoregressive.model <- function(data, fixed.formula = value ~ 0) {
  result <- lme(
    fixed = fixed.formula,
    random = list(participant.id = ~ 1),
    correlation = corExp(form = ~ time | participant.id, nugget = TRUE),
    data = data
  )
  return(result)
}

get.generalised.autoregressive.model.parameters <- function(fit) {
  sd.random.intercept <- as.numeric(VarCorr(fit)["(Intercept)", "StdDev"])
  sd.error <- as.numeric(VarCorr(fit)["Residual", "StdDev"])
  range <- unname(coef(fit$modelStruct$corStruct, unconstrained = FALSE))[1]
  nugget <- unname(coef(fit$modelStruct$corStruct, unconstrained = FALSE))[2]
  temporal.correlation <- exp(-1 / range)
  result <- c(
    fixed.effects(fit),
    "sd.random.intercept" = sd.random.intercept,
    "temporal.correlation" = temporal.correlation,
    "sd.shock" = sd.error * sqrt(1 - nugget) * sqrt(1 - temporal.correlation ^ 2),
    "sd.uncorrelated.error" = sd.error * sqrt(nugget)
  )
  return(result)
}
