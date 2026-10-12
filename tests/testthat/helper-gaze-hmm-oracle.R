# Dense log-space forward-backward for a general HMM. It is a test oracle
# only: the package uses the scaled low-rank Replay recursion
# (gaze_replay_forward_backward()), which these tests check against this
# straightforward implementation.
gaze_hmm_forward_backward <- function(log_emission, initial, transition) {
  n_time <- nrow(log_emission)
  n_states <- ncol(log_emission)
  if (length(initial) != n_states || any(dim(transition) != c(n_states, n_states)) ||
      any(initial < 0) || any(transition < 0) ||
      abs(sum(initial) - 1) > 1e-10 ||
      any(abs(rowSums(transition) - 1) > 1e-10)) {
    stop("Initial and transition probabilities do not match the emissions.")
  }
  log_initial <- log(initial)
  log_transition <- log(transition)
  forward <- matrix(-Inf, n_time, n_states)
  forward[1, ] <- log_initial + log_emission[1, ]
  if (n_time > 1L) {
    for (time in 2:n_time) {
      for (state in seq_len(n_states)) {
        forward[time, state] <- log_emission[time, state] +
          gaze_log_sum_exp(forward[time - 1L, ] + log_transition[, state])
      }
    }
  }
  log_likelihood <- gaze_log_sum_exp(forward[n_time, ])

  backward <- matrix(0, n_time, n_states)
  if (n_time > 1L) {
    for (time in seq.int(n_time - 1L, 1L)) {
      for (state in seq_len(n_states)) {
        backward[time, state] <- gaze_log_sum_exp(
          log_transition[state, ] + log_emission[time + 1L, ] +
            backward[time + 1L, ]
        )
      }
    }
  }
  log_posterior <- forward + backward - log_likelihood
  posterior <- exp(log_posterior)
  posterior <- posterior / rowSums(posterior)

  transition_posterior <- array(0, dim = c(max(n_time - 1L, 0L), n_states, n_states))
  if (n_time > 1L) {
    for (time in seq_len(n_time - 1L)) {
      log_xi <- outer(
        forward[time, ],
        log_emission[time + 1L, ] + backward[time + 1L, ],
        FUN = "+"
      ) + log_transition - log_likelihood
      xi <- exp(log_xi)
      transition_posterior[time, , ] <- xi / sum(xi)
    }
  }

  list(
    log_likelihood = log_likelihood,
    forward = forward,
    backward = backward,
    posterior = posterior,
    transition_posterior = transition_posterior
  )
}
