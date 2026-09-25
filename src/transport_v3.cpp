#include <RcppArmadillo.h>
#include <string>

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]]

namespace {

struct ObjectiveV3 {
  double scientific;
  double optimization;
  double spatial;
  double chronology;
  double reference_selection;
  double source_selection;
  double smoothing;
  double correspondence_information;
  arma::vec reference_selected_unit;
  arma::vec source_selected_unit;
  arma::mat gradient;
};

double js_divergence(const arma::vec& observed, const arma::vec& target) {
  arma::vec midpoint = (observed + target) / 2.0;
  double result = 0.0;
  for (arma::uword i = 0; i < observed.n_elem; ++i) {
    if (observed[i] > 0.0) {
      result += 0.5 * observed[i] * std::log(observed[i] / midpoint[i]);
    }
    if (target[i] > 0.0) {
      result += 0.5 * target[i] * std::log(target[i] / midpoint[i]);
    }
  }
  return result;
}

arma::vec js_gradient(const arma::vec& observed, const arma::vec& target) {
  arma::vec midpoint = (observed + target) / 2.0;
  arma::vec result(observed.n_elem);
  const double floor = std::numeric_limits<double>::min();
  for (arma::uword i = 0; i < observed.n_elem; ++i) {
    result[i] = 0.5 * std::log(
      std::max(observed[i], floor) / std::max(midpoint[i], floor)
    );
  }
  return result;
}

double matrix_median(const arma::mat& values) {
  arma::vec sorted = arma::sort(arma::vectorise(values));
  arma::uword n = sorted.n_elem;
  if (n % 2 == 1) return sorted[n / 2];
  return (sorted[n / 2 - 1] + sorted[n / 2]) / 2.0;
}

ObjectiveV3 objective_v3(
    const arma::mat& coupling,
    const arma::vec& reference_mass,
    const arma::vec& source_mass,
    const arma::mat& reference_relation,
    const arma::mat& source_relation,
    const arma::mat& spatial_cost,
    double temporal_weight,
    double selection_weight,
    double entropy,
    bool need_gradient,
    bool revised = false) {
  const double coverage = arma::accu(coupling);
  arma::mat correspondence = coupling / coverage;
  arma::vec reference_selected = arma::sum(correspondence, 1);
  arma::vec source_selected = arma::sum(correspondence, 0).t();

  double spatial = arma::accu(correspondence % spatial_cost);
  double reference_edge_mass = arma::as_scalar(
    reference_selected.t() * reference_relation * reference_selected
  );
  double source_edge_mass = arma::as_scalar(
    source_selected.t() * source_relation * source_selected
  );
  arma::mat forward = reference_relation * correspondence * source_relation.t();
  double agreement = arma::accu(correspondence % forward);
  double denominator = reference_edge_mass + source_edge_mass;
  double chronology = denominator <= 1e-15 ? 0.0 :
    1.0 - 2.0 * agreement / denominator;
  chronology = std::min(1.0, std::max(0.0, chronology));

  double reference_selection = js_divergence(reference_selected, reference_mass);
  double source_selection = js_divergence(source_selected, source_mass);
  double mutual_information = 0.0;
  const double floor = std::numeric_limits<double>::min();
  for (arma::uword i = 0; i < correspondence.n_rows; ++i) {
    for (arma::uword j = 0; j < correspondence.n_cols; ++j) {
      double value = correspondence(i, j);
      if (value > 0.0) {
        double independent = reference_selected[i] * source_selected[j];
        double ratio = value / independent;
        if (revised && (!(independent > 0.0) || !std::isfinite(ratio))) {
          // Revision 2026.10: the product of two tiny marginals can underflow
          // although the cell itself is positive. Evaluate the log ratio in
          // the log domain instead of returning an infinite objective.
          mutual_information += value * (
            std::log(value) - std::log(reference_selected[i]) -
              std::log(source_selected[j])
          );
        } else {
          mutual_information += value * std::log(ratio);
        }
      }
    }
  }
  double scientific = coverage * (
    spatial + temporal_weight * chronology +
      selection_weight * (reference_selection + source_selection)
  );
  double smoothing = entropy * coverage * mutual_information;

  ObjectiveV3 result;
  result.scientific = scientific;
  result.optimization = scientific + smoothing;
  result.spatial = spatial;
  result.chronology = chronology;
  result.reference_selection = reference_selection;
  result.source_selection = source_selection;
  result.smoothing = smoothing;
  result.correspondence_information = mutual_information;
  result.reference_selected_unit = reference_selected;
  result.source_selected_unit = source_selected;
  if (!need_gradient) return result;

  arma::mat gradient = spatial_cost;
  if (denominator > 1e-15) {
    arma::mat agreement_gradient = forward +
      reference_relation.t() * correspondence * source_relation;
    arma::vec reference_edge_gradient =
      (reference_relation + reference_relation.t()) * reference_selected;
    arma::vec source_edge_gradient =
      (source_relation + source_relation.t()) * source_selected;
    arma::mat denominator_gradient =
      arma::repmat(reference_edge_gradient, 1, correspondence.n_cols) +
      arma::repmat(source_edge_gradient.t(), correspondence.n_rows, 1);
    arma::mat chronology_gradient = -2.0 * (
      agreement_gradient * denominator - agreement * denominator_gradient
    ) / (denominator * denominator);
    gradient += temporal_weight * chronology_gradient;
  }
  arma::vec reference_js = js_gradient(reference_selected, reference_mass);
  arma::vec source_js = js_gradient(source_selected, source_mass);
  gradient += selection_weight * (
    arma::repmat(reference_js, 1, correspondence.n_cols) +
    arma::repmat(source_js.t(), correspondence.n_rows, 1)
  );
  arma::mat mi_gradient(correspondence.n_rows, correspondence.n_cols);
  for (arma::uword i = 0; i < correspondence.n_rows; ++i) {
    for (arma::uword j = 0; j < correspondence.n_cols; ++j) {
      mi_gradient(i, j) = std::log(std::max(correspondence(i, j), floor)) -
        std::log(std::max(reference_selected[i], floor)) -
        std::log(std::max(source_selected[j], floor));
    }
  }
  gradient += entropy * mi_gradient;
  result.gradient = gradient;
  return result;
}

bool project_standard(
    const arma::mat& kernel,
    const arma::vec& row_target,
    const arma::vec& column_target,
    int maxit,
    double tolerance,
    arma::mat& plan,
    double& error,
    int& iterations) {
  arma::mat scaled = kernel;
  scaled(scaled.n_rows - 1, scaled.n_cols - 1) = 0.0;
  double scale = scaled.max();
  if (!std::isfinite(scale) || scale <= 0.0) return false;
  scaled /= scale;
  arma::uvec positive = arma::find(scaled > 0.0);
  if (positive.n_elem == 0 || !scaled.elem(positive).is_finite() ||
      arma::any(scaled.elem(positive) <= 0.0)) return false;

  arma::vec u(scaled.n_rows, arma::fill::ones);
  arma::vec v(scaled.n_cols, arma::fill::ones);
  error = std::numeric_limits<double>::infinity();
  bool converged = false;
  for (int iteration = 1; iteration <= maxit; ++iteration) {
    arma::vec kernel_v = scaled * v;
    if (!kernel_v.is_finite() || arma::any(kernel_v <= 0.0)) return false;
    u = row_target / kernel_v;
    arma::vec kernel_u = scaled.t() * u;
    if (!kernel_u.is_finite() || arma::any(kernel_u <= 0.0)) return false;
    v = column_target / kernel_u;
    if (!u.is_finite() || !v.is_finite()) return false;
    if (iteration == 1 || iteration % 10 == 0 || iteration == maxit) {
      arma::vec row_current = u % (scaled * v);
      arma::vec column_current = v % (scaled.t() * u);
      error = std::max(
        arma::abs(row_current - row_target).max(),
        arma::abs(column_current - column_target).max()
      );
      if (std::isfinite(error) && error <= tolerance) {
        converged = true;
        iterations = iteration;
        break;
      }
    }
    iterations = iteration;
  }
  plan = scaled % (u * v.t());
  plan(plan.n_rows - 1, plan.n_cols - 1) = 0.0;
  return converged;
}

// Damped Newton ascent on the dual of the masked KL projection (revision
// 2026.10 finisher after standard Sinkhorn exhausts its iterations). Mirrors
// masked_newton_projection() in R/gaze_weave_transport_projection.R.
bool project_newton(
    const arma::mat& kernel,
    const arma::vec& row_target,
    const arma::vec& column_target,
    double tolerance,
    arma::mat& plan,
    double& error,
    int& iterations,
    int maxit = 100) {
  const arma::uword m = kernel.n_rows;
  const arma::uword n = kernel.n_cols;
  const double negative_infinity = -std::numeric_limits<double>::infinity();
  error = std::numeric_limits<double>::infinity();
  // Work with the log kernel so that entries spanning hundreds of orders of
  // magnitude neither overflow nor underflow.
  arma::mat log_kernel(m, n);
  double top = negative_infinity;
  for (arma::uword i = 0; i < m; ++i) {
    for (arma::uword j = 0; j < n; ++j) {
      const double value = kernel(i, j);
      if (i == m - 1 && j == n - 1) {
        log_kernel(i, j) = negative_infinity;
        continue;
      }
      if (std::isnan(value) || value < 0.0 || !std::isfinite(value)) {
        return false;
      }
      log_kernel(i, j) = value > 0.0 ? std::log(value) : negative_infinity;
      if (log_kernel(i, j) > top) top = log_kernel(i, j);
    }
  }
  if (!std::isfinite(top)) return false;
  log_kernel -= top;
  auto log_sum_exp = [&](const arma::vec& values) {
    const double largest = values.max();
    if (!std::isfinite(largest)) return largest;
    return largest + std::log(arma::accu(arma::exp(values - largest)));
  };
  arma::vec row_potential(m, arma::fill::zeros);
  arma::vec column_potential(n, arma::fill::zeros);
  // Three log-domain Sinkhorn sweeps place the potentials near the solution.
  for (int pass = 0; pass < 3; ++pass) {
    for (arma::uword i = 0; i < m; ++i) {
      arma::vec values = log_kernel.row(i).t() + column_potential;
      row_potential[i] = std::log(row_target[i]) - log_sum_exp(values);
    }
    for (arma::uword j = 0; j < n; ++j) {
      arma::vec values = log_kernel.col(j) + row_potential;
      column_potential[j] = std::log(column_target[j]) - log_sum_exp(values);
    }
  }
  if (!row_potential.is_finite() || !column_potential.is_finite()) {
    return false;
  }
  auto plan_at = [&](const arma::vec& rows, const arma::vec& columns) {
    arma::mat potential = log_kernel + arma::repmat(rows, 1, n) +
      arma::repmat(columns.t(), m, 1);
    return arma::mat(arma::exp(potential));
  };
  auto dual_at = [&](const arma::vec& rows, const arma::vec& columns,
                     const arma::mat& current) {
    return arma::dot(row_target, rows) + arma::dot(column_target, columns) -
      arma::accu(current);
  };
  arma::mat current = plan_at(row_potential, column_potential);
  const arma::uword free = m + n - 1;
  int iteration = 0;
  for (iteration = 1; iteration <= maxit; ++iteration) {
    arma::vec row_sums = arma::sum(current, 1);
    arma::vec column_sums = arma::sum(current, 0).t();
    arma::vec gradient = arma::join_cols(
      row_target - row_sums, column_target - column_sums
    );
    error = arma::abs(gradient).max();
    if (!std::isfinite(error)) return false;
    if (error <= tolerance) break;
    arma::mat hessian(m + n, m + n, arma::fill::zeros);
    hessian.submat(0, 0, m - 1, m - 1) = arma::diagmat(row_sums);
    hessian.submat(0, m, m - 1, m + n - 1) = current;
    hessian.submat(m, 0, m + n - 1, m - 1) = current.t();
    hessian.submat(m, m, m + n - 1, m + n - 1) = arma::diagmat(column_sums);
    arma::mat reduced = hessian.submat(0, 0, free - 1, free - 1);
    const double ridge = 1e-14 * reduced.diag().max();
    reduced.diag() += ridge;
    arma::vec reduced_direction;
    bool solved = arma::solve(
      reduced_direction, reduced, gradient.subvec(0, free - 1),
      arma::solve_opts::no_approx
    );
    if (!solved || !reduced_direction.is_finite()) return false;
    arma::vec direction = arma::join_cols(reduced_direction, arma::vec({0.0}));
    arma::vec row_direction = direction.subvec(0, m - 1);
    arma::vec column_direction = direction.subvec(m, m + n - 1);
    const double current_dual = dual_at(
      row_potential, column_potential, current
    );
    const double slope = arma::dot(gradient, direction);
    double step = 1.0;
    bool accepted = false;
    arma::vec trial_rows;
    arma::vec trial_columns;
    arma::mat trial_plan;
    while (step >= 1e-12) {
      trial_rows = row_potential + step * row_direction;
      trial_columns = column_potential + step * column_direction;
      trial_plan = plan_at(trial_rows, trial_columns);
      const double trial_dual = dual_at(trial_rows, trial_columns, trial_plan);
      // Armijo ascent on the dual, or (once the dual increase falls below
      // rounding) a decrease of the marginal error.
      const double trial_error = std::max(
        arma::abs(row_target - arma::sum(trial_plan, 1)).max(),
        arma::abs(column_target - arma::sum(trial_plan, 0).t()).max()
      );
      if (std::isfinite(trial_dual) &&
          (trial_dual >= current_dual + 1e-4 * step * slope ||
           trial_error < error)) {
        accepted = true;
        break;
      }
      step /= 2.0;
    }
    if (!accepted) break;
    row_potential = trial_rows;
    column_potential = trial_columns;
    current = trial_plan;
  }
  iterations = std::min(iteration, maxit);
  error = std::max(
    arma::abs(row_target - arma::sum(current, 1)).max(),
    arma::abs(column_target - arma::sum(current, 0).t()).max()
  );
  current(m - 1, n - 1) = 0.0;
  plan = current;
  return std::isfinite(error) && error <= tolerance;
}

// Revision 2026.10 projection: standard Sinkhorn, then the Newton finisher.
bool project_revised(
    const arma::mat& kernel,
    const arma::vec& row_target,
    const arma::vec& column_target,
    int maxit,
    double tolerance,
    arma::mat& plan,
    double& error,
    int& iterations,
    bool& used_newton) {
  used_newton = false;
  if (project_standard(kernel, row_target, column_target, maxit, tolerance,
                       plan, error, iterations)) {
    return true;
  }
  used_newton = true;
  return project_newton(
    kernel, row_target, column_target, tolerance, plan, error, iterations
  );
}

Rcpp::List named_objective(const ObjectiveV3& objective, double coverage) {
  Rcpp::NumericVector components = Rcpp::NumericVector::create(
    Rcpp::Named("spatial") = coverage * objective.spatial,
    Rcpp::Named("chronology") = coverage * objective.chronology,
    Rcpp::Named("reference_selection") =
      coverage * objective.reference_selection,
    Rcpp::Named("source_selection") = coverage * objective.source_selection,
    Rcpp::Named("warp_penalty") = 0.0,
    Rcpp::Named("correspondence_smoothing") = objective.smoothing
  );
  Rcpp::NumericVector conditional = Rcpp::NumericVector::create(
    Rcpp::Named("spatial") = objective.spatial,
    Rcpp::Named("chronology") = objective.chronology,
    Rcpp::Named("reference_selection") = objective.reference_selection,
    Rcpp::Named("source_selection") = objective.source_selection,
    Rcpp::Named("correspondence_information") =
      objective.correspondence_information
  );
  return Rcpp::List::create(
    Rcpp::Named("scientific") = objective.scientific,
    Rcpp::Named("optimization") = objective.optimization,
    Rcpp::Named("coverage") = coverage,
    Rcpp::Named("selected_reference_mass") =
      coverage * objective.reference_selected_unit,
    Rcpp::Named("selected_source_mass") =
      coverage * objective.source_selected_unit,
    Rcpp::Named("components") = components,
    Rcpp::Named("conditional") = conditional,
    Rcpp::Named("zero_coverage") = false
  );
}

} // namespace



namespace {

struct NodeFitV3 {
  std::string start;
  arma::mat augmented;
  Rcpp::List history;
  bool all_converged;
  bool numerical_ok;
  std::string status;
  double final_change;
  double final_objective_change;
  double final_step;
  double final_projection_error;
  ObjectiveV3 final;
};

} // namespace

// Solver revisions: 0 reproduces the frozen 2026.08 solver exactly; 1 is the
// corrected 2026.10 solver. The R reference oracle
// (transport_v3_reference_stage_revised() and
// solve_transport_v3_mass_reference()) implements the same 2026.10 rules:
//
// * Projection: standard Sinkhorn, then a damped dual Newton finisher when
//   Sinkhorn exhausts its iterations.
// * Convergence: the first projected trial of an iteration taken at an
//   unbacktracked step (at least step_size / 8) has projected update residual
//   max|P_trial - P| / step <= tolerance and a converged projection.
// * Noise-limited iterations: both revisions accept a trial that rises by at
//   most 1e-12 relative, so a line search "fails" only when every trial rose
//   by more than that, which at 1e-8 projection accuracy is noise-dominated.
//   A failed line search therefore does not establish stationarity. When no
//   trial is accepted, or the accepted decrease is below
//   100 * projection_tolerance * max(1, |f|), the stage evaluates the
//   mirror-descent fixed-point residual at the fixed step step_size with a
//   projection 100 times tighter, and checks the plan's feasibility. A
//   residual <= tolerance converges the stage ("stationary_residual"). If no
//   unbacktracked progress is possible, or the residual is below what the
//   projection can resolve (10 * projection_tolerance / step_size), the stage
//   ends "stalled_projection_limited": not converged, not a failure.
// * Backtracking continues past 21 trials while the clamped mirror exponent
//   is still large.
// * Starts: every coverage node is solved from the independent start and,
//   when available, from the adjacent-coverage continuation start. The
//   converged fit with the lower regularized objective is kept (ties keep
//   the independent start), so the result does not depend on how far the
//   previous node was optimized.
// [[Rcpp::export]]
Rcpp::List transport_v3_profile_native_cpp(
    const arma::vec& reference_mass,
    const arma::vec& source_mass,
    const arma::mat& reference_relation,
    const arma::mat& source_relation,
    const arma::mat& spatial_cost,
    const arma::vec& coverage_values,
    const arma::vec& entropy_schedule,
    double temporal_weight,
    double selection_weight,
    int maxit,
    double step_size,
    double tolerance,
    int projection_maxit,
    double projection_tolerance,
    int revision = 0) {
  const bool revised = revision >= 1;
  const arma::uword nr = reference_mass.n_elem;
  const arma::uword ns = source_mass.n_elem;
  Rcpp::List fits(coverage_values.n_elem);
  arma::mat previous_augmented;
  bool have_previous_coverage = false;

  for (arma::uword coverage_index = 0;
       coverage_index < coverage_values.n_elem; ++coverage_index) {
    double coverage = coverage_values[coverage_index];
    arma::mat independent(nr + 1, ns + 1, arma::fill::zeros);
    independent.submat(0, 0, nr - 1, ns - 1) =
      coverage * (reference_mass * source_mass.t());
    independent.submat(0, ns, nr - 1, ns) =
      (1.0 - coverage) * reference_mass;
    independent.submat(nr, 0, nr, ns - 1) =
      (1.0 - coverage) * source_mass.t();
    arma::vec row_target = arma::join_cols(
      reference_mass, arma::vec({1.0 - coverage})
    );
    arma::vec column_target = arma::join_cols(
      source_mass, arma::vec({1.0 - coverage})
    );
    bool have_continuation = false;
    arma::mat continuation;
    if (have_previous_coverage) {
      arma::mat warm_plan;
      double warm_error = std::numeric_limits<double>::infinity();
      int warm_iterations = 0;
      bool warm_newton = false;
      bool warm_ok = revised ?
        project_revised(
          previous_augmented, row_target, column_target,
          projection_maxit, projection_tolerance,
          warm_plan, warm_error, warm_iterations, warm_newton
        ) :
        project_standard(
          previous_augmented,
          row_target,
          column_target,
          projection_maxit,
          projection_tolerance,
          warm_plan,
          warm_error,
          warm_iterations
        );
      if (warm_ok) {
        have_continuation = true;
        continuation = warm_plan;
      }
    }

    auto solve_node = [&](const std::string& start_name,
                          arma::mat augmented) -> NodeFitV3 {
      bool all_converged = true;
      bool numerical_ok = true;
      bool any_stalled = false;
      bool any_maxit = false;
      double final_change = std::numeric_limits<double>::infinity();
      double final_objective_change = std::numeric_limits<double>::infinity();
      double final_step = std::numeric_limits<double>::quiet_NaN();
      double final_projection_error = std::numeric_limits<double>::infinity();
      Rcpp::List history(entropy_schedule.n_elem);

      for (arma::uword stage = 0; stage < entropy_schedule.n_elem; ++stage) {
        double entropy = entropy_schedule[stage];
        double step = step_size;
        bool stage_converged = false;
        std::string termination = "maxit";
        bool stalled = false;
        bool used_newton = false;
        int iteration = 0;
        int projection_iterations = 0;
        for (iteration = 1; iteration <= maxit; ++iteration) {
          arma::mat coupling = augmented.submat(0, 0, nr - 1, ns - 1);
          ObjectiveV3 current = objective_v3(
            coupling, reference_mass, source_mass,
            reference_relation, source_relation, spatial_cost,
            temporal_weight, selection_weight, entropy, true, revised
          );
          arma::mat gradient =
            current.gradient - matrix_median(current.gradient);
          const double gradient_scale = arma::abs(gradient).max();
          bool accepted = false;
          double trial_step = step;
          arma::mat proposal = augmented;
          ObjectiveV3 proposal_objective = current;
          // Revision 2026.10 bookkeeping.
          bool have_certificate = false;
          double certificate_change = std::numeric_limits<double>::infinity();
          double certificate_step = std::numeric_limits<double>::quiet_NaN();
          double certificate_error = std::numeric_limits<double>::infinity();
          double best_increase = std::numeric_limits<double>::infinity();
          for (int backtrack = 0; ; ++backtrack) {
            if (!revised && backtrack > 20) break;
            if (revised && backtrack > 20 &&
                (trial_step * gradient_scale <= 1e-4 || backtrack > 80)) {
              break;
            }
            arma::mat kernel = augmented;
            arma::mat exponent =
              arma::clamp(-trial_step * gradient, -50.0, 50.0);
            kernel.submat(0, 0, nr - 1, ns - 1) =
              arma::clamp(coupling, 1e-300,
                          std::numeric_limits<double>::max()) %
              arma::exp(exponent);
            double projection_error = std::numeric_limits<double>::infinity();
            int projection_iteration = 0;
            bool trial_newton = false;
            bool projected = revised ?
              project_revised(
                kernel, row_target, column_target,
                projection_maxit, projection_tolerance,
                proposal, projection_error, projection_iteration, trial_newton
              ) :
              project_standard(
                kernel, row_target, column_target,
                projection_maxit, projection_tolerance,
                proposal, projection_error, projection_iteration
              );
            used_newton = trial_newton;
            final_projection_error = projection_error;
            projection_iterations = projection_iteration;
            if (!projected) {
              trial_step /= 2.0;
              continue;
            }
            arma::mat proposal_coupling =
              proposal.submat(0, 0, nr - 1, ns - 1);
            proposal_objective = objective_v3(
              proposal_coupling, reference_mass, source_mass,
              reference_relation, source_relation, spatial_cost,
              temporal_weight, selection_weight, entropy, false, revised
            );
            if (revised) {
              double trial_change = arma::abs(proposal - augmented).max();
              if (!have_certificate && trial_step >= step_size / 8.0) {
                have_certificate = true;
                certificate_change = trial_change;
                certificate_step = trial_step;
                certificate_error = projection_error;
              }
              double increase = proposal_objective.optimization -
                current.optimization;
              if (std::isfinite(increase) && increase < best_increase) {
                best_increase = increase;
              }
            }
            if (std::isfinite(proposal_objective.optimization) &&
                proposal_objective.optimization <= current.optimization +
                  1e-12 * std::max(1.0, std::abs(current.optimization))) {
              accepted = true;
              break;
            }
            trial_step /= 2.0;
          }
          const bool certified = revised && have_certificate &&
            certificate_change / certificate_step <= tolerance &&
            certificate_error <= projection_tolerance;
          if (revised) {
            const double scale = std::max(1.0, std::abs(current.optimization));
            // Both revisions accept a trial that rises by at most 1e-12
            // relative. At 1e-8 projection accuracy a rise or fall below
            // 100 * projection_tolerance * scale is noise: it neither shows
            // descent nor establishes stationarity. Stationarity is decided
            // only by the fixed-step residual check below.
            const bool noise_limited = !accepted || (
              !certified &&
                current.optimization - proposal_objective.optimization <=
                  100.0 * projection_tolerance * scale
            );
            if (noise_limited) {
              // Mirror-descent fixed-point residual at the fixed reference
              // step step_size, with a projection 100 times tighter than the
              // solver's, plus a feasibility check of the current plan.
              // Certified: converged. Otherwise the stage stalls when no
              // unbacktracked progress is possible or when the residual is
              // already below what the projection accuracy can resolve
              // (10 * projection_tolerance / step_size); else it continues.
              arma::mat kernel = augmented;
              arma::mat exponent =
                arma::clamp(-step_size * gradient, -50.0, 50.0);
              kernel.submat(0, 0, nr - 1, ns - 1) =
                arma::clamp(coupling, 1e-300,
                            std::numeric_limits<double>::max()) %
                arma::exp(exponent);
              arma::mat reference_plan;
              double reference_error = std::numeric_limits<double>::infinity();
              int reference_iterations = 0;
              bool reference_newton = false;
              const bool reference_projected = project_revised(
                kernel, row_target, column_target,
                projection_maxit, projection_tolerance * 1e-2,
                reference_plan, reference_error, reference_iterations,
                reference_newton
              );
              const double feasibility = std::max(
                arma::abs(arma::sum(augmented, 1) - row_target).max(),
                arma::abs(arma::sum(augmented, 0).t() - column_target).max()
              );
              const double reference_change = reference_projected ?
                arma::abs(reference_plan - augmented).max() :
                std::numeric_limits<double>::infinity();
              const double residual = reference_change / step_size;
              const bool feasible = feasibility <= projection_tolerance;
              const bool resolvable =
                residual > 10.0 * projection_tolerance / step_size;
              const bool progressing = accepted &&
                trial_step >= step_size / 8.0;
              if (reference_projected && feasible && residual <= tolerance) {
                final_change = reference_change;
                final_step = step_size;
                final_projection_error = reference_error;
                final_objective_change = std::abs(best_increase) / scale;
                stage_converged = true;
                termination = "stationary_residual";
                break;
              }
              if (!progressing || !reference_projected || !feasible ||
                  !resolvable) {
                final_change = reference_change;
                final_step = step_size;
                final_projection_error = reference_error;
                final_objective_change = std::abs(best_increase) / scale;
                stalled = true;
                termination = "stalled_projection_limited";
                break;
              }
            }
          } else if (!accepted) {
            numerical_ok = false;
            break;
          }
          final_change = arma::abs(proposal - augmented).max();
          final_objective_change = std::abs(
            proposal_objective.optimization - current.optimization
          ) / std::max(1.0, std::abs(current.optimization));
          final_step = trial_step;
          augmented = proposal;
          step = std::min(trial_step * 1.1, step_size);
          if (revised) {
            if (certified) {
              stage_converged = true;
              termination = "residual";
              if (certificate_step != trial_step) {
                // Report the residual that certified convergence.
                final_change = certificate_change;
                final_step = certificate_step;
              }
              break;
            }
          } else if (final_objective_change <= tolerance &&
                     final_projection_error <= projection_tolerance) {
            stage_converged = true;
            break;
          }
        }
        all_converged = all_converged && stage_converged;
        any_stalled = any_stalled || stalled;
        any_maxit = any_maxit || (!stage_converged && !stalled);
        Rcpp::List entry = Rcpp::List::create(
          Rcpp::Named("entropy") = entropy,
          Rcpp::Named("iterations") = iteration,
          Rcpp::Named("converged") = stage_converged,
          Rcpp::Named("coupling_change") = final_change,
          Rcpp::Named("relative_objective_change") = final_objective_change,
          Rcpp::Named("projected_update_residual") = final_change /
            std::max(final_step, std::numeric_limits<double>::epsilon()),
          Rcpp::Named("projection_error") = final_projection_error,
          Rcpp::Named("projection_converged") =
            final_projection_error <= projection_tolerance,
          Rcpp::Named("projection_method") =
            used_newton ? "native_newton" : "native_standard",
          Rcpp::Named("projection_iterations") = projection_iterations,
          Rcpp::Named("projection_fallback") = false
        );
        if (revised) {
          entry.push_back(termination, "termination");
          entry.push_back(final_step, "step");
          entry.push_back(start_name, "start");
        }
        history[stage] = entry;
      }
      arma::mat coupling = augmented.submat(0, 0, nr - 1, ns - 1);
      NodeFitV3 fit;
      fit.start = start_name;
      fit.augmented = augmented;
      fit.history = history;
      fit.all_converged = all_converged;
      fit.status = !numerical_ok ? "numerical_failure" :
        all_converged ? "converged" :
        any_maxit ? "not_converged" : "stalled_projection_limited";
      (void) any_stalled;
      fit.numerical_ok = numerical_ok;
      fit.final_change = final_change;
      fit.final_objective_change = final_objective_change;
      fit.final_step = final_step;
      fit.final_projection_error = final_projection_error;
      fit.final = objective_v3(
        coupling, reference_mass, source_mass,
        reference_relation, source_relation, spatial_cost,
        temporal_weight, selection_weight,
        entropy_schedule[entropy_schedule.n_elem - 1], false, revised
      );
      return fit;
    };

    std::vector<NodeFitV3> candidates;
    if (revised) {
      candidates.push_back(solve_node("independent_native", independent));
      if (have_continuation) {
        candidates.push_back(
          solve_node("adjacent_coverage_native", continuation)
        );
      }
    } else {
      candidates.push_back(solve_node(
        "independent_native",
        have_continuation ? continuation : independent
      ));
    }
    // Keep the fit with the lowest finite regularized objective (ties keep
    // the independent start); its convergence status is reported as is.
    arma::uword chosen = 0;
    for (arma::uword index = 1; index < candidates.size(); ++index) {
      const double value = candidates[index].final.optimization;
      const double incumbent = candidates[chosen].final.optimization;
      if (std::isfinite(value) &&
          (!std::isfinite(incumbent) || value < incumbent)) {
        chosen = index;
      }
    }
    const NodeFitV3& fit = candidates[chosen];
    Rcpp::NumericVector start_objectives(candidates.size());
    Rcpp::NumericVector start_scientific(candidates.size());
    Rcpp::CharacterVector start_names(candidates.size());
    for (arma::uword index = 0; index < candidates.size(); ++index) {
      start_objectives[index] = candidates[index].final.optimization;
      start_scientific[index] = candidates[index].final.scientific;
      start_names[index] = candidates[index].start;
    }
    start_objectives.attr("names") = start_names;
    start_scientific.attr("names") = start_names;
    const double start_spread = candidates.size() < 2 ? 0.0 :
      std::abs(candidates[0].final.optimization -
               candidates[1].final.optimization);
    const double start_scientific_spread = candidates.size() < 2 ? 0.0 :
      std::abs(candidates[0].final.scientific -
               candidates[1].final.scientific);
    arma::mat coupling = fit.augmented.submat(0, 0, nr - 1, ns - 1);
    Rcpp::List node = Rcpp::List::create(
      Rcpp::Named("start") = fit.start,
      Rcpp::Named("coverage") = coverage,
      Rcpp::Named("coupling") = coupling,
      Rcpp::Named("augmented") = fit.augmented,
      Rcpp::Named("objective") = named_objective(fit.final, coverage),
      Rcpp::Named("history") = fit.history,
      Rcpp::Named("converged") = fit.all_converged && fit.numerical_ok,
      Rcpp::Named("projection_converged") =
        fit.final_projection_error <= projection_tolerance,
      Rcpp::Named("stages_converged") = fit.all_converged,
      Rcpp::Named("native_ok") = fit.numerical_ok,
      Rcpp::Named("final_change") = fit.final_change,
      Rcpp::Named("final_objective_change") = fit.final_objective_change,
      Rcpp::Named("projected_update_residual") = fit.final_change /
        std::max(fit.final_step, std::numeric_limits<double>::epsilon()),
      Rcpp::Named("start_objectives") = start_objectives,
      Rcpp::Named("start_scientific_objectives") = start_scientific,
      Rcpp::Named("start_spread") = start_spread,
      Rcpp::Named("start_scientific_spread") = start_scientific_spread
    );
    if (revised) node.push_back(fit.status, "status");
    fits[coverage_index] = node;
    previous_augmented = fit.augmented;
    have_previous_coverage = true;
  }
  return Rcpp::List::create(
    Rcpp::Named("fits") = fits,
    Rcpp::Named("native_available") = true,
    Rcpp::Named("backend") = "rcpparmadillo"
  );
}
