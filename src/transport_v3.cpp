#include <RcppArmadillo.h>

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
    bool need_gradient) {
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
        mutual_information += value * std::log(value / independent);
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
    double projection_tolerance) {
  const arma::uword nr = reference_mass.n_elem;
  const arma::uword ns = source_mass.n_elem;
  Rcpp::List fits(coverage_values.n_elem);
  arma::mat previous_augmented;
  bool have_previous_coverage = false;

  for (arma::uword coverage_index = 0;
       coverage_index < coverage_values.n_elem; ++coverage_index) {
    double coverage = coverage_values[coverage_index];
    arma::mat augmented(nr + 1, ns + 1, arma::fill::zeros);
    augmented.submat(0, 0, nr - 1, ns - 1) =
      coverage * (reference_mass * source_mass.t());
    augmented.submat(0, ns, nr - 1, ns) =
      (1.0 - coverage) * reference_mass;
    augmented.submat(nr, 0, nr, ns - 1) =
      (1.0 - coverage) * source_mass.t();
    arma::vec row_target = arma::join_cols(
      reference_mass, arma::vec({1.0 - coverage})
    );
    arma::vec column_target = arma::join_cols(
      source_mass, arma::vec({1.0 - coverage})
    );
    if (have_previous_coverage) {
      arma::mat warm_plan;
      double warm_error = std::numeric_limits<double>::infinity();
      int warm_iterations = 0;
      bool warm_ok = project_standard(
        previous_augmented,
        row_target,
        column_target,
        projection_maxit,
        projection_tolerance,
        warm_plan,
        warm_error,
        warm_iterations
      );
      if (warm_ok) augmented = warm_plan;
    }

    bool all_converged = true;
    bool numerical_ok = true;
    double final_change = std::numeric_limits<double>::infinity();
    double final_objective_change = std::numeric_limits<double>::infinity();
    double final_step = std::numeric_limits<double>::quiet_NaN();
    double final_projection_error = std::numeric_limits<double>::infinity();
    Rcpp::List history(entropy_schedule.n_elem);

    for (arma::uword stage = 0; stage < entropy_schedule.n_elem; ++stage) {
      double entropy = entropy_schedule[stage];
      double step = step_size;
      bool stage_converged = false;
      int iteration = 0;
      int projection_iterations = 0;
      for (iteration = 1; iteration <= maxit; ++iteration) {
        arma::mat coupling = augmented.submat(0, 0, nr - 1, ns - 1);
        ObjectiveV3 current = objective_v3(
          coupling, reference_mass, source_mass,
          reference_relation, source_relation, spatial_cost,
          temporal_weight, selection_weight, entropy, true
        );
        arma::mat gradient = current.gradient - matrix_median(current.gradient);
        bool accepted = false;
        double trial_step = step;
        arma::mat proposal = augmented;
        ObjectiveV3 proposal_objective = current;
        for (int backtrack = 0; backtrack <= 20; ++backtrack) {
          arma::mat kernel = augmented;
          arma::mat exponent = arma::clamp(-trial_step * gradient, -50.0, 50.0);
          kernel.submat(0, 0, nr - 1, ns - 1) =
            arma::clamp(coupling, 1e-300,
                        std::numeric_limits<double>::max()) % arma::exp(exponent);
          double projection_error = std::numeric_limits<double>::infinity();
          int projection_iteration = 0;
          bool projected = project_standard(
            kernel, row_target, column_target,
            projection_maxit, projection_tolerance,
            proposal, projection_error, projection_iteration
          );
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
            temporal_weight, selection_weight, entropy, false
          );
          if (std::isfinite(proposal_objective.optimization) &&
              proposal_objective.optimization <= current.optimization +
                1e-12 * std::max(1.0, std::abs(current.optimization))) {
            accepted = true;
            break;
          }
          trial_step /= 2.0;
        }
        if (!accepted) {
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
        if (final_objective_change <= tolerance &&
            final_projection_error <= projection_tolerance) {
          stage_converged = true;
          break;
        }
      }
      all_converged = all_converged && stage_converged;
      history[stage] = Rcpp::List::create(
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
        Rcpp::Named("projection_method") = "native_standard",
        Rcpp::Named("projection_iterations") = projection_iterations,
        Rcpp::Named("projection_fallback") = false
      );
    }
    arma::mat coupling = augmented.submat(0, 0, nr - 1, ns - 1);
    ObjectiveV3 final = objective_v3(
      coupling, reference_mass, source_mass,
      reference_relation, source_relation, spatial_cost,
      temporal_weight, selection_weight,
      entropy_schedule[entropy_schedule.n_elem - 1], false
    );
    fits[coverage_index] = Rcpp::List::create(
      Rcpp::Named("start") = "independent_native",
      Rcpp::Named("coverage") = coverage,
      Rcpp::Named("coupling") = coupling,
      Rcpp::Named("augmented") = augmented,
      Rcpp::Named("objective") = named_objective(final, coverage),
      Rcpp::Named("history") = history,
      Rcpp::Named("converged") = all_converged && numerical_ok,
      Rcpp::Named("projection_converged") =
        final_projection_error <= projection_tolerance,
      Rcpp::Named("stages_converged") = all_converged,
      Rcpp::Named("native_ok") = numerical_ok,
      Rcpp::Named("final_change") = final_change,
      Rcpp::Named("final_objective_change") = final_objective_change,
      Rcpp::Named("projected_update_residual") = final_change /
        std::max(final_step, std::numeric_limits<double>::epsilon()),
      Rcpp::Named("start_objectives") = Rcpp::NumericVector::create(
        Rcpp::Named("independent_native") = final.optimization
      ),
      Rcpp::Named("start_scientific_objectives") = Rcpp::NumericVector::create(
        Rcpp::Named("independent_native") = final.scientific
      ),
      Rcpp::Named("start_spread") = 0.0,
      Rcpp::Named("start_scientific_spread") = 0.0
    );
    previous_augmented = augmented;
    have_previous_coverage = true;
  }
  return Rcpp::List::create(
    Rcpp::Named("fits") = fits,
    Rcpp::Named("native_available") = true,
    Rcpp::Named("backend") = "rcpparmadillo"
  );
}
