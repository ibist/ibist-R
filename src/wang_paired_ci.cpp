#include <Rcpp.h>

#include <algorithm>
#include <cmath>
#include <limits>
#include <vector>

using Rcpp::IntegerVector;
using Rcpp::List;
using Rcpp::NumericVector;

namespace {

constexpr double kEps = 1e-10;

struct State {
  int n10;
  int t;
  int n01;
  double score;
  double log_coef;
};

double tail_prob_at_p0(const std::vector<State>& states,
                       const std::vector<int>& indices,
                       double delta,
                       double p0) {
  const double p10 = (1.0 + delta - p0) / 2.0;
  const double p01 = (1.0 - p0 - delta) / 2.0;

  if (p10 <= 0.0 || p0 <= 0.0 || p01 <= 0.0) return 0.0;

  const double lp10 = std::log(p10);
  const double lp0 = std::log(p0);
  const double lp01 = std::log(p01);
  double out = 0.0;

  for (int idx : indices) {
    const State& s = states[idx];
    out += std::exp(s.log_coef +
                    static_cast<double>(s.n10) * lp10 +
                    static_cast<double>(s.t) * lp0 +
                    static_cast<double>(s.n01) * lp01);
  }

  return out;
}

double grid_max_tail(const std::vector<State>& states,
                     const std::vector<int>& indices,
                     double delta,
                     int grid_one,
                     int grid_two) {
  const double lo = kEps;
  const double hi = std::min(1.0 - delta - kEps, 1.0 + delta - kEps);
  if (!(hi > lo)) return 0.0;

  int grid = std::max(2, grid_one);
  double best = -1.0;
  double best_p0 = lo;

  for (int i = 0; i < grid; ++i) {
    const double p0 = lo + (hi - lo) * static_cast<double>(i) /
                             static_cast<double>(grid - 1);
    const double prob = tail_prob_at_p0(states, indices, delta, p0);
    if (prob > best) {
      best = prob;
      best_p0 = p0;
    }
  }

  if (grid_two <= 1) return best;

  const double step = (hi - lo) / static_cast<double>(grid);
  const double lo2 = std::max(lo, best_p0 - step) + kEps;
  const double hi2 = std::min(hi, best_p0 + step) - kEps;
  if (!(hi2 > lo2)) return best;

  grid = std::max(2, grid_two);
  for (int i = 0; i < grid; ++i) {
    const double p0 = lo2 + (hi2 - lo2) * static_cast<double>(i) /
                              static_cast<double>(grid - 1);
    best = std::max(best, tail_prob_at_p0(states, indices, delta, p0));
  }

  return best;
}

double solve_lower_limit(const std::vector<State>& states,
                         const std::vector<int>& indices,
                         double alpha,
                         double upper,
                         double precision,
                         int grid_one,
                         int grid_two) {
  double lo = -1.0 + kEps;
  double hi = std::min(1.0 - 100.0 * kEps, upper);
  const double tol = std::max(precision, 1e-12);

  while (std::abs(hi - lo) >= tol) {
    const double mid = (hi + lo) / 2.0;
    const double prob = grid_max_tail(states, indices, mid, grid_one, grid_two);
    if (prob >= alpha) {
      hi = mid;
    } else {
      lo = mid;
    }
  }

  return lo;
}

std::vector<State> make_states(int n) {
  std::vector<double> lfact(n + 1);
  for (int i = 0; i <= n; ++i) lfact[i] = std::lgamma(static_cast<double>(i) + 1.0);

  std::vector<State> states;
  states.reserve((n + 1) * (n + 2) / 2);

  for (int n10 = 0; n10 <= n; ++n10) {
    for (int t = 0; t <= n - n10; ++t) {
      const int n01 = n - n10 - t;
      states.push_back(State{
        n10,
        t,
        n01,
        (static_cast<double>(n10) - static_cast<double>(n01)) /
          static_cast<double>(n),
        lfact[n] - lfact[n10] - lfact[t] - lfact[n01]
      });
    }
  }

  std::sort(states.begin(), states.end(), [](const State& a, const State& b) {
    if (a.score != b.score) return a.score > b.score;
    if (a.n10 != b.n10) return a.n10 > b.n10;
    return a.t > b.t;
  });

  return states;
}

double round_to_precision(double x, double precision) {
  if (!(precision > 0.0)) return x;
  const double scale = std::pow(10.0, std::ceil(std::log10(1.0 / precision)));
  return std::round(x * scale) / scale;
}

double wang_lower_one_sided(int target_n10,
                            int target_t,
                            int target_n01,
                            double conf_level,
                            double precision,
                            int grid_one,
                            int grid_two) {
  const int n = target_n10 + target_t + target_n01;
  if (n == 0) return -1.0;
  if (target_n10 == 0 && target_t == 0) return -1.0;

  std::vector<State> states = make_states(n);
  const int m = static_cast<int>(states.size());

  std::vector<std::vector<int> > state_index(n + 1, std::vector<int>(n + 1, -1));
  int target = -1;
  for (int i = 0; i < m; ++i) {
    state_index[states[i].n10][states[i].t] = i;
    if (states[i].n10 == target_n10 &&
        states[i].t == target_t &&
        states[i].n01 == target_n01) {
      target = i;
    }
  }

  if (target < 0) Rcpp::stop("Target table is outside the paired table space.");

  const double alpha = 1.0 - conf_level;
  std::vector<int> selected;
  selected.reserve(m);
  selected.push_back(0);

  std::vector<unsigned char> remaining(m, 1);
  remaining[0] = 0;

  double current_limit = solve_lower_limit(
    states, selected, alpha, 1.0 - 100.0 * kEps, 1e-5, 500, 0
  );
  current_limit = round_to_precision(current_limit, precision);
  if (target == 0) return current_limit;

  while (static_cast<int>(selected.size()) < m) {
    std::vector<int> eligible;
    eligible.reserve(n + 1);

    for (int idx = 0; idx < m; ++idx) {
      if (!remaining[idx]) continue;

      const State& s = states[idx];
      bool blocked = false;

      if (s.n10 + s.t + 1 <= n) {
        const int pred = state_index[s.n10][s.t + 1];
        if (pred >= 0 && remaining[pred]) blocked = true;
      }

      if (!blocked && s.t >= 1 && s.n10 + s.t <= n) {
        const int pred = state_index[s.n10 + 1][s.t - 1];
        if (pred >= 0 && remaining[pred]) blocked = true;
      }

      if (!blocked) eligible.push_back(idx);
    }

    if (eligible.empty()) break;

    double best_limit = -std::numeric_limits<double>::infinity();
    std::vector<double> candidate_limits(eligible.size());

    for (std::size_t i = 0; i < eligible.size(); ++i) {
      std::vector<int> trial = selected;
      trial.push_back(eligible[i]);
      candidate_limits[i] = solve_lower_limit(
        states, trial, alpha, current_limit, precision, grid_one, grid_two
      );
      best_limit = std::max(best_limit, candidate_limits[i]);
    }

    std::vector<int> group;
    for (std::size_t i = 0; i < eligible.size(); ++i) {
      if (candidate_limits[i] >= best_limit - kEps) group.push_back(eligible[i]);
    }

    std::vector<int> trial = selected;
    trial.insert(trial.end(), group.begin(), group.end());
    current_limit = solve_lower_limit(
      states, trial, alpha, current_limit, precision, grid_one, grid_two
    );
    current_limit = round_to_precision(current_limit, precision);

    bool found_target = false;
    for (int idx : group) {
      remaining[idx] = 0;
      selected.push_back(idx);
      if (idx == target) found_target = true;
    }

    if (found_target) return current_limit;
  }

  return -1.0;
}

NumericVector wang_lower_all_one_sided(int n,
                                       double conf_level,
                                       double precision,
                                       int grid_one,
                                       int grid_two) {
  if (n == 0) {
    NumericVector empty(0);
    return empty;
  }

  std::vector<State> states = make_states(n);
  const int m = static_cast<int>(states.size());

  std::vector<std::vector<int> > state_index(n + 1, std::vector<int>(n + 1, -1));
  for (int i = 0; i < m; ++i) state_index[states[i].n10][states[i].t] = i;

  const double alpha = 1.0 - conf_level;
  std::vector<int> selected;
  selected.reserve(m);
  selected.push_back(0);

  std::vector<unsigned char> remaining(m, 1);
  remaining[0] = 0;

  NumericVector lower_by_code((n + 1) * (n + 2), NumericVector::get_na());

  double current_limit = solve_lower_limit(
    states, selected, alpha, 1.0 - 100.0 * kEps, 1e-5, 500, 0
  );
  current_limit = round_to_precision(current_limit, precision);
  lower_by_code[states[0].n10 * (n + 2) + states[0].t] = current_limit;

  while (static_cast<int>(selected.size()) < m) {
    std::vector<int> eligible;
    eligible.reserve(n + 1);

    for (int idx = 0; idx < m; ++idx) {
      if (!remaining[idx]) continue;

      const State& s = states[idx];
      bool blocked = false;

      if (s.n10 + s.t + 1 <= n) {
        const int pred = state_index[s.n10][s.t + 1];
        if (pred >= 0 && remaining[pred]) blocked = true;
      }

      if (!blocked && s.t >= 1 && s.n10 + s.t <= n) {
        const int pred = state_index[s.n10 + 1][s.t - 1];
        if (pred >= 0 && remaining[pred]) blocked = true;
      }

      if (!blocked) eligible.push_back(idx);
    }

    if (eligible.empty()) break;

    double best_limit = -std::numeric_limits<double>::infinity();
    std::vector<double> candidate_limits(eligible.size());

    for (std::size_t i = 0; i < eligible.size(); ++i) {
      std::vector<int> trial = selected;
      trial.push_back(eligible[i]);
      candidate_limits[i] = solve_lower_limit(
        states, trial, alpha, current_limit, precision, grid_one, grid_two
      );
      best_limit = std::max(best_limit, candidate_limits[i]);
    }

    std::vector<int> group;
    for (std::size_t i = 0; i < eligible.size(); ++i) {
      if (candidate_limits[i] >= best_limit - kEps) group.push_back(eligible[i]);
    }

    std::vector<int> trial = selected;
    trial.insert(trial.end(), group.begin(), group.end());
    current_limit = solve_lower_limit(
      states, trial, alpha, current_limit, precision, grid_one, grid_two
    );
    current_limit = round_to_precision(current_limit, precision);

    for (int idx : group) {
      remaining[idx] = 0;
      selected.push_back(idx);
      lower_by_code[states[idx].n10 * (n + 2) + states[idx].t] = current_limit;
    }
  }

  for (int n10 = 0; n10 <= n; ++n10) {
    for (int t = 0; t <= n - n10; ++t) {
      const int code = n10 * (n + 2) + t;
      if (NumericVector::is_na(lower_by_code[code])) lower_by_code[code] = -1.0;
    }
  }

  return lower_by_code;
}

}  // namespace

// [[Rcpp::export]]
NumericVector wang_paired_ci_cpp(int n10,
                                 int t,
                                 int n01,
                                 double conf_level = 0.95,
                                 std::string ci_type = "Two.sided",
                                 double precision = 0.00001,
                                 int grid_one = 30,
                                 int grid_two = 20) {
  if (n10 < 0 || t < 0 || n01 < 0) Rcpp::stop("Counts must be non-negative.");
  if (!(conf_level > 0.0 && conf_level < 1.0)) {
    Rcpp::stop("conf_level must be between 0 and 1.");
  }
  if (!(precision > 0.0)) Rcpp::stop("precision must be positive.");
  if (grid_one < 2 || grid_two < 2) Rcpp::stop("Grid sizes must be at least 2.");

  const int n = n10 + t + n01;
  const double estimate = n == 0 ? 0.0 :
    (static_cast<double>(n10) - static_cast<double>(n01)) / static_cast<double>(n);

  double lower = -1.0;
  double upper = 1.0;

  if (ci_type == "Lower") {
    lower = wang_lower_one_sided(n10, t, n01, conf_level, precision, grid_one, grid_two);
  } else if (ci_type == "Upper") {
    upper = -wang_lower_one_sided(n01, t, n10, conf_level, precision, grid_one, grid_two);
  } else if (ci_type == "Two.sided") {
    const double one_sided_level = 1.0 - (1.0 - conf_level) / 2.0;
    lower = wang_lower_one_sided(n10, t, n01, one_sided_level, precision, grid_one, grid_two);
    upper = -wang_lower_one_sided(n01, t, n10, one_sided_level, precision, grid_one, grid_two);
  } else {
    Rcpp::stop("ci_type must be 'Lower', 'Upper', or 'Two.sided'.");
  }

  NumericVector out = NumericVector::create(
    Rcpp::Named("estimate") = estimate,
    Rcpp::Named("lower") = lower,
    Rcpp::Named("upper") = upper
  );
  return out;
}

// [[Rcpp::export]]
NumericVector wang_paired_lower_all_cpp(int n,
                                        double conf_level = 0.95,
                                        double precision = 0.00001,
                                        int grid_one = 30,
                                        int grid_two = 20) {
  if (n < 0) Rcpp::stop("n must be non-negative.");
  if (!(conf_level > 0.0 && conf_level < 1.0)) {
    Rcpp::stop("conf_level must be between 0 and 1.");
  }
  if (!(precision > 0.0)) Rcpp::stop("precision must be positive.");
  if (grid_one < 2 || grid_two < 2) Rcpp::stop("Grid sizes must be at least 2.");

  return wang_lower_all_one_sided(n, conf_level, precision, grid_one, grid_two);
}
