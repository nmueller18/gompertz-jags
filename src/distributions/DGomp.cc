#include "DGomp.h"
#include <cmath>      // For exp, log
#include <stdexcept>  // For exceptions
#include <rng/RNG.h>  // Provides random functions
#include <util/nainf.h> // Provides na and inf functions, etc.
#include <limits>

using std::vector;

namespace Gompertz {
DGomp::DGomp() : jags::RScalarDist("dgomp", 2, jags::DIST_POSITIVE) {}

// Probability Density Function (PDF)
double DGomp::d(double x, jags::PDFType type, std::vector<double const *> const &params, bool give_log) const {
  double b = *params[0]; // Scale parameter
  double a = *params[1]; // Shape parameter

  if (x < 0) {
    return give_log ? -std::numeric_limits<double>::infinity() : 0.0; // PDF is 0 for x < 0
  }

  // Compute log-density to avoid underflow
  double log_exp_bx = b * x;
  double exp_bx = std::exp(log_exp_bx);
  double log_density = std::log(a) + log_exp_bx - (a / b) * (exp_bx - 1);

  return give_log ? log_density : std::exp(log_density);
}

// Cumulative Distribution Function (CDF)
double DGomp::p(double x, std::vector<double const *> const &params, bool lower, bool give_log) const {
  double b = *params[0]; // Scale parameter
  double a = *params[1]; // Shape parameter

  if (x < 0) return lower ? 0.0 : 1.0;

  double cdf = 1 - std::exp(- (a / b) * (std::exp(b * x) - 1));

  if (!lower) cdf = 1.0 - cdf;
  return give_log ? std::log(cdf) : cdf;
}

// Quantile Function (Inverse CDF)
double DGomp::q(double p, std::vector<double const *> const &params, bool lower, bool log_p) const {
  double b = *params[0]; // Scale parameter
  double a = *params[1]; // Shape parameter

  if (log_p) p = std::exp(p);
  if (!lower) p = 1.0 - p;

  if (p <= 0.0) return 0.0;  // q(0) = 0
  if (p >= 1.0) return std::numeric_limits<double>::infinity();  // q(1) = ∞

  return (1.0 / b) * std::log(1 - (b / a) * std::log(1 - p));
}

// Random Number Generation (RNG)
double DGomp::r(std::vector<double const *> const &params, jags::RNG *rng) const {
  double b = *params[0]; // Scale parameter
  double a = *params[1]; // Shape parameter

  double U = rng->uniform(); // Generate a uniform random number
  return (1.0 / b) * std::log(1 - (b / a) * std::log(1 - U));
}

// Parameter validation
bool DGomp::checkParameterValue(vector<double const *> const &params) const {
  double b = *params[0];
  double a = *params[1];
  return (b > 0.0 && a > 0.0);  // Ensure both parameters are positive
}

} // namespace Gompertz
