#include "DGomp.h"
#include <cmath>      // For exp, log
#include <stdexcept>  // For exceptions
#include <rng/RNG.h>  // Provides random functions
#include <util/nainf.h> // Provides na and inf functions, etc.

using std::vector;

namespace Gompertz {

DGomp::DGomp() : ScalarDist("dgomp", 3, jags::DIST_POSITIVE) {}

// PDF (logDensity) for the Gompertz distribution
double DGomp::logDensity(double x, jags::PDFType type,
                         std::vector<double const *> const &params,
                         double const *lower, double const *upper) const {
    double b = *params[0]; // Scale parameter
    double a = *params[1]; // Shape parameter
    double max_age = *params[2];

    // Handle x < 0, as the Gompertz distribution is only valid for x >= 0
    if (x < 0 || x > max_age) {
        return JAGS_NEGINF;  // Log(0) is -infinity, so we return JAGS_NEGINF
    }

    // Compute log(density) directly to avoid numerical issues
    double log_exp_bx = b * x;  // b * x is the argument of the first exp()
    double exp_bx = std::exp(log_exp_bx);  // exp(b * x)
    double log_density = std::log(a) + log_exp_bx - (a / b) * (exp_bx - 1);

    // Return log(density) if valid; otherwise, return JAGS_NEGINF
    return std::isfinite(log_density) ? log_density : JAGS_NEGINF;
}



// RNG (random generation)
double DGomp::randomSample(std::vector<double const *> const &params,
                           double const *lower, double const *upper,
                           jags::RNG *rng) const {
  double b = *params[0];
  double a = *params[1];
  double max_age = *params[2];

  // Truncate to effective maximum age
  double exp_neg_b_max_age = std::exp(-b * max_age);
  double F_upper = 1 - std::exp(-a * (1 - exp_neg_b_max_age));

  // Generate a uniform random variable in (0, F_upper)
  double u = rng->uniform() * F_upper;

  // Inverse CDF transformation
  double x = -std::log(1 - std::log(1 - u) / a) / b;

  return x;
}


// Typical value
double DGomp::typicalValue(std::vector<double const *> const &params, double const *lower, double const *upper) const {
    double b = *params[0];
    double a = *params[1];
    double max_age = *params[2];

    // The mode of the Gompertz distribution is always 0 (typical value)
    double value = 0.0;

    // Ensure the mode respects the bounds
    if (lower && value < *lower) {
        value = *lower;
    }
    if (upper && value > *upper) {
        value = *upper;
    }
    if (value > max_age) {
        value = max_age;
    }

    return value;
}


// Parameter validity check
bool DGomp::checkParameterValue(std::vector<double const *> const &params) const {
    double b = *params[0]; // Dereference pointers
    double a = *params[1];
    double max_age = *params[2]; // Dereference the new parameter
    return (a > 0.0 && b > 0.0 && max_age > 0.0); // All parameters must be positive
}

} // namespace Gompertz
