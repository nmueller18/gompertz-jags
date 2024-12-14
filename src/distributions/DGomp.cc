#include "DGomp.h"
#include <cmath>      // For exp, log
#include <stdexcept>  // For exceptions
#include <rng/RNG.h>  // Provides random functions
#include <util/nainf.h> // Provides na and inf functions, etc.

using std::vector;

namespace Gompertz {

DGomp::DGomp() : ScalarDist("dgomp", 3, jags::DIST_POSITIVE) {}

// PDF (density)
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

    // Compute the density for x >= 0
    double exp_bx = std::exp(b * x);
    double density = a * exp_bx * std::exp(-a/b * (exp_bx - 1));

    // Return log(density) if the density is positive; otherwise, return JAGS_NEGINF
    return (density == 0) ? JAGS_NEGINF : std::log(density);
}



// RNG (random generation)
double DGomp::randomSample(std::vector<double const *> const &params,
                           double const *lower, double const *upper,
                           jags::RNG *rng) const {
    double b = *params[0];
    double a = *params[1];
    double max_age = *params[2];

    // Generate a valid sample respecting bounds
    double x;
    do {
        double u = rng->uniform(); // Generate a uniform random variable in (0,1)
        x = std::log(1 - (std::log(1 - u) / a)) / b; // Transform uniform random variable
    } while ((lower && x < *lower) || (upper && x > *upper) || x > max_age); // Check bounds

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
