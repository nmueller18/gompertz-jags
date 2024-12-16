#ifndef DGOMP_H
#define DGOMP_H

#include <distribution/ScalarDist.h> // JAGS base class for distributions

namespace Gompertz {

class DGomp : public jags::ScalarDist {
public:
    DGomp(); // Constructor

    // Corrected method signatures with proper overrides
    double logDensity(double x, jags::PDFType type,
                      std::vector<double const *> const &params, double const *lower, double const *upper) const;

    double randomSample(std::vector<double const *> const &params, double const *lower, double const *upper, jags::RNG *rng) const;

    double typicalValue(std::vector<double const *> const &params, double const *lower, double const *upper) const;

    bool checkParameterValue(std::vector<double const *> const &params) const;
};

} // namespace Gompertz

#endif // GOMPERTZ_H
