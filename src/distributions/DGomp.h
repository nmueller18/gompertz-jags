#ifndef DGOMP_H
#define DGOMP_H

#include <distribution/RScalarDist.h>  // JAGS base class for RScalarDist

namespace Gompertz {

class DGomp : public jags::RScalarDist {
public:
  // Constructor
  DGomp();

  // Probability density function (PDF)
  double d(double x, jags::PDFType type, std::vector<double const *> const &params, bool give_log) const;

  // Cumulative distribution function (CDF)
  double p(double x, std::vector<double const *> const &params, bool lower, bool give_log) const;

  // Quantile function (Inverse CDF)
  double q(double p, std::vector<double const *> const &params, bool lower, bool log_p) const;

  // Random number generation (RNG)
  double r(std::vector<double const *> const &params, jags::RNG *rng) const;

  // Parameter validation
  bool checkParameterValue(std::vector<double const *> const &params) const;

};

} // namespace Gompertz

#endif // GOMPERTZ_H
