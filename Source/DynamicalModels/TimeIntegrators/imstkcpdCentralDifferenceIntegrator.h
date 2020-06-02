#ifndef CPDCENTRALDIFFERENCEINTEGRATOR_H
#define CPDCENTRALDIFFERENCEINTEGRATOR_H

#include "imstkcpdMath.h"
#include "imstkcpdBaseTimeIntegrator.h"

namespace cpd
{
  class CentralDifferenceIntegrator :public BaseTimeIntegrator {

  public:
    CentralDifferenceIntegrator() {}
    Type getType() const override { return Type::CD; }
    void integrate(Vec3d& p, Vec3d& v, const Vec3d& a) override;

  private:

  };

  using CentralDifferenceIntegratorPtr = std::shared_ptr<CentralDifferenceIntegrator>;
}

#endif // !1CPDSEMIIMPLICITEULERINTEGRATOR_H

