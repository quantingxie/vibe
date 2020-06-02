#ifndef CPDSEMIIMPLICITEULERINTEGRATOR_H
#define CPDSEMIIMPLICITEULERINTEGRATOR_H

#include "imstkcpdMath.h"
#include "imstkcpdBaseTimeIntegrator.h"

namespace cpd
{
  class SemiImplicitEulerIntegrator :public BaseTimeIntegrator {

  public:
    SemiImplicitEulerIntegrator() {}
    Type getType() const override { return Type::SIE; }
    void integrate(Vec3d& p, Vec3d& v, const Vec3d& a) override;

  private:

  };

  using SemiImplicitEulerIntegratorPtr = std::shared_ptr<SemiImplicitEulerIntegrator>;
}

#endif // !1CPDSEMIIMPLICITEULERINTEGRATOR_H

