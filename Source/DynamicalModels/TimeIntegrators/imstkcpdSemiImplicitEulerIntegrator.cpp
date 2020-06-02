#include "imstkcpdSemiImplicitEulerIntegrator.h"

namespace cpd
{
  void SemiImplicitEulerIntegrator::integrate(Vec3d & p, Vec3d & v, const Vec3d & a)
  {
    v = v + m_timeStepSize*a;
    p = p + m_timeStepSize*v;
  }

}