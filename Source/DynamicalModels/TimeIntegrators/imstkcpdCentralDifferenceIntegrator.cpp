#include "imstkcpdCentralDifferenceIntegrator.h"

namespace cpd
{
  void CentralDifferenceIntegrator::integrate(Vec3d & p, Vec3d & v, const Vec3d & a)
  {
 //   v = v + m_timeStepSize*a;
    p = 2 * p - v + m_timeStepSize*m_timeStepSize*a;
  }

}