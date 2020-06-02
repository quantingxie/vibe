#ifndef CPDBASETIMEINTEGRATOR_H
#define CPDBASETIMEINTEGRATOR_H

#include <memory>

#include "imstkcpdMath.h"
#include "imstkcpdConstants.h"

namespace cpd
{

  class BaseTimeIntegrator {
  public:
    enum class Type
    {
      SIE, //Semi-Implicit Euler
      CD, //Central Difference
    };

  public:
    BaseTimeIntegrator() {}

    virtual Type getType() const = 0;
    void setTimeStepSize(float p_timeStepSize) { m_timeStepSize = p_timeStepSize; }
    float getTimeStepSize() const { return m_timeStepSize; }

    virtual void integrate(Vec3d& p, Vec3d& v, const Vec3d& a) = 0;

  protected:
    float m_timeStepSize = TIMESTEP;

  };

  using BaseTimeIntegratorPtr = std::shared_ptr<BaseTimeIntegrator>;
}
#endif // !CPDBASETIMEINTEGRATOR_H