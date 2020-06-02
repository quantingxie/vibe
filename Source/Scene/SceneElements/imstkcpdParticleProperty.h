#ifndef CPDPARTICLEPROPERTY_H
#define CPDPARTICLEPROPERTY_H

#include <memory>

#include "imstkcpdConstants.h"

namespace cpd
{
  enum class PropertyType
  {
    RADIUS,
    MASS,
    PROXIMITY,
    STATICFRICTION,
    DYNAMICFRICTION,
    STOPTHRESHOLD
  };

  class ParticleProperty
  {
  public:
    ParticleProperty();
    ParticleProperty(float p_radius, float p_mass = PARTICLE_MASS);

    void setProperty(PropertyType p_type, float p_value);
    float getProperty(PropertyType p_type) const;

  private:
    float m_radius;
    float m_mass;
    float m_staticFrictionCoefficient;
    float m_dynamicFrictionCoefficient;
    float m_stopThreshold;
    float m_proximity;
  };

  using ParticlePropertyPtr = std::shared_ptr<ParticleProperty>;
}
#endif // !CPDPARTICLEPROPERTY_H
