#include "imstkcpdParticleProperty.h"

namespace cpd
{
  ParticleProperty::ParticleProperty()
  {
    m_radius = PARTICLE_RADIUS;
    m_mass = PARTICLE_MASS;
    m_proximity = PARTICLE_PROXIMITY;
    m_staticFrictionCoefficient = STATIC_FRICTION;
    m_dynamicFrictionCoefficient = DYNAMIC_FRICTION;
    m_stopThreshold = STOP_THRESHOLD;
  }

  ParticleProperty::ParticleProperty(float p_radius, float p_mass)
  {
    m_radius = p_radius;
    m_mass = p_mass;

    m_staticFrictionCoefficient = STATIC_FRICTION;
    m_dynamicFrictionCoefficient = DYNAMIC_FRICTION;
    m_stopThreshold = STOP_THRESHOLD;
  }

  void ParticleProperty::setProperty(PropertyType p_type, float p_value)
  {
    switch (p_type)
    {
    case PropertyType::RADIUS:
      m_radius = p_value;
      break;
    case PropertyType::MASS:
      m_mass = p_value;
      break;
    case PropertyType::PROXIMITY:
      m_proximity = p_value;
      break;
    case PropertyType::STATICFRICTION:
      m_staticFrictionCoefficient = p_value;
      break;
    case PropertyType::DYNAMICFRICTION:
      m_dynamicFrictionCoefficient = p_value;
      break;
    case PropertyType::STOPTHRESHOLD:
      m_stopThreshold = p_value;
      break;
    default:
      break;
    }
  }

  float ParticleProperty::getProperty(PropertyType p_type) const
  {
    float property_value = -1.0;

    switch (p_type)
    {
    case PropertyType::RADIUS:
      property_value = m_radius;
      break;
    case PropertyType::MASS:
      property_value = m_mass;
      break;
    case PropertyType::PROXIMITY:
      property_value = m_proximity;
      break;
    case PropertyType::STATICFRICTION:
      property_value = m_staticFrictionCoefficient;
      break;
    case PropertyType::DYNAMICFRICTION:
      property_value = m_dynamicFrictionCoefficient;
      break;
    case PropertyType::STOPTHRESHOLD:
      property_value = m_stopThreshold;
      break;
    default:
      break;
    }
    return property_value;
  }

}