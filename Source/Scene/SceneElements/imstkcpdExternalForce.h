#ifndef CPDEXTERNALFORCE_H
#define CPDEXTERNALFORCE_H

#include <memory>
#include <vector>

#include "imstkcpdMath.h"
#include "imstkcpdParticleObject.h"

namespace cpd
{
  class ExternalForce
  {
  public:
    ExternalForce() {}

    void setForce(const Vec3d& p_force) { m_force = p_force; }
    void setForce(double x, double y, double z) { m_force = Vec3d(x, y, z); }
    const Vec3d& getForce() const { return m_force; }

    bool isActive() const { return m_isActive; }
    void activiate();
    void deactiviate();

    bool isForAll() const { return m_isForAll; }
    void addAffectedObject(ParticleObjectPtr p_object);
    const std::vector<ParticleObjectPtr>& getAffectedObjects() const { return m_affectedObjects; }

  private:
    Vec3d m_force = Vec3d(0.0, 0.0, 0.0);
    std::vector<ParticleObjectPtr> m_affectedObjects;

    bool m_isActive = true;
    bool m_isForAll = true;
  };

  using ExternalForcePtr = std::shared_ptr<ExternalForce>;
}


#endif // !CPDEXTERNALFORCE_H
