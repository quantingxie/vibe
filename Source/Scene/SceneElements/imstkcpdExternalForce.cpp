#include "imstkcpdExternalForce.h"

namespace cpd
{
  void ExternalForce::activiate()
  {
    if (!m_isActive)
      m_isActive = true;
  }

  void ExternalForce::deactiviate()
  {
    if (m_isActive)
      m_isActive = false;
  }

  void ExternalForce::addAffectedObject(ParticleObjectPtr p_object)
  {
    if (m_isForAll)
      m_isForAll = false;

    m_affectedObjects.push_back(p_object);
  }

}