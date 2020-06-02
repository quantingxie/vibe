#include "imstkcpdCollisionConstraint.h"

namespace cpd
{

  template<int SIZE, int SIZEC, int DIM>
  void CollisionConstraintBase<SIZE, SIZEC, DIM>::updateDeltaDisplacement()
  {
    ParticleObjectPtr tempPOptr = m_object;
    bool possibleSwitch = true;

    for (unsigned i = 0; i < 4; i++)
    {
      unsigned idx = m_particleIDs[i];

      if (possibleSwitch && (i >= m_sizeCol[0]))
      {
        tempPOptr = m_objectOther;
        possibleSwitch = false;
      }

      if (tempPOptr->isForceUpdated(idx) && m_movable[i])
        m_deltaDisplacement[i] = tempPOptr->getConstraintForce(idx);

    }

  }

}