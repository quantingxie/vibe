#include "imstkcpdExternalForceSolver.h"

namespace cpd
{
  void ExternalForceSolver::updateAffectedObject()
  {
    for (auto& f : m_forces)
    {
      for (auto& obj : f->getAffectedObjects())
      {
        m_objects.insert(obj);
      }
    }

    for (auto& f : m_distributedForces)
    {
      for (auto& obj : f->getAffectedObjects())
      {
        m_objects.insert(obj);
      }
    }

  }

  void ExternalForceSolver::exertForces()
  {
    for (auto& object : m_objects)
    {

      auto& positions = object->getPositions();
      auto& prevPositions = object->getTemporaryPositions(); // temporaryPositions holds prevPositions at this step
      auto& velocities = object->getVelocities();
      auto& masses = object->getMassList();
      auto& accelerations = object->getAccelerations();
      auto& tempPostions = object->getTemporaryPositions(); // tempPostions and prevPositions points to the same data

      size_t n = positions.size();

      StdVectorOfVec3d para2;
      para2.resize(n);

      switch (m_timeIntegrator->getType())
      {
      case BaseTimeIntegrator::Type::SIE:
        para2 = velocities;
        break;
      case BaseTimeIntegrator::Type::CD:
        para2 = prevPositions;
        break;
      default:
        break;
      }

      Vec3d p, acceleration;
      unsigned i = 0;
      for (auto& position : positions)
      {
        p = position;
        if (masses[i] > -EPS)
        {
          acceleration = accelerations[i];
          m_timeIntegrator->integrate(p, para2[i], acceleration);
        }
        tempPostions[i] = p;
        i++;
      }

      //object->tempApplyTempPositions();
    }
  }

}