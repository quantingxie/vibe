#include "imstkcpdOneDEnergyConstraint.h"

namespace cpd
{
  int OneDEnergyConstraint::m_idGenerator = 0;

  void OneDEnergyConstraint::initConstraint(ParticleObjectPtr p_object, unsigned p_idx1, unsigned p_idx2, double p_E, double p_mu)
  {
    m_ID = m_idGenerator++;

    std::array<unsigned, 2> idx;
    idx[0] = p_idx1;
    idx[1] = p_idx2;

    EnergyConstraint::initConstraint(p_object, idx, p_E, p_mu);
  }

  void OneDEnergyConstraint::initVariables()
  {
    m_elasticityMatrix(0, 0) = 1.0 / (m_YoungsModulus*m_dt * m_dt);
  }

  void OneDEnergyConstraint::initDerivatives()
  {
    computeDerivative();
    m_updateDerivative = false;
  }

  void OneDEnergyConstraint::computeDerivative()
  {
    if (m_updateDerivative)
    {
      const Vec3d& x0 = m_object->getPosition(m_particleIDs[0]);
      const Vec3d& x1 = m_object->getPosition(m_particleIDs[1]);

      double length = x1.x() - x0.x();
      double sqrtLength = sqrt(length);

      m_derivatives[0][0] = -1.0 / sqrtLength;
      m_derivatives[1][0] = 1.0 / sqrtLength;

      for (unsigned i = 0; i < 2; i++)
      {
        m_B(0, i) = m_derivatives[i][0];
      }
    }
  }

  void OneDEnergyConstraint::computeForce()
  {
    Vec3d force;
    for (unsigned i = 0; i < m_size; i++)
    {
      force[0] = m_derivatives[i][0] * m_lamdaList[0];
      m_object->addToConstraintForce(m_particleIDs[i], m_ID, m_invMass[i] * force);
    }
  }

  void OneDEnergyConstraint::updateConstraint()
  {
  }

  void OneDEnergyConstraint::reset()
  {
  }

  void OneDEnergyConstraint::integrate(bool p_all)
  {
  }

  void OneDEnergyConstraint::initShapeFunction()
  {
  }

}