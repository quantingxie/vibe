#include "imstkcpdAreaConstraint.h"

namespace cpd
{
  int AreaConstraint::m_idGenerator = 0;

  void AreaConstraint::initConstraint(ParticleObjectPtr p_object, unsigned p_idx1, unsigned p_idx2, unsigned p_idx3, const double k)
  {
    m_ID = m_idGenerator++;
    m_stiffness = k;

    std::array<unsigned, 3> idx;
    idx[0] = p_idx1;
    idx[1] = p_idx2;
    idx[2] = p_idx3;

    ConstraintTemplate::initConstraint(p_object, idx, m_stiffness);
  }

  void AreaConstraint::initVariables()
  {
    m_elasticityMatrix(0, 0) = m_stiffness / (m_dt * m_dt);
  }

  void AreaConstraint::initDerivatives()
  {
    computeDerivative();
    m_volume = (m_delta[0].cross(m_delta[1])).norm();
  }

  void AreaConstraint::computeDerivative()
  {
    const Vec3d& d0 = m_object->getTemporaryPosition(m_particleIDs[0]);
    const Vec3d& d1 = m_object->getTemporaryPosition(m_particleIDs[1]);
    const Vec3d& d2 = m_object->getTemporaryPosition(m_particleIDs[2]);

    m_delta[0] = d1 - d0;
    m_delta[1] = d2 - d0;

    Vec3d n = (m_delta[0].cross(m_delta[1])).normalized();

    Vec3d derivatives[3];
    derivatives[0] = -(d2 - d1).cross(n);
    derivatives[1] = -(d0 - d2).cross(n);
    derivatives[2] = -(d1 - d0).cross(n);

    for (unsigned i = 0; i < 3; i++)
    {
      for (unsigned j = 0; j < 3; j++)
      {
        m_derivatives[i][j] = derivatives[i][j];
      }
    }

  }

  bool AreaConstraint::updateDenominator()
  {
    Vec3d e1 = m_delta[0] + m_deltaDisplacement[1] - m_deltaDisplacement[0];
    Vec3d e2 = m_delta[1] + m_deltaDisplacement[2] - m_deltaDisplacement[0];

    Vec3d n = (e1.cross(e2)).normalized();
    Vec3d temp1, temp2;
    for (unsigned i = 0; i < 3; i++)
    {
      temp1[i] = m_invMass[1] * m_derivatives[1][i] - m_invMass[0] * m_derivatives[0][i];
      temp2[i] = m_invMass[2] * m_derivatives[2][i] - m_invMass[0] * m_derivatives[0][i];
    }
    double temp = n.dot(temp1.cross(e2) + e1.cross(temp2)) + m_elasticityMatrix(0, 0);
    if (abs(temp) < EPS*EPS)
    {
      m_deno(0, 0) = 1.0 / m_elasticityMatrix(0, 0);
      m_lamdaList.fill(0.0);
    }
    else
    {
      m_deno(0, 0) = 1.0 / temp;
    }

    return true;
  }

  void AreaConstraint::computeForce()
  {
    Vec3d force;
    double coeff;
    for (unsigned i = 0; i < m_size; i++)
    {
      coeff = m_lamdaList[0];
      for (unsigned j = 0; j < 3; j++)
      {
        force[j] = coeff * m_derivatives[i][j];
      }
      m_object->addToConstraintForce(m_particleIDs[i], m_ID, m_invMass[i] * force);
    }
  }

  void AreaConstraint::updateConstraint()
  {
    Vec3d e1 = m_delta[0] + m_deltaDisplacement[1] - m_deltaDisplacement[0];
    Vec3d e2 = m_delta[1] + m_deltaDisplacement[2] - m_deltaDisplacement[0];

    m_Constraint[0] = (e1.cross(e2)).norm() - m_volume;
  }

  void AreaConstraint::reset()
  {
  }

  void AreaConstraint::integrate(bool p_all)
  {
  }

  void AreaConstraint::initShapeFunction()
  {
  }

}