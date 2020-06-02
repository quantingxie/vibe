#include "imstkcpdEdgeEdgeCollisionConstraint.h"

namespace cpd
{
  void EdgeEdgeCollisionConstraint::initConstraint(ParticleObjectPtr p_object1, const size_t p_idx1, const size_t p_idx2, ParticleObjectPtr p_object2, const size_t p_idx3, const size_t p_idx4, double p_stiffness)
  {
    m_stiffness = p_stiffness;

    m_particleIDs[0] = p_idx1;
    m_particleIDs[1] = p_idx1;
    m_particleIDs[2] = p_idx1;
    m_particleIDs[3] = p_idx1;

    for (unsigned i = 0; i < 4; i++)
    {
      m_invMass[i] = m_object->getInvMass(m_particleIDs[i]);
      m_deltaDisplacement[i].setZero();
      if (m_invMass[i] > EPS)
        m_movable[i] = true;
      else
        m_movable[i] = false;
    }
    m_volume = p_object1->getParticleProximity() + p_object2->getParticleProximity();
    m_elasticityMatrix[0] = 0.0;
    m_deno[0] = 0.0;
    clearLamda();
    computeDerivative();
  }

  void EdgeEdgeCollisionConstraint::updateConstraint()
  {
    m_delta[0] = m_delta[0] + m_deltaDisplacement[1] - m_deltaDisplacement[0];
    m_delta[1] = m_delta[1] + m_deltaDisplacement[3] - m_deltaDisplacement[2];
    m_delta[2] = m_delta[2] + m_deltaDisplacement[0] - m_deltaDisplacement[2];

    double a = m_delta[1].dot(m_delta[0]);
    double b = m_delta[0].dot(m_delta[0]);
    double c = m_delta[2].dot(m_delta[0]);
    double d = m_delta[1].dot(m_delta[1]);
    double e = a;
    double f = m_delta[2].dot(m_delta[1]);

    double det = a * e - d * b;
    m_st[0] = 0.5;
    m_st[1] = 0.5;
    if (fabs(det) > EPS*EPS)
    {
      m_st[0] = (c*e - b * f) / det;
      m_st[1] = (c*d - a * f) / det;

      if (m_st[0] < 0)
        m_st[0] = 0.0;
      else if (m_st[0] > 1.0)
        m_st[0] = 1.0;

      if (m_st[1] < 0)
        m_st[1] = 0.0;
      else if (m_st[1] > 1.0)
        m_st[1] = 1.0;

    }

    for (unsigned i = 0; i < 4; i++)
      m_object->setParticleCollided(m_particleIDs[i]);

    Vec3d P = m_object->getTemporaryPosition(m_particleIDs[0]) + m_object->getTemporaryPosition(m_particleIDs[0]) + m_st[1] * m_delta[0];
    Vec3d Q = m_object->getTemporaryPosition(m_particleIDs[2]) + m_object->getTemporaryPosition(m_particleIDs[2]) + m_st[0] * m_delta[1];

    m_normal = Q - P;

    m_Constraint[0] = m_normal.norm() - m_volume;

  }

  void EdgeEdgeCollisionConstraint::computeDerivative()
  {
    m_delta[0] = m_object->getTemporaryPosition(m_particleIDs[1]) - m_object->getTemporaryPosition(m_particleIDs[0]);
    m_delta[1] = m_object->getTemporaryPosition(m_particleIDs[3]) - m_object->getTemporaryPosition(m_particleIDs[2]);
    m_delta[2] = m_object->getTemporaryPosition(m_particleIDs[0]) - m_object->getTemporaryPosition(m_particleIDs[2]);

    for (unsigned i = 0; i < 4; i++)
    {
      m_derivatives[i].fill(0.0);
    }
  }

  bool EdgeEdgeCollisionConstraint::updateDenominator()
  {
    if (m_Constraint[0] > EPS)
      return false;


    double coeff[4];

    coeff[0] = -(1.0 - m_st[1]);
    coeff[1] = -m_st[1];
    coeff[2] = (1.0 - m_st[0]);
    coeff[3] = m_st[0];

    for (unsigned i = 0; i < 4; i++)
      m_object->setParticleCollided(m_particleIDs[i]);

    m_normal.normalize();
    for (unsigned i = 0; i < 4; i++)
    {
      for (unsigned j = 0; j < 3; j++)
      {
        m_derivatives[i][j] = coeff[i] * m_normal[j];
      }
    }

    double cmc = 0;
    for (unsigned i = 0; i < 4; i++)
    {
      cmc += m_invMass[i] * coeff[i] * coeff[i];
    }

    if (abs(cmc) > EPS)
      m_deno[0] = 1.0 / cmc;
    else
      m_deno[0] = 0.0;

    return true;
  }

  void EdgeEdgeCollisionConstraint::computeForce()
  {
    Vec3d force;
    double coeff;
    for (unsigned i = 0; i < 4; i++)
    {
      coeff = m_lamdaList[0];
      for (unsigned j = 0; j < 3; j++)
      {
        force[j] = coeff * m_derivatives[i][j];
      }
      m_object->addToConstraintForce(m_particleIDs[i], m_ID, m_invMass[i] * force);
    }
  }

  void EdgeEdgeCollisionConstraint::integrate(bool p_all)
  {
  }

  void EdgeEdgeCollisionConstraint::initShapeFunction()
  {
  }

}