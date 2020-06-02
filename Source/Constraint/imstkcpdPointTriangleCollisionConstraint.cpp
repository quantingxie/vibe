#include "imstkcpdPointTriangleCollisionConstraint.h"

namespace cpd
{
  void PointTriangleCollisionConstraint::initConstraint(ParticleObjectPtr p_object1, const size_t p_idx1, ParticleObjectPtr p_object2, const size_t p_idx2, const size_t p_idx3, const size_t p_idx4, double p_stiffness)
  {
    m_object = p_object1;
    m_objectOther = p_object2;

    m_stiffness = p_stiffness;
    m_elasticityMatrix(0,0) = 0;
    
    m_particleIDs[0] = p_idx1;
    m_particleIDs[1] = p_idx2;
    m_particleIDs[2] = p_idx3;
    m_particleIDs[3] = p_idx4;

    for (unsigned i = 0; i < 4; i++)
    {
      if (i < m_sizeCol[0])
        m_invMass[i] = m_object->getInvMass(m_particleIDs[i]);
      else
        m_invMass[i] = m_objectOther->getInvMass(m_particleIDs[i]);

      m_deltaDisplacement[i].setZero();
      if (m_invMass[i] > EPS)
        m_movable[i] = true;
      else
        m_movable[i] = false;
    }
    m_volume = m_object->getParticleProximity() + m_objectOther->getParticleProximity();
    m_elasticityMatrix[0] = 0.0;
    m_deno[0] = 0.0;
    clearLamda();
    computeDerivative();
  }

  void PointTriangleCollisionConstraint::updateConstraint()
  {
    m_delta[0] = m_delta[0] + m_deltaDisplacement[0] - m_deltaDisplacement[1];
    m_delta[1] = m_delta[1] + m_deltaDisplacement[2] - m_deltaDisplacement[1];
    m_delta[2] = m_delta[2] + m_deltaDisplacement[3] - m_deltaDisplacement[1];
    //m_delta[1] = m_delta[1] + m_deltaDisplacement[3] - m_deltaDisplacement[1];
    //m_delta[2] = m_delta[2] + m_deltaDisplacement[2] - m_deltaDisplacement[1];

    m_Constraint[0] = m_delta[0].dot((m_delta[1].cross(m_delta[2])).normalized()) - m_volume;
  }

  void PointTriangleCollisionConstraint::computeDerivative()
  {
    m_delta[0] = m_object->getTemporaryPosition(m_particleIDs[0]) - m_objectOther->getTemporaryPosition(m_particleIDs[1]);
    m_delta[1] = m_objectOther->getTemporaryPosition(m_particleIDs[2]) - m_objectOther->getTemporaryPosition(m_particleIDs[1]);
    m_delta[2] = m_objectOther->getTemporaryPosition(m_particleIDs[3]) - m_objectOther->getTemporaryPosition(m_particleIDs[1]);
    //m_delta[1] = m_objectOther->getTemporaryPosition(m_particleIDs[3]) - m_objectOther->getTemporaryPosition(m_particleIDs[1]);
    //m_delta[2] = m_objectOther->getTemporaryPosition(m_particleIDs[2]) - m_objectOther->getTemporaryPosition(m_particleIDs[1]);

    for (unsigned i = 0; i < 4; i++)
    {
      m_derivatives[i].fill(0.0);
    }
  }

  bool PointTriangleCollisionConstraint::updateDenominator()
  {
    //std::cout << "C" << m_Constraint[0] << ", eps = " << EPS << std::endl;
    if (m_Constraint[0] > EPS)
      return false;

    double coeff[4];

    Vec3d n = m_delta[1].cross(m_delta[2]);
    coeff[1] = n.dot(m_delta[1].cross(m_delta[0])) / (n.dot(n));
    coeff[2] = n.dot(m_delta[0].cross(m_delta[2])) / (n.dot(n));
    coeff[3] = (1.0 - coeff[1] - coeff[2]);
    //coeff[1] = -n.dot(m_delta[1].cross(m_delta[0])) / (n.dot(n));
    //coeff[2] = -n.dot(m_delta[0].cross(m_delta[2])) / (n.dot(n));
    //coeff[3] = -(1.0 - coeff[1] - coeff[2]);
    coeff[0] = 1.0;

    if (coeff[1] < 0 || coeff[2] < 0 || coeff[3] < 0)
    {
    //if (coeff[1] > 0 || coeff[2] > 0 || coeff[3] > 0)
    //{
      return false;
    }

    for (unsigned i = 0; i < 4; i++)
    {
      //m_object->setParticleCollided(m_particleIDs[i]);
      if (i < m_sizeCol[0])
        m_object->setParticleCollided(m_particleIDs[i]);
      else
        m_objectOther->setParticleCollided(m_particleIDs[i]);
    }

    n.normalize();
    for (unsigned i = 0; i < 4; i++)
    {
      for (unsigned j = 0; j < 3; j++)
      {
        m_derivatives[i][j] = coeff[i] * n[j];
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

  void PointTriangleCollisionConstraint::computeForce()
  {
    Vec3d force;
    double coeff;
    for (unsigned i = 0; i < 4; i++)
    {
        coeff = m_lamdaList[0] * 0.5;
        for (unsigned j = 0; j < 3; j++)
        {
            force[j] = coeff * m_derivatives[i][j];
        }
    }
#pragma omp critical
    for (unsigned i = 0; i < 4; i++)
    {
      //m_object->addToConstraintForce(m_particleIDs[i], m_invMass[i] * force);
      if (i < m_sizeCol[0])
        m_object->addToConstraintForce(m_particleIDs[i], m_ID, m_invMass[i] * force);
      else
        m_objectOther->addToConstraintForce(m_particleIDs[i], m_ID, m_invMass[i] * force);
    }
  }

  void PointTriangleCollisionConstraint::integrate(bool p_all)
  {
  }

  void PointTriangleCollisionConstraint::initShapeFunction()
  {
  }

}
