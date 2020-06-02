#include "imstkcpdCSTEnergyConstraint.h"

namespace cpd
{
  int CSTEnergyConstraint::m_idGenerator = 0;

  void CSTEnergyConstraint::initConstraint(ParticleObjectPtr p_object, unsigned p_idx1, unsigned p_idx2, unsigned p_idx3, double p_E, double p_mu)
  {
    m_ID = m_idGenerator++;

    std::array<unsigned, 3> idx;
    idx[0] = p_idx1;
    idx[1] = p_idx2;
    idx[2] = p_idx3;

    EnergyConstraint::initConstraint(p_object, idx, p_E, p_mu);
  }

  void CSTEnergyConstraint::initVariables()
  {
    m_elasticityMatrix.setIdentity();
    m_elasticityMatrix(0, 1) = m_PossionRatio;
    m_elasticityMatrix(1, 0) = m_PossionRatio;
    m_elasticityMatrix(2, 2) = (1 - m_PossionRatio) / 2.0;

    m_elasticityMatrix *= m_YoungsModulus / (1.0 - m_PossionRatio * m_PossionRatio);

    m_elasticityMatrix = (m_elasticityMatrix.inverse()).eval() / (m_dt * m_dt);
  }

  void CSTEnergyConstraint::initDerivatives()
  {
    computeDerivative();
    m_updateDerivative = false;
  }

  void CSTEnergyConstraint::computeDerivative()
  {
    if (m_updateDerivative)
    {
      const Vec3d& x0 = m_object->getPosition(m_particleIDs[0]);
      const Vec3d& x1 = m_object->getPosition(m_particleIDs[1]);
      const Vec3d& x2 = m_object->getPosition(m_particleIDs[2]);

      double xx1 = x0.x();
      double yy1 = x0.y();
      double xx2 = x1.x();
      double yy2 = x1.y();
      double xx3 = x2.x();
      double yy3 = x2.y();

      Matrix3d AA;
      AA(0, 0) = 1.0;
      AA(1, 0) = 1.0;
      AA(2, 0) = 1.0;

      AA(0, 1) = xx1;
      AA(1, 1) = xx2;
      AA(2, 1) = xx3;

      AA(0, 2) = yy1;
      AA(1, 2) = yy2;
      AA(2, 2) = yy3;

      m_volume = 0.5*AA.determinant();

      double sqrtArea = sqrt(m_volume);

      double a1, a2, a3, b1, b2, b3, c1, c2, c3;

      a1 = xx2 * yy3 - xx3 * yy2;
      a2 = xx3 * yy1 - xx1 * yy3;
      a3 = xx1 * yy2 - xx2 * yy1;

      b1 = yy2 - yy3;
      b2 = yy3 - yy1;
      b3 = yy1 - yy2;

      c1 = xx3 - xx2;
      c2 = xx1 - xx3;
      c3 = xx2 - xx1;

      m_derivatives[0][0] = 0.5*b1 / sqrtArea;
      m_derivatives[1][0] = 0.5*b2 / sqrtArea;
      m_derivatives[2][0] = 0.5*b3 / sqrtArea;

      m_derivatives[0][1] = 0.5*c1 / sqrtArea;
      m_derivatives[1][1] = 0.5*c2 / sqrtArea;
      m_derivatives[2][1] = 0.5*c3 / sqrtArea;

      m_B.setZero();
      for (unsigned i = 0; i < 3; i++)
      {
        m_B(0, 2 * i + 0) = m_derivatives[i][0];
        m_B(1, 2 * i + 1) = m_derivatives[i][1];

        m_B(2, 2 * i + 0) = m_derivatives[i][1];
        m_B(2, 2 * i + 1) = m_derivatives[i][0];
      }
    }
  }

  void CSTEnergyConstraint::computeForce()
  {
    Vec3d force(0.0, 0.0, 0.0);
    for (unsigned i = 0; i < m_size; i++)
    {
      force[0] = (m_derivatives[i][0] * m_lamdaList[0] + m_derivatives[i][1] * m_lamdaList[2]);
      force[1] = (m_derivatives[i][1] * m_lamdaList[1] + m_derivatives[i][0] * m_lamdaList[2]);
#pragma omp critical
      m_object->addToConstraintForce(m_particleIDs[i], m_ID, m_invMass[i] * force);
    }
  }

  void CSTEnergyConstraint::updateConstraint()
  {
    m_Constraint.setZero();

    Vec3d pt;

    for (unsigned i = 0; i < m_size; i++)
    {
      pt = m_unconstrainedDisplacements[i] + m_deltaDisplacement[i];

      m_Constraint[0] += m_derivatives[i][0] * pt[0];
      m_Constraint[1] += m_derivatives[i][1] * pt[1];

      m_Constraint[2] += m_derivatives[i][1] * pt[0] + m_derivatives[i][0] * pt[1];
    }
  }

  double CSTEnergyConstraint::getResidual()
  {
    evaluate();
    auto c = m_Constraint + m_elasticityMatrix * m_lamdaList;
    return c[0]* c[0] + c[1]* c[1] + c[2]* c[2];
  }

  void CSTEnergyConstraint::reset()
  {
  }

  void CSTEnergyConstraint::integrate(bool p_all)
  {
  }

  void CSTEnergyConstraint::initShapeFunction()
  {
  }

};