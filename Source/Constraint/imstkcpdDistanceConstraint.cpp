#include "imstkcpdDistanceConstraint.h"

namespace cpd
{
  int DistanceConstraint::m_idGenerator = 0;

  void DistanceConstraint::initConstraint(ParticleObjectPtr p_object, unsigned p_idx1, unsigned p_idx2, const double k)
  {
    m_ID = m_idGenerator++;
    m_stiffness = k;

    std::array<unsigned, 2> idx;
    idx[0] = p_idx1;
    idx[1] = p_idx2;

    ConstraintTemplate::initConstraint(p_object, idx, m_stiffness);
  }

  void DistanceConstraint::initVariables()
  {
    m_elasticityMatrix(0, 0) = 1.0 / (m_stiffness*m_dt * m_dt);
  }

  void DistanceConstraint::initDerivatives()
  {
    computeDerivative();

    m_volume = m_delta.norm();
  }

  bool DistanceConstraint::updateDenominator()
  {
    Vec3d n = (m_delta + m_deltaDisplacement[0] - m_deltaDisplacement[1]).normalized();

    double cmc = (m_invMass[0] + m_invMass[1])*(n[0] * m_derivatives[0][0] + n[1] * m_derivatives[0][1]
      + n[2] * m_derivatives[0][2]);
    double temp = cmc + m_elasticityMatrix(0, 0);

    if (abs(temp) < EPS)
    {
      if (abs(m_elasticityMatrix(0.0)) > EPS)
      {
        m_deno(0, 0) = 1.0 / m_elasticityMatrix(0, 0);
        m_lamdaList.fill(0.0);
      }
      else
      {
        m_deno(0, 0) = 0.0;
        m_lamdaList.fill(0.0);
      }
    }
    else
    {
      m_deno(0, 0) = 1.0 / temp;
    }

    return true;
  }

  void DistanceConstraint::computeDerivative()
  {
    m_delta = m_object->getTemporaryPosition(m_particleIDs[0]) - m_object->getTemporaryPosition(m_particleIDs[1]);
    Vec3d n = m_delta.normalized();

    for (unsigned i = 0; i < 3; i++)
    {
      m_derivatives[0][i] = n[i];
      m_derivatives[1][i] = -n[i];
    }
  }

  void DistanceConstraint::computeForce()
  {
    Vec3d force[2];
    double coeff;
    for (unsigned i = 0; i < m_size; i++)
    {
        coeff = m_lamdaList[0];
        for (unsigned j = 0; j < 3; j++)
        {
            force[i][j] = coeff * m_derivatives[i][j];
        }
    }
#pragma omp critical
	{
		for (unsigned i = 0; i < m_size; i++)
		{
			m_object->addToConstraintForce(m_particleIDs[i], m_ID, m_invMass[i] * force[i]);
		}
	}
  }

  void DistanceConstraint::updateConstraint()
  {
    const Vec3d& n = m_delta + m_deltaDisplacement[0] - m_deltaDisplacement[1];
    m_Constraint[0] = n.norm() - m_volume;
  }

  void DistanceConstraint::reset()
  {
  }

  double DistanceConstraint::getResidual()
  {
    evaluate();
    return m_Constraint[0] + m_elasticityMatrix[0] * m_lamdaList[0];
  }

  void DistanceConstraint::integrate(bool p_all)
  {
  }

  void DistanceConstraint::initShapeFunction()
  {

  }

}