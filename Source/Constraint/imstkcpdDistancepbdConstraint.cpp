#include "imstkcpdDistancepbdConstraint.h"

namespace cpd
{
  int DistancepbdConstraint::m_idGenerator = 0;

  void DistancepbdConstraint::initConstraint(ParticleObjectPtr p_object, unsigned p_idx1, unsigned p_idx2, const double k)
  {
    m_ID = m_idGenerator++;
    m_stiffness = k;

    std::array<unsigned, 2> idx;
    idx[0] = p_idx1;
    idx[1] = p_idx2;

    ConstraintTemplate::initConstraint(p_object, idx, m_stiffness);

  }

  void DistancepbdConstraint::resolve()
  {
    m_delta = m_object->getTemporaryPosition(m_particleIDs[0]) - m_object->getTemporaryPosition(m_particleIDs[1]);

    m_Constraint[0] = m_delta.norm() - m_volume;

    updateDenominator();

    if (m_isCompliant)
      m_Constraint += m_elasticityMatrix * m_lamdaList;

    m_deltaLamdaList = -m_deno * m_Constraint;
    m_lamdaList += m_deltaLamdaList;

    computeForce();
  }

  bool DistancepbdConstraint::updateDenominator()
  {
    Vec3d n = m_delta.normalized();

    for (unsigned i = 0; i < 3; i++)
    {
      m_derivatives[0][i] = n[i];
      m_derivatives[1][i] = -n[i];
    }

    double cmc = (m_invMass[0] + m_invMass[1])*(n.dot(n));

    if (m_isCompliant)
      cmc += m_elasticityMatrix(0, 0);

    if (abs(cmc) < EPS)
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
      m_deno(0, 0) = 1.0 / cmc;
    }

    return true;
  }

  void DistancepbdConstraint::computeForce()
  {
    Vec3d force;
    double coeff;

    for (unsigned i = 0; i < m_size; i++)
    {
      coeff = m_deltaLamdaList[0];
      for (unsigned j = 0; j < 3; j++)
      {
        force[j] = coeff * m_derivatives[i][j];
      }

      if (!m_isCompliant)
        force *= m_stiffness;

      m_object->addToTemporaryPosition(m_particleIDs[i], m_invMass[i] * force);
    }
  }

  void DistancepbdConstraint::setPBDstiffness(float p_stiffness)
  {
    m_stiffness = p_stiffness;
    m_isCompliant = false;
  }

  ConstraintBase::ConstraintType DistancepbdConstraint::getType() const
  {
    if (m_isCompliant)
      return ConstraintType(ConstraintBase::Type::Distance, ConstraintBase::SubType::XPBD);
    else
      return ConstraintType(ConstraintBase::Type::Distance, ConstraintBase::SubType::PBD);
  }

}