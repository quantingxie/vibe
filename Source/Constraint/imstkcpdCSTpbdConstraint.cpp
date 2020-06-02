#include "imstkcpdCSTpbdConstraint.h"

namespace cpd
{
  int CSTpbdConstraint::m_idGenerator = 0;

  void CSTpbdConstraint::initConstraint(ParticleObjectPtr p_object, unsigned p_idx1, unsigned p_idx2, unsigned p_idx3, double p_E, double p_mu)
  {
    m_ID = m_idGenerator++;

    std::array<unsigned, 3> idx;
    idx[0] = p_idx1;
    idx[1] = p_idx2;
    idx[2] = p_idx3;

    EnergyConstraint::initConstraint(p_object, idx, p_E, p_mu);
  }

  void CSTpbdConstraint::resolve()
  {
    updateConstraint();

    if (m_isCompliant)
      m_Constraint += m_elasticityMatrix * m_lamdaList;

    m_deltaLamdaList = -m_deno * m_Constraint;

    m_lamdaList += m_deltaLamdaList;

    computeForce();

  }

  void CSTpbdConstraint::updateConstraint()
  {
    m_Constraint.setZero();

    Vec3d pt;

    for (unsigned i = 0; i < 3; i++)
    {
      pt = m_object->getTempDisplacement(m_particleIDs[i]);

      m_Constraint[0] += m_derivatives[i][0] * pt[0];
      m_Constraint[1] += m_derivatives[i][1] * pt[1];

      m_Constraint[2] += m_derivatives[i][1] * pt[0] + m_derivatives[i][0] * pt[1];
    }
  }

  void CSTpbdConstraint::computeForce()
  {
    Vec3d force(0.0, 0.0, 0.0);
    for (unsigned i = 0; i < m_size; i++)
    {
      force[0] = (m_derivatives[i][0] * m_deltaLamdaList[0] + m_derivatives[i][1] * m_deltaLamdaList[2]);
      force[1] = (m_derivatives[i][1] * m_deltaLamdaList[1] + m_derivatives[i][0] * m_deltaLamdaList[2]);

      if (!m_isCompliant)
        force *= m_stiffness;

      m_object->addToTemporaryPosition(m_particleIDs[i], m_invMass[i] * force);
    }
  }

  void CSTpbdConstraint::setPBDstiffness(float p_stiffness)
  {
    m_stiffness = p_stiffness;
    m_isCompliant = false;
  }


  ConstraintBase::ConstraintType CSTpbdConstraint::getType() const
  {
    if (m_isCompliant)
      return ConstraintType(ConstraintBase::Type::CST2D, ConstraintBase::SubType::XPBD);
    else
      return ConstraintType(ConstraintBase::Type::CST2D, ConstraintBase::SubType::PBD);
  }

};