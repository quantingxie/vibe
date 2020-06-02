#include "imstkcpdCSVpbdConstraint.h"

namespace cpd
{
  int CSVpbdConstraint::m_idGenerator = 0;

  ConstraintBase::ConstraintType CSVpbdConstraint::getType() const
  {
    if (m_isCompliant)
      return ConstraintType(ConstraintBase::Type::CST3D, ConstraintBase::SubType::XPBD);
    else
      return ConstraintType(ConstraintBase::Type::CST3D, ConstraintBase::SubType::PBD);
  }

  void CSVpbdConstraint::initConstraint(ParticleObjectPtr p_object, unsigned p_idx1, unsigned p_idx2, unsigned p_idx3, unsigned p_idx4, double p_E, double p_mu)
  {
    m_ID = m_idGenerator++;

    std::array<unsigned, 4> idx;
    idx[0] = p_idx1;
    idx[1] = p_idx2;
    idx[2] = p_idx3;
    idx[3] = p_idx4;

    EnergyConstraint::initConstraint(p_object, idx, p_E, p_mu);
  }

  void CSVpbdConstraint::resolve()
  {
    updateConstraint();

    if (m_isCompliant)
      m_Constraint += m_elasticityMatrix * m_lamdaList;

    m_lamdaList += -m_deno * m_Constraint;

    computeForce();
  }

  void CSVpbdConstraint::computeForce()
  {
    Vec3d force(0.0, 0.0, 0.0);
    for (int i = 0; i < m_size; i++)
    {
      force[0] = (m_derivatives[i][0] * m_lamdaList[0] + m_derivatives[i][1] * m_lamdaList[3] + m_derivatives[i][2] * m_lamdaList[5]);
      force[1] = (m_derivatives[i][1] * m_lamdaList[1] + m_derivatives[i][0] * m_lamdaList[3] + m_derivatives[i][2] * m_lamdaList[4]);
      force[2] = (m_derivatives[i][2] * m_lamdaList[2] + m_derivatives[i][1] * m_lamdaList[4] + m_derivatives[i][0] * m_lamdaList[5]);

      if (!m_isCompliant)
        force *= m_stiffness;
      m_object->addToTemporaryPosition(i, force);

    }
  }

  void CSVpbdConstraint::updateConstraint()
  {
    m_Constraint.setZero();

    Vec3d pt;

    for (int i = 0; i < 4; i++)
    {
      pt = m_object->getTempDisplacement(m_particleIDs[i]);

      m_Constraint[0] += m_derivatives[i][0] * pt[0];
      m_Constraint[1] += m_derivatives[i][1] * pt[1];
      m_Constraint[2] += m_derivatives[i][2] * pt[2];

      m_Constraint[3] += m_derivatives[i][1] * pt[0] + m_derivatives[i][0] * pt[1];
      m_Constraint[4] += m_derivatives[i][2] * pt[1] + m_derivatives[i][1] * pt[2];
      m_Constraint[5] += m_derivatives[i][2] * pt[0] + m_derivatives[i][0] * pt[2];
    }
  }

  void CSVpbdConstraint::setPBDstiffness(float p_stiffness)
  {
    m_stiffness = p_stiffness;
    m_isCompliant = false;
  }

}