#include "imstkcpdEnergyConstraint.h"

namespace cpd
{
  template <int SIZE, int SIZEC, int DIM>
  void EnergyConstraint<SIZE, SIZEC, DIM>::initConstraint(ParticleObjectPtr p_object, std::array<unsigned, SIZE> p_idx, double p_E, double p_mu)
  {
    m_YoungsModulus = p_E;
    m_PossionRatio = p_mu;

    ConstraintTemplate::initConstraint(p_object, p_idx, p_E, p_mu);

    VecXd invMassVec;
    invMassVec.resize(SIZE*DIM);
    for (unsigned i = 0; i < SIZE; i++)
    {
      for (unsigned j = 0; j < DIM; j++)
      {
        invMassVec[i*DIM + j] = m_invMass[i];
      }
    }

    m_invM = invMassVec.asDiagonal();
    m_CMC = m_B * m_invM*m_B.transpose();
    m_deno = (m_CMC + m_elasticityMatrix).inverse();
  }

  template<int SIZE, int SIZEC, int DIM>
  void EnergyConstraint<SIZE, SIZEC, DIM>::getStrain(std::array<double, SIZEC>& p_straint) const
  {
    for (unsigned i = 0; i < SIZEC; i++)
    {
      p_straint[i] = 0.0;
    }

    for (unsigned i = 0; i < SIZE; i++)
    {
      const std::array<double, DIM>& derivative = m_derivatives[i];
      auto& pt = m_object->getTempDisplacement(m_particleIDs[i]);
      double p[3] = { pt[0], pt[1], pt[2] };

      p_straint[0] += derivative[0] * p[0];
      p_straint[1] += derivative[1] * p[1];

      p_straint[2] += derivative[1] * p[0] + derivative[0] * p[1];
    }

    for (unsigned i = 0; i < SIZEC; i++)
    {
      p_straint[i] /= sqrt(m_volume);
    }

  }

  template<int SIZE, int SIZEC, int DIM>
  VecXd EnergyConstraint<SIZE, SIZEC, DIM>::getBI(size_t i)
  {
    VecXd p;
    p.resize(SIZEC);
    for (size_t k = 0; k < SIZEC; k++)
    {
      p[k] = m_B(k, i);
    }
    return p;
  }

  template<int SIZE, int SIZEC, int DIM>
  void EnergyConstraint<SIZE, SIZEC, DIM>::updateUnconstrainedDisplacements()
  {
    for (unsigned i = 0; i < SIZE; i++)
    {
      m_unconstrainedDisplacements[i] = m_object->getTempDisplacement(m_particleIDs[i]);
    }
  }

  template<int SIZE, int SIZEC, int DIM>
  void EnergyConstraint<SIZE, SIZEC, DIM>::updateConstraintNew()
  {
    updateConstraint();
    m_Constraint += m_elasticityMatrix * m_lamdaList;
  }

  template<int SIZE, int SIZEC, int DIM>
  void EnergyConstraint<SIZE, SIZEC, DIM>::updateStiffnessMatrix()
  {

  }
  template<int SIZE, int SIZEC, int DIM>
  void EnergyConstraint<SIZE, SIZEC, DIM>::updateNN()
  {
  }
}