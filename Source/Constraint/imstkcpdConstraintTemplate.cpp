#include "imstkcpdConstraintTemplate.h"

namespace cpd
{
  template<int SIZE, int SIZEC, int DIM>
  void ConstraintTemplate<SIZE, SIZEC, DIM>::initConstraint(ParticleObjectPtr p_object, std::array<unsigned, SIZE> p_idx, double p_para1, double p_para2)
  {
    m_dt = p_object->getTimeStepSize();

    m_object = p_object;

    m_stiffness = p_para1;
    for (unsigned i = 0; i < SIZE; i++)
    {
      m_particleIDs[i] = p_idx[i];

      m_invMass[i] = m_object->getInvMass(p_idx[i]);

      m_deltaDisplacement[i].setZero();
      if (m_invMass[i] > -EPS)
        m_movable[i] = true;
      else
        m_movable[i] = false;
    }

    initVariables();
    initDerivatives();
  }

  template<int SIZE, int SIZEC, int DIM>
  void ConstraintTemplate<SIZE, SIZEC, DIM>::evaluate()
  {
    updateDeltaDisplacement();
    updateConstraint();
  }

  template<int SIZE, int SIZEC, int DIM>
  void ConstraintTemplate<SIZE, SIZEC, DIM>::resolve()
  {
    evaluate();
    if (updateDenominator())
    {
      m_lamdaList -= m_deno * (m_Constraint + m_elasticityMatrix * m_lamdaList);
      computeForce();
    }
  }

  template<int SIZE, int SIZEC, int DIM>
  void ConstraintTemplate<SIZE, SIZEC, DIM>::updateDeltaDisplacement()
  {
    //#pragma omp parallel for num_threads(SIZE)
    for (int i = 0; i < SIZE; i++)
    {
      unsigned idx = m_particleIDs[i];
      if (m_object->isForceUpdated(idx) && m_movable[i])
      {
        m_deltaDisplacement[i] = m_object->getConstraintForce(idx);
      }
      else
      {
        m_deltaDisplacement[i].setZero();
      }
    }
  }

  template <int SIZE, int SIZEC, int DIM>
  void ConstraintTemplate<SIZE, SIZEC, DIM>::clearLamda()
  {
    m_lamdaList.setZero();
  }

  template<int SIZE, int SIZEC, int DIM>
  void ConstraintTemplate<SIZE, SIZEC, DIM>::clearDeltaDisplacement()
  {
    for (auto& d : m_deltaDisplacement)
    {
      d.fill(0.0);
    }
  }

  template<int SIZE, int SIZEC, int DIM>
  void ConstraintTemplate<SIZE, SIZEC, DIM>::writeResults(std::string p_fileName)
  {
    std::array<double, SIZEC> straint;
    getStrain(straint);

    std::ofstream mfile;
    mfile.open(p_fileName, std::ofstream::out | std::ofstream::app);
    if (!mfile.is_open())
    {
      std::cout << "Unable to create or open file.";
      return;
    }

    mfile << m_ID << ' ';
    for (double v : straint)
    {
      mfile << v << ' ';
    }
    mfile << ';' << std::endl;
    mfile.close();
  }

  template<int SIZE, int SIZEC, int DIM>
  VecXd ConstraintTemplate<SIZE, SIZEC, DIM>::evaluateNew(const VecXd& p_lambda)
  {
    // constraintForces can be computed once at the beginning.
    for (size_t i = 0; i < m_sizeC; i++)
    {
      m_lamdaList(i, 0) = p_lambda[i];
    }
    
    //clearConstraintForces();
    computeForce();

    for (int i = 0; i < SIZE; i++)
    {
      if ((m_movable[i]))
      {
        m_deltaDisplacement[i] = m_object->getConstraintForce(m_particleIDs[i]);
      }
    }

    updateConstraintNew();

    return m_Constraint;

  }

  template<int SIZE, int SIZEC, int DIM>
  void ConstraintTemplate<SIZE, SIZEC, DIM>::setLambda(VecXd & p_lambda)
  {
    for (size_t i = 0; i < m_sizeC; i++)
    {
      m_lamdaList(i, 0) = p_lambda[i];
    }
  }

  template<int SIZE, int SIZEC, int DIM>
  VecXd ConstraintTemplate<SIZE, SIZEC, DIM>::getDeno()
  {
    VecXd vec;
    vec.resize(m_sizeC*m_sizeC);
    for (size_t i = 0; i < m_sizeC; i++)
      for (size_t j = 0; j < m_sizeC; j++)
        vec[i*m_sizeC + j] = m_deno(i, j);
    return vec;
  }

  template<int SIZE, int SIZEC, int DIM>
  void ConstraintTemplate<SIZE, SIZEC, DIM>::clearConstraintForces()
  {
    for (unsigned i = 0; i < m_size; i++)
    {
      m_object->clearConstraintForce(m_particleIDs[i]);
    }
  }

  template<int SIZE, int SIZEC, int DIM>
  VecXd ConstraintTemplate<SIZE, SIZEC, DIM>::getBI(size_t i)
  {
    return VecXd();
  }

  template<int SIZE, int SIZEC, int DIM>
  VecXi ConstraintTemplate<SIZE, SIZEC, DIM>::getpParticleIDs()
  {
    VecXi p;
    p.resize(SIZE);
    for (size_t i = 0; i < SIZE; i++)
    {
      p[i] = m_particleIDs[i];
    }
    return p;
  }

  template<int SIZE, int SIZEC, int DIM>
  VecXd ConstraintTemplate<SIZE, SIZEC, DIM>::getElasticityMatrix(size_t i)
  {
    VecXd p;
    p.resize(SIZEC);
    for (size_t k = 0; k < SIZEC; k++)
    {
      p[k] = m_elasticityMatrix(k, i);
    }
    return p;
  }

}