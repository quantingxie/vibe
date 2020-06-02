#include "imstkcpdCSVEnergyConstraint.h"

namespace cpd
{
  int CSVEnergyConstraint::m_idGenerator = 0;

  void CSVEnergyConstraint::initConstraint(ParticleObjectPtr p_object, unsigned p_idx1, unsigned p_idx2, unsigned p_idx3, unsigned p_idx4, double p_E, double p_mu)
  {
    m_ID = m_idGenerator++;

    initShapeFunction();
    m_quadrature.setType(Quadrature::Type::TETRA4);

    m_abcd.resize(4);
    for (auto& v : m_abcd)
      v.resize(4);

    std::array<unsigned, 4> idx;
    idx[0] = p_idx1;
    idx[1] = p_idx2;
    idx[2] = p_idx3;
    idx[3] = p_idx4;

    EnergyConstraint::initConstraint(p_object, idx, p_E, p_mu);
  }

  void CSVEnergyConstraint::initVariables()
  {
    m_elasticityMatrix.setIdentity();

    m_elasticityMatrix(0, 0) = 1.0 - m_PossionRatio;
    m_elasticityMatrix(0, 1) = m_PossionRatio;
    m_elasticityMatrix(0, 2) = m_PossionRatio;

    m_elasticityMatrix(1, 0) = m_PossionRatio;
    m_elasticityMatrix(1, 1) = 1.0 - m_PossionRatio;
    m_elasticityMatrix(1, 2) = m_PossionRatio;

    m_elasticityMatrix(2, 0) = m_PossionRatio;
    m_elasticityMatrix(2, 1) = m_PossionRatio;
    m_elasticityMatrix(2, 2) = 1.0 - m_PossionRatio;

    m_elasticityMatrix(3, 3) = (0.5 - m_PossionRatio);
    m_elasticityMatrix(4, 4) = (0.5 - m_PossionRatio);
    m_elasticityMatrix(5, 5) = (0.5 - m_PossionRatio);

    m_elasticityMatrix *= m_YoungsModulus / ((1.0 + m_PossionRatio)*(1.0 - 2.0 * m_PossionRatio));

    m_elasticityMatrix = (m_elasticityMatrix.inverse()).eval() / (m_dt * m_dt);
  }

  void CSVEnergyConstraint::initDerivatives()
  {
    computeDerivative();
    m_updateDerivative = false;
  }

  void CSVEnergyConstraint::computeDerivative()
  {
    
    if (m_updateDerivative)
    {
      const Vec3d& x0 = m_object->getPosition(m_particleIDs[0]);
      const Vec3d& x1 = m_object->getPosition(m_particleIDs[1]);
      const Vec3d& x2 = m_object->getPosition(m_particleIDs[2]);
      const Vec3d& x3 = m_object->getPosition(m_particleIDs[3]);

      double xx1 = x0.x();
      double yy1 = x0.y();
      double zz1 = x0.z();
      double xx2 = x1.x();
      double yy2 = x1.y();
      double zz2 = x1.z();
      double xx3 = x2.x();
      double yy3 = x2.y();
      double zz3 = x2.z();
      double xx4 = x3.x();
      double yy4 = x3.y();
      double zz4 = x3.z();

      Matrix4d VV;
      VV(0, 0) = 1.0;
      VV(0, 1) = 1.0;
      VV(0, 2) = 1.0;
      VV(0, 3) = 1.0;

      VV(1, 0) = xx1;
      VV(1, 1) = xx2;
      VV(1, 2) = xx3;
      VV(1, 3) = xx4;

      VV(2, 0) = yy1;
      VV(2, 1) = yy2;
      VV(2, 2) = yy3;
      VV(2, 3) = yy4;

      VV(3, 0) = zz1;
      VV(3, 1) = zz2;
      VV(3, 2) = zz3;
      VV(3, 3) = zz4;

      double volume = (1.0 / 6.0)*VV.determinant();

      double sqrtVolume = sqrt(volume);
      double sqrtVolume6 = sqrtVolume * 6;

      double a1, a2, a3, a4, b1, b2, b3, b4, c1, c2, c3, c4, d1, d2, d3, d4;

      a1 = xx2 * yy3*zz4 - xx2 * yy4*zz3 - xx3 * yy2*zz4 + xx3 * yy4*zz2 + xx4 * yy2*zz3 - xx4 * yy3*zz2;
      b1 = yy3 * zz2 - yy2 * zz3 + yy2 * zz4 - yy4 * zz2 - yy3 * zz4 + yy4 * zz3;
      c1 = xx2 * zz3 - xx3 * zz2 - xx2 * zz4 + xx4 * zz2 + xx3 * zz4 - xx4 * zz3;
      d1 = xx3 * yy2 - xx2 * yy3 + xx2 * yy4 - xx4 * yy2 - xx3 * yy4 + xx4 * yy3;

      a2 = xx1 * yy4*zz3 - xx1 * yy3*zz4 + xx3 * yy1*zz4 - xx3 * yy4*zz1 - xx4 * yy1*zz3 + xx4 * yy3*zz1;
      b2 = yy1 * zz3 - yy3 * zz1 - yy1 * zz4 + yy4 * zz1 + yy3 * zz4 - yy4 * zz3;
      c2 = xx3 * zz1 - xx1 * zz3 + xx1 * zz4 - xx4 * zz1 - xx3 * zz4 + xx4 * zz3;
      d2 = xx1 * yy3 - xx3 * yy1 - xx1 * yy4 + xx4 * yy1 + xx3 * yy4 - xx4 * yy3;

      a3 = xx1 * yy2*zz4 - xx1 * yy4*zz2 - xx2 * yy1*zz4 + xx2 * yy4*zz1 + xx4 * yy1*zz2 - xx4 * yy2*zz1;
      b3 = yy2 * zz1 - yy1 * zz2 + yy1 * zz4 - yy4 * zz1 - yy2 * zz4 + yy4 * zz2;
      c3 = xx1 * zz2 - xx2 * zz1 - xx1 * zz4 + xx4 * zz1 + xx2 * zz4 - xx4 * zz2;
      d3 = xx2 * yy1 - xx1 * yy2 + xx1 * yy4 - xx4 * yy1 - xx2 * yy4 + xx4 * yy2;

      a4 = xx1 * yy3*zz2 - xx1 * yy2*zz3 + xx2 * yy1*zz3 - xx2 * yy3*zz1 - xx3 * yy1*zz2 + xx3 * yy2*zz1;
      b4 = yy1 * zz2 - yy2 * zz1 - yy1 * zz3 + yy3 * zz1 + yy2 * zz3 - yy3 * zz2;
      c4 = xx2 * zz1 - xx1 * zz2 + xx1 * zz3 - xx3 * zz1 - xx2 * zz3 + xx3 * zz2;
      d4 = xx1 * yy2 - xx2 * yy1 - xx1 * yy3 + xx3 * yy1 + xx2 * yy3 - xx3 * yy2;

      m_abcd[0][0] = a1;
      m_abcd[0][1] = b1;
      m_abcd[0][2] = c1;
      m_abcd[0][3] = d1;

      m_abcd[1][0] = a2;
      m_abcd[1][1] = b2;
      m_abcd[1][2] = c2;
      m_abcd[1][3] = d2;

      m_abcd[2][0] = a3;
      m_abcd[2][1] = b3;
      m_abcd[2][2] = c3;
      m_abcd[2][3] = d3;

      m_abcd[3][0] = a4;
      m_abcd[3][1] = b4;
      m_abcd[3][2] = c4;
      m_abcd[3][3] = d4;


      m_derivatives[0][0] = b1 / sqrtVolume6;
      m_derivatives[1][0] = b2 / sqrtVolume6;
      m_derivatives[2][0] = b3 / sqrtVolume6;
      m_derivatives[3][0] = b4 / sqrtVolume6;

      m_derivatives[0][1] = c1 / sqrtVolume6;
      m_derivatives[1][1] = c2 / sqrtVolume6;
      m_derivatives[2][1] = c3 / sqrtVolume6;
      m_derivatives[3][1] = c4 / sqrtVolume6;

      m_derivatives[0][2] = d1 / sqrtVolume6;
      m_derivatives[1][2] = d2 / sqrtVolume6;
      m_derivatives[2][2] = d3 / sqrtVolume6;
      m_derivatives[3][2] = d4 / sqrtVolume6;

      m_B.setZero();

      for (unsigned i = 0; i < 4; i++)
      {
        m_B(0, 3 * i + 0) = m_derivatives[i][0];
        m_B(1, 3 * i + 1) = m_derivatives[i][1];
        m_B(2, 3 * i + 2) = m_derivatives[i][2];

        m_B(3, 3 * i + 0) = m_derivatives[i][1];
        m_B(3, 3 * i + 1) = m_derivatives[i][0];
        m_B(4, 3 * i + 1) = m_derivatives[i][2];
        m_B(4, 3 * i + 2) = m_derivatives[i][1];
        m_B(5, 3 * i + 2) = m_derivatives[i][0];
        m_B(5, 3 * i + 0) = m_derivatives[i][2];
      }


      integrate(true);
      updateNN();
    }
  }

  void CSVEnergyConstraint::computeForce()
  {
	  Vec3d force[4];
	  //#pragma omp parallel for num_threads(m_size)
	  //std::cout << "step 4, inside computeForce" << std::endl;
	  //std::cout << omp_get_thread_num() << '/' << omp_get_num_threads() << std::endl;
	  for (int i = 0; i < m_size; i++)
	  {

		  //      std::cout << omp_get_thread_num() << '/' << omp_get_num_threads() << std::endl;

		  force[i][0] = m_derivatives[i][0] * m_lamdaList[0] + m_derivatives[i][1] * m_lamdaList[3] + m_derivatives[i][2] * m_lamdaList[5];
		  force[i][1] = m_derivatives[i][1] * m_lamdaList[1] + m_derivatives[i][0] * m_lamdaList[3] + m_derivatives[i][2] * m_lamdaList[4];
		  force[i][2] = m_derivatives[i][2] * m_lamdaList[2] + m_derivatives[i][1] * m_lamdaList[4] + m_derivatives[i][0] * m_lamdaList[5];
	  }
#pragma omp critical
	  {
		  for (int i = 0; i < m_size; i++)
		  {
			  //std::cout << omp_get_thread_num() << '/' << omp_get_num_threads() << std::endl;
			  m_object->addToConstraintForce(m_particleIDs[i], m_ID, m_invMass[i] * force[i]);
		  }
	  }
	  //std::cout << "step 5, inside computeForce, after critical" << std::endl;
  }

  void CSVEnergyConstraint::updateConstraint()
  {
    m_Constraint.setZero();

    Vec3d pt;

    for (int i = 0; i < m_size; i++)
    {

      pt = m_unconstrainedDisplacements[i] + m_deltaDisplacement[i];


      m_Constraint[0] += m_derivatives[i][0] * pt[0];
      m_Constraint[1] += m_derivatives[i][1] * pt[1];
      m_Constraint[2] += m_derivatives[i][2] * pt[2];

      m_Constraint[3] += m_derivatives[i][1] * pt[0] + m_derivatives[i][0] * pt[1];
      m_Constraint[4] += m_derivatives[i][2] * pt[1] + m_derivatives[i][1] * pt[2];
      m_Constraint[5] += m_derivatives[i][2] * pt[0] + m_derivatives[i][0] * pt[2];

    }


  }

  void CSVEnergyConstraint::reset()
  {
  }

  void CSVEnergyConstraint::updateDeno()
  {
    updateDeltaDisplacement();

    //Eigen::Matrix<double, 6, 24> intN;
    //intN =
    //m_CMC = m_B * m_invM*m_B.transpose();
    //m_deno = (m_CMC + m_elasticityMatrix).inverse();
  }

  void CSVEnergyConstraint::integrate(bool p_all)
  {
    StdVectorOfVec3d points;
    for (size_t i = 0; i < 4; i++)
    {
      points.push_back(m_object->getPosition(m_particleIDs[i]));
    }
    VecXd integral(4);
 
    m_quadrature.integrate(points, m_shapeFunction, integral);

    for (unsigned i = 0; i < 4; i++)
    {
      m_N[i] = integral[i] / 6.0;
    }

  }

  void CSVEnergyConstraint::initShapeFunction()
  {
    m_shapeFunction = [this](VecXd &N, const VecXd &u)
    {
      for (size_t i = 0; i < 4; i++)
      {
        N[i] = m_abcd[i][0] + m_abcd[i][1] * u[0] + m_abcd[i][2] * u[1] + m_abcd[i][3] * u[2];
      }
    };
  }

  void CSVEnergyConstraint::updateNN()
  {
    for (size_t i = 0; i < m_size; i++)
    {
      for (size_t j = 0; j < m_size; j++)
      {
        m_NN[i][j] = m_N[i] * m_N[j];
      }
    }
  }

}