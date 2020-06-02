#include <cassert>
#include "imstkcpdQuadrature.h"

namespace cpd
{

  std::vector<Quadrature::QuadratureInfo> Quadrature::s_typesInfo;
  Quadrature::_init Quadrature::s_initializer;

  Quadrature::Quadrature(Type p_type)
  {
    setType(p_type);
  }

  Quadrature::~Quadrature()
  {
  }

  void Quadrature::setType(Type p_type)
  {
    for (auto& info : s_typesInfo)
    {
      if (info.m_type == p_type)
        m_quadratureInfo = info;
    }
  }

  void Quadrature::integrate(const StdVectorOfVec3d& p_points, VectorFunction & p_func, VecXd& p_integral)
  {
    assert(p_points.size() == m_quadratureInfo.m_points[0][0].size());

    int numNode = p_points.size();
    VecXd tempIntegral(p_integral.size());
    VecXd singleIntegral(p_integral.size());
    Vec3d integralPoint;

    // Go over all weights and then all points associated with that weight, evaluate matrix at
    // those points and add their contribution
    p_integral.setZero();
    for (int i = 0; i < m_quadratureInfo.m_sizeOfWeight; i++)
    {
      tempIntegral.setZero();
      for (int j = 0; j < m_quadratureInfo.m_sizeOfWiPoints[i]; j++)
      {
        integralPoint.setZero();
        for (int k = 0; k < numNode; k++)
        {
          integralPoint += m_quadratureInfo.m_points[i][j][k] * p_points[k];
        }
        p_func(singleIntegral, integralPoint);
        tempIntegral += singleIntegral;
      }
      p_integral += m_quadratureInfo.m_weights[i] * tempIntegral;
    }
  }

  Quadrature::QuadratureInfo::QuadratureInfo(Type p_type, std::vector<std::vector<std::vector<double>>> p_points, std::vector<double> p_weights)
  {
    m_type = p_type;

    m_points = p_points;
    m_weights = p_weights;
    m_sizeOfWeight = m_points.size();

    for (unsigned i = 0; i < m_sizeOfWeight; i++)
    {
      m_sizeOfWiPoints.push_back(m_points[i].size());
    }
  }

  Quadrature::_init::_init()
  {
    std::cout << "init() called" << std::endl;

    std::vector<std::vector<std::vector<double>>> points;
    std::vector<double> weights;

    weights = { 1.0 };
    points = { { {0.333333333333333,0.333333333333333,0.333333333333333} } };
    s_typesInfo.push_back(QuadratureInfo(Type::TRIA1, points, weights));

    weights = { 0.333333333333333 };
    points = { permute33(0.666666666666667, 0.166666666666667) };
    s_typesInfo.push_back(QuadratureInfo(Type::TRIA3, points, weights));

    weights = { 1.0 };
    points = { { {0.25,0.25,0.25,0.25} } };
    s_typesInfo.push_back(QuadratureInfo(Type::TETRA1, points, weights));

    weights = { 0.25 };
    points = { permute44(0.585410196624969,0.138196601125011) };
    s_typesInfo.push_back(QuadratureInfo(Type::TETRA4, points, weights));

    weights = { -0.8, 0.45 };
    points = { { { 0.25,0.25,0.25,0.25 } },permute44(0.5,00.166666666666667) };
    s_typesInfo.push_back(QuadratureInfo(Type::TETRA5, points, weights));
  }


  std::vector<std::vector<double>> Quadrature::_init::permute33(double a, double b)
  {
    return { {a,b,b},{b,a,b},{b,b,a} };
  }

  std::vector<std::vector<double>> Quadrature::_init::permute36(double a, double b, double c)
  {
    return { {a,b,c},{a,c,b},{b,c,a},{b,a,c},{c,b,a},{c,a,b} };
  }

  std::vector<std::vector<double>> Quadrature::_init::permute44(double a, double b)
  {
    return { {a,b,b,b},{b,a,b,b},{b,b,a,b},{b,b,b,a} };
  }

  std::vector<std::vector<double>> Quadrature::_init::permute46(double a, double b)
  {
    return { {a,a,b,b},{a,b,a,b},{a,b,b,a},{b,b,a,a},{b,a,b,a},{b,a,a,b} };
  }

}