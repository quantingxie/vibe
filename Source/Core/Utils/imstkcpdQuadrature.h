#ifndef  CPDQUADRATURE_H
#define  CPDQUADRATURE_H

#include "imstkcpdMath.h"

namespace cpd
{
  using VectorFunction = std::function < void(VecXd&, const VecXd&) >;
  using MatrixFunction = std::function < void(SparseMatrixd&, const VecXd&) >;


  class Quadrature
  {
  public:
    enum class Type
    {
      TRIA1,
      TRIA3,
      TETRA1,
      TETRA4,
      TETRA5
    };

    class QuadratureInfo
    {
    public:
      QuadratureInfo() {}
      QuadratureInfo(Type p_type, std::vector<std::vector<std::vector<double>>> p_points, std::vector<double> p_weights);
    public:
      Type m_type;
      int m_sizeOfWeight;
      std::vector<int> m_sizeOfWiPoints;
      std::vector<std::vector<std::vector<double>>> m_points;
      std::vector<double> m_weights;
    };


  public:
    Quadrature() {}
    Quadrature(Type p_type);
    ~Quadrature();

    void setType(Type p_type);
    Type getType() { return m_quadratureInfo.m_type; }

    void integrate(const StdVectorOfVec3d& p_points, VectorFunction& p_func, VecXd& p_integral);

  private:

  private:
    QuadratureInfo m_quadratureInfo;

    static std::vector<QuadratureInfo> s_typesInfo;

  public:
    class _init
    {
    public:
      _init();
    private:
      std::vector<std::vector<double>> permute33(double a, double b);
      std::vector<std::vector<double>> permute36(double a, double b, double c);
      std::vector<std::vector<double>> permute44(double a, double b);
      std::vector<std::vector<double>> permute46(double a, double b);
    };
  private:
    static _init s_initializer;
  };

}

#endif // ! CPDQUADRATURE_H

