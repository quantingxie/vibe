#ifndef  CPDMATH_H
#define  CPDMATH_H

#include <Eigen/Core>
#include <Eigen/Geometry>
#include <Eigen/Sparse>
#include <Eigen/StdVector>

#include <iostream>

#include "imstkcpdMathConstants.h"

namespace cpd
{
  using Vec2f = Eigen::Vector2f;
  using Vec2d = Eigen::Vector2d;
  using Vec3f = Eigen::Vector3f;
  using Vec3d = Eigen::Vector3d;
  typedef Eigen::Matrix< double, 6, 1> Vec6d; 
  using VecXd = Eigen::VectorXd;
  using VecXi = Eigen::VectorXi;
  using StdVectorOfVec3f = std::vector<Vec3f, Eigen::aligned_allocator<Vec3f>>;
  using StdVectorOfVec3d = std::vector<Vec3d, Eigen::aligned_allocator<Vec3d>>;
  using Matrix2f = Eigen::Matrix<float, 2, 2>;
  using Matrix3f = Eigen::Matrix<float, 3, 3>;
  using Matrix4f = Eigen::Matrix<float, 4, 4>;
  using Matrix2d = Eigen::Matrix<double, 2, 2>;
  using Matrix3d = Eigen::Matrix<double, 3, 3>;
  using Matrix4d = Eigen::Matrix<double, 4, 4>;
  using Matrix8d = Eigen::Matrix<double, 8, 8>;
  using Matrix6d = Eigen::Matrix<double, 6, 6>;
  using Matrix12d = Eigen::Matrix<double, 12, 12>;

  // A dynamic size sparse matrix of doubles
  using SparseMatrixf = Eigen::SparseMatrix < float >;
  // A dynamic size sparse matrix of doubles
  using SparseMatrixd = Eigen::SparseMatrix < double, /*Eigen::ColMajor*/ Eigen::RowMajor>;

  void printVector(const Vec3d& p_vector);
  double dotStdVectord(const std::vector<double>& p_vector1, const std::vector<double>& p_vector2);

#define CPD_UP_VECTOR Vec3d(0.0, 0.0, 1.0);
#define CPD_DOWN_VECTOR Vec3d(0.0, 0.0, -1.0);
#define CPD_RIGHT_VECTOR Vec3d(1.0, 0.0, 0.0);
#define CPD_LEFT_VECTOR Vec3d(-1.0, 0.0, 0.0);
#define CPD_FORWARD_VECTOR Vec3d(0.0, 1.0, 0.0);
#define CPD_BACKWARD_VECTOR Vec3d(0.0, -1.0, 0.0);
#define CPD_WORLD_ORIGIN Vec3d(0.0, 0.0, 0.0);

  //const double PI = 3.141592653;

  //const double MAX_D = 1e16;
  //const double MIN_D = 1e-16;
  //const double MAX_F = 1e16;
  //const double MIN_F = 1e-16;
}

#endif // ! CPDMATH_H

