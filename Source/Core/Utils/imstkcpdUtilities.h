#ifndef CPDUTILITIES_H
#define CPDUTILITIES_H

#include "imstkcpdMath.h"
#include "imstkcpdConstants.h"
#include <fstream>
#include <iomanip>

//#include "mex.h"
#include <Eigen/Sparse>
#include <type_traits>
#include <limits>
//#include "mat.h"

namespace cpd
{
	const double max_d = 1e16;
	const double min_d = 1e-16;


  class BoundingBox {
  public:
    BoundingBox() :m_min(Vec3d(max_d, max_d, max_d)), m_max(Vec3d(-max_d, -max_d, -max_d)) {}
    BoundingBox(const Vec3d& p_min, const Vec3d& p_max, bool p_prechecked = false);
    BoundingBox(const Vec3d& p_point): m_min(p_point), m_max(p_point){}

    /// Check if min and max nodes are correctly selected
    bool check();
    /// Set correct max and min
    void computeMaxMin();
    void setValue(const Vec3d& p_min, const Vec3d& p_max);
    void extend(float p_prox);

    const Vec3d& getMax() const { return m_max; }
    const Vec3d& getMin() const { return m_min; }

    void print();

  private:
    Vec3d m_min;
    Vec3d m_max;
  };

  /// Test if interval [a,b] and [c,d] intersect
  bool isOverlap(const double& a, const double& b, const double& c, const double& d);

  bool testAABBToAABB(const BoundingBox& aabb1, const BoundingBox& aabb2);
  bool testLineToLineAABB(const Vec3d& v1, const Vec3d& v2, const Vec3d& v3, const Vec3d& v4, const double& prox1, const double& prox2);
  bool testPointToTriAABB(const Vec3d& v1, const Vec3d& v2, const Vec3d& v3, const Vec3d& v4, const double& prox1, const double& prox2);

  void writeVectorMatlabPlot(Eigen::VectorXd& u, const char* fileName);
  void writeVectorVectorMatlabPlot(Eigen::VectorXd& u, Eigen::VectorXd& v, const char* fileName);
  //void writeVectorMatlab(Eigen::VectorXd& u, const char* fileName);
  //void writeVectorMatlab(Eigen::VectorXd& u, const char* fileName);




  //using namespace Eigen;

  //typedef SparseMatrix<double, ColMajor, std::make_signed<mwIndex>::type> MatlabSparse;


  //Map<MatlabSparse >
  //  matlab_to_eigen_sparse(const mxArray * mat);
  ////{
  ////  mxAssert(mxGetClassID(mat) == mxDOUBLE_CLASS,
  ////    "Type of the input matrix isn't double");
  ////  mwSize     m = mxGetM(mat);
  ////  mwSize     n = mxGetN(mat);
  ////  mwSize    nz = mxGetNzmax(mat);
  ////  /*Theoretically fails in very very large matrices*/
  ////  mxAssert(nz <= std::numeric_limits< std::make_signed<mwIndex>::type>::max(),
  ////    "Unsupported Data size."
  ////  );
  ////  double  * pr = mxGetPr(mat);
  ////  MatlabSparse::StorageIndex* ir = reinterpret_cast<MatlabSparse::StorageIndex*>(mxGetIr(mat));
  ////  MatlabSparse::StorageIndex* jc = reinterpret_cast<MatlabSparse::StorageIndex*>(mxGetJc(mat));
  ////  Map<MatlabSparse> result(m, n, nz, jc, ir, pr);
  ////  return result;
  ////}

  //mxArray*
  //  eigen_to_matlab_sparse(const Ref<const MatlabSparse, StandardCompressedFormat>& mat);
  ///*{
  //  mxArray * result = mxCreateSparse(mat.rows(), mat.cols(), mat.nonZeros(), mxREAL);
  //  const MatlabSparse::StorageIndex* ir = mat.innerIndexPtr();
  //  const MatlabSparse::StorageIndex* jc = mat.outerIndexPtr();
  //  const double* pr = mat.valuePtr();

  //  mwIndex * ir2 = mxGetIr(result);
  //  mwIndex * jc2 = mxGetJc(result);
  //  double  * pr2 = mxGetPr(result);

  //  for (mwIndex i = 0; i < mat.nonZeros(); i++) {
  //    pr2[i] = pr[i];
  //    ir2[i] = ir[i];
  //  }
  //  for (mwIndex i = 0; i < mat.cols() + 1; i++) {
  //    jc2[i] = jc[i];
  //  }
  //  return result;
  //}*/

}

#endif // !CPDUTILITIES_H

