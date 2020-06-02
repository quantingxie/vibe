#include "imstkcpdUtilities.h"

namespace cpd
{
  BoundingBox::BoundingBox(const Vec3d& p_min, const Vec3d& p_max, bool p_prechecked) :m_min(p_min), m_max(p_max)
  {
    if (!p_prechecked)
        if (!check())
            computeMaxMin();
  }

  bool BoundingBox::check()
  {
    return ((m_min.x() < m_max.x()) && (m_min.y() < m_max.y()) && (m_min.z() < m_max.z()));
  }

  void BoundingBox::computeMaxMin()
  {
    float min_x, min_y, min_z, max_x, max_y, max_z;

    min_x = std::min(m_min.x(), m_max.x());
    min_y = std::min(m_min.y(), m_max.y());
    min_z = std::min(m_min.z(), m_max.z());

    max_x = std::max(m_min.x(), m_max.x());
    max_y = std::max(m_min.y(), m_max.y());
    max_z = std::max(m_min.z(), m_max.z());

    m_min = Vec3d(min_x, min_y, min_z);
    m_max = Vec3d(max_x, max_y, max_z);
  }

  void BoundingBox::setValue(const Vec3d& p_min, const Vec3d& p_max)
  {
    m_min = p_min;
    m_max = p_max;

    if (!check())
      computeMaxMin();

  }

  void BoundingBox::extend(float p_prox)
  {
    if (p_prox < EPS)
      p_prox = abs(p_prox);
    Vec3d ext(p_prox, p_prox, p_prox);
    m_min -= ext;
    m_max += ext;
  }

  void BoundingBox::print()
  {
    printVector(m_min);
    printVector(m_max);
  }


  bool isOverlap(const double& a, const double& b, const double& c, const double& d)
  {
    return ((a <= d && a >= c) || (c <= b && c >= a)) ? true : false;
  }

  bool testAABBToAABB(const BoundingBox& aabb1, const BoundingBox& aabb2)
  {
    const Vec3d& min1 = aabb1.getMin();
    const Vec3d& max1 = aabb1.getMax();
    const Vec3d& min2 = aabb2.getMin();
    const Vec3d& max2 = aabb2.getMax();

    return (isOverlap(min1.x(), max1.x(), min2.x(), max2.x()) &&
      isOverlap(min1.y(), max1.y(), min2.y(), max2.y()) &&
      isOverlap(min1.z(), max1.z(), min2.z(), max2.z()));
  }

  bool testLineToLineAABB(const Vec3d& v1, const Vec3d& v2, const Vec3d& v3, const Vec3d& v4, const double& prox1, const double& prox2)
  {
    BoundingBox aabb1(v1, v2);
    BoundingBox aabb2(v3, v4);

    aabb1.extend(prox1);
    aabb2.extend(prox2);

    return testAABBToAABB(aabb1, aabb2);
  }

  bool testPointToTriAABB(const Vec3d& v1, const Vec3d& v2, const Vec3d& v3, const Vec3d& v4, const double& prox1, const double& prox2)
  {
	  auto& trix = std::minmax({ v2[0],v3[0], v4[0] });
	  auto& triy = std::minmax({ v2[1],v3[1], v4[1] });
	  auto& triz = std::minmax({ v2[2],v3[2], v4[2] });

	  float x = v1[0];
	  float y = v1[1];
	  float z = v1[2];

	  BoundingBox aabb1(Vec3d(x - prox1, y - prox1, z - prox1), Vec3d(x + prox1, y + prox1, z + prox1), true);
	  BoundingBox aabb2(Vec3d(trix.first - prox2, triy.first - prox2, triz.first - prox2),
		  Vec3d(trix.second + prox2, triy.second + prox2, triz.second + prox2), true);

	  return testAABBToAABB(aabb1, aabb2);
  }

  void writeVectorMatlabPlot(Eigen::VectorXd & u, const char * fileName)
  {
    std::ofstream mfile(fileName);

    if (!mfile.is_open())
    {
      std::cout << "Unable to create or open file.";
      return;
    }

    std::string tempstr = std::string(fileName);
    mfile << "u"+ tempstr.substr(0,tempstr.size()-2) +"=[\n";
    for (auto i = 0; i < u.size(); ++i)
    {
      mfile << std::setprecision(std::numeric_limits<long double>::digits10 + 1) << u[i] << "\n";
    }
    mfile << "];\n";

    //mfile << "plot(u" + tempstr.substr(0, tempstr.size() - 2) + ", 'r-', 'LineWidth', 2, 'MarkerSize', 10);\n";
    mfile.close();
  }

  void writeVectorVectorMatlabPlot(Eigen::VectorXd & u, Eigen::VectorXd & v, const char * fileName)
  {
    if (u.size() != v.size())
    {
      std::cout << "Vector size mismatch!";
      return;
    }

    std::ofstream mfile(fileName);

    if (!mfile.is_open())
    {
      std::cout << "Unable to create or open file!";
      return;
    }

    mfile << "u=[\n";
    for (auto i = 0; i < u.size(); ++i)
    {
      mfile << u[i] << "\n";
    }
    mfile << "];\n";

    mfile << "v=[\n";
    for (auto i = 0; i < v.size(); ++i)
    {
      mfile << v[i] << "\n";
    }
    mfile << "];\n";

    mfile << "plot(u, v, 'r-', 'LineWidth', 2, 'MarkerSize', 10);\n";
    mfile.close();
  }

  //using namespace Eigen;

  //typedef SparseMatrix<double, ColMajor, std::make_signed<mwIndex>::type> MatlabSparse;

  //Map<MatlabSparse >
  //  matlab_to_eigen_sparse(const mxArray * mat)
  //{
  //  mxAssert(mxGetClassID(mat) == mxDOUBLE_CLASS,
  //    "Type of the input matrix isn't double");
  //  mwSize     m = mxGetM(mat);
  //  mwSize     n = mxGetN(mat);
  //  mwSize    nz = mxGetNzmax(mat);
  //  /*Theoretically fails in very very large matrices*/
  //  mxAssert(nz <= std::numeric_limits< std::make_signed<mwIndex>::type>::max(),
  //    "Unsupported Data size."
  //  );
  //  double  * pr = mxGetPr(mat);
  //  MatlabSparse::StorageIndex* ir = reinterpret_cast<MatlabSparse::StorageIndex*>(mxGetIr(mat));
  //  MatlabSparse::StorageIndex* jc = reinterpret_cast<MatlabSparse::StorageIndex*>(mxGetJc(mat));
  //  Map<MatlabSparse> result(m, n, nz, jc, ir, pr);
  //  return result;
  //}

  //mxArray*
  //  eigen_to_matlab_sparse(const Ref<const MatlabSparse, StandardCompressedFormat>& mat)
  //{
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

  //  MATFile *pmat;
  //  const char *myFile = "matimport.mat";
  //  pmat = matOpen(myFile, "w");

  //  const char *myDouble = "array5";
  //  auto status = matPutVariable(pmat, myDouble, result);

  //  matClose(pmat);
  //  mxDestroyArray(result);

  //  return result;
  //}

}