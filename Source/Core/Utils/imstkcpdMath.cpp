#include "imstkcpdMath.h"

namespace cpd
{
  void printVector(const Vec3d& p_vector)
  {
    std::cout << '(' << p_vector.x() << ',' << p_vector.y()
      << ',' << p_vector.z() << ')' << std::endl;
  }

  double dotStdVectord(const std::vector<double>& p_vector1, const std::vector<double>& p_vector2)
  {
    size_t n = p_vector1.size();

    assert(n == p_vector2.size());

    double dot = 0.0;
    for (size_t i = 0; i < n; i++)
    {
      dot += p_vector1[i] * p_vector2[i];
    }

    return dot;
  }

}