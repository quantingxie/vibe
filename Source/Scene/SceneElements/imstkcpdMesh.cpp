#include "imstkcpdMesh.h"

namespace cpd
{
  const int Mesh::getNumVertices() const
  {
    return m_vertexPositions.size();
  }

  const StdVectorOfVec3d & Mesh::getVertexPositions() const
  {
    return m_vertexPositions;
  }

  const Vec3d& Mesh::getVertexPosition(const size_t & p_vertID) const
  {
    return this->getVertexPositions().at(p_vertID);
  }

  void Mesh::setVertexPosition(int p_ID, const Vec3d& p_pos)
  {
    m_vertexPositions[p_ID] = p_pos;
  }

  void Mesh::scale(double p_scale)
  {
      for (auto& v : m_vertexPositions)
      {
          v *= p_scale;
      }
  }

  const std::vector<Mesh::TriangleArray>& SurfaceMesh::getTrianglesVertices() const
  {
    return m_trianglesVertices;
  }

  const std::vector<Mesh::TetraArray>& TetrahedronMesh::getElementVertices() const
  {
    return m_tetrahedronVertices;
  }

  const Mesh::TetraArray & TetrahedronMesh::getElementVertices(const size_t & tetId) const
  {
    return m_tetrahedronVertices[tetId];
  }

  size_t TetrahedronMesh::getNumTetrahedron() const
  {
    return m_tetrahedronVertices.size();
  }

  void VolumeMesh::setWeightsArray(std::vector<size_t> p_verticeTetraIDs, std::vector<std::array<double, 4>> p_verticeWeights)
  {
    m_verticeTetraIDs = p_verticeTetraIDs;
    m_verticeWeights = p_verticeWeights;    
  }

  void VolumeMesh::setWeights(size_t p_verticeTetraID, std::array<double, 4> p_verticeWeight)
  {
    m_verticeTetraIDs.push_back(p_verticeTetraID);
    m_verticeWeights.push_back(p_verticeWeight);
  }

}