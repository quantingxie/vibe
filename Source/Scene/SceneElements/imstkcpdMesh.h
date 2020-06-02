#ifndef CPDMESH_H
#define CPDMESH_H

#include <memory>
#include <array>

#include "imstkcpdMath.h"

namespace cpd
{
  class Mesh
  {
  public:
    using TriangleArray = std::array<size_t, 3>;
    using QuadArray = std::array<size_t, 4>;
    using TetraArray = std::array<size_t, 4>;

  public:
    Mesh() = default;
    Mesh(StdVectorOfVec3d p_vertexPositions/*, std::vector<TriangleArray> p_trianglesVertices*/)
      :m_vertexPositions(p_vertexPositions)/*, m_trianglesVertices(p_trianglesVertices)*/ {}

    const int getNumVertices() const;
    const StdVectorOfVec3d&	getVertexPositions() const;
    const Vec3d& getVertexPosition(const size_t& p_vertID) const;
    void setVertexPosition(int p_ID, const Vec3d& p_pos);
    void scale(double p_scale);

  protected:
    StdVectorOfVec3d m_vertexPositions;
  };
  using MeshPtr = std::shared_ptr<Mesh>;

  class SurfaceMesh :public Mesh
  {
  public:
    SurfaceMesh(): Mesh(){}
    SurfaceMesh(StdVectorOfVec3d p_vertexPositions, std::vector<TriangleArray> p_trianglesVertices) :Mesh(p_vertexPositions), m_trianglesVertices(p_trianglesVertices) {}
    const std::vector<TriangleArray>& getTrianglesVertices() const;
  private:
    std::vector<TriangleArray> m_trianglesVertices;
  };
  using SurfaceMeshPtr = std::shared_ptr<SurfaceMesh>;

  class VolumeMesh :public Mesh
  {
  public:
    VolumeMesh() : Mesh() {}
    VolumeMesh(StdVectorOfVec3d p_vertexPositions) :Mesh(p_vertexPositions) {}
    virtual const std::vector<TetraArray>& getElementVertices() const = 0;
    virtual const TetraArray& getElementVertices(const size_t& tetId) const = 0;
    void setWeightsArray(std::vector<size_t> p_verticeTetraIDs, std::vector<std::array<double, 4>> p_verticeWeights);
    void setWeights(size_t p_verticeTetraId, std::array<double, 4> p_verticeWeight);
    void setAttachedSurfMesh(SurfaceMeshPtr p_surfMesh) { m_attachedSurfaceMesh = p_surfMesh; }
    SurfaceMeshPtr getAttachedSurfMesh() { return m_attachedSurfaceMesh; }
  private:
    SurfaceMeshPtr m_attachedSurfaceMesh = nullptr;
    std::vector<std::array<double, 4>> m_verticeWeights;
    std::vector<size_t> m_verticeTetraIDs;
  };
  using VolumeMeshPtr = std::shared_ptr<VolumeMesh>;

  class TetrahedronMesh :public VolumeMesh
  {
  public:
    TetrahedronMesh() : VolumeMesh() {}
    TetrahedronMesh(StdVectorOfVec3d p_vertexPositions, std::vector<TetraArray> p_tetrahedraVertices) :VolumeMesh(p_vertexPositions), m_tetrahedronVertices(p_tetrahedraVertices) {}
    
    virtual const std::vector<TetraArray>& getElementVertices() const override;
    virtual const TetraArray& getElementVertices(const size_t& tetId) const override;

    //const std::vector<TetraArray>& getTetrahedronVertices() const;
    //const TetraArray& getTetrahedronVertices(const size_t& tetId) const;

    size_t getNumTetrahedron() const;

  private:
    std::vector<TetraArray> m_tetrahedronVertices;
  };
  using TetrahedronMeshPtr = std::shared_ptr<TetrahedronMesh>;

}

#endif // !CPDMESH_H

