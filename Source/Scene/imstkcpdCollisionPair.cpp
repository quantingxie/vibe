#include "imstkcpdCollisionPair.h"
#include "imstkcpdEdgeEdgeCollisionConstraint.h"
#include "imstkcpdPointTriangleCollisionConstraint.h"

namespace cpd
{
  void CollisionPair::initCollisionConstraint()
  {

  }

  void CollisionPair::resetConstraints()
  {
    for (auto& collisionConstraint : m_collisionConstraints)
      collisionConstraint->reset();
  }

  void CollisionPair::resetLambda()
  {
    for (auto& c : m_collisionConstraints)
    {
      // TODO: collisionConstraint is not child of constraint yet!!
    }
  }

  void CollisionPair::resolve()
  {
      int n = m_collisionConstraints.size();
#pragma omp parallel for/* num_threads(32)*//*schedule(guided, 40)*/
      for (int idx = 0; idx < n; idx++)
      {
          m_collisionConstraints[idx]->resolve();
      }
    m_collisionConstraints.clear();
  }

  void CollisionPair::detectAndResolve()
  {
    if (broadPhaseCollisionDetection())
    {
      narrowPhaseCollisionDetetion();
      resolve();
    }
  }

  bool CollisionPair::broadPhaseCollisionDetection()
  {
    m_object1->updateBoundingBox();
    m_object2->updateBoundingBox();
    auto& aabb1 = m_object1->getBoundingBox();
    auto& aabb2 = m_object2->getBoundingBox();

    return testAABBToAABB(aabb1, aabb2);
  }

  void CollisionPair::narrowPhaseCollisionDetetion()
  {
    auto proximity1 = 1 * m_object1->getParticleProximity();
    auto proximity2 = 1 * m_object2->getParticleProximity();

    auto mesh1 = m_object1->getSurfaceMesh(); // Should avoid using mesh!!!!
    //if (mesh1 == nullptr)
    //  mesh1 = m_object1->getVolumeMesh();
    auto mesh2 = m_object2->getSurfaceMesh();

    const auto nV = m_object1->getParticleCount();
    std::vector<std::vector<bool>> E(nV, std::vector<bool>(nV, 1));

    const auto nV2 = mesh2->getNumVertices();

    // point-triangle
    auto mesh2elements = mesh2->getTrianglesVertices();
    int ne2 = mesh2elements.size();
#pragma omp parallel for/* num_threads(32)*/
    for (int i = 0; i < nV; ++i)
    {
      //#pragma omp critical
      //      std::cout << omp_get_thread_num() << '/' << omp_get_num_threads() << std::endl;
      const Vec3d& p = m_object1->getTemporaryPosition(i);
//#pragma omp parallel for num_threads(2)
      for (int j = 0; j < ne2; ++j)
      {
        SurfaceMesh::TriangleArray& e = mesh2elements[j];
        const Vec3d& p0 = m_object2->getTemporaryPosition(e[0]);
        const Vec3d& p1 = m_object2->getTemporaryPosition(e[1]);
        const Vec3d& p2 = m_object2->getTemporaryPosition(e[2]);

        if (testPointToTriAABB(p, p0, p1, p2, proximity1, proximity2))
        {
          auto c = std::make_unique<PointTriangleCollisionConstraint>();
          c->initConstraint(m_object1, i, m_object2, e[0], e[1], e[2]);
#pragma omp critical
          m_collisionConstraints.push_back(std::move(c));
        }
      }
    }

    if (nV2 > 1)
    {
      //std::cout << "No E-E collision here." << std::endl;
      return;
    }
    else
    {
      std::cout << "nV2 = " << nV2 << std::endl;
    }
    // EE
    if (!mesh1)
    {
      std::vector<SurfaceMesh::TriangleArray> elements2 = mesh2->getTrianglesVertices();
      const Vec3d& P = m_object1->getPosition(0);
      const Vec3d& Q = m_object1->getPosition(1);
      std::vector<std::vector<bool>> E2(nV2, std::vector<bool>(nV2, 1));
//#pragma omp parallel for num_threads(32)
      for (int j = 0; j < elements2.size(); ++j)
      {
        SurfaceMesh::TriangleArray& e = elements2[j];
        const Vec3d& p0 = mesh2->getVertexPosition(e[0]);
        const Vec3d& p1 = mesh2->getVertexPosition(e[1]);
        const Vec3d& p2 = mesh2->getVertexPosition(e[2]);
        if (E2[e[0]][e[1]] && E2[e[1]][e[0]])
        {
          if (testLineToLineAABB(P, Q, p0, p1, proximity1, proximity2))
          {
            auto c = std::make_unique<EdgeEdgeCollisionConstraint>();
            c->initConstraint(m_object1, 0, 1, m_object2, e[0], e[1]);
#pragma omp critical
            m_collisionConstraints.push_back(std::move(c));
            E2[e[0]][e[1]] = 0;
          }
        }
        if (E2[e[1]][e[2]] && E2[e[2]][e[1]])
        {
          if (testLineToLineAABB(P, Q, p1, p2, proximity1, proximity2))
          {
            auto c = std::make_unique<EdgeEdgeCollisionConstraint>();
            c->initConstraint(m_object1, 0, 1, m_object2, e[1], e[2]);
#pragma omp critical
            m_collisionConstraints.push_back(std::move(c));
            E2[e[1]][e[2]] = 0;
          }
        }
        if (E2[e[2]][e[0]] && E2[e[0]][e[2]])
        {
          if (testLineToLineAABB(P, Q, p2, p0, proximity1, proximity2))
          {
            auto c = std::make_unique<EdgeEdgeCollisionConstraint>();
            c->initConstraint(m_object1, 0, 1, m_object2, e[2], e[1]);
#pragma omp critical
            m_collisionConstraints.push_back(std::move(c));
            E2[e[2]][e[0]] = 0;
          }
        }

      }
      return;
    }


    std::vector<SurfaceMesh::TriangleArray> elements = mesh1->getTrianglesVertices();
//#pragma omp parallel for num_threads(32)
    for (int k = 0; k < elements.size(); ++k)
    {
      std::vector<std::vector<bool>> E2(nV2, std::vector<bool>(nV2, 1));

      SurfaceMesh::TriangleArray& tri = elements[k];
      size_t i1 = tri[0];
      size_t i2 = tri[1];

      if (E[i1][i2] && E[i2][i1])
      {
        const Vec3d& P = mesh1->getVertexPosition(i1);
        const Vec3d& Q = mesh1->getVertexPosition(i2);
        std::vector<SurfaceMesh::TriangleArray> elements2 = mesh2->getTrianglesVertices();
        for (size_t j = 0; j < elements2.size(); ++j)
        {
          SurfaceMesh::TriangleArray& e = elements2[j];
          const Vec3d& p0 = mesh2->getVertexPosition(e[0]);
          const Vec3d& p1 = mesh2->getVertexPosition(e[1]);
          const Vec3d& p2 = mesh2->getVertexPosition(e[2]);
          if (E2[e[0]][e[1]] && E2[e[1]][e[0]])
          {
            if (testLineToLineAABB(P, Q, p0, p1, proximity1, proximity2))
            {
              auto c = std::make_unique<EdgeEdgeCollisionConstraint>();
              c->initConstraint(m_object1, i1, i2, m_object2, e[0], e[1]);
#pragma omp critical
              m_collisionConstraints.push_back(std::move(c));
              E2[e[0]][e[1]] = 0;
            }
          }
          if (E2[e[1]][e[2]] && E2[e[2]][e[1]])
          {
            if (testLineToLineAABB(P, Q, p1, p2, proximity1, proximity2))
            {
              auto c = std::make_unique<EdgeEdgeCollisionConstraint>();
              c->initConstraint(m_object1, i1, i2, m_object2, e[1], e[2]);
#pragma omp critical
              m_collisionConstraints.push_back(std::move(c));
              E2[e[1]][e[2]] = 0;
            }
          }
          if (E2[e[2]][e[0]] && E2[e[0]][e[2]])
          {
            if (testLineToLineAABB(P, Q, p2, p0, proximity1, proximity2))
            {
              auto c = std::make_unique<EdgeEdgeCollisionConstraint>();
              c->initConstraint(m_object1, i1, i2, m_object2, e[2], e[1]);
#pragma omp critical
              m_collisionConstraints.push_back(std::move(c));
              E2[e[2]][e[0]] = 0;
            }
          }
        }
        E[i1][i2] = 0;
      }

      i1 = tri[1];
      i2 = tri[2];
      if (E[i1][i2] && E[i2][i1])
      {
        const Vec3d& P = mesh1->getVertexPosition(i1);
        const Vec3d& Q = mesh1->getVertexPosition(i2);
        std::vector<SurfaceMesh::TriangleArray> elements2 = mesh2->getTrianglesVertices();

        for (size_t j = 0; j < elements2.size(); ++j)
        {
          SurfaceMesh::TriangleArray& e = elements2[j];
          const Vec3d& p0 = mesh2->getVertexPosition(e[0]);
          const Vec3d& p1 = mesh2->getVertexPosition(e[1]);
          const Vec3d& p2 = mesh2->getVertexPosition(e[2]);
          if (E2[e[0]][e[1]] && E2[e[1]][e[0]])
          {
            if (testLineToLineAABB(P, Q, p0, p1, proximity1, proximity2))
            {
              auto c = std::make_unique<EdgeEdgeCollisionConstraint>();
              c->initConstraint(m_object1, i1, i2, m_object2, e[0], e[1]);
#pragma omp critical
              m_collisionConstraints.push_back(std::move(c));
              E2[e[0]][e[1]] = 0;
            }
          }
          if (E2[e[1]][e[2]] && E2[e[2]][e[1]])
          {
            if (testLineToLineAABB(P, Q, p1, p2, proximity1, proximity2))
            {
              auto c = std::make_unique<EdgeEdgeCollisionConstraint>();
              c->initConstraint(m_object1, i1, i2, m_object2, e[1], e[2]);
#pragma omp critical
              m_collisionConstraints.push_back(std::move(c));
              E2[e[1]][e[2]] = 0;
            }
          }
          if (E2[e[2]][e[0]] && E2[e[0]][e[2]])
          {
            if (testLineToLineAABB(P, Q, p2, p0, proximity1, proximity2))
            {
              auto c = std::make_unique<EdgeEdgeCollisionConstraint>();
              c->initConstraint(m_object1, i1, i2, m_object2, e[2], e[0]);
#pragma omp critical
              m_collisionConstraints.push_back(std::move(c));
              E2[e[2]][e[0]] = 0;
            }
          }
        }
        E[i1][i2] = 0;
      }

      i1 = tri[2];
      i2 = tri[0];
      if (E[i1][i2] && E[i2][i1])
      {
        const Vec3d& P = mesh1->getVertexPosition(i1);
        const Vec3d& Q = mesh1->getVertexPosition(i2);
        std::vector<SurfaceMesh::TriangleArray> elements2 = mesh2->getTrianglesVertices();
        for (size_t j = 0; j < elements2.size(); ++j)
        {
          SurfaceMesh::TriangleArray& e = elements2[j];
          const Vec3d& p0 = mesh2->getVertexPosition(e[0]);
          const Vec3d& p1 = mesh2->getVertexPosition(e[1]);
          const Vec3d& p2 = mesh2->getVertexPosition(e[2]);
          if (E2[e[0]][e[1]] && E2[e[1]][e[0]])
          {
            if (testLineToLineAABB(P, Q, p0, p1, proximity1, proximity2))
            {
              auto c = std::make_unique<EdgeEdgeCollisionConstraint>();
              c->initConstraint(m_object1, i1, i2, m_object2, e[0], e[1]);
#pragma omp critical
              m_collisionConstraints.push_back(std::move(c));
              E2[e[0]][e[1]] = 0;
            }
          }
          if (E2[e[1]][e[2]] && E2[e[2]][e[1]])
          {
            if (testLineToLineAABB(P, Q, p1, p2, proximity1, proximity2))
            {
              auto c = std::make_unique<EdgeEdgeCollisionConstraint>();
              c->initConstraint(m_object1, i1, i2, m_object2, e[1], e[2]);
#pragma omp critical
              m_collisionConstraints.push_back(std::move(c));
              E2[e[1]][e[2]] = 0;
            }
          }
          if (E2[e[2]][e[0]] && E2[e[0]][e[2]])
          {
            if (testLineToLineAABB(P, Q, p2, p0, proximity1, proximity2))
            {
              auto c = std::make_unique<EdgeEdgeCollisionConstraint>();
              c->initConstraint(m_object1, i1, i2, m_object2, e[2], e[0]);
#pragma omp critical
              m_collisionConstraints.push_back(std::move(c));
              E2[e[2]][e[0]] = 0;
            }
          }
        }
        E[i1][i2] = 0;
      }
    }
  }

}
