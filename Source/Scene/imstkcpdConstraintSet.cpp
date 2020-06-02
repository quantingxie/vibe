#include "imstkcpdConstraintSet.h"

namespace cpd
{
  void ConstraintSet::setTimeStepSize(double p_dt)
  {
    for (auto& c : m_constraints)
    {
      c->setTimeStepSize(p_dt);
    }
  }

  void ConstraintSet::addConstraint(ConstraintBasePtr p_constraint)
  {
    auto p = std::find(m_constraints.begin(), m_constraints.end(), p_constraint);
    if (p == m_constraints.end())
    {
      m_constraints.push_back(std::move(p_constraint));
      m_orders.push_back(m_constraints.size() - 1);
    }
  }

  void ConstraintSet::clearConstraints()
  {
    m_constraints.clear();
    m_orders.clear();
  }

  //void ConstraintSet::clearAffectedParticles()
  //{
  //  m_affectedParticles.clear();
  //}

  void ConstraintSet::resetLambda()
  {
    for (auto& c : m_constraints)
    {
      c->clearLamda();
    }
  }

  void ConstraintSet::shuffle(std::mt19937& p_g)
  {
    std::shuffle(m_orders.begin(), m_orders.end(), p_g);
  }

  void ConstraintSet::reverse()
  {
    std::reverse(m_orders.begin(), m_orders.end());
  }

  bool ConstraintSet::initializeDistanceConstraints(ConstraintBase::SubType p_subType)
  {
    double stiffness = m_object->getProperties(0);
    bool isCompliant = (p_subType == ConstraintBase::SubType::XPBD);

    if (m_object->getSurfaceMesh() == nullptr)
    { // if no mesh, connect particles one by one

      unsigned np = m_object->getParticleCount() - 1;
      for (size_t k = 0; k < np; ++k)
      {
        initDistanceConstraint(k, k + 1, p_subType, stiffness, isCompliant);
      }
    }
    else
    {
      size_t pairs[3][2] = { {0,1},{1,2},{2,0} };
      auto surMesh = std::static_pointer_cast<SurfaceMesh>(m_object->getSurfaceMesh());
      int nV = surMesh->getNumVertices();
      std::vector<std::vector<bool>> E(nV, std::vector<bool>(nV, 1));
      auto elements = surMesh->getTrianglesVertices();

      for (size_t i = 0; i < elements.size(); i++)
      {
        auto& tri = elements[i];

        for (size_t j = 0; j < 3; j++)
        {
          size_t idx1 = tri[pairs[j][0]];
          size_t idx2 = tri[pairs[j][1]];

          if (E[idx1][idx2] && E[idx2][idx1])
          {
            initDistanceConstraint(idx1, idx2, p_subType, stiffness, isCompliant);
            E[idx1][idx2] = 0;
          }
        }
      }
    }
    //checkTypeAndID();
    return true;
  }

  void ConstraintSet::initDistanceConstraint(size_t p_idx1, size_t p_idx2, ConstraintBase::SubType p_subType, double p_stiffness, bool p_isCompliant)
  {
    if ((m_object->getInvMass(p_idx1) < EPS) && (m_object->getInvMass(p_idx2) < EPS))
      return;

    auto c = std::make_unique<DistanceConstraint>();

    if (p_subType != ConstraintBase::SubType::XCPD)
    {

      auto b = std::make_unique<DistancepbdConstraint>(p_isCompliant);
      c = static_cast<DistancepbdConstraintPtr>(std::move(b));
    }

    c->initConstraint(m_object, p_idx1, p_idx2, p_stiffness);

    m_constraints.push_back(std::move(c));
    m_orders.push_back(m_constraints.size() - 1);
  }

  bool ConstraintSet::initializeDihedralConstraints()
  {
    return true;
  }

  bool ConstraintSet::initializeAreaConstraints(ConstraintBase::SubType p_subType)
  {
    auto triMesh = std::static_pointer_cast<SurfaceMesh>(m_object->getSurfaceMesh());
    std::vector<SurfaceMesh::TriangleArray> elements = triMesh->getTrianglesVertices();

    double stiffness = m_object->getProperties(0);
    bool isCompliant = (p_subType == ConstraintBase::SubType::XPBD);

    for (size_t k = 0; k < elements.size(); ++k)
    {
      auto& tri = elements[k];

      bool allfixed = true;
      for (unsigned i = 0; i < 3; i++)
      {
        if (m_object->getInvMass(tri[i]) > EPS)
        {
          allfixed = false;
          break;
        }
      }
      if (allfixed)
        continue;

      auto c = std::make_unique<AreaConstraint>();
      c->initConstraint(m_object, tri[0], tri[1], tri[2], stiffness);
      m_constraints.push_back(std::move(c));
      m_orders.push_back(m_constraints.size() - 1);
    }
    return true;
  }

  bool ConstraintSet::initializeCSTEnergyConstraints(ConstraintBase::SubType p_subType)
  {
    auto triMesh = std::static_pointer_cast<SurfaceMesh>(m_object->getSurfaceMesh());
    std::vector<SurfaceMesh::TriangleArray> elements = triMesh->getTrianglesVertices();

    double youngsModulus = m_object->getProperties(0);
    double poissonRatio = m_object->getProperties(1);

    bool isCompliant = (p_subType == ConstraintBase::SubType::XPBD);

    //#pragma omp parallel for num_threads(16)
    for (size_t k = 0; k < elements.size(); ++k)
    {
      auto& tri = elements[k];

      bool allfixed = true;
      for (unsigned i = 0; i < 3; i++)
      {
        if (m_object->getInvMass(tri[i]) > EPS)
        {
          allfixed = false;
          break;
        }
      }
      if (allfixed)
        continue;

      auto c = std::make_unique<CSTEnergyConstraint>();

      if (p_subType != ConstraintBase::SubType::XCPD)
      {
        auto b = std::make_unique<CSTpbdConstraint>(isCompliant);
        c = static_cast<CSTpbdConstraintPtr>(std::move(b));
      }

      c->initConstraint(m_object, tri[0], tri[1], tri[2], youngsModulus, poissonRatio);
      m_constraints.push_back(std::move(c));
      m_orders.push_back(m_constraints.size() - 1);
    }
    //checkTypeAndID();
    return true;
  }

  bool ConstraintSet::initializeOneDEnergyConstraints()
  {

    for (size_t k = 0; k < m_object->getParticleCount() - 1; ++k)
    {
      auto c = std::make_unique<OneDEnergyConstraint>();
      c->initConstraint(m_object, k, k + 1, 1e6, 0.3);
      m_constraints.push_back(std::move(c));
      m_orders.push_back(m_constraints.size() - 1);
    }
    return true;
  }

  bool ConstraintSet::initializeCSVEnergyConstraints(ConstraintBase::SubType p_subType)
  {
    auto tetraMesh = std::static_pointer_cast<TetrahedronMesh>(m_object->getVolumeMesh());
    std::vector<TetrahedronMesh::TetraArray> elements = tetraMesh->getElementVertices();

    double youngsModulus = m_object->getProperties(0);
    double poissonRatio = m_object->getProperties(1);

    bool isCompliant = (p_subType == ConstraintBase::SubType::XPBD);

    for (size_t k = 0; k < elements.size(); ++k)
    {
      auto& tetra = elements[k];

      bool allfixed = true;
      for (unsigned i = 0; i < 4; i++)
      {
        if (m_object->getInvMass(tetra[i]) > EPS*EPS)
        {
          allfixed = false;
          break;
        }
      }
      if (allfixed)
        continue;

      auto c = std::make_unique<CSVEnergyConstraint>();

      if (p_subType != ConstraintBase::SubType::XCPD)
      {
        auto b = std::make_unique<CSVpbdConstraint>(isCompliant);
        c = static_cast<CSVpbdConstraintPtr>(std::move(b));
      }

      c->initConstraint(m_object, tetra[0], tetra[1], tetra[2], tetra[3], youngsModulus, poissonRatio);
      m_constraints.push_back(std::move(c));
      m_orders.push_back(m_constraints.size() - 1);
    }
    //checkTypeAndID();
    return true;
  }

  void ConstraintSet::initializeSystemVariables()
  {
    auto& tempPos = m_object->getTemporaryPositions();

    size_t nbrConstraint = m_constraints.size();
    size_t csz = m_constraints[0]->getConstraintSize();

    m_f.resize(nbrConstraint*csz);

    for (size_t i = 0; i < nbrConstraint; i++)
    {
      for (size_t j = 0; j < csz; j++)
      {
      }
    }
  }

  void ConstraintSet::checkTypeAndID() const
  {
    for (auto& c : m_constraints)
    {
      std::string type;
      auto tp = c->getType().getSubType();
      switch (tp)
      {
      case ConstraintBase::SubType::XCPD:
        type = "xcpd";
        break;
      case ConstraintBase::SubType::XPBD:
        type = "xpbd";
        break;
      case ConstraintBase::SubType::PBD:
        type = "pbd";
        break;
      default:
        break;
      }
      std::cout << "subtype = " << type << ", ID = " << c->getID() << std::endl;
    }
  }

}