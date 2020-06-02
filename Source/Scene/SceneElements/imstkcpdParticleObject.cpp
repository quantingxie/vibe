#include "imstkcpdParticleObject.h"
#include "imstkcpdParticleObject.h"
#include "imstkcpdParticleObject.h"
#include "imstkcpdParticleObject.h"

namespace cpd {

  ParticleObject::ParticleObject()
  {
    m_particleProperty = std::make_shared<ParticleProperty>();
  }

  ParticleObject::ParticleObject(StdVectorOfVec3d& p_positions)
  {
    m_positions = p_positions;
    m_particleInitTpye = ParticlesInitTpye::PARTICLES;
    m_particleProperty = std::make_shared<ParticleProperty>();

	initialize();
  }

  ParticleObject::ParticleObject(VolumeMeshPtr p_volumeMesh)
  {
    m_volumeMesh = p_volumeMesh;
    m_particleInitTpye = ParticlesInitTpye::VOLUMETRICMESH;
    m_particleProperty = std::make_shared<ParticleProperty>();

	initialize();
  }

  ParticleObject::ParticleObject(SurfaceMeshPtr p_surfaceMesh)
  {
    m_surfaceMesh = p_surfaceMesh;
    m_particleInitTpye = ParticlesInitTpye::SURFACEMESH;
    m_particleProperty = std::make_shared<ParticleProperty>();

	initialize();
  }

  void ParticleObject::initialize()
  {
    if (m_particleInitTpye == cpd::ParticleObject::ParticlesInitTpye::SURFACEMESH)
    {
      auto& vertices = m_surfaceMesh->getVertexPositions();
      for (auto& v : vertices) {
        m_positions.push_back(v);
      }
    }
    else if (m_particleInitTpye == cpd::ParticleObject::ParticlesInitTpye::VOLUMETRICMESH)
    {
      auto& vertices = m_volumeMesh->getVertexPositions();
      for (auto& v : vertices) {
        m_positions.push_back(v);
      }
    }

	m_numNodes = m_positions.size();

	m_temporaryPositions.resize(m_numNodes);
	m_velocities.resize(m_numNodes, Vec3d(0, 0, 0));
	m_accelerations.resize(m_numNodes, Vec3d(0, 0, 0));
	m_constraintForces.resize(m_numNodes, Vec3d(0, 0, 0));
	m_displacements.resize(m_numNodes, Vec3d(0, 0, 0));
	m_constraintForceVectors.resize(m_numNodes);
	m_invMass.resize(m_numNodes);
	m_collided.resize(m_numNodes, false);
	m_forceUpdated.resize(m_numNodes, false);
	m_excludeCollision.resize(m_numNodes, false);

	for (unsigned i = 0; i < m_numNodes; i++)
	{
		m_temporaryPositions[i] = m_positions[i];

		if (m_positions[i][1] > -2.0) //temporary, for dragon example
			m_velocities[i][2] = (m_positions[i][1] + 2.0) / 5.0;

		m_constraintForceVectors[i].clear();
	}

	m_initialized = true;

  }

  void ParticleObject::initializeDerivativeVariables()
  {
	  if (!m_initPositionGiven)
	  {
		  m_initPositions.resize(m_numNodes);
		  for (unsigned i = 0; i < m_numNodes; i++)
		  {
			  m_initPositions[i] = m_positions[i];

		  }
	  }

	  if (m_massList.size() == 0)
	  {
		  std::cout << "MassList is not provided. Use default particle mass instead." << std::endl;
		  for (unsigned i = 0; i < m_numNodes; i++)
		  {
			  m_invMass[i] = m_particleProperty->getProperty(PropertyType::MASS);
		  }
	  }
	  else if (m_massList.size() == m_numNodes)
	  {
		  for (unsigned i = 0; i < m_numNodes; i++)
		  {
			  m_invMass[i] = ((m_massList[i] < EPS) || (abs(m_massList[i] + 1.0) < EPS)) ? 0.0 : 1 / m_massList[i];
		  }
	  }
	  else
	  {
		  std::cout << "MassList size WRONG!!" << std::endl;
	  }

	  for (auto& i : m_fixedParticleIndex)
	  {
		  m_invMass[i] = 0;
	  }

	  for (auto& i : m_excludeIDs)
	  {
		  m_excludeCollision[i] = true;
	  }
  }

  void ParticleObject::setInitPositions(StdVectorOfVec3d p_particles)
  {
    m_initPositionGiven = true;
    m_initPositions = p_particles;
  }


  void ParticleObject::setParticleProperty(ParticlePropertyPtr p_property)
  {
    m_particleProperty = p_property;
  }

  void ParticleObject::setProperties(double p_property1, double p_property2)
  {
    m_properties[0] = p_property1;
    m_properties[1] = p_property2;
  }

  void ParticleObject::updateMesh()
  {
    if (m_volumeMesh != nullptr)
    {
      for (unsigned i = 0; i < m_volumeMesh->getNumVertices(); i++)
      {
        m_volumeMesh->setVertexPosition(i, m_positions[i]);
      }
    }
    else if (m_surfaceMesh != nullptr)
    {
      for (unsigned i = 0; i < m_surfaceMesh->getNumVertices(); i++)
      {
        m_surfaceMesh->setVertexPosition(i, m_positions[i]);
      }
    }
  }

  void ParticleObject::reset()
  {
    for (unsigned i = 0; i < m_numNodes; i++) {
      m_positions[i] = m_initPositions[i];
      m_velocities[i].setZero();
      m_temporaryPositions[i] = m_initPositions[i];
    }
    updateMesh();

    m_toBeReset = false;
  }

  void ParticleObject::writePositions(Mode p_mode)
  {
    Eigen::VectorXd v1[3], v2[3];
    std::string name1, name2;

    for (unsigned i = 0; i < 3; i++)
    {
      v1[i].resize(m_numNodes);

      if (p_mode == Mode::Both)
        v2[i].resize(m_numNodes);
    }

    std::string prefix[] = { "x","y","z" };
    for (unsigned i = 0; i < 3; i++)
    {
      switch (p_mode)
      {
      case ParticleObject::Mode::POS:
        name1 = prefix[i] + "Position.m";
        break;
      case ParticleObject::Mode::DISP:
        name1 = prefix[i] + "Displacement.m";
        break;
      case ParticleObject::Mode::Both:
        name1 = prefix[i] + "Position.m";
        name2 = prefix[i] + "Displacement.m";
        break;
      default:
        break;
      }

      auto& var1 = v1[i];
      auto& var2 = v2[i];
      for (unsigned j = 0; j < m_numNodes; j++)
      {
        switch (p_mode)
        {
        case ParticleObject::Mode::POS:
          var1[j] = m_positions[j][i];
          break;
        case ParticleObject::Mode::DISP:
          var1[j] = getDisplacement(j)[i];
          break;
        case ParticleObject::Mode::Both:
          var1[j] = m_positions[j][i];
          var2[j] = getDisplacement(j)[i];
          break;
        default:
          break;
        }
      }

      writeVectorMatlabPlot(var1, name1.c_str());

      if (p_mode == ParticleObject::Mode::Both)
      {
        writeVectorMatlabPlot(var2, name2.c_str());
      }

    }
  }

  void ParticleObject::preSimulation(double p_dt)
  {
    m_dt = p_dt;

	if(!m_initialized)
		initialize();
	initializeDerivativeVariables();

    std::cout << "Number of particles = " << m_numNodes << std::endl;
  }

  void ParticleObject::postSimulation()
  {

  }

  void ParticleObject::simulate()
  {

  }

  void ParticleObject::clearConstraintForces()
  {
    for (unsigned i = 0; i < m_numNodes; i++) {
      m_constraintForces[i].setZero();
      m_constraintForceVectors[i].clear();
      m_forceUpdated[i] = false;
    }
  }

  void ParticleObject::addToConstraintForce(unsigned p_idx, Vec3d p_force)
  {
    m_constraintForces[p_idx] += p_force;
    m_forceUpdated[p_idx] = true;
  }

  void ParticleObject::addToConstraintForce(unsigned p_idx, unsigned p_cid, Vec3d p_force)
  {
    Vec3d pref;
    auto& pf = m_constraintForceVectors[p_idx].find(p_cid);
    if (pf == m_constraintForceVectors[p_idx].end()) {
      pref.setZero();
    }
    else
    {
      pref = m_constraintForceVectors[p_idx][p_cid];
    }

    m_constraintForces[p_idx] += p_force - pref;
    m_constraintForceVectors[p_idx][p_cid] = p_force;
    m_forceUpdated[p_idx] = true;
  }

  void ParticleObject::updateTempDisplacement()
  {
    for (int i = 0; i < m_numNodes; i++)
    {
      m_displacements[i] = m_temporaryPositions[i] - m_initPositions[i];
    }
  }

  double ParticleObject::testConvergence()
  {
    double totalNorm = 0.0;
    for (auto& f : m_constraintForces)
    {
      totalNorm += f.norm();
    }
    return totalNorm;
  }

  void ParticleObject::updateParticleAcceleration()
  {
    for (auto& i : m_fixedParticleIndex)
    {
      m_massList[i] = -1.0;
    }
    for (unsigned i = 0; i < m_numNodes; i++) {
      if (m_massList[i] > EPS)
        m_accelerations[i] += m_externalForce / m_massList[i] + m_distributedForce;
    }
  }

  void ParticleObject::addExternalForce(const Vec3d& p_force, unsigned p_idx, double p_mass)
  {
	  if (!m_initialized)
	  {
		  std::cout << "Object not initialized!" << std::endl;
		  return;
	  }

	  m_accelerations[p_idx] += p_force / p_mass;
  }

  void ParticleObject::addExternalForce(const Vec3d& p_force, std::vector<unsigned> p_idx, std::vector<double> p_mass)
  {
	  if (!m_initialized)
	  {
		  std::cout << "Object not initialized!" << std::endl;
		  return;
	  }

	  if (p_idx.size() != p_mass.size())
	  {
		  std::cout << "Number of mass and node do not match!" << std::endl;
		  return;
	  }
	  for (unsigned i = 0; i < p_idx.size(); i++)
	  {
		  m_accelerations[p_idx[i]] += p_force / p_mass[i];
	  }
  }

  void ParticleObject::updateTempPositions()
  {
    //    int thread = 16;
    //#pragma omp parallel for num_threads(thread)
    for (int i = 0; i < m_numNodes; i++)
    {
      if (m_massList[i] < EPS)
        continue;

      m_temporaryPositions[i] += m_constraintForces[i];
    }
  }

  void ParticleObject::updatePositions(bool p_updateVelocity)
  {
    for (unsigned i = 0; i < m_numNodes; i++)
    {
      if (p_updateVelocity)
      {
          double damp = 0.995;// 1.0;
          if (m_velocities[i].norm() > 0.1)
              damp = 0.997; //1.0;//1.0;


        m_velocities[i] = damp * (m_temporaryPositions[i] - m_positions[i]) / m_dt;

        if (m_collided[i])
        {
          m_velocity *= -0.0;
          m_velocity[2] = -0.0*m_velocity[1];
          m_collided[i] = false;
        }
      }

      auto prePos = m_positions[i];
      m_positions[i] = m_temporaryPositions[i];
      m_temporaryPositions[i] = prePos; // temporaryPositions will hold prevPositions for now
    }
  }


  void ParticleObject::updateVelocity()
  {

  }

  Vec3d ParticleObject::getDisplacement(unsigned p_id)
  {
    return m_positions[p_id] - m_initPositions[p_id];
  }

  Vec3d ParticleObject::getTempDisplacement(unsigned p_id)
  {
    return m_displacements[p_id];
  }

  VecXd ParticleObject::getTempDisplacement()
  {
    VecXd disp;
    disp.resize(m_displacements.size() * 3);

    for (unsigned i = 0; i < m_displacements.size(); i++)
    {
      for (unsigned j = 0; j < 3; j++)
      {
        disp[i * 3 + j] = m_displacements[i][j];
      }
    }

    return disp;
  }

  void ParticleObject::addConstraintType(ConstraintBase::ConstraintType p_type)
  {
    bool isPBD = (p_type.getSubType() != ConstraintBase::SubType::XCPD);

    if (m_constriantTypes.size() > 0)
    {
      if (isPBD != m_isPBD)
      {
        std::cout << "ConstraintTypes inconsistent: CPD and PBD constraints can't be mixed!";
        return;
      }
    }


    m_constriantTypes.push_back(p_type);

    m_isPBD = isPBD;
  }

  void ParticleObject::updateBoundingBox()
  {
	  Vec3d min(max_d, max_d, max_d);
	  Vec3d max(-max_d, -max_d, -max_d);

	  int id = 0;
	  for (auto& tempPosition : m_temporaryPositions)
	  {
		  if (!m_excludeCollision[id])
		  {
			  for (int i = 0; i < 3; i++)
			  {
				  min[i] = std::min(min[i], tempPosition[i]);
				  max[i] = std::max(max[i], tempPosition[i]);
			  }
		  }
		  id++;
	  }
	  m_boundingBox.setValue(min, max);

	  double prox = m_particleProperty->getProperty(PropertyType::PROXIMITY);
	  m_boundingBox.extend(prox);
  }

  void ParticleObject::setSurfaceMesh(StdVectorOfVec3d p_vertexPositions, std::vector<std::array<size_t, 3>> p_trianglesVertices)
  {
    m_surfaceMesh = std::make_shared<SurfaceMesh>(p_vertexPositions, p_trianglesVertices);
  }

  SurfaceMeshPtr ParticleObject::getSurfaceMesh()
  {
    if (m_surfaceMesh)
      return m_surfaceMesh;

    if (m_volumeMesh)
    {
      if (m_volumeMesh->getAttachedSurfMesh())
      {
        m_surfaceMesh = m_volumeMesh->getAttachedSurfMesh();
        return m_surfaceMesh;
      }
    }
    return nullptr;
  }

}