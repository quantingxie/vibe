/*=========================================================================

   Library: iMSTK

   Copyright (c) Kitware, Inc. & Center for Modeling, Simulation,
   & Imaging in Medicine, Rensselaer Polytechnic Institute.

   Licensed under the Apache License, Version 2.0 (the "License");
   you may not use this file except in compliance with the License.
   You may obtain a copy of the License at

      http://www.apache.org/licenses/LICENSE-2.0.txt

   Unless required by applicable law or agreed to in writing, software
   distributed under the License is distributed on an "AS IS" BASIS,
   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
   See the License for the specific language governing permissions and
   limitations under the License.

=========================================================================*/

#include "imstkCPDObject.h"
#include "imstkGeometryMap.h"

#include <g3log/g3log.hpp>

namespace imstk
{
	void CPDObject::convertToMesh(std::shared_ptr<PointSet> p_pointSet)
	{
		StdVectorOfVec3d vertices;

		for (auto& p : m_particleObject->getPositions())
		{
			vertices.push_back(p);
		}
		p_pointSet->setInitialVertexPositions(vertices);
		p_pointSet->setVertexPositions(vertices);
	}

	void CPDObject::convertToMesh(std::shared_ptr<TetrahedralMesh> p_tetraMesh, std::shared_ptr<SurfaceMesh> p_surfMesh)
	{
		auto& volmesh = m_particleObject->getVolumeMesh();

		p_tetraMesh->setInitialVertexPositions(volmesh->getVertexPositions());
		p_tetraMesh->setVertexPositions(volmesh->getVertexPositions());
		p_tetraMesh->setTetrahedraVertices(volmesh->getElementVertices());

		if (!p_surfMesh)
			p_surfMesh = std::make_shared<SurfaceMesh>();

		if (!m_particleObject->getSurfaceMesh())
		{
			p_tetraMesh->extractSurfaceMesh(p_surfMesh, true);
			m_particleObject->setSurfaceMesh(p_surfMesh->getInitialVertexPositions(), p_surfMesh->getTrianglesVertices());
		}
		else
		{
			p_surfMesh->setInitialVertexPositions(m_particleObject->getSurfaceMesh()->getVertexPositions());
			p_surfMesh->setVertexPositions(m_particleObject->getSurfaceMesh()->getVertexPositions());
			p_surfMesh->setTrianglesVertices(m_particleObject->getSurfaceMesh()->getTrianglesVertices());
		}
	}

	void CPDObject::convertToMesh(std::shared_ptr<SurfaceMesh> p_surfMesh)
	{
		auto& surfmesh = m_particleObject->getSurfaceMesh();

		p_surfMesh->setInitialVertexPositions(surfmesh->getVertexPositions());
		p_surfMesh->setVertexPositions(surfmesh->getVertexPositions());
		p_surfMesh->setTrianglesVertices(surfmesh->getTrianglesVertices());
	}


	bool CPDObject::initialize()
	{
		auto type = m_physicsGeometry->getType();

		if (type == Geometry::Type::SurfaceMesh)
		{
			auto mesh = std::static_pointer_cast<SurfaceMesh>(m_physicsGeometry);

			cpd::SurfaceMeshPtr cpdSurfaceMesh = std::make_shared<cpd::SurfaceMesh>(mesh->getVertexPositions(), mesh->getTrianglesVertices());
			m_particleObject->setSurfaceMesh(cpdSurfaceMesh);
			m_particleObject->setParticleInitType(cpd::ParticleObject::ParticlesInitTpye::SURFACEMESH);
		}
		else if (type == Geometry::Type::TetrahedralMesh)
		{
			auto mesh = std::static_pointer_cast<TetrahedralMesh>(m_physicsGeometry);

			cpd::TetrahedronMeshPtr cpdTetraMesh = std::make_shared<cpd::TetrahedronMesh>(mesh->getVertexPositions(), mesh->getTetrahedraVertices());
			m_particleObject->setVolumeMesh(cpdTetraMesh);
			m_particleObject->setParticleInitType(cpd::ParticleObject::ParticlesInitTpye::VOLUMETRICMESH);
		}
		else if (type == Geometry::Type::PointSet)
		{
			//auto mesh = std::static_pointer_cast<PointSet>(m_physicsGeometry);
			//std::vector<cpd::ParticlePtr> particles;
			//for (int i = 0; i < mesh->getNumVertices(); i++)
			//{
			//  particles.push_back(std::make_shared<cpd::Particle>(mesh->getVertexPosition(i)));
			//}
			//m_particleObject->setParticles(particles); // TODO:
			//m_particleObject->setParticleInitType(cpd::ParticleObject::ParticlesInitTpye::PARTICLES);
		}


		m_CPDModel = std::dynamic_pointer_cast<CPDModel>(m_dynamicalModel);
		m_CPDModel->setParticleObject(m_particleObject);
		m_CPDModel->setPhysicsGeometry(m_physicsGeometry);

		if (m_CPDModel)
		{
			return DynamicObject::initialize();
		}
		else
		{
			LOG(WARNING) << "Dynamics pointer cast failure in CPDObject::initialize()";
			return false;
		}
	}

} //imstk
