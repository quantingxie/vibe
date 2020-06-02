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
#include "imstkCPDModel.h"
#include "imstkTetrahedralMesh.h"
#include "imstkSurfaceMesh.h"

#include <g3log/g3log.hpp>

namespace imstk {

	CPDModel::CPDModel() : DynamicalModel(DynamicalModelType::CPD) {
		m_initialState = std::make_shared<CPDState>();
		m_previousState = std::make_shared<CPDState>();
		m_currentState = std::make_shared<CPDState>();
	}

	bool CPDModel::initialize()
	{
		if ((!m_particleObject) || (!m_physicsGeometry))
			return false;

		bool option[3] = { 1, 1, 1 };
		m_initialState->initialize(m_particleObject, option);
		m_previousState->initialize(m_particleObject, option);
		m_currentState->initialize(m_particleObject, option);

		return true;
	}

	void CPDModel::updatePhysicsGeometry() {

		if (m_physicsGeometry->isMesh()) {
			auto& mesh = std::dynamic_pointer_cast<PointSet>(m_physicsGeometry);

			StdVectorOfVec3d pos = m_particleObject->getPositions();
			m_particleObject->updateMesh(); // duplicated due to method used in narrowCCD
			//for (unsigned i = 0; i < mesh->getNumVertices(); i++) {
			//	mesh->setVertexPosition(i, pos[i]);
			//}
			mesh->setVertexPositions(pos);
		}
	}

	void CPDModel::resetToInitialState()
	{
		m_currentState->setState(m_initialState);
		m_previousState->setState(m_initialState);

		m_particleObject->setToReset();
	}

} // imstk
