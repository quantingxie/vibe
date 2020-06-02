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

#include "imstkCPDState.h"

namespace imstk
{

	void CPDState::initialize(const size_t numNodes, const bool(&options)[3])
	{
		if (options[0])
			m_pos.resize(numNodes, Vec3d(0, 0, 0));

		if (options[1])
			m_vel.resize(numNodes, Vec3d(0, 0, 0));

		if (options[2])
			m_acc.resize(numNodes, Vec3d(0, 0, 0));
	}

	void CPDState::initialize(const std::shared_ptr<PointSet>& m, const bool(&options)[3])
	{
		this->initialize(m->getNumVertices(), options);
	}

	void CPDState::initialize(const cpd::ParticleObjectPtr p_particleObject, const bool(&options)[3]) {
		initialize(p_particleObject->getParticleCount(), options);
		updateStateFromParticleObject(p_particleObject);
	}

	void CPDState::setState(std::shared_ptr<CPDState> p_state)
	{
		m_pos = p_state->getPositions();
		m_vel = p_state->getVelocities();
		m_acc = p_state->getAccelerations();
	}

	void CPDState::updateStateFromParticleObject(const cpd::ParticleObjectPtr p_particleObject) {

		auto& pos = p_particleObject->getPositions();
		for (unsigned i = 0; i < m_pos.size(); i++) {
			setVertexPosition(i, pos[i]);
		}
	}

} // imstk
