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

#ifndef imstkCPDState_h
#define imstkCPDState_h

#include <Eigen/Dense>
#include <vector>

#include "imstkPointSet.h"
#include "imstkMath.h"

#include "imstkcpdParticleObject.h"

namespace imstk
{

	///
	/// \class CPDState
	///
	/// \brief State of the body governed by CPD mathematical model
	///
	class CPDState
	{
	public:
		///
		/// \brief Default constructor/destructor
		///
		CPDState() = default;
		~CPDState() = default;

		///
		/// \brief Initialize the CPD state
		///
		void initialize(const size_t numNodes, const bool(&options)[3]);
		void initialize(const std::shared_ptr<PointSet>& m, const bool(&options)[3]);

		///
		/// \brief Initialize the CPD state from CPD ParticleObject
		///
		void initialize(const cpd::ParticleObjectPtr p_particleObject, const bool(&options)[3]);

		///
		/// \brief Get/Set nodal position given the index
		///
		void setVertexPosition(const size_t idx, const Vec3d pos) { m_pos.at(idx) = pos; };
		Vec3d& getVertexPosition(const size_t idx) { return m_pos.at(idx); };

		///
		/// \brief Returns the vector of current nodal positions
		///
		StdVectorOfVec3d& getPositions() { return m_pos; };
		void setPositions(const StdVectorOfVec3d& p) { m_pos = p; };

		///
		/// \brief Returns the vector of current nodal velocities
		///
		StdVectorOfVec3d& getVelocities() { return m_vel; };

		///
		/// \brief Returns the vector of current nodal accelerations
		///
		StdVectorOfVec3d& getAccelerations() { return m_acc; };

		void setState(std::shared_ptr<CPDState> p_state);

		///
		/// \brief Update the CPD state from the CPD particle Object
		///
		void updateStateFromParticleObject(const cpd::ParticleObjectPtr p_particleObject);

	private:
		StdVectorOfVec3d m_pos; ///> Nodal positions
		StdVectorOfVec3d m_vel; ///> Nodal velocities
		StdVectorOfVec3d m_acc; ///> Nodal acelerations
	};

} // imstk

#endif // imstkCPDState_h

