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

#ifndef imstkCPDModel_h
#define imstkCPDModel_h

#include <vector>
#include <Eigen/Dense>

#include "imstkcpdParticleObject.h"
#include "imstkDynamicalModel.h"
#include "imstkCPDState.h"
#include "imstkMath.h"

namespace imstk
{

	///
	/// \class CPDModel
	///
	/// \brief This class implements constrained particle dynamics mathematical model
	///
	class CPDModel : public DynamicalModel<CPDState>
	{
	public:
		///
		/// \brief Constructor
		///
		CPDModel();

		///
		/// \brief Destructor
		///
		~CPDModel() = default;

		///
		/// \brief Initialize the states
		///
		bool initialize() override;

		///
		/// \brief
		///
		void updateBodyStates(const Vectord& q, const StateUpdateType updateType = StateUpdateType::displacement) override {};

		///
		/// \brief
		///
		void setPhysicsGeometry(std::shared_ptr<Geometry> p_physicsGeometry) { m_physicsGeometry = p_physicsGeometry; }

		void setParticleObject(cpd::ParticleObjectPtr p_particleObject) { m_particleObject = p_particleObject; }

		///
		/// \brief
		///
		void updatePhysicsGeometry();

		void resetToInitialState() override;

		void setTimeStep(const double timeStep) override {}
		double getTimeStep() const override { return 0.01; }


	protected:
		std::shared_ptr<Geometry> m_physicsGeometry;
		cpd::ParticleObjectPtr m_particleObject;
	};

} // imstk

#endif // imstkCPDModel_h
