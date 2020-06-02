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

#ifndef imstkCPDObject_h
#define imstkCPDObject_h

#include "imstkDynamicObject.h"
#include "imstkDynamicalModel.h"
#include "imstkCPDModel.h"
#include "imstkSurfaceMesh.h"
#include "imstkTetrahedralMesh.h"

#include <stdarg.h>

namespace imstk
{

	class Geometry;
	class GeometryMap;

	///
	/// \class CPDObject
	///
	/// \brief Base class for scene objects that move and/or deform under constrained
	/// particle dynamics formulation
	///
	class CPDObject : public DynamicObject
	{
	public:
		///
		/// \brief Constructor
		///
		CPDObject(std::string name) : DynamicObject(name) { m_type = SceneObject::Type::CPD; }

		CPDObject(cpd::ParticleObjectPtr p_particleObject, std::string name) : DynamicObject(name), m_particleObject(p_particleObject)
		{
			m_type = SceneObject::Type::CPD;
		}

		///
		/// \brief Destructor
		///
		virtual ~CPDObject() = default;

		void setParticleObject(cpd::ParticleObjectPtr p_particleObject) { m_particleObject = p_particleObject; }

		void convertToMesh(std::shared_ptr<PointSet> p_pointSet);
		void convertToMesh(std::shared_ptr<TetrahedralMesh> p_tetraMesh, std::shared_ptr<SurfaceMesh> p_surfMesh = nullptr);
		void convertToMesh(std::shared_ptr<SurfaceMesh> p_surfMesh);

		bool initialize() override;

	protected:

		std::shared_ptr<CPDModel> m_CPDModel; ///> CPD mathematical model
		cpd::ParticleObjectPtr m_particleObject; ///> CPD object
	};

} // imstk

#endif // imstkCPDObject_h

