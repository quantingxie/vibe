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

#include "imstkSimulationManager.h"
#include "imstkCPDModel.h"
#include "imstkCPDObject.h"
#include "imstkAPIUtilities.h"
#include "imstkcpdAPIUtilities.h"
#include "imstkScene.h"
#include <omp.h>
#include <Map\imstkTetraTriangleMap.h>
#include <Map\imstkOneToOneMap.h>
#include "imstkcpdMeshIO.h"

using namespace imstk;

Color tiffanyBlue = Color(0.039, 0.729, 0.71);
Color liverRed = Color(0.92, 0.21, 0.19);

///
/// \brief This example demonstrates the 2D Beam simulation
/// using Generalized Constrained Dynamics
///
int
main()
{
	auto simManager = std::make_shared<SimulationManager>();
	cpd::ScenePtr CPDScene = std::make_shared<cpd::Scene>();
	auto scene = simManager->createNewScene("GCD_Cloth");
	scene->setCPDScene(CPDScene);

	// a. Create surface mesh
	auto surfMesh = std::make_shared<imstk::SurfaceMesh>();

	const double lengthY = 10.0;
	const double lengthX = 10.0;
	const int nRowX = 31;
	const int nRowY = 31;
	const double rho = 800.0;
	const double g = -9.8;
	const double thickness = 1.0;
	const double area = lengthY * lengthX;
	const double element_area = 0.5 * area / ((nRowX - 1) * (nRowY - 1));
	const double nodeMass = element_area * rho * thickness / 3.0;

	// a.1 Vertices coordinates
	cpd::StdVectorOfVec3d vertList;
	vertList.resize(nRowX * nRowY);
	const double dy = lengthY / (double)(nRowY - 1);
	const double dx = lengthX / (double)(nRowX - 1);
	for (int i = 0; i < nRowY; ++i)
	{
		for (int j = 0; j < nRowX; j++)
		{
			vertList[i * nRowX + j] = Vec3d((double)dx * j, (double)dy * i, 1.0);
		}
	}

	// a.2 Connectivity
	std::vector<imstk::SurfaceMesh::TriangleArray> triangles;
	std::vector<double> massList;
	massList.resize(nRowX * nRowY);
	for (unsigned i = 0; i < massList.size(); i++)
	{
		massList[i] = 0.0;
	}

	for (std::size_t i = 0; i < nRowY - 1; ++i)
	{
		imstk::SurfaceMesh::TriangleArray tri[2];
		for (std::size_t j = 0; j < nRowX - 1; j++)
		{
			tri[0] = { { i * nRowX + j, (i + 1) * nRowX + j, i * nRowX + j + 1 } };
			tri[1] = { { (i + 1) * nRowX + j + 1, i * nRowX + j + 1, (i + 1) * nRowX + j } };
			triangles.push_back(tri[0]);
			triangles.push_back(tri[1]);

			for (unsigned m = 0; m < 3; m++)
			{
				massList[tri[0][m]] += nodeMass;
				massList[tri[1][m]] += nodeMass;
			}
		}
	}

	// a.3 Initialize surfMesh
	surfMesh->setInitialVertexPositions(vertList);
	surfMesh->setVertexPositions(vertList);
	surfMesh->setTrianglesVertices(triangles);

	// b. Create SceneObject
	// b.1 Create cpdParticleObject
	cpd::ParticleObjectPtr particleObject = std::make_shared<cpd::ParticleObject>();

	std::vector<unsigned> fixed;
	fixed.push_back(0);
	fixed.push_back((nRowY - 1) * nRowX);
	fixed.push_back(nRowX - 1);
	fixed.push_back(nRowY * nRowX - 1);

	particleObject->setFixedPoints(fixed);
	particleObject->setProperties(1E5);
	//particleObject->addConstraintType(cpd::ConstraintBase::ConstraintType::Area);
	particleObject->addConstraintType(cpd::ConstraintBase::ConstraintType(cpd::ConstraintBase::Type::Distance, cpd::ConstraintBase::SubType::XCPD));

	cpd::ParticlePropertyPtr property = std::make_shared<cpd::ParticleProperty>();
	property->setProperty(cpd::PropertyType::MASS, 20.0 / nRowX);
	particleObject->setParticleProperty(property);
	particleObject->setMassList(massList);

	CPDScene->addObject(particleObject);

	// b.2 Create CPDObject and its visual/dynamical models
	auto deformableObj = std::make_shared<CPDObject>("Cloth");
	auto surfMeshModel = std::make_shared<VisualModel>(surfMesh);
	auto cpdModel = std::make_shared<CPDModel>();

	auto material = std::make_shared<RenderMaterial>();
	material->setDisplayMode(RenderMaterial::DisplayMode::WIREFRAME_SURFACE);
	material->setColor(Color::Green);
	material->setBackFaceCulling(false);
	surfMeshModel->setRenderMaterial(material);

	deformableObj->setDynamicalModel(cpdModel);
	deformableObj->addVisualModel(surfMeshModel);
	deformableObj->setPhysicsGeometry(surfMesh);
	deformableObj->setParticleObject(particleObject);

	scene->addSceneObject(deformableObj);

	// c. Create force and time integrator
	/* this part should be from a file reader*/
	size_t nbrForces = 1;
	std::vector<std::array<double, 3>> forces;
	std::vector<bool> isDistributed;
	forces.resize(nbrForces);
	isDistributed.resize(nbrForces);
	forces[0][0] = 0.0;
	forces[0][1] = 0.0;
	forces[0][2] = 9.8;
	isDistributed[0] = true;
	cpd::BaseTimeIntegrator::Type type = cpd::BaseTimeIntegrator::Type::CD;
	double timestep = 0.005;
	/* this part should be from a file reader*/

	createForces(CPDScene, forces, isDistributed);
	createTimeIntegrator(CPDScene, timestep, type);

	// Light and camera
	auto whiteLight = std::make_shared<imstk::DirectionalLight>("whiteLight1");
	whiteLight->setFocalPoint(imstk::Vec3d(-5, 4, 5));
	whiteLight->setIntensity(10);
	scene->addLight(whiteLight);

	scene->getCamera()->setFocalPoint(-5, 5, 6.5);
	scene->getCamera()->setPosition(20, 5, 6.5);
	scene->getCamera()->setViewUp(0, 0, -1);

	simManager->setActiveScene(scene);
	simManager->start(SimulationStatus::paused);

    return 0;
}
