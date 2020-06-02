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
	auto scene = simManager->createNewScene("CPD Liver Collision");

	cpd::ScenePtr CPDScene = std::make_shared<cpd::Scene>();
	scene->setCPDScene(CPDScene);

	/* this part should be from a file reader*/
	size_t nbrForces = 1;
	std::vector<std::array<double, 3>> forces;
	std::vector<bool> isDistributed;
	forces.resize(nbrForces);
	isDistributed.resize(nbrForces);
	forces[0][0] = 0.0;
	forces[0][1] = -9.8;
	forces[0][2] = 0.0;
	isDistributed[0] = true;
	cpd::BaseTimeIntegrator::Type type = cpd::BaseTimeIntegrator::Type::SIE;
	double timestep = 0.01;
	/* this part should be from a file reader*/

	// Create forces and time integrator
	createForces(CPDScene, forces, isDistributed);
	createTimeIntegrator(CPDScene, timestep, type);

	// Create liver cpdParticleObject from meshIO
	auto& liverCPDobject = cpd::createObjectFromMeshIO(CPDScene, "D:/GCD_ExampleFiles/liver");

	// Create and initialize liver imstk::CPDObject from cpdParticleObject
	auto liver_CSV = std::make_shared<CPDObject>("liver");
	liver_CSV->setParticleObject(liverCPDobject);
	auto liver_tetraMesh = std::make_shared<imstk::TetrahedralMesh>();
	auto liver_surfMesh = std::make_shared<imstk::SurfaceMesh>();
	liver_CSV->convertToMesh(liver_tetraMesh, liver_surfMesh);

	// Set up mapping between tet and surface mesh, by loading precomputed map
	auto TetTriMap = std::make_shared<imstk::TetraTriangleMap>();
	TetTriMap->setMaster(liver_tetraMesh);
	TetTriMap->setSlave(liver_surfMesh);
	//TetTriMap->compute();
	TetTriMap->load("D:/GCD_ExampleFiles/liver.wet");

	// Set up material and visual model
	auto material = std::make_shared<RenderMaterial>();
	material->setDisplayMode(RenderMaterial::DisplayMode::SURFACE);
	material->setColor(liverRed);
	auto liverSurfMeshModel = std::make_shared<VisualModel>(liver_surfMesh);
	liverSurfMeshModel->setRenderMaterial(material);

	// Create CPD dynamical Model
	auto liverCpdModel = std::make_shared<CPDModel>();

	// Set Dynamical Model, Visual Model and Map
	liver_CSV->setDynamicalModel(liverCpdModel);
	liver_CSV->addVisualModel(liverSurfMeshModel);
	liver_CSV->setPhysicsGeometry(liver_tetraMesh);
	liver_CSV->setPhysicsToVisualMap(TetTriMap);


	// Create plane cpdParticleObject
	auto& planeCPDobject = cpd::createPlane1(50, CPDScene, "");

	// Create and initialize plane imstk::CPDObject from cpdParticleObject
	auto planeObj = std::make_shared<CPDObject>("Plane");
	planeObj->setParticleObject(planeCPDobject);
	auto planeMesh = std::make_shared<imstk::SurfaceMesh>();
	planeObj->convertToMesh(planeMesh);

	// Set up material and visual model
	auto materialplane = std::make_shared<RenderMaterial>();
	materialplane->setDisplayMode(RenderMaterial::DisplayMode::SURFACE);
	materialplane->setColor(tiffanyBlue);
	materialplane->backfaceCullingOff();
	auto planesurfMeshModel = std::make_shared<VisualModel>(planeMesh);
	planesurfMeshModel->setRenderMaterial(materialplane);

	// Create CPD dynamical Model
	auto cpdModelPlane = std::make_shared<CPDModel>();

	// Set Dynamical Model, Visual Model and Map
	planeObj->setDynamicalModel(cpdModelPlane);
	planeObj->addVisualModel(planesurfMeshModel);
	planeObj->setPhysicsGeometry(planeMesh);

	// Add collision pair
	CPDScene->addCollisionPair(liverCPDobject, planeCPDobject);

	// Add sceneObjects
	scene->addSceneObject(liver_CSV);
	scene->addSceneObject(planeObj);

	// Light
	auto whiteLight = std::make_shared<imstk::PointLight>("whiteLight");
	whiteLight->setPosition(0, 8, -15);
	whiteLight->setFocalPoint(imstk::Vec3d(-1, -4, -1));
	whiteLight->setIntensity(300);
	scene->addLight(whiteLight);

	// Camera
	scene->getCamera()->setFocalPoint(-2, -2, -10);
	scene->getCamera()->setPosition(2, 8, -32);

	simManager->setActiveScene(scene);
	simManager->start(SimulationStatus::paused);

    return 0;
}
