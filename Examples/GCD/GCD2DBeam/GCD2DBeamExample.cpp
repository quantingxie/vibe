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
	auto scene = simManager->createNewScene("GCD_2D_Beam");

	cpd::ScenePtr CPDScene = std::make_shared<cpd::Scene>();
	scene->setCPDScene(CPDScene);

	std::string resourceDir = "F:\Research\Resources\APIUtilities";
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
	cpd::BaseTimeIntegrator::Type type = cpd::BaseTimeIntegrator::Type::CD;
	double timestep = 0.01;
	/* this part should be from a file reader*/

	auto& beam = cpd::create2DBeam(CPDScene, resourceDir);

	createForces(CPDScene, forces, isDistributed);
	createTimeIntegrator(CPDScene, timestep, type);

	auto deformableObj = std::make_shared<CPDObject>("2D Beam");
	deformableObj->setParticleObject(beam);
	auto surfMesh = std::make_shared<imstk::SurfaceMesh>();
	deformableObj->convertToMesh(surfMesh);

	// set up visual model
	auto material = std::make_shared<RenderMaterial>();
	material->setDisplayMode(RenderMaterial::DisplayMode::WIREFRAME_SURFACE);
	material->setColor(tiffanyBlue);
	auto SurfMeshModel = std::make_shared<VisualModel>(surfMesh);
	SurfMeshModel->setRenderMaterial(material);

	auto cpdModel = std::make_shared<CPDModel>();
	deformableObj->setDynamicalModel(cpdModel);
	deformableObj->addVisualModel(SurfMeshModel);
	deformableObj->setPhysicsGeometry(surfMesh);

	scene->addSceneObject(deformableObj);

	// Light (white)
	auto whiteLight1 = std::make_shared<imstk::DirectionalLight>("whiteLight1");
	whiteLight1->setFocalPoint(imstk::Vec3d(-1, 4, -1));
	whiteLight1->setIntensity(10);
	auto whiteLight2 = std::make_shared<imstk::DirectionalLight>("whiteLight2");
	whiteLight2->setFocalPoint(imstk::Vec3d(1, 4, 1));
	whiteLight2->setIntensity(10);

	// Add in scene
	scene->addLight(whiteLight1);
	scene->addLight(whiteLight2);

	scene->getCamera()->setFocalPoint(10, 0.5, -15);
	scene->getCamera()->setPosition(10, 0.5, 20);

	simManager->setActiveScene(scene);
	simManager->start(SimulationStatus::paused);

    return 0;
}
