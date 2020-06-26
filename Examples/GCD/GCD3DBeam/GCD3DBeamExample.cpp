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

void testCPD3DBeam();
void testCPD3DSlab();    // created by Jose James for the slab deformation using GCD

///
/// \brief This example demonstrates the 2D Beam simulation
/// using Generalized Constrained Dynamics
///
int
main()
{
	//testCPD3DBeam();
	testCPD3DSlab();    // created by Jose James
	
    return 0;
}

void testCPD3DBeam()
{
	auto simManager = std::make_shared<SimulationManager>();
	auto scene = simManager->createNewScene("GCD 3D Beam");

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
	double timestep = 0.001;
	/* this part should be from a file reader*/

	auto& beam = cpd::create3DBeam(CPDScene, resourceDir);

	cpd::createForces(CPDScene, forces, isDistributed);
	cpd::createTimeIntegrator(CPDScene, timestep, type);

	auto deformableObj = std::make_shared<CPDObject>("3D Beam");
	deformableObj->setParticleObject(beam);
	auto tetraMesh = std::make_shared<imstk::TetrahedralMesh>();
	deformableObj->convertToMesh(tetraMesh);

	// set up visual model
	auto material = std::make_shared<RenderMaterial>();
	material->setDisplayMode(RenderMaterial::DisplayMode::WIREFRAME_SURFACE);
	material->setColor(Color::Green);
	auto TetraMeshModel = std::make_shared<VisualModel>(tetraMesh);
	TetraMeshModel->setRenderMaterial(material);

	// Object & Model
	auto cpdModel = std::make_shared<CPDModel>();
	deformableObj->setDynamicalModel(cpdModel);
	//deformableObj->setVisualGeometry(tetraMesh);
	deformableObj->addVisualModel(TetraMeshModel);
	deformableObj->setPhysicsGeometry(tetraMesh);

	//auto material = std::make_shared<RenderMaterial>();
	//material->setBackFaceCulling(false);
	////material->setDiffuseColor(Color::White);
	//material->setDisplayMode(RenderMaterial::DisplayMode::WIREFRAME_SURFACE);
	//material->setSphereGlyphSize(.5);
	////tetraMesh->setRenderMaterial(material);

	scene->addSceneObject(deformableObj);

	// Light (white)
	auto whiteLight1 = std::make_shared<imstk::DirectionalLight>("whiteLight1");
	whiteLight1->setFocalPoint(imstk::Vec3d(-1, -4, -1));
	whiteLight1->setIntensity(1);
	auto whiteLight2 = std::make_shared<imstk::DirectionalLight>("whiteLight2");
	whiteLight2->setFocalPoint(imstk::Vec3d(1, 4, 1));
	whiteLight2->setIntensity(1);

	// Add in scene
	scene->addLight(whiteLight1);
	scene->addLight(whiteLight2);

	scene->getCamera()->setFocalPoint(1, 10, 50);
	scene->getCamera()->setPosition(-90, 40, 50);

	simManager->setActiveScene(scene);
	simManager->start(SimulationStatus::paused);

}

void testCPD3DSlab()
{

	auto sdk = std::make_shared<SimulationManager>();
	auto scene = sdk->createNewScene("GCD_3D_Slab");

	cpd::ScenePtr CPDScene = std::make_shared<cpd::Scene>();
	scene->setCPDScene(CPDScene);


	//std::string resourceDir = "F:\Research\Resources\APIUtilities";
	std::string resourceDir;
	/* this part should be from a file reader*/
	//size_t nbrForces = 1;
	size_t nbrForces = 1;
	std::vector<std::array<double, 3>> forces;
	std::vector<bool> isDistributed;
	forces.resize(nbrForces);
	isDistributed.resize(nbrForces);
	//forces[0][0] = 0.0;
	//forces[0][1] = -9.8;
	//forces[0][2] = 0.0;

	forces[0][0] = 0;
	forces[0][1] = -100;
	forces[0][2] = 0;
	isDistributed[0] = true;
	//isDistributed[0] = false;

	cpd::BaseTimeIntegrator::Type type = cpd::BaseTimeIntegrator::Type::CD;

	//double timestep = 0.001;
	double timestep = 0.001;

	/* this part should be from a file reader*/

	//auto& beam = cpd::create3DBeam(CPDScene, resourceDir);
	auto& beam = cpd::create3DSlab(CPDScene, resourceDir);

	cpd::createForces(CPDScene, forces, isDistributed);
	cpd::createTimeIntegrator(CPDScene, timestep, type);

	auto deformableObj = std::make_shared<CPDObject>("3D Slab");
	deformableObj->setParticleObject(beam);
	auto tetraMesh = std::make_shared<imstk::TetrahedralMesh>();
	deformableObj->convertToMesh(tetraMesh);
	int numTet = tetraMesh->getNumTetrahedra();
	std::cout << "Number of tetrahedral elements: " << numTet << "\n";

	// set up visual model
	auto material = std::make_shared<RenderMaterial>();
	material->setDisplayMode(RenderMaterial::DisplayMode::WIREFRAME_SURFACE);
	material->setColor(Color::Green);
	auto TetraMeshModel = std::make_shared<VisualModel>(tetraMesh);
	TetraMeshModel->setRenderMaterial(material);

	// Object & Model
	auto cpdModel = std::make_shared<CPDModel>();
	deformableObj->setDynamicalModel(cpdModel);
	//deformableObj->setVisualGeometry(tetraMesh);
	deformableObj->addVisualModel(TetraMeshModel);
	deformableObj->setPhysicsGeometry(tetraMesh);

	//auto material = std::make_shared<RenderMaterial>();
	//material->setBackFaceCulling(false);
	////material->setDiffuseColor(Color::White);
	//material->setDisplayMode(RenderMaterial::DisplayMode::WIREFRAME_SURFACE);
	//material->setSphereGlyphSize(.5);
	////tetraMesh->setRenderMaterial(material);

	scene->addSceneObject(deformableObj);

	// Light (white)
	auto whiteLight1 = std::make_shared<imstk::DirectionalLight>("whiteLight1");
	whiteLight1->setFocalPoint(imstk::Vec3d(-1, -4, -1));
	whiteLight1->setIntensity(1);
	auto whiteLight2 = std::make_shared<imstk::DirectionalLight>("whiteLight2");
	whiteLight2->setFocalPoint(imstk::Vec3d(1, 4, 1));
	whiteLight2->setIntensity(1);

	// Add in scene
	scene->addLight(whiteLight1);
	scene->addLight(whiteLight2);

	scene->getCamera()->setFocalPoint(10, -2, 5);
	scene->getCamera()->setPosition(-7.5, 5, 5.0);

	sdk->setActiveScene(scene);
	sdk->start();
}