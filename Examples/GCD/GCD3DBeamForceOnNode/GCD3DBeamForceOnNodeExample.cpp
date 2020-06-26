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
int testCPD3DBeamFonN();    
int testGCD3DSlabFonN();    // created by Jose James for the deformation by external force on selected nodes using GCD.


int
main()
{
	//testCPD3DBeamFonN();
	testGCD3DSlabFonN();

}

int testCPD3DBeamFonN()
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

	std::vector<unsigned> fixed;
	for (unsigned i = 0; i < 3 * 12; i++)
	{
		fixed.push_back(i);
		fixed.push_back(269 - i);
	}
	beam->setFixedPoints(fixed);
	beam->setProperties(8e7);

	//bool forceOnSingleNode = false;
	bool forceOnSingleNode = true;
	if (forceOnSingleNode)
	{
		std::cout << "forceOnSingleNode : " << forceOnSingleNode << endl;
		// (1) force, nodeID, nodeMass(=1.0 by default)
		//beam->addExternalForce(Vec3d(-300.0, 1000.0, 0.0), 134);
		//beam->addExternalForce(Vec3d(0, -1000.0, 800.0), 72);

		beam->addExternalForce(Vec3d(0.0, -1000.0, 0.0), 201);
		//beam->addExternalForce(Vec3d(0.0, -1000.0, 0.0), 202);
	}
	else
	{
		std::cout << "forceOnSingleNode : " << forceOnSingleNode << endl;
		// (2) force, nodeID list, nodeMass list
		//std::vector<unsigned> idx{ 132,133,134,168,169,170 };
		//std::vector<double> mass(6, 1.0);
		//beam->addExternalForce(Vec3d(-100.0, 500.0, 0.0), idx, mass);

		std::vector<unsigned> idx{ 132,133,134,168,169,170 };
		std::vector<double> mass(6, 1.0);
		beam->addExternalForce(Vec3d(-100.0, 500.0, 0.0), idx, mass);


	}

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

	return 0;

}

int testGCD3DSlabFonN()
{
	auto simManager = std::make_shared<SimulationManager>();
	auto scene = simManager->createNewScene("GCD 3D Slab");

	cpd::ScenePtr CPDScene = std::make_shared<cpd::Scene>();
	scene->setCPDScene(CPDScene);


	std::string resourceDir;
	//std::string resourceDir = "F:\Research\Resources\APIUtilities";
	/* this part should be from a file reader*/
	size_t nbrForces = 1;
	std::vector<std::array<double, 3>> forces;
	std::vector<bool> isDistributed;
	forces.resize(nbrForces);
	isDistributed.resize(nbrForces);
	forces[0][0] = 0.0;
	forces[0][1] = 0.0;
	forces[0][2] = 0.0;
	isDistributed[0] = true;
	cpd::BaseTimeIntegrator::Type type = cpd::BaseTimeIntegrator::Type::CD;
	//double timestep = 0.0001;
	double timestep = 0.001;

	/* this part should be from a file reader*/

	//auto& beam = cpd::create3DBeam(CPDScene, resourceDir);
	auto& beam = cpd::create3DSlab(CPDScene, resourceDir);

	cpd::createForces(CPDScene, forces, isDistributed);
	cpd::createTimeIntegrator(CPDScene, timestep, type);

	std::vector<unsigned> fixed;
	for (unsigned i = 0; i < 3 * 12; i++)
	{
		//fixed.push_back(i);
		//fixed.push_back(269 - i);
	}
	//beam->setFixedPoints(fixed);
	beam->setProperties(8e7);

	//bool forceOnSingleNode = false;
	bool forceOnSingleNode = true;
	if (forceOnSingleNode)
	{
		std::cout << "forceOnSingleNode : " << forceOnSingleNode << endl;
		// (1) force, nodeID, nodeMass(=1.0 by default)
		//beam->addExternalForce(Vec3d(-300.0, 1000.0, 0.0), 134);
		//beam->addExternalForce(Vec3d(0, -1000.0, 800.0), 72);

		beam->addExternalForce(Vec3d(0.0, -1000.0, 0.0), 230);
		//beam->addExternalForce(Vec3d(0.0, -2000.0, 0.0), 315);
		//beam->addExternalForce(Vec3d(0.0, -2000.0, 0.0), 335);
		beam->addExternalForce(Vec3d(0.0, -1000.0, 0.0), 630);
	}
	else
	{
		std::cout << "forceOnSingleNode : " << forceOnSingleNode << endl;
		// (2) force, nodeID list, nodeMass list
		//std::vector<unsigned> idx{ 132,133,134,168,169,170 };
		//std::vector<double> mass(6, 1.0);
		//beam->addExternalForce(Vec3d(-100.0, 500.0, 0.0), idx, mass);

		std::vector<unsigned> idx{ 132,133,134,168,169,170 };
		std::vector<double> mass(6, 1.0);
		beam->addExternalForce(Vec3d(-100.0, 500.0, 0.0), idx, mass);


	}

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
	//simManager->start(SimulationStatus::paused);
	simManager->start(SimulationStatus::running);

	return 0;

}