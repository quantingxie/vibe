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
	auto scene = simManager->createNewScene("CPD Dragon");

	cpd::ScenePtr CPDScene = std::make_shared<cpd::Scene>();
	scene->setCPDScene(CPDScene);

	bool convertMesh = false;
	if (convertMesh)
	{
		//auto TetMesh = std::dynamic_pointer_cast<TetrahedralMesh>(MeshIO::read("D:/Research/IMSTK_GCDDebug/install/data/asianDragon/asianDragon.veg"/*"D:/GCD_ExampleFiles/dragon77k/dragon-coarse.tet.veg"*/));
		auto TetMesh = std::dynamic_pointer_cast<TetrahedralMesh>(MeshIO::read(iMSTK_DATA_ROOT"/asianDragon/asianDragon.veg"/*"D:/GCD_ExampleFiles/dragon77k/dragon-coarse.tet.veg"*/));

		//auto TetMesh = std::dynamic_pointer_cast<TetrahedralMesh>(MeshIO::read("F:/VIBE/Resources/vibe/3D_New/Stomach/tetramodels/1/stomach_new_SB.tet"/*"D:/GCD_ExampleFiles/dragon77k/dragon-coarse.tet.veg"*/));

		auto& vertexPositions = TetMesh->getVertexPositions();
		auto& tetrahedraVertices = TetMesh->getTetrahedraVertices();
		auto tetra = std::make_shared<cpd::TetrahedronMesh>(vertexPositions, tetrahedraVertices);

		auto& tetIO = std::make_shared<cpd::TETMeshIO>();
		tetIO->write(tetra, "");

		auto& surfIO = std::make_shared<cpd::SURFMeshIO>();
		//auto& surf = surfIO->read("D:/GCD_ExampleFiles/dragon77k/dragon77K.srf");
		auto surfmesh = std::make_shared<imstk::SurfaceMesh>();
		//surfmesh->initialize(surf->getVertexPositions(), surf->getTrianglesVertices());
		TetMesh->extractSurfaceMesh(surfmesh, true);
		auto& vertexPositions_s = surfmesh->getVertexPositions();
		auto& triangleVertices_s = surfmesh->getTrianglesVertices();
		auto surf = std::make_shared<cpd::SurfaceMesh>(vertexPositions_s, triangleVertices_s);
		surfIO->write(surf, "");

		auto TetTriMap1 = std::make_shared<imstk::TetraTriangleMap>();
		TetTriMap1->setMaster(TetMesh);
		TetTriMap1->setSlave(surfmesh);
		TetTriMap1->compute();

		std::cout << "Done" << std::endl;
	}

	/* this part should be from a file reader*/
	size_t nbrForces = 1;
	std::vector<std::array<double, 3>> forces;
	std::vector<bool> isDistributed;
	forces.resize(nbrForces);
	isDistributed.resize(nbrForces);
	forces[0][0] = 0.0;
	forces[0][1] = -0.8;
	forces[0][2] = 0.0;
	isDistributed[0] = true;
	cpd::BaseTimeIntegrator::Type type = cpd::BaseTimeIntegrator::Type::SIE;
	double timestep = 0.02;
	/* this part should be from a file reader*/

	// Create forces and time integrator
	createForces(CPDScene, forces, isDistributed);
	createTimeIntegrator(CPDScene, timestep, type);

	// Create liver cpdParticleObject from meshIO
	//std::string path2files("D:/GCD_ExampleFiles/dragon");
	std::string path2files(iMSTK_DATA_ROOT"/GCD_ExampleFiles/dragon");
	auto& dragonCPDobject = cpd::createObjectFromMeshIO(CPDScene, path2files);
	std::vector<unsigned> fixed;
	auto& vmesh = dragonCPDobject->getVolumeMesh();
	auto& vertices = vmesh->getVertexPositions();
	for (int i = 0; i < vertices.size(); ++i)
	{
		if (vertices[i][1] < -2.0)
			fixed.push_back(i);
	}
	dragonCPDobject->setFixedPoints(fixed);
	const double youngs = 5E11; // read
	const double mu = 0.4; // read
	dragonCPDobject->setProperties(youngs, mu);
	std::vector<double> massList;
	size_t n = vmesh->getNumVertices();
	massList.resize(n);
	double nodeMass = 0.00000001;
	for (int i = 0; i < n; i++)
	{
		massList[i] = nodeMass;
	}
	dragonCPDobject->setMassList(massList);

	// Create and initialize liver imstk::CPDObject from cpdParticleObject
	auto dragon_CSV = std::make_shared<CPDObject>("dragon");
	dragon_CSV->setParticleObject(dragonCPDobject);
	auto dragon_tetraMesh = std::make_shared<imstk::TetrahedralMesh>();
	auto dragon_surfMesh = std::make_shared<imstk::SurfaceMesh>();
	dragon_CSV->convertToMesh(dragon_tetraMesh, dragon_surfMesh);

	// Set up mapping between tet and surface mesh, by loading precomputed map
	auto TetTriMap = std::make_shared<imstk::TetraTriangleMap>();
	TetTriMap->setMaster(dragon_tetraMesh);
	TetTriMap->setSlave(dragon_surfMesh);
	//TetTriMap->compute();
	TetTriMap->load(path2files + ".wet");

	// Set up material and visual model
	auto material = std::make_shared<RenderMaterial>();
	material->setDisplayMode(RenderMaterial::DisplayMode::SURFACE);
	material->setColor(Color::Yellow);
	material->setMetalness(0.8);
	auto dragonSurfMeshModel = std::make_shared<VisualModel>(dragon_surfMesh);
	dragonSurfMeshModel->setRenderMaterial(material);

	// Create CPD dynamical Model
	auto dragonCpdModel = std::make_shared<CPDModel>();

	// Set Dynamical Model, Visual Model and Map
	dragon_CSV->setDynamicalModel(dragonCpdModel);
	dragon_CSV->addVisualModel(dragonSurfMeshModel);
	dragon_CSV->setPhysicsGeometry(dragon_tetraMesh);
	dragon_CSV->setPhysicsToVisualMap(TetTriMap);

	/*
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
	CPDScene->addCollisionPair(dragonCPDobject, planeCPDobject);
	*/

	// Add sceneObjects
	scene->addSceneObject(dragon_CSV);
	//scene->addSceneObject(planeObj);

	// Light
	auto whiteLight = std::make_shared<imstk::PointLight>("whiteLight");
	whiteLight->setPosition(-100, 80, 50);
	whiteLight->setFocalPoint(imstk::Vec3d(-1, -4, -1));
	whiteLight->setIntensity(900);
	//scene->addLight(whiteLight);

	// Camera
	scene->getCamera()->setFocalPoint(0, -1, -10);
	scene->getCamera()->setPosition(0, 3, 12);

	simManager->setActiveScene(scene);
	simManager->start(SimulationStatus::paused);

    return 0;
}
