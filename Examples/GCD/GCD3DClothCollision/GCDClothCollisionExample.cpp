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
void clothCollisionGCD();
void clothCollisionGCDTest();
///
/// \brief This example demonstrates the 2D Beam simulation
/// using Generalized Constrained Dynamics
///
int
main()
{
	//clothCollisionGCD();
	clothCollisionGCDTest();
	
	return 0;
}

void clothCollisionGCD()
{
	auto simManager = std::make_shared<SimulationManager>();
	cpd::ScenePtr CPDScene = std::make_shared<cpd::Scene>();
	auto scene = simManager->createNewScene("GCD_Cloth");
	scene->setCPDScene(CPDScene);

	// a. Create surface mesh
	auto surfMesh = std::make_shared<imstk::SurfaceMesh>();

	const double lengthY = 10.0;
	const double lengthX = 10.0;
	const int nRowX = 51;
	const int nRowY = 51;
	const double rho = 300.0;
	const double g = 9.8;
	const double thickness = 1.0;
	const double area = lengthY * lengthX;
	const double element_area = 0.5 * area / ((nRowX - 1) * (nRowY - 1));
	const double nodeMass = element_area * rho * thickness / 3.0;
	const double youngs = 8E4;

	double shift = -lengthX / 4.0;

	// a.1 Vertices coordinates
	cpd::StdVectorOfVec3d vertList;
	vertList.resize(nRowX * nRowY);
	const double dy = lengthY / (double)(nRowY - 1);
	const double dx = lengthX / (double)(nRowX - 1);
	for (int i = 0; i < nRowY; ++i)
	{
		for (int j = 0; j < nRowX; j++)
		{
			vertList[i * nRowX + j] = Vec3d((double)dx * j + shift, (double)dy * i - lengthY / 2.0, 3.0);
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
	surfMesh->flipNormals();

	// b. Create SceneObject
	// b.1 Create cpdParticleObject
	cpd::ParticleObjectPtr particleObject = std::make_shared<cpd::ParticleObject>();

	std::vector<unsigned> fixed;
	unsigned fixedsize = 2;
	for (unsigned i = 0; i < fixedsize; i++)
	{
		for (unsigned j = 0; j < fixedsize - i; j++)
		{
			fixed.push_back(i * nRowX + j);
			fixed.push_back(i * nRowX + nRowX - 1 - j);

			fixed.push_back((nRowY - 1 - i) * nRowX + j);
			fixed.push_back((nRowY - 1 - i) * nRowX + nRowX - 1 - j);
		}
	}
	/*fixed.push_back(0);
	fixed.push_back((nRowY - 1) * nRowX);
	fixed.push_back(nRowX - 1);
	fixed.push_back(nRowY * nRowX - 1);*/

	unsigned regionsize = 15;
	std::vector<unsigned> excluded;
	for (unsigned i = 0; i < regionsize; i++)
	{
		for (unsigned j = 0; j < regionsize - i; j++)
		{
			excluded.push_back(i * nRowX + j);
			excluded.push_back(i * nRowX + nRowX - 1 - j);

			excluded.push_back((nRowY - 1 - i) * nRowX + j);
			excluded.push_back((nRowY - 1 - i) * nRowX + nRowX - 1 - j);
		}
	}
	particleObject->setFixedPoints(fixed);
	particleObject->setExculdeCollision(excluded);
	particleObject->setProperties(youngs);
	//particleObject->addConstraintType(cpd::ConstraintBase::ConstraintType::Area);
	particleObject->addConstraintType(cpd::ConstraintBase::ConstraintType(cpd::ConstraintBase::Type::Distance, cpd::ConstraintBase::SubType::XCPD));

	cpd::ParticlePropertyPtr property = std::make_shared<cpd::ParticleProperty>();
	property->setProperty(cpd::PropertyType::PROXIMITY, 0.1);
	particleObject->setParticleProperty(property);
	particleObject->setMassList(massList);

	CPDScene->addObject(particleObject);

	// b.2 Create CPDObject and its visual/dynamical models
	auto deformableObj = std::make_shared<CPDObject>("Cloth");
	auto surfMeshModel = std::make_shared<VisualModel>(surfMesh);
	auto cpdModel = std::make_shared<CPDModel>();

	auto material = std::make_shared<RenderMaterial>();
	material->setDisplayMode(RenderMaterial::DisplayMode::SURFACE);
	material->setColor(tiffanyBlue);
	material->setBackFaceCulling(false);
	surfMeshModel->setRenderMaterial(material);

	deformableObj->setDynamicalModel(cpdModel);
	deformableObj->addVisualModel(surfMeshModel);
	deformableObj->setPhysicsGeometry(surfMesh);
	deformableObj->setParticleObject(particleObject);

	scene->addSceneObject(deformableObj);

	// c. Create second cpdParticleObject
	auto& sphereCPDobject = cpd::createObjectFromSurfMesh(CPDScene, "F:/VIBE/Resources/GCD_ExampleFiles/sphere");
	sphereCPDobject->setParticleProperty(property);

	// c.1 Create and initialize plane imstk::CPDObject from cpdParticleObject
	auto sphereObj = std::make_shared<CPDObject>("sphere");
	sphereObj->setParticleObject(sphereCPDobject);
	auto sphereMesh = std::make_shared<imstk::SurfaceMesh>();
	sphereObj->convertToMesh(sphereMesh);

	// c.2 Set up material and visual model
	auto materialsphere = std::make_shared<RenderMaterial>();
	materialsphere->setDisplayMode(RenderMaterial::DisplayMode::SURFACE);
	materialsphere->setColor(liverRed);
	materialsphere->backfaceCullingOff();
	auto spheresurfMeshModel = std::make_shared<VisualModel>(sphereMesh);
	spheresurfMeshModel->setRenderMaterial(materialsphere);

	// c.3 Create CPD dynamical Model
	auto cpdModelPlane = std::make_shared<CPDModel>();

	// c.4 Set Dynamical Model, Visual Model and Map
	sphereObj->setDynamicalModel(cpdModelPlane);
	sphereObj->addVisualModel(spheresurfMeshModel);
	sphereObj->setPhysicsGeometry(sphereMesh);

	scene->addSceneObject(sphereObj);

	// Add collision pair
	CPDScene->addCollisionPair(particleObject, sphereCPDobject);



	// d. Create force and time integrator
	/* this part should be from a file reader*/
	size_t nbrForces = 1;
	std::vector<std::array<double, 3>> forces;
	std::vector<bool> isDistributed;
	forces.resize(nbrForces);
	isDistributed.resize(nbrForces);
	forces[0][0] = 0.0;
	forces[0][1] = 0.0;
	forces[0][2] = -9.8;
	isDistributed[0] = true;
	cpd::BaseTimeIntegrator::Type type = cpd::BaseTimeIntegrator::Type::SIE;
	double timestep = 0.005;
	/* this part should be from a file reader*/

	createForces(CPDScene, forces, isDistributed);
	createTimeIntegrator(CPDScene, timestep, type);

	// Light and camera
	auto whiteLight = std::make_shared<imstk::PointLight>("whiteLight1");
	whiteLight->setPosition(0, 4, 0);
	//whiteLight->setFocalPoint(imstk::Vec3d(0, 4, 0));
	whiteLight->setIntensity(100);
	//scene->addLight(whiteLight);

	scene->getCamera()->setFocalPoint(-shift, -1, -0.4);
	scene->getCamera()->setPosition(-shift, 1.2 * lengthY, 8);
	scene->getCamera()->setViewUp(0, 0, 1);

	simManager->setActiveScene(scene);
	simManager->start(SimulationStatus::paused);
	
}

void clothCollisionGCDTest()
{
	
	auto simManager = std::make_shared<SimulationManager>();
	cpd::ScenePtr CPDScene = std::make_shared<cpd::Scene>();
	auto scene = simManager->createNewScene("GCD_Cloth");
	scene->setCPDScene(CPDScene);

	// a. Create surface mesh
	auto surfMesh = std::make_shared<imstk::SurfaceMesh>();

	const double lengthY = 10.0;
	const double lengthX = 10.0;
	const int nRowX = 51;
	const int nRowY = 51;
	const double rho = 300.0;
	const double g = 9.8;
	const double thickness = 1.0;
	const double area = lengthY * lengthX;
	const double element_area = 0.5 * area / ((nRowX - 1) * (nRowY - 1));
	const double nodeMass = element_area * rho * thickness / 3.0;
	const double youngs = 8E4;

	double shift = -lengthX / 4.0;

	// a.1 Vertices coordinates
	cpd::StdVectorOfVec3d vertList;
	vertList.resize(nRowX * nRowY);
	const double dy = lengthY / (double)(nRowY - 1);
	const double dx = lengthX / (double)(nRowX - 1);
	for (int i = 0; i < nRowY; ++i)
	{
		for (int j = 0; j < nRowX; j++)
		{
			vertList[i * nRowX + j] = Vec3d((double)dx * j + shift, (double)dy * i - lengthY / 2.0, 3.0);
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
	surfMesh->flipNormals();

	// b. Create SceneObject
	// b.1 Create cpdParticleObject
	cpd::ParticleObjectPtr particleObject = std::make_shared<cpd::ParticleObject>();

	std::vector<unsigned> fixed;
	unsigned fixedsize = 2;
	for (unsigned i = 0; i < fixedsize; i++)
	{
		for (unsigned j = 0; j < fixedsize - i; j++)
		{
			fixed.push_back(i * nRowX + j);
			fixed.push_back(i * nRowX + nRowX - 1 - j);

			fixed.push_back((nRowY - 1 - i) * nRowX + j);
			fixed.push_back((nRowY - 1 - i) * nRowX + nRowX - 1 - j);
		}
	}
	/*fixed.push_back(0);
	fixed.push_back((nRowY - 1) * nRowX);
	fixed.push_back(nRowX - 1);
	fixed.push_back(nRowY * nRowX - 1);*/

	unsigned regionsize = 15;
	std::vector<unsigned> excluded;
	for (unsigned i = 0; i < regionsize; i++)
	{
		for (unsigned j = 0; j < regionsize - i; j++)
		{
			excluded.push_back(i * nRowX + j);
			excluded.push_back(i * nRowX + nRowX - 1 - j);

			excluded.push_back((nRowY - 1 - i) * nRowX + j);
			excluded.push_back((nRowY - 1 - i) * nRowX + nRowX - 1 - j);
		}
	}
	particleObject->setFixedPoints(fixed);
	particleObject->setExculdeCollision(excluded);
	particleObject->setProperties(youngs);
	//particleObject->addConstraintType(cpd::ConstraintBase::ConstraintType::Area);
	particleObject->addConstraintType(cpd::ConstraintBase::ConstraintType(cpd::ConstraintBase::Type::Distance, cpd::ConstraintBase::SubType::XCPD));

	cpd::ParticlePropertyPtr property = std::make_shared<cpd::ParticleProperty>();
	property->setProperty(cpd::PropertyType::PROXIMITY, 0.1);
	particleObject->setParticleProperty(property);
	particleObject->setMassList(massList);

	CPDScene->addObject(particleObject);

	// b.2 Create CPDObject and its visual/dynamical models
	auto deformableObj = std::make_shared<CPDObject>("Cloth");
	auto surfMeshModel = std::make_shared<VisualModel>(surfMesh);
	auto cpdModel = std::make_shared<CPDModel>();

	auto material = std::make_shared<RenderMaterial>();
	material->setDisplayMode(RenderMaterial::DisplayMode::SURFACE);
	material->setColor(tiffanyBlue);
	material->setBackFaceCulling(false);
	surfMeshModel->setRenderMaterial(material);

	deformableObj->setDynamicalModel(cpdModel);
	deformableObj->addVisualModel(surfMeshModel);
	deformableObj->setPhysicsGeometry(surfMesh);
	deformableObj->setParticleObject(particleObject);

	scene->addSceneObject(deformableObj);

	// c. Create second cpdParticleObject
	auto& sphereCPDobject = cpd::createObjectFromSurfMesh(CPDScene, "F:/VIBE/Resources/GCD_ExampleFiles/sphere");
	sphereCPDobject->setParticleProperty(property);
	
	// c.1 Create and initialize plane imstk::CPDObject from cpdParticleObject
	auto sphereObj = std::make_shared<CPDObject>("sphere");
	sphereObj->setParticleObject(sphereCPDobject);
	auto sphereMesh = std::make_shared<imstk::SurfaceMesh>();
	sphereObj->convertToMesh(sphereMesh);
	sphereMesh->setTranslation(4,0,0);

	// c.2 Set up material and visual model
	auto materialsphere = std::make_shared<RenderMaterial>();
	materialsphere->setDisplayMode(RenderMaterial::DisplayMode::SURFACE);
	materialsphere->setColor(liverRed);
	materialsphere->backfaceCullingOff();
	auto spheresurfMeshModel = std::make_shared<VisualModel>(sphereMesh);
	spheresurfMeshModel->setRenderMaterial(materialsphere);
	

	// c.3 Create CPD dynamical Model
	auto cpdModelPlane = std::make_shared<CPDModel>();

	// c.4 Set Dynamical Model, Visual Model and Map
	sphereObj->setDynamicalModel(cpdModelPlane);
	sphereObj->addVisualModel(spheresurfMeshModel);
	sphereObj->setPhysicsGeometry(sphereMesh);

	
	scene->addSceneObject(sphereObj);

	// Add collision pair
	CPDScene->addCollisionPair(particleObject, sphereCPDobject);
	

	// d. Create force and time integrator
	/* this part should be from a file reader*/
	size_t nbrForces = 1;
	std::vector<std::array<double, 3>> forces;
	std::vector<bool> isDistributed;
	forces.resize(nbrForces);
	isDistributed.resize(nbrForces);
	forces[0][0] = 0.0;
	forces[0][1] = 0.0;
	forces[0][2] = -9.8;
	isDistributed[0] = true;
	cpd::BaseTimeIntegrator::Type type = cpd::BaseTimeIntegrator::Type::SIE;
	double timestep = 0.005;
	/* this part should be from a file reader*/

	createForces(CPDScene, forces, isDistributed);
	createTimeIntegrator(CPDScene, timestep, type);

	// Light and camera
	auto whiteLight = std::make_shared<imstk::PointLight>("whiteLight1");
	whiteLight->setPosition(0, 4, 0);
	//whiteLight->setFocalPoint(imstk::Vec3d(0, 4, 0));
	whiteLight->setIntensity(100);
	//scene->addLight(whiteLight);

	scene->getCamera()->setFocalPoint(-shift, -1, -0.4);
	scene->getCamera()->setPosition(-shift, 1.2 * lengthY, 8);
	scene->getCamera()->setViewUp(0, 0, 1);

	simManager->setActiveScene(scene);
	simManager->start(SimulationStatus::paused);
	
}