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

void testCPDVolume();
void testCPDVolumeVibeMesh(); // created by Jose James
void testLoadVIBECPD3DModels(); // created by Jose James
void loadStomach3DModel(); // created by Jose James
void loadESGTool(); // created by Jose James
void testESGToolParts(); // created by Jose James

///
/// \brief This example demonstrates the 2D Beam simulation
/// using Generalized Constrained Dynamics
///
int
main()
{
	//testCPDVolume();
	//testCPDVolumeVibeMesh(); // created by Jose James
	testLoadVIBECPD3DModels();  // created by Jose James
	
    return 0;
}

void testCPDVolume()
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
	auto& liverCPDobject = cpd::createObjectFromMeshIO(CPDScene, "F:/VIBE/Resources/GCD_ExampleFiles/liver");

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
	TetTriMap->load("F:/VIBE/Resources/GCD_ExampleFiles/liver.wet");

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
}

void testCPDVolumeVibeMesh()
{
	auto sdk = std::make_shared<SimulationManager>();
	auto scene = sdk->createNewScene("CPD_VolumeVibeMesh");

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

	//auto& meshobject = cpd::createObjectFromMeshIO(CPDScene, "D:/GCD_ExampleFiles/liver");
	auto& meshobject = cpd::createObjectFromMeshIO(CPDScene, "F:/VIBE/Resources/GCD_ExampleFiles/liver");
	//auto& meshobject = cpd::createObjectFromMeshIO(CPDScene, "E:/Jose James/Work/VIBE/3DModels/stomach/Stomach_wip_01");

	createForces(CPDScene, forces, isDistributed);
	createTimeIntegrator(CPDScene, timestep, type);

	auto deformableObj = std::make_shared<CPDObject>("meshes");
	deformableObj->setParticleObject(meshobject);
	auto tetraMesh = std::make_shared<imstk::TetrahedralMesh>();
	auto surfMesh = std::make_shared<imstk::SurfaceMesh>();
	deformableObj->convertToMesh(tetraMesh, surfMesh);
	//surfMesh->flipNormals();


	auto one2oneMap = std::make_shared<imstk::TetraTriangleMap>();
	one2oneMap->setMaster(tetraMesh);
	one2oneMap->setSlave(surfMesh);
	//one2oneMap->compute();
	one2oneMap->load("F:/VIBE/Resources\GCD_ExampleFiles/liver.wet");
	//one2oneMap->load("E:/Jose James/Work/VIBE/3DModels/stomach/Stomach_wip_01");

	// set up visual model
	auto material = std::make_shared<RenderMaterial>();
	material->setDisplayMode(RenderMaterial::DisplayMode::SURFACE);
	material->setColor(Color::Red);
	auto surfMeshModel = std::make_shared<VisualModel>(surfMesh);
	surfMeshModel->setRenderMaterial(material);

	// Object & Model
	auto cpdModel = std::make_shared<CPDModel>();
	deformableObj->setDynamicalModel(cpdModel);
	//deformableObj->setVisualGeometry(tetraMesh);
	deformableObj->addVisualModel(surfMeshModel);
	deformableObj->setPhysicsGeometry(tetraMesh);
	deformableObj->setPhysicsToVisualMap(one2oneMap);


	auto& plane = cpd::createPlane1(50, CPDScene, "F:/VIBE/Resources/GCD_ExampleFiles/ss.tet");
	//auto& plane = cpd::createPlane1(50, CPDScene, "D:/CPDFiles/ss.tet");

	auto planeObj = std::make_shared<CPDObject>("Plane");
	planeObj->setParticleObject(plane);
	auto planeMesh = std::make_shared<imstk::SurfaceMesh>();
	planeObj->convertToMesh(planeMesh);
	//plane->setTemporaryPosition()

	// set up visual model
	auto materialplane = std::make_shared<RenderMaterial>();
	materialplane->setDisplayMode(RenderMaterial::DisplayMode::SURFACE);
	materialplane->setColor(Color::LightGray);
	materialplane->backfaceCullingOff();
	auto planesurfMeshModel = std::make_shared<VisualModel>(planeMesh);
	planesurfMeshModel->setRenderMaterial(materialplane);

	// Object & Model
	auto cpdModelPlane = std::make_shared<CPDModel>();
	planeObj->setDynamicalModel(cpdModelPlane);
	planeObj->addVisualModel(planesurfMeshModel);
	planeObj->setPhysicsGeometry(planeMesh);

	//----------------------------
	auto sphereMat = std::make_shared<RenderMaterial>();
	auto texture = std::make_shared<Texture>("F:/VIBE/Resources/Textures/green.png", Texture::Type::DIFFUSE);
	sphereMat->setEmissivity(1.0);
	sphereMat->addTexture(texture);

	auto sphere1 = std::make_shared<imstk::Sphere>();
	//auto sphere1 = std::make_shared<CPDObject>("Sphere");
	//sphereY->setTranslation(Vec3d(0, 0.25, 0));
	sphere1->setRadius(1);
	sphere1->setPosition(10, 0, 0);
	auto sphere1Model = std::make_shared<VisualModel>(sphere1);
	sphere1Model->setRenderMaterial(sphereMat);
	auto sphere1Object = std::make_shared<SceneObject>("Sphere1");
	sphere1Object->addVisualModel(sphere1Model);
	scene->addSceneObject(sphere1Object);
	//----------------------------

	auto sphere2 = std::make_shared<imstk::Sphere>();
	//visualSphere11->setTranslation(Vec3d(-14.4, -1.4, 7.9)*scalexzh);
	sphere2->setRadius(2);
	sphere2->setPosition(-10, 0, 0);
	auto sphere2Model = std::make_shared<VisualModel>(sphere2);
	sphere2Model->setRenderMaterial(sphereMat);
	auto sphere2Object = std::make_shared<CollidingObject>("Sphere2");
	sphere2Object->setCollidingGeometry(sphere2);
	sphere2Object->addVisualModel(sphere2Model);
	scene->addSceneObject(sphere2Object);
	//----------------------------

	//std::vector<unsigned> fixed;

	//fixed.push_back(10);
	//fixed.push_back(11);
	//fixed.push_back(12);
	//fixed.push_back(13);
	//fixed.push_back(14);

	//std::cout << "N_fixed = " << fixed.size() << std::endl;

	//meshobject->setFixedPoints(fixed);
	std::cout << "N = " << meshobject->getParticleCount() << std::endl;

	CPDScene->addCollisionPair(meshobject, plane);
	//CPDScene->addCollisionPair(sphere1, plane);


	scene->addSceneObject(deformableObj);
	scene->addSceneObject(planeObj);

	// Light (white)
	auto whiteLight1 = std::make_shared<imstk::PointLight>("whiteLight1");
	whiteLight1->setPosition(10.0, 5, 0.0);
	whiteLight1->setFocalPoint(imstk::Vec3d(-1, -4, -1));
	whiteLight1->setIntensity(20);
	auto whiteLight2 = std::make_shared<imstk::PointLight>("whiteLight2");
	whiteLight2->setPosition(10.0, -5, 0.0);
	whiteLight2->setFocalPoint(imstk::Vec3d(-1, -4, -1));
	whiteLight2->setIntensity(20);

	// Add in scene
	//scene->addLight(whiteLight1);
	//scene->addLight(whiteLight2);

	scene->getCamera()->setFocalPoint(20, -20, 20);
	scene->getCamera()->setPosition(-15, 4, -20);

	sdk->setActiveScene(scene);
	sdk->start();
}

void testLoadVIBECPD3DModels()
{
	//loadStomach3DModel();
	//loadESGTool();
	testESGToolParts();

}

void loadStomach3DModel()
{
	// SDK and Scene
	auto sdk = std::make_shared<SimulationManager>();
	auto scene = sdk->createNewScene("Stomach");

	// Add IBL Probe
	auto globalIBLProbe = std::make_shared<IBLProbe>(
		iMSTK_DATA_ROOT "/IBL/roomIrradiance.dds",
		iMSTK_DATA_ROOT "/IBL/roomRadiance.dds",
		iMSTK_DATA_ROOT "/IBL/roomBRDF.png");
	scene->setGlobalIBLProbe(globalIBLProbe);

	//Create mesh and scale it down:
	//auto stomachObject = MeshIO::read(iMSTK_DATA_ROOT "/stomach/Stomach_wip_01.OBJ");
	auto stomachObject = MeshIO::read("F:/VIBE/Resources/vibe/3D/Stomach_wip_01.OBJ");

	stomachObject->scale(0.05, Geometry::TransformType::ConcatenateToTransform);	//The stomach mesh apparently is large and needs to be substantially scaled down
	//stomachMesh->rotate(Vec3d(0.0, 1.0, 0.0), PI / 2, Geometry::TransformType::ConcatenateToTransform);
	//stomachObject->getNumVertices();
	std::cout << "Num vetices : " << stomachObject->getNumVertices();
	std::cout << "Num tets : " << stomachObject->getTotalNumberGeometries();
	//Create mesh render material:
	auto stomachMeshRenderMaterial = std::make_shared<RenderMaterial>();
	stomachMeshRenderMaterial->setDisplayMode(RenderMaterial::DisplayMode::SURFACE);
	stomachMeshRenderMaterial->setBackFaceCulling(false);

	//Importing and applying textures:
	//auto stomachBaseTexture = std::make_shared<Texture>(iMSTK_DATA_ROOT "/stomach/Stomach_Textures/Stomach_mat_BaseColor.png", Texture::Type::DIFFUSE);
	//auto stomachNormalTexture = std::make_shared<Texture>(iMSTK_DATA_ROOT "/stomach/Stomach_Textures/Stomach_mat_Normal.png", Texture::Type::NORMAL);
	//auto stomachRoughnessTexture = std::make_shared<Texture>(iMSTK_DATA_ROOT "/stomach/Stomach_Textures/Stomach_mat_Roughness.png", Texture::Type::ROUGHNESS);
	//auto stomachMetalnessTexture = std::make_shared<Texture>(iMSTK_DATA_ROOT "/stomach/Stomach_Textures/Stomach_mat_Specular.png", Texture::Type::METALNESS);
	//auto stomachOcclusionTexture = std::make_shared<Texture>(iMSTK_DATA_ROOT "/stomach/Stomach_Textures/Stomach_mat_Occlusion.png", Texture::Type::AMBIENT_OCCLUSION);

	auto stomachBaseTexture = std::make_shared<Texture>("F:/VIBE/Resources/vibe/Stomach_Textures/Stomach_mat_BaseColor.png", Texture::Type::DIFFUSE);
	auto stomachNormalTexture = std::make_shared<Texture>("F:/VIBE/Resources/vibe/Stomach_Textures/Stomach_mat_Normal.png", Texture::Type::NORMAL);
	auto stomachRoughnessTexture = std::make_shared<Texture>("F:/VIBE/Resources/vibe/Stomach_Textures/Stomach_mat_Roughness.png", Texture::Type::ROUGHNESS);
	auto stomachMetalnessTexture = std::make_shared<Texture>("F:/VIBE/Resources/vibe/Stomach_Textures/Stomach_mat_Specular.png", Texture::Type::METALNESS);
	auto stomachOcclusionTexture = std::make_shared<Texture>("F:/VIBE/Resources/vibe/Stomach_Textures/Stomach_mat_Occlusion.png", Texture::Type::AMBIENT_OCCLUSION);
	stomachMeshRenderMaterial->addTexture(stomachBaseTexture);
	//stomachMeshRenderMaterial->setRoughness(0.132f);
	stomachMeshRenderMaterial->setRoughness(0.080f);
	stomachMeshRenderMaterial->setMetalness(0.0f);
	//stomachMeshRenderMaterial->setEmissivity(0.2f);

	//Create visual object:
	auto stomachMeshObject = std::make_shared<VisualObject>("stomachMesh");
	//Create visual model:
	auto stomachMeshVisualModel = std::make_shared<VisualModel>(stomachObject);

	//Add mesh to object:
	stomachMeshObject->addVisualModel(stomachMeshVisualModel);

	//Add render material onto model:
	stomachMeshVisualModel->setRenderMaterial(stomachMeshRenderMaterial);

	scene->addSceneObject(stomachMeshObject);

	// Position camera
	auto cam = scene->getCamera();
	//cam->setPosition(0, 0.25, 1);
	//cam->setFocalPoint(0, 0.25, 0);
	//cam->setPosition(-15., 0.0, 15.0);
	//cam->setFocalPoint(0, -5, 5);
	//scene->getCamera()->setFocalPoint(10, -2, 5);
	//scene->getCamera()->setPosition(-7.5, 5, 5.0);
	cam->setPosition(0, 0, 20);
	cam->setFocalPoint(0, 0, 0);
	//cam->setFieldOfView(92); 
	//cam->setViewUp(0, 0, 0);

	//Lights
	//auto pointLight = std::make_shared<PointLight>("PointLight");
	//pointLight->setIntensity(0.5f);
	//pointLight->setPosition(0.467f, 0.321f, -0.188f);
	//pointLight->setIntensity(0.1);
	//pointLight->setPosition(5, 0.2, 0.5);
	//pointLight->setPosition(0, 0, 10);
	//pointLight->setColor(pointLightColor);
	//scene->addLight(pointLight);

	// Light (white)
	auto light1 = std::make_shared<DirectionalLight>("light1");
	light1->setIntensity(5);
	light1->setColor(Color(1.0, 0.95, 0.8));
	light1->setFocalPoint(imstk::Vec3d(0, 0, 0));
	//directionalLight->setCastsShadow(false);
	//directionalLight->setShadowRange(1.5);
	//scene->addLight(light1);

	// Light (white)
	auto whiteLight1 = std::make_shared<imstk::DirectionalLight>("whiteLight1");
	whiteLight1->setFocalPoint(imstk::Vec3d(-1, -4, -1));
	whiteLight1->setIntensity(5);
	//whiteLight1->setColor(Color(1.0, 0.95, 0.8));
	whiteLight1->setColor(Color(1.0, 0.8, 0.8));
	scene->addLight(whiteLight1);

	auto whiteLight2 = std::make_shared<imstk::DirectionalLight>("whiteLight2");
	whiteLight2->setFocalPoint(imstk::Vec3d(1, 4, 1));
	whiteLight2->setIntensity(5);
	//whiteLight2->setColor(Color(1.0, 0.95, 0.8));
	whiteLight2->setColor(Color(1.0, 0.8, 0.8));
	scene->addLight(whiteLight2);
	
	auto spotLight1 = std::make_shared<imstk::SpotLight>("spotLight1");
	spotLight1->setPosition(Vec3d(0, 0, 7));
	spotLight1->setIntensity(50);
	//spotLight1->setCastsShadow(false);
	scene->addLight(spotLight1);

	auto spotLight2 = std::make_shared<imstk::SpotLight>("spotLight2");
	spotLight2->setPosition(Vec3d(0, 0, -7));
	spotLight2->setIntensity(50);
	//spotLight2->setCastsShadow(false);
	scene->addLight(spotLight2);
	
	// Run
	sdk->setActiveScene(scene);
	//sdk->getViewer()->setBackgroundColors(Vec3d(0, 0, 0));

   #ifdef iMSTK_USE_Vulkan
	auto viewer = std::dynamic_pointer_cast<VulkanViewer>(sdk->getViewer());
	viewer->setResolution(1000, 800);
	viewer->disableVSync();
	//viewer->enableFullscreen();
   #endif

//sdk->startSimulation(SimulationStatus::RUNNING);
	sdk->start();
	
	//return 0;
}

void loadESGTool()
{
	//Importing and applying textures:
	//auto esgToolNormalTexture = std::make_shared<Texture>(iMSTK_DATA_ROOT "/VIBE/Renders/overstitch01_01.jpg", Texture::Type::NORMAL);
	//auto esgToolBaseTexture = std::make_shared<Texture>(iMSTK_DATA_ROOT "/VIBE/Textures/Stomach_mat_BaseColor.png", Texture::Type::DIFFUSE);
	//==========================

	// SDK and Scene
	auto sdk = std::make_shared<SimulationManager>();
	auto scene = sdk->createNewScene("esgToolMesh");

	// Add IBL Probe
	auto globalIBLProbe = std::make_shared<IBLProbe>(
		iMSTK_DATA_ROOT "/IBL/roomIrradiance.dds",
		iMSTK_DATA_ROOT "/IBL/roomRadiance.dds",
		iMSTK_DATA_ROOT "/IBL/roomBRDF.png");
	scene->setGlobalIBLProbe(globalIBLProbe);

	//Create mesh and scale it down:
	//auto stomachObject = MeshIO::read(iMSTK_DATA_ROOT "/stomach/Stomach_wip_01.OBJ");
	auto esgToolObject = MeshIO::read("F:/VIBE/Resources/vibe/3D/Medical_Instrument_wip01.OBJ");

	esgToolObject->scale(0.05, Geometry::TransformType::ConcatenateToTransform);
	//esgToolObject->rotate(Vec3d(0.0, 1.0, 0.0), PI / 2, Geometry::TransformType::ConcatenateToTransform);
	//esgToolMesh->setTranslation(10, 0, 0);

	//Create mesh render material:
	auto esgToolMeshRenderMaterial = std::make_shared<RenderMaterial>();
	esgToolMeshRenderMaterial->setDisplayMode(RenderMaterial::DisplayMode::SURFACE);
	esgToolMeshRenderMaterial->setBackFaceCulling(false);

	//Importing and applying textures:
	//auto stomachBaseTexture = std::make_shared<Texture>(iMSTK_DATA_ROOT "/stomach/Stomach_Textures/Stomach_mat_BaseColor.png", Texture::Type::DIFFUSE);
	//auto stomachNormalTexture = std::make_shared<Texture>(iMSTK_DATA_ROOT "/stomach/Stomach_Textures/Stomach_mat_Normal.png", Texture::Type::NORMAL);
	//auto stomachRoughnessTexture = std::make_shared<Texture>(iMSTK_DATA_ROOT "/stomach/Stomach_Textures/Stomach_mat_Roughness.png", Texture::Type::ROUGHNESS);
	//auto stomachMetalnessTexture = std::make_shared<Texture>(iMSTK_DATA_ROOT "/stomach/Stomach_Textures/Stomach_mat_Specular.png", Texture::Type::METALNESS);
	//auto stomachOcclusionTexture = std::make_shared<Texture>(iMSTK_DATA_ROOT "/stomach/Stomach_Textures/Stomach_mat_Occlusion.png", Texture::Type::AMBIENT_OCCLUSION);

	auto esgToolBaseTexture = std::make_shared<Texture>("F:/VIBE/Resources/vibe/Overstitch_Textures/Overstitch_low_base_mat_BaseColor.png", Texture::Type::DIFFUSE);
	auto esgToolNormalTexture = std::make_shared<Texture>("F:/VIBE/Resources/vibe/Overstitch_Textures/Overstitch_low_base_mat_Normal.png", Texture::Type::NORMAL);
	auto esgToolRoughnessTexture = std::make_shared<Texture>("F:/VIBE/Resources/vibe/Overstitch_Textures/Overstitch_low_base_mat_Roughness.png", Texture::Type::ROUGHNESS);
	auto esgToolMetalnessTexture = std::make_shared<Texture>("F:/VIBE/Resources/vibe/Overstitch_Textures/Overstitch_low_base_mat_Metallic.png", Texture::Type::METALNESS);
	auto esgToolOcclusionTexture = std::make_shared<Texture>("F:/VIBE/Resources/vibe/Overstitch_Textures/Overstitch_low_base_mat_AO.png", Texture::Type::AMBIENT_OCCLUSION);
	//esgToolMeshRenderMaterial->addTexture(esgToolBaseTexture);
	//esgToolMeshRenderMaterial->addTexture(esgToolNormalTexture);
	//esgToolMeshRenderMaterial->addTexture(esgToolRoughnessTexture);
	//esgToolMeshRenderMaterial->addTexture(esgToolMetalnessTexture);
	//esgToolMeshRenderMaterial->addTexture(esgToolOcclusionTexture);

	//stomachMeshRenderMaterial->setRoughness(0.132f);
	//esgToolMeshRenderMaterial->setRoughness(0.080f);
	//esgToolMeshRenderMaterial->setMetalness(0.0f);
	//stomachMeshRenderMaterial->setEmissivity(0.2f);

	//Create visual object:
	auto esgToolMeshObject = std::make_shared<VisualObject>("esgToolMesh");
	//Create visual model:
	auto esgToolMeshVisualModel = std::make_shared<VisualModel>(esgToolObject);

	//Add mesh to object:
	esgToolMeshObject->addVisualModel(esgToolMeshVisualModel);

	//Add render material onto model:
	esgToolMeshVisualModel->setRenderMaterial(esgToolMeshRenderMaterial);

	scene->addSceneObject(esgToolMeshObject);

	// Position camera
	auto cam = scene->getCamera();
	//cam->setPosition(0, 0.25, 1);
	//cam->setFocalPoint(0, 0.25, 0);
	//cam->setPosition(-15., 0.0, 15.0);
	//cam->setFocalPoint(0, -5, 5);
	//scene->getCamera()->setFocalPoint(10, -2, 5);
	//scene->getCamera()->setPosition(-7.5, 5, 5.0);
	cam->setPosition(0, 0, 20);
	cam->setFocalPoint(0, 0, 0);
	//cam->setFieldOfView(92); 
	//cam->setViewUp(0, 0, 0);

	//Lights
	//auto pointLight = std::make_shared<PointLight>("PointLight");
	//pointLight->setIntensity(0.5f);
	//pointLight->setPosition(0.467f, 0.321f, -0.188f);
	//pointLight->setIntensity(0.1);
	//pointLight->setPosition(5, 0.2, 0.5);
	//pointLight->setPosition(0, 0, 10);
	//pointLight->setColor(pointLightColor);
	//scene->addLight(pointLight);

	// Light (white)
	auto light1 = std::make_shared<DirectionalLight>("light1");
	light1->setIntensity(5);
	//light1->setColor(Color(1.0, 0.95, 0.8));
	light1->setFocalPoint(imstk::Vec3d(0, 0, 0));
	//directionalLight->setCastsShadow(false);
	//directionalLight->setShadowRange(1.5);
	//scene->addLight(light1);

	// Light (white)
	auto whiteLight1 = std::make_shared<imstk::DirectionalLight>("whiteLight1");
	whiteLight1->setFocalPoint(imstk::Vec3d(-1, -4, -1));
	whiteLight1->setIntensity(5);
	//whiteLight1->setColor(Color(1.0, 0.95, 0.8));
	//whiteLight1->setColor(Color(1.0, 0.8, 0.8));
	scene->addLight(whiteLight1);

	auto whiteLight2 = std::make_shared<imstk::DirectionalLight>("whiteLight2");
	whiteLight2->setFocalPoint(imstk::Vec3d(1, 4, 1));
	whiteLight2->setIntensity(5);
	//whiteLight2->setColor(Color(1.0, 0.95, 0.8));
	//whiteLight2->setColor(Color(1.0, 0.8, 0.8));
	scene->addLight(whiteLight2);


	auto spotLight1 = std::make_shared<imstk::SpotLight>("spotLight1");
	spotLight1->setPosition(Vec3d(0, 0, 7));
	spotLight1->setIntensity(50);
	//spotLight1->setCastsShadow(false);
	//scene->addLight(spotLight1);

	auto spotLight2 = std::make_shared<imstk::SpotLight>("spotLight2");
	spotLight2->setPosition(Vec3d(0, 0, -7));
	spotLight2->setIntensity(50);
	//spotLight2->setCastsShadow(false);
	//scene->addLight(spotLight2);


	// Run
	sdk->setActiveScene(scene);
	//sdk->getViewer()->setBackgroundColors(Vec3d(0, 0, 0));

   #ifdef iMSTK_USE_Vulkan
	auto viewer = std::dynamic_pointer_cast<VulkanViewer>(sdk->getViewer());
	viewer->setResolution(1000, 800);
	viewer->disableVSync();
	//viewer->enableFullscreen();
    #endif

//sdk->startSimulation(SimulationStatus::RUNNING);
	sdk->start();
	
	//return 0;
}

void testESGToolParts()
{
	//Importing and applying textures:
	//auto esgToolNormalTexture = std::make_shared<Texture>(iMSTK_DATA_ROOT "/VIBE/Renders/overstitch01_01.jpg", Texture::Type::NORMAL);
	//auto esgToolBaseTexture = std::make_shared<Texture>(iMSTK_DATA_ROOT "/VIBE/Textures/Stomach_mat_BaseColor.png", Texture::Type::DIFFUSE);
	//==========================

	// SDK and Scene
	auto sdk = std::make_shared<SimulationManager>();
	auto scene = sdk->createNewScene("esgToolMesh");

	// Add IBL Probe
	auto globalIBLProbe = std::make_shared<IBLProbe>(
		iMSTK_DATA_ROOT "/IBL/roomIrradiance.dds",
		iMSTK_DATA_ROOT "/IBL/roomRadiance.dds",
		iMSTK_DATA_ROOT "/IBL/roomBRDF.png");
	scene->setGlobalIBLProbe(globalIBLProbe);

	//Create mesh and scale it down:
	//auto stomachObject = MeshIO::read(iMSTK_DATA_ROOT "/stomach/Stomach_wip_01.OBJ");
	auto esgToolObject = MeshIO::read("F:/VIBE/Resources/vibe/3D_New/tool_2/overstitch.obj");

	esgToolObject->scale(0.05, Geometry::TransformType::ConcatenateToTransform);
	//esgToolObject->rotate(Vec3d(0.0, 1.0, 0.0), PI / 2, Geometry::TransformType::ConcatenateToTransform);
	//esgToolMesh->setTranslation(10, 0, 0);

	//Create mesh render material:
	auto esgToolMeshRenderMaterial = std::make_shared<RenderMaterial>();
	esgToolMeshRenderMaterial->setDisplayMode(RenderMaterial::DisplayMode::SURFACE);
	esgToolMeshRenderMaterial->setBackFaceCulling(false);

	//Importing and applying textures:
	//auto stomachBaseTexture = std::make_shared<Texture>(iMSTK_DATA_ROOT "/stomach/Stomach_Textures/Stomach_mat_BaseColor.png", Texture::Type::DIFFUSE);
	//auto stomachNormalTexture = std::make_shared<Texture>(iMSTK_DATA_ROOT "/stomach/Stomach_Textures/Stomach_mat_Normal.png", Texture::Type::NORMAL);
	//auto stomachRoughnessTexture = std::make_shared<Texture>(iMSTK_DATA_ROOT "/stomach/Stomach_Textures/Stomach_mat_Roughness.png", Texture::Type::ROUGHNESS);
	//auto stomachMetalnessTexture = std::make_shared<Texture>(iMSTK_DATA_ROOT "/stomach/Stomach_Textures/Stomach_mat_Specular.png", Texture::Type::METALNESS);
	//auto stomachOcclusionTexture = std::make_shared<Texture>(iMSTK_DATA_ROOT "/stomach/Stomach_Textures/Stomach_mat_Occlusion.png", Texture::Type::AMBIENT_OCCLUSION);

	auto esgToolBaseTexture = std::make_shared<Texture>("F:/VIBE/Resources/vibe/Overstitch_Textures/Overstitch_low_base_mat_BaseColor.png", Texture::Type::DIFFUSE);
	auto esgToolNormalTexture = std::make_shared<Texture>("F:/VIBE/Resources/vibe/Overstitch_Textures/Overstitch_low_base_mat_Normal.png", Texture::Type::NORMAL);
	auto esgToolRoughnessTexture = std::make_shared<Texture>("F:/VIBE/Resources/vibe/Overstitch_Textures/Overstitch_low_base_mat_Roughness.png", Texture::Type::ROUGHNESS);
	auto esgToolMetalnessTexture = std::make_shared<Texture>("F:/VIBE/Resources/vibe/Overstitch_Textures/Overstitch_low_base_mat_Metallic.png", Texture::Type::METALNESS);
	auto esgToolOcclusionTexture = std::make_shared<Texture>("F:/VIBE/Resources/vibe/Overstitch_Textures/Overstitch_low_base_mat_AO.png", Texture::Type::AMBIENT_OCCLUSION);
	//esgToolMeshRenderMaterial->addTexture(esgToolBaseTexture);
	//esgToolMeshRenderMaterial->addTexture(esgToolNormalTexture);
	//esgToolMeshRenderMaterial->addTexture(esgToolRoughnessTexture);
	//esgToolMeshRenderMaterial->addTexture(esgToolMetalnessTexture);
	//esgToolMeshRenderMaterial->addTexture(esgToolOcclusionTexture);

	//stomachMeshRenderMaterial->setRoughness(0.132f);
	//esgToolMeshRenderMaterial->setRoughness(0.080f);
	//esgToolMeshRenderMaterial->setMetalness(0.0f);
	//stomachMeshRenderMaterial->setEmissivity(0.2f);

	//Create visual object:
	auto esgToolMeshObject = std::make_shared<VisualObject>("esgToolMesh");
	//Create visual model:
	auto esgToolMeshVisualModel = std::make_shared<VisualModel>(esgToolObject);

	//Add mesh to object:
	esgToolMeshObject->addVisualModel(esgToolMeshVisualModel);

	//Add render material onto model:
	esgToolMeshVisualModel->setRenderMaterial(esgToolMeshRenderMaterial);

	scene->addSceneObject(esgToolMeshObject);

	// Position camera
	auto cam = scene->getCamera();
	//cam->setPosition(0, 0.25, 1);
	//cam->setFocalPoint(0, 0.25, 0);
	//cam->setPosition(-15., 0.0, 15.0);
	//cam->setFocalPoint(0, -5, 5);
	//scene->getCamera()->setFocalPoint(10, -2, 5);
	//scene->getCamera()->setPosition(-7.5, 5, 5.0);
	cam->setPosition(0, 0, 20);
	cam->setFocalPoint(0, 0, 0);
	//cam->setFieldOfView(92); 
	//cam->setViewUp(0, 0, 0);

	//Lights
	//auto pointLight = std::make_shared<PointLight>("PointLight");
	//pointLight->setIntensity(0.5f);
	//pointLight->setPosition(0.467f, 0.321f, -0.188f);
	//pointLight->setIntensity(0.1);
	//pointLight->setPosition(5, 0.2, 0.5);
	//pointLight->setPosition(0, 0, 10);
	//pointLight->setColor(pointLightColor);
	//scene->addLight(pointLight);

	// Light (white)
	auto light1 = std::make_shared<DirectionalLight>("light1");
	light1->setIntensity(5);
	//light1->setColor(Color(1.0, 0.95, 0.8));
	light1->setFocalPoint(imstk::Vec3d(0, 0, 0));
	//directionalLight->setCastsShadow(false);
	//directionalLight->setShadowRange(1.5);
	//scene->addLight(light1);

	// Light (white)
	auto whiteLight1 = std::make_shared<imstk::DirectionalLight>("whiteLight1");
	whiteLight1->setFocalPoint(imstk::Vec3d(-1, -4, -1));
	whiteLight1->setIntensity(5);
	//whiteLight1->setColor(Color(1.0, 0.95, 0.8));
	//whiteLight1->setColor(Color(1.0, 0.8, 0.8));
	scene->addLight(whiteLight1);

	auto whiteLight2 = std::make_shared<imstk::DirectionalLight>("whiteLight2");
	whiteLight2->setFocalPoint(imstk::Vec3d(1, 4, 1));
	whiteLight2->setIntensity(5);
	//whiteLight2->setColor(Color(1.0, 0.95, 0.8));
	//whiteLight2->setColor(Color(1.0, 0.8, 0.8));
	scene->addLight(whiteLight2);


	auto spotLight1 = std::make_shared<imstk::SpotLight>("spotLight1");
	spotLight1->setPosition(Vec3d(0, 0, 7));
	spotLight1->setIntensity(50);
	//spotLight1->setCastsShadow(false);
	//scene->addLight(spotLight1);

	auto spotLight2 = std::make_shared<imstk::SpotLight>("spotLight2");
	spotLight2->setPosition(Vec3d(0, 0, -7));
	spotLight2->setIntensity(50);
	//spotLight2->setCastsShadow(false);
	//scene->addLight(spotLight2);


	// Run
	sdk->setActiveScene(scene);
	//sdk->getViewer()->setBackgroundColors(Vec3d(0, 0, 0));

   #ifdef iMSTK_USE_Vulkan
	auto viewer = std::dynamic_pointer_cast<VulkanViewer>(sdk->getViewer());
	viewer->setResolution(1000, 800);
	viewer->disableVSync();
	//viewer->enableFullscreen();
   #endif

//sdk->startSimulation(SimulationStatus::RUNNING);
	sdk->start();
	
	//return 0;
}