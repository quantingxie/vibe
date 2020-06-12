#include "imstkcpdAPIUtilities.h"
#include "imstkcpdMeshIO.h"

namespace cpd
{
  void createForces(const ScenePtr p_scene, const std::vector<std::array<double, 3>>& forces,
    const std::vector<bool>& isDistributed)
  {
    size_t nbrForces = forces.size();

    for (size_t i = 0; i < nbrForces; i++)
    {
      ExternalForcePtr force = std::make_shared<ExternalForce>();
      force->setForce(forces[i][0], forces[i][1], forces[i][2]);
      p_scene->addDistributedForce(force);
    }

    return;
  }

  void createTimeIntegrator(const ScenePtr p_scene, const double timestep,
    const BaseTimeIntegrator::Type type)
  {
    switch (type)
    {
    case cpd::BaseTimeIntegrator::Type::SIE:
    {
      auto& integrator = std::make_shared<SemiImplicitEulerIntegrator>();
      p_scene->getExternalForceSolver()->setTimeIntegrator(integrator);
      break;
    }
    case cpd::BaseTimeIntegrator::Type::CD:
    {
      auto& integrator = std::make_shared<CentralDifferenceIntegrator>();
      p_scene->getExternalForceSolver()->setTimeIntegrator(integrator);
      break;
    }
    default:
      break;
    }

    p_scene->setTimeStepSize(timestep);
  }

  ParticleObjectPtr createPlane(double p_width, const ScenePtr p_scene, const std::string & p_resourceDir)
  {
    const double lengthX = 30.0; // read
    const double lengthY = p_width; // read
    const double startX = -lengthX / 2.0; // read
    const double startY = -lengthY / 2.0; // read
    const int nRowX = 3; // read
    const int nRowY = 2; // read
    const double rho = 800.0; // read
    const double thickness = 1.0; // read
    std::vector<unsigned> fixed;
    for (int i = 0; i < nRowY*nRowX; ++i) // read
    {
      fixed.push_back(i);
    }
    const double youngs = 0.5E7; // read
    const double mu = 0.4; // read
    ConstraintBase::ConstraintType constraintType = ConstraintBase::ConstraintType(ConstraintBase::Type::CST2D,
      ConstraintBase::SubType::XCPD); // read

    const double area = lengthY * lengthX;
    const double element_area = 0.5*area / ((nRowX - 1)*(nRowY - 1));
    const double nodeMass = element_area * rho*thickness / 3.0;
    const double dx = lengthX / (double)(nRowX - 1);
    const double dy = lengthY / (double)(nRowY - 1);

    StdVectorOfVec3d vertList;
    vertList.resize(nRowX*nRowY);
    for (int i = 0; i < nRowY; ++i)
      for (int j = 0; j < nRowX; j++)
        vertList[i*nRowX + j] = Vec3d(startY + (double)dy*i, -6.0 /*+ 0.1*(double)dy*i*/, startX + (double)dx*j);

    std::vector<double> massList;
    massList.resize(nRowX*nRowY);
    for (unsigned i = 0; i < massList.size(); i++)
      massList[i] = 0.0;

    using TriangleArray = std::array<size_t, 3>;
    std::vector<TriangleArray> triangles;
    for (std::size_t i = 0; i < nRowY - 1; ++i)
    {
      for (std::size_t j = 0; j < nRowX - 1; j++)
      {
        TriangleArray tri[2];
        tri[0] = { { i *nRowX + j , i*nRowX + j + 1, (i + 1)*nRowX + j } };
        tri[1] = { { (i + 1)*nRowX + j + 1, (i + 1)*nRowX + j, i*nRowX + j + 1} };
        triangles.push_back(tri[0]);
        triangles.push_back(tri[1]);

        for (unsigned m = 0; m < 3; m++)
        {
          massList[tri[0][m]] += nodeMass;
          massList[tri[1][m]] += nodeMass;
        }
      }
    }

    SurfaceMeshPtr surfaceMesh = std::make_shared<SurfaceMesh>(vertList, triangles);

    ParticleObjectPtr particleObject = std::make_shared<ParticleObject>(surfaceMesh);

    particleObject->setProperties(youngs, mu);
    particleObject->addConstraintType(constraintType);

    ParticlePropertyPtr property = std::make_shared<ParticleProperty>();
    particleObject->setParticleProperty(property);

    particleObject->setMassList(massList);

    particleObject->setFixedPoints(fixed);

    p_scene->addObject(particleObject);

    return particleObject;
  }

  ParticleObjectPtr createPlane1(double p_width, const ScenePtr p_scene, const std::string & p_resourceDir)
  {
    const double lengthX = /*30.0*/(p_width > 40.0) ? 50.0 : p_width; // read
    const double lengthY = p_width; // read
    const double startX = -lengthX / 2.0; // read
    const double startY = -lengthY / 2.0; // read   
    const int nRowX = (p_width > 20.0) ? 4 : 5; // read
    const int nRowY = (p_width > 20.0) ? 4 : 5; // read
    const double rho = 25.5000; // read
    const double thickness = 1.0; // read
    std::vector<unsigned> fixed;

    for (int i = 0; i < nRowY*nRowX*(p_width > 10.0); ++i) // read
    {
      fixed.push_back(i);
    }

    //fixed.push_back(0);
    //fixed.push_back((nRowY - 1) *nRowX);
    //fixed.push_back(nRowX - 1);
    //fixed.push_back(nRowY *nRowX - 1);

    const double youngs = 0.8E5; // read
    const double mu = 0.4; // read
    ConstraintBase::ConstraintType constraintType = ConstraintBase::ConstraintType(ConstraintBase::Type::Distance, /*ConstraintBase::Type::CST2D,*/
      ConstraintBase::SubType::XCPD); // read

    const double area = lengthY * lengthX;
    const double element_area = 0.5*area / ((nRowX - 1)*(nRowY - 1));
    const double nodeMass = element_area * rho*thickness / 3.0;
    const double dx = lengthX / (double)(nRowX - 1);
    const double dy = lengthY / (double)(nRowY - 1);

    const double height = (p_width > 20.0) ? -00.0 : -20.0; // read
    StdVectorOfVec3d vertList;
    vertList.resize(nRowX*nRowY);
    for (int i = 0; i < nRowY; ++i)
      for (int j = 0; j < nRowX; j++)
        vertList[i*nRowX + j] = Vec3d(startY + (double)dy*i, -8.0 - height + 0*0.02*(double)dy*abs(i - 3 * nRowY / 4), startX + (double)dx*j - 5.0);

    std::vector<double> massList;
    massList.resize(nRowX*nRowY);
    for (unsigned i = 0; i < massList.size(); i++)
      massList[i] = nodeMass/*0.0*/;

    using TriangleArray = std::array<size_t, 3>;
    std::vector<TriangleArray> triangles;
    for (std::size_t i = 0; i < nRowY - 1; ++i)
    {
      for (std::size_t j = 0; j < nRowX - 1; j++)
      {
        TriangleArray tri[2];
        tri[0] = { { i *nRowX + j , i*nRowX + j + 1, (i + 1)*nRowX + j, } };
        tri[1] = { { (i + 1)*nRowX + j,   i*nRowX + j + 1 ,(i + 1)*nRowX + j + 1} };
        triangles.push_back(tri[0]);
        triangles.push_back(tri[1]);

        //for (unsigned m = 0; m < 3; m++)
        //{
        //  massList[tri[0][m]] += nodeMass;
        //  massList[tri[1][m]] += nodeMass;
        //}
      }
    }

    SurfaceMeshPtr surfaceMesh = std::make_shared<SurfaceMesh>(vertList, triangles);

    ParticleObjectPtr particleObject = std::make_shared<ParticleObject>(surfaceMesh);

    particleObject->setProperties(youngs, mu);
    particleObject->addConstraintType(constraintType);

    ParticlePropertyPtr property = std::make_shared<ParticleProperty>();
    particleObject->setParticleProperty(property);

    particleObject->setMassList(massList);

    particleObject->setFixedPoints(fixed);

    p_scene->addObject(particleObject);

    return particleObject;
  }


  ParticleObjectPtr createCylinder(const ScenePtr p_scene, const std::string & p_resourceDir)
  {
    const double length = 30.0; // read
    const double radius = 1.0; // read
    const Vec3d position(0.0, 6, 0.0); // read
    const Vec3d normal(1.0, 0.0, 0.0); // read   

    std::vector<unsigned> fixed;
    //fixed.push_back(0);
    //fixed.push_back(1);

    std::vector<double> massList;
    massList.push_back(1.0);
    massList.push_back(1.0);

    ConstraintBase::ConstraintType constraintType = ConstraintBase::ConstraintType(ConstraintBase::Type::Distance,
      ConstraintBase::SubType::XCPD); // read

    /*CylinderPtr cylinder = std::make_shared<Cylinder>(radius, length, position, normal);

    ParticleObjectPtr particleObject = std::make_shared<ParticleObject>(cylinder);

    particleObject->setProperties(1e7);
    particleObject->addConstraintType(constraintType);

    ParticlePropertyPtr property = std::make_shared<ParticleProperty>();
    particleObject->setParticleProperty(property);

    particleObject->setMassList(massList);

    particleObject->setFixedPoints(fixed);

    p_scene->addObject(particleObject);

    return particleObject;*/
    return std::make_shared<ParticleObject>();
  }

  ParticleObjectPtr create2DBeam(ScenePtr p_scene, const std::string& p_resourceDir)
  {
    const double lengthX = 20.0; // read
    const double lengthY = 10.0; // read
    const int factor = 10;
    const int nRowX = 2 * factor + 1; // read

    const int nRowY = factor + 1; // read
    const double rho = 800.0; // read
    const double thickness = 1.0; // read
    std::vector<unsigned> fixed;
    for (int i = 0; i < nRowY; ++i) // read
    {
      fixed.push_back(i *nRowX);
    }
    const double youngs = 5E6; // read
    const double mu = 0.3; // read
    ConstraintBase::ConstraintType constraintType = ConstraintBase::ConstraintType(ConstraintBase::Type::CST2D,
      ConstraintBase::SubType::XCPD); // read

    const double area = lengthY * lengthX;
    const double element_area = 0.5*area / ((nRowX - 1)*(nRowY - 1));
    const double nodeMass = element_area * rho*thickness / 3.0;
    const double dx = lengthX / (double)(nRowX - 1);
    const double dy = lengthY / (double)(nRowY - 1);

    StdVectorOfVec3d vertList;
    vertList.resize(nRowX*nRowY);
    for (int i = 0; i < nRowY; ++i)
      for (int j = 0; j < nRowX; j++)
        vertList[i*nRowX + j] = Vec3d((double)dx*j, (double)dy*i, 1.0);

    std::vector<double> massList;
    massList.resize(nRowX*nRowY);
    for (unsigned i = 0; i < massList.size(); i++)
      massList[i] = 0.0;

    using TriangleArray = std::array<size_t, 3>;
    std::vector<TriangleArray> triangles;
    for (std::size_t i = 0; i < nRowY - 1; ++i)
    {
      for (std::size_t j = 0; j < nRowX - 1; j++)
      {
        TriangleArray tri[2];
        tri[0] = { { i *nRowX + j , i*nRowX + j + 1, (i + 1)*nRowX + j + 1 } };
        tri[1] = { { (i + 1)*nRowX + j + 1, (i + 1)*nRowX + j, i*nRowX + j } };
        triangles.push_back(tri[0]);
        triangles.push_back(tri[1]);

        for (unsigned m = 0; m < 3; m++)
        {
          massList[tri[0][m]] += nodeMass;
          massList[tri[1][m]] += nodeMass;
        }
      }
    }

    SurfaceMeshPtr surfaceMesh = std::make_shared<SurfaceMesh>(vertList, triangles);

    ParticleObjectPtr particleObject = std::make_shared<ParticleObject>(surfaceMesh);

    particleObject->setProperties(youngs, mu);
    particleObject->addConstraintType(constraintType);

    ParticlePropertyPtr property = std::make_shared<ParticleProperty>();
    particleObject->setParticleProperty(property);

    particleObject->setMassList(massList);

    particleObject->setFixedPoints(fixed);

    p_scene->addObject(particleObject);

    return particleObject;
  }

  ParticleObjectPtr create3DBeam(const ScenePtr p_scene, const std::string & p_resourceDir)
  {
    /* this part should be from a file reader*/
    const double length = 20;
    const double width = 30;
    const double height = 100.0;

    const int factor = 3;

    const int nRows1 = 1 * factor;
    const int nRows2 = 2 * factor;
    const int nRows3 = 5 * factor;

    const double rho = 7800/*7800*/;

    std::vector<unsigned> fixed;
    for (unsigned i = 0; i < nRows1*nRows2; i++)
      fixed.push_back(i);

    const double youngs = 8.0E8; // read
    const double mu = 0.4; // read
    ConstraintBase::ConstraintType constraintType = ConstraintBase::ConstraintType(ConstraintBase::Type::CST3D,
      ConstraintBase::SubType::XCPD); // read

    /* this part should be from a file reader*/

    const double volume = length * width*height;
    const double element_volume = volume / (6 * (nRows1 - 1)*(nRows2 - 1)*(nRows3 - 1));
    const double nodeMass = element_volume * rho / 4.0;

    const double dx = length / (double)(nRows1 - 1);
    const double dy = width / (double)(nRows2 - 1);
    const double dz = height / (double)(nRows3 - 1);

    StdVectorOfVec3d vertList;
    vertList.resize(nRows1*nRows2*nRows3);
    for (int i = 0; i < nRows1; i++)
      for (int j = 0; j < nRows2; j++)
        for (int k = 0; k < nRows3; k++)
          vertList[k*nRows1 *nRows2 + j * nRows1 + i] = Vec3d((double)dx*i, (double)dy*j, (double)dz*k);

    std::vector<double> massList;
    massList.resize(vertList.size());
    for (unsigned i = 0; i < massList.size(); i++)
      massList[i] = 0.0;

    using TetraArray = std::array<size_t, 4>;
    std::vector<TetraArray> tetras;
    for (std::size_t i = 0; i < nRows1 - 1; i++)
    {
      for (std::size_t j = 0; j < nRows2 - 1; j++)
      {
        for (std::size_t k = 0; k < nRows3 - 1; k++)
        {
          unsigned idx[8] =
          { k*nRows1 *nRows2 + j * nRows1 + i, k*nRows1 *nRows2 + j * nRows1 + i + 1,
            k*nRows1 *nRows2 + (j + 1)*nRows1 + i, k*nRows1 *nRows2 + (j + 1)*nRows1 + i + 1,
            (k + 1)*nRows1 *nRows2 + j * nRows1 + i, (k + 1)*nRows1 *nRows2 + j * nRows1 + i + 1,
            (k + 1)*nRows1 *nRows2 + (j + 1)*nRows1 + i, (k + 1)*nRows1 *nRows2 + (j + 1)*nRows1 + i + 1
          };

          TetraArray tri[6];
          tri[0] = { idx[0],idx[1],idx[2],idx[6] };
          tri[1] = { idx[0],idx[4],idx[1],idx[6] };
          tri[2] = { idx[4],idx[5],idx[1],idx[6] };
          tri[3] = { idx[2],idx[1],idx[3],idx[6] };
          tri[4] = { idx[1],idx[5],idx[7],idx[6] };
          tri[5] = { idx[1],idx[7],idx[3],idx[6] };

          for (std::size_t t = 0; t < 6; t++)
            tetras.push_back(tri[t]);

          for (unsigned m = 0; m < 4; m++)
            for (unsigned n = 0; n < 6; n++)
              massList[tri[n][m]] += nodeMass;

        }
      }
    }

    TetrahedronMeshPtr tetraMesh = std::make_shared<TetrahedronMesh>(vertList, tetras);

    ParticleObjectPtr particleObject = std::make_shared<ParticleObject>(tetraMesh);

    particleObject->setProperties(youngs, mu);
    particleObject->addConstraintType(constraintType);

    ParticlePropertyPtr property = std::make_shared<ParticleProperty>();
    particleObject->setParticleProperty(property);

    particleObject->setMassList(massList);

    particleObject->setFixedPoints(fixed);

    p_scene->addObject(particleObject);

    return particleObject;
  }

  ParticleObjectPtr createChain(const ScenePtr p_scene, const std::string & p_resourceDir)
  {
    /* this part should be from a file reader*/
    const double lengthX = 2.5;
    const int nRowX = 21;
    const double rho = 8.0;
    std::vector<unsigned> fixed;
    for (int i = 0; i < 1; ++i)
    {
      fixed.push_back(i);
    }
    const double stiffness = 1.0E4;
    ConstraintBase::ConstraintType constraintType = ConstraintBase::ConstraintType(
      ConstraintBase::Type::Distance, ConstraintBase::SubType::XPBD);
    /* this part should be from a file reader*/

    const double thickness = 1.0;
    const double area = lengthX;
    const double nodeMass = area * rho*thickness / (nRowX - 1);
    const double dx = lengthX / (double)(nRowX - 1);

    StdVectorOfVec3d vertList;
    std::vector<double> massList;

    vertList.resize(nRowX);
    massList.resize(nRowX);

    for (int i = 0; i < nRowX; i++)
    {
      vertList[i] = Vec3d((double)dx*i, 0.0, 0.0);
      massList[i] = nodeMass;
    }
    massList[nRowX - 1] += 2000.0;

    ParticleObjectPtr particleObject = std::make_shared<ParticleObject>(vertList);

    particleObject->addConstraintType(constraintType);
    particleObject->setProperties(stiffness);

    ParticlePropertyPtr property = std::make_shared<ParticleProperty>();
    particleObject->setParticleProperty(property);

    particleObject->setMassList(massList);
    particleObject->setFixedPoints(fixed);

    p_scene->addObject(particleObject);

    return particleObject;


    return std::make_shared<ParticleObject>();
  }

  ParticleObjectPtr createHarmonicOscillator(const ScenePtr p_scene, const std::string & p_resourceDir)
  {
    /* this part should be from a file reader*/
    const double lengthX = 1.0;
    const int nRowX = 2;

    std::vector<unsigned> fixed;
    for (int i = 0; i < 1; ++i)
    {
      fixed.push_back(i);
    }

    ConstraintBase::ConstraintType constraintType = ConstraintBase::ConstraintType(
      ConstraintBase::Type::Distance, ConstraintBase::SubType::XPBD);
    /* this part should be from a file reader*/

    const double nodeMass = 1.0;
    const double dx = 1.5;

    StdVectorOfVec3d vertList;
     std::vector<double> massList;

    vertList.resize(nRowX);
    massList.resize(nRowX);

    for (int i = 0; i < nRowX; i++)
    {
      vertList[i] = Vec3d((double)dx*i, 0.0, 0.0);
      massList[i] = nodeMass;
    }

    ParticleObjectPtr particleObject = std::make_shared<ParticleObject>(vertList);

    particleObject->addConstraintType(constraintType);
    particleObject->setProperties(1000);

    ParticlePropertyPtr property = std::make_shared<ParticleProperty>();
    particleObject->setParticleProperty(property);

    particleObject->setMassList(massList);
    particleObject->setFixedPoints(fixed);

    StdVectorOfVec3d vertInitList;
    vertInitList.resize(2);
    vertInitList[0] = Vec3d(0.0, 0.0, 0.0);
    vertInitList[1] = Vec3d(1.0, 0.0, 0.0);
    particleObject->setInitPositions(vertInitList);

    p_scene->addObject(particleObject);

    return particleObject;

  }

  ParticleObjectPtr createObjectFromMeshIO(const ScenePtr p_scene, const std::string & p_resourceDir)
  {
    std::string root = std::string("F:/iMSTK/build15Release/Innerbuild/ExternalData/Data");

    /* this part should be from a file reader*/
    const double rho = 7800.0;
    std::vector<unsigned> fixed;
    for (int i = 0; i < 0; ++i)
    {
      fixed.push_back(i);
    }
    const double youngs = 5E4; // read
    const double mu = 0.4; // read
    ConstraintBase::ConstraintType constraintType = ConstraintBase::ConstraintType(ConstraintBase::Type::CST3D,
      ConstraintBase::SubType::XCPD); // read
    /* this part should be from a file reader*/

    auto tetraMesh = TETMeshIO::read(p_resourceDir + ".tet");
    auto surfMesh = SURFMeshIO::read(p_resourceDir + ".srf");

    double scale = 1.0;
    tetraMesh->scale(scale);
    surfMesh->scale(scale);

    tetraMesh->setAttachedSurfMesh(surfMesh);

    //for (int i = 0; i < tetraMesh->getNumVertices(); ++i)
    //{
    //  fixed.push_back(i);
    //}

    std::vector<double> massList;
    size_t n = tetraMesh->getNumVertices();
    massList.resize(n);
    double nodeMass = 0.1;
    for (int i = 0; i < n; i++)
    {
      massList[i] = nodeMass;
    }

    ParticleObjectPtr particleObject = std::make_shared<ParticleObject>(tetraMesh);

    particleObject->addConstraintType(constraintType);
    particleObject->setProperties(youngs, mu);

    ParticlePropertyPtr property = std::make_shared<ParticleProperty>();
    particleObject->setParticleProperty(property);

    particleObject->setMassList(massList);
    particleObject->setFixedPoints(fixed);

    p_scene->addObject(particleObject);

    return particleObject;
  }

  ParticleObjectPtr createObjectFromTetraMesh(const ScenePtr p_scene, const TetrahedronMeshPtr p_tetraMesh)
  {
    return ParticleObjectPtr();
  }

  ParticleObjectPtr createObjectFromSurfMesh(const ScenePtr p_scene, const std::string & p_resourceDir)
  {
    /* this part should be from a file reader*/
    const double rho = 7800.0;
    std::vector<unsigned> fixed;
    for (int i = 0; i < 0; ++i)
    {
      fixed.push_back(i);
    }
    const double youngs = 1E6; // read
    const double mu = 0.4; // read
    ConstraintBase::ConstraintType constraintType = ConstraintBase::ConstraintType(ConstraintBase::Type::Distance,
      ConstraintBase::SubType::XCPD); // read
    /* this part should be from a file reader*/
    
    auto surfMesh = SURFMeshIO::read(p_resourceDir + ".srf");

    for (int i = 0; i < surfMesh->getNumVertices(); ++i)
    {
      fixed.push_back(i);
    }

    std::vector<double> massList;
    size_t n = surfMesh->getNumVertices();
    massList.resize(n);
    double nodeMass = 1;
    for (int i = 0; i < n; i++)
    {
      massList[i] = nodeMass;
    }

    ParticleObjectPtr particleObject = std::make_shared<ParticleObject>(surfMesh);

    particleObject->addConstraintType(constraintType);
    particleObject->setProperties(youngs, mu);

    ParticlePropertyPtr property = std::make_shared<ParticleProperty>();
    particleObject->setParticleProperty(property);

    particleObject->setMassList(massList);
    particleObject->setFixedPoints(fixed);

    p_scene->addObject(particleObject);

    return particleObject;
  }

  void writeTetra(TetrahedronMeshPtr p_mesh)
  {
    TETMeshIO::write(p_mesh, "");
  }

  void writeSurf(SurfaceMeshPtr p_mesh)
  {
    SURFMeshIO::write(p_mesh, "");
  }

  // added by Jose James
  ParticleObjectPtr create2DPlane(ScenePtr p_scene, const std::string& p_resourceDir)
  {
	  //std::cout << " Inside create2DPlane \n";
	  const double lengthX = 20.0; // read
	  const double lengthY = 10.0; // read
	  const int factor = 10;
	  const int nRowX = 2 * factor + 1; // read
	  const int nRowY = factor + 1; // read
	  const double rho = 800.0; // read
	  const double thickness = 1.0; // read
	  std::vector<unsigned> fixed;

	  // Fixed all the points
	  /*for (int i = 0; i < nRowY*nRowX; ++i) // read
	  {
		 fixed.push_back(i);
	  }*/

	  // Fixed half the points
	  /*for (int i = 0; i < nRowY; ++i)
		  for (int j = 0; j < nRowX; ++i)
		  {
			  fixed.push_back(i *nRowX + j);
		  }*/

		  // Fixed the left side points
	  for (int i = 0; i < nRowY; ++i) // read
	  {
		  fixed.push_back(i *nRowX);
	  }

	  // Fixed the right side points
	  for (int i = 1; i <= nRowY; ++i) // read
	  {
		  fixed.push_back(i*nRowX - 1);
	  }

	  /*
	  // Fixed the top side points
	  for (int i = 0; i < nRowX; ++i) // read
	  {
		  fixed.push_back(i);
	  }*/

	  /*
	  // Fixed the bottom side points
	  for (int i = 0; i < nRowX; ++i) // read
	  {
		  fixed.push_back(i);
	  }*/


	  const double youngs = 5E6; // read
	  const double mu = 0.3; // read
	  ConstraintBase::ConstraintType constraintType = ConstraintBase::ConstraintType(ConstraintBase::Type::CST2D,
		  ConstraintBase::SubType::XCPD); // read

	  const double area = lengthY * lengthX;
	  const double element_area = 0.5*area / ((nRowX - 1)*(nRowY - 1));
	  const double nodeMass = element_area * rho*thickness / 3.0;
	  //const double nodeMass = 1000.0;
	  std::cout << "nodeMass: " << nodeMass << "\n";
	  const double dx = lengthX / (double)(nRowX - 1);
	  const double dy = lengthY / (double)(nRowY - 1);

	  StdVectorOfVec3d vertList;
	  vertList.resize(nRowX*nRowY);
	  for (int i = 0; i < nRowY; ++i)
		  for (int j = 0; j < nRowX; j++)
		  {
			  vertList[i*nRowX + j] = Vec3d((double)dx*j, (double)dy*i, 1.0); // vertical plane
			  //vertList[i*nRowX + j] = Vec3d((double)dx*j, 1.0, (double)dy*i*(-1)); // horizontal plane
		  }

	  std::vector<double> massList;
	  massList.resize(nRowX*nRowY);
	  for (unsigned i = 0; i < massList.size(); i++)
	  {
		  massList[i] = 0.0;
		  //std::cout << "massList: " << i << "\n";
	  }

	  using TriangleArray = std::array<size_t, 3>;
	  std::vector<TriangleArray> triangles;

	  for (std::size_t i = 0; i < nRowY - 1; ++i)
	  {
		  for (std::size_t j = 0; j < nRowX - 1; j++)
		  {
			  TriangleArray tri[2];
			  tri[0] = { { i *nRowX + j , i*nRowX + j + 1, (i + 1)*nRowX + j + 1 } };
			  tri[1] = { { (i + 1)*nRowX + j + 1, (i + 1)*nRowX + j, i*nRowX + j } };

			  //tri[0] = { { i *nRowX + j , i*nRowX + j + 1, (i + 1)*nRowX + j, } };
			  //tri[1] = { { (i + 1)*nRowX + j,   i*nRowX + j + 1 ,(i + 1)*nRowX + j + 1} };

			  //tri[0] = { { i *nRowX + j , i*nRowX + j + 1, (i + 1)*nRowX + j + 1 } };
			  //tri[1] = { { (i + 1)*nRowX + j + 1, (i + 1)*nRowX + j, i*nRowX + j } };

			  //tri[0] = { { i *nRowX + j , i*nRowX + j + 1, (i + 1)*nRowX + j } };
			  //tri[1] = { { i *nRowX + j+1,  (i + 1)*nRowX + j + 1, (i + 1)*nRowX + j } };

			  //tri[0] = { { i *nRowX + j , (i + 1)*nRowX + j + 1, i*nRowX + j + 1} };
			  //tri[1] = { { (i + 1)*nRowX + j + 1, i*nRowX + j, (i + 1)*nRowX + j } };

			  //tri[0] = { vertList[i*nRowX + j] , vertList[i*nRowX + j+1] , vertList[(i+1)*nRowX + j+1] };

			  triangles.push_back(tri[0]);
			  triangles.push_back(tri[1]);

			  for (unsigned m = 0; m < 3; m++)
			  {
				  massList[tri[0][m]] += nodeMass;
				  massList[tri[1][m]] += nodeMass;
			  }
		  }
	  }

	  SurfaceMeshPtr surfaceMesh = std::make_shared<SurfaceMesh>(vertList, triangles);

	  ParticleObjectPtr particleObject = std::make_shared<ParticleObject>(surfaceMesh);

	  particleObject->setProperties(youngs, mu);
	  particleObject->addConstraintType(constraintType);

	  ParticlePropertyPtr property = std::make_shared<ParticleProperty>();
	  particleObject->setParticleProperty(property);

	  particleObject->setMassList(massList);

	  particleObject->setFixedPoints(fixed);

	  p_scene->addObject(particleObject);

	  return particleObject;
  }

}