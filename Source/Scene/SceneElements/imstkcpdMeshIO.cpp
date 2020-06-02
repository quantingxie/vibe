#include "imstkcpdMeshIO.h"
#include <fstream>

namespace cpd
{
  std::shared_ptr<Mesh> MeshIO::read(const std::string & p_filePath)
  {
    MeshFileType meshType = MeshIO::getFileType(p_filePath);

    switch (meshType)
    {
    case cpd::MeshIO::MeshFileType::TET:
      return TETMeshIO::read(p_filePath);
      break;
    case cpd::MeshIO::MeshFileType::VEG:
      break;
    case cpd::MeshIO::MeshFileType::UNDEFINED:
      break;
    default:
      break;
    }

    return std::shared_ptr<Mesh>();
  }

  bool MeshIO::write(const std::shared_ptr<Mesh> p_Mesh, const std::string & p_filePath)
  {
    return false;
  }

  const MeshIO::MeshFileType MeshIO::getFileType(const std::string & p_filePath)
  {
    MeshFileType meshFileType = MeshFileType::UNDEFINED;

    std::string fileExtension = p_filePath.substr(p_filePath.find_last_of(".") + 1);
    if (fileExtension.empty())
    {
      std::cout << "MeshIO::getFileType error: invalid file name";
      return meshFileType;
    }


    if (fileExtension == "tet" || fileExtension == "TET")
    {
      meshFileType = MeshFileType::TET;
    }
    else if (fileExtension == "veg" || fileExtension == "VEG")
    {
      meshFileType = MeshFileType::VEG;
    }

    return meshFileType;
  }


  std::shared_ptr<TetrahedronMesh> TETMeshIO::read(const std::string & p_filePath)
  {
    StdVectorOfVec3d vertexPositions;
    std::vector<Mesh::TetraArray> tetrahedraVertices;
    std::vector<std::array<double, 4>> verticeWeights;
    std::vector<size_t> verticeTetraIDs;

    std::ifstream infile(p_filePath);
    if (!infile.is_open())
    {
      std::cout << "Unable to create or open file.";
      return nullptr;
    }

    char type;
    const int size = 256;
    char str[size];
    double x, y, z;
    size_t a, b, c, d;
    while (!infile.eof())
    {
      infile >> type;
      switch (type)
      {
      case '#':
        //std::cout << "comments" << std::endl;
        infile.getline(str, size);
        break;
      case 'v':
        //std::cout << "vertex" << std::endl;
        infile >> x;
        infile >> y;
        infile >> z;        
        vertexPositions.push_back(Vec3d(x, y, z));
        break;
      case 't':
        //std::cout << "tetra" << std::endl;
        infile >> a;
        infile >> b;
        infile >> c;
        infile >> d;
        tetrahedraVertices.push_back({ a, b, c, d });
        break;
      case 'l':
        //std::cout << "link" << std::endl;
        infile >> a;
        infile >> x;
        infile >> y;
        infile >> z;
        verticeWeights.push_back({ x, y, z, 1 - x - y - z });
        verticeTetraIDs.push_back(a);
        break;
      default:        
        break;
      }
    }
    infile.close();

    auto tetra = std::make_shared<TetrahedronMesh>(vertexPositions, tetrahedraVertices);
    tetra->setWeightsArray(verticeTetraIDs, verticeWeights);
    return tetra;
  }

  bool TETMeshIO::write(const std::shared_ptr<TetrahedronMesh> p_Mesh, const std::string & p_filePath)
  {
    std::ofstream meshfile("dragon77K.tet");
    if (!meshfile.is_open())
    {
      std::cout << "Unable to create or open file.";
      return false;
    }
    meshfile << "### Asian Dragon mesh" << std::endl;
    meshfile << "### Contains " << p_Mesh->getNumVertices() << " vertices" << std::endl;
    meshfile << "### Contains " << p_Mesh->getNumTetrahedron() << " tetrahedra" << std::endl;

    auto& vertices = p_Mesh->getVertexPositions();
    for (auto& v : vertices)
    {
      meshfile << "v " << v.x() << ' ' << v.y() << ' ' << v.z() << std::endl;
    }

    auto& tetras = p_Mesh->getElementVertices();
    for (auto& tet : tetras)
    {
      meshfile << "t " << tet[0] << ' ' << tet[1] << ' ' << tet[2] << ' ' << tet[3] << std::endl;
    }

    meshfile.close();
    return true;
  }

  std::shared_ptr<SurfaceMesh> SURFMeshIO::read(const std::string & p_filePath)
  {
    StdVectorOfVec3d vertexPositions;
    std::vector<std::array<size_t, 3>> trianglesvertices;

    std::ifstream infile(p_filePath);
    if (!infile.is_open())
    {
      std::cout << "Unable to create or open file.";
      return nullptr;
    }

    char type;
    const int size = 256;
    char str[size];
    double x, y, z;
    size_t a, b, c;
    size_t offset = 0;
    double scale = 1.0/1.0;
    bool write2file = false;
    while (!infile.eof())
    {
      infile >> type;
      switch (type)
      {
      case '#':
        //std::cout << "comments" << std::endl;
        infile.getline(str, size);
        break;
      case 'v':
        //std::cout << "vertex" << std::endl;
        infile >> x;
        infile >> y;
        infile >> z;
        vertexPositions.push_back(scale*Vec3d(x, y, z));
        break;
      case 't':
      case 'f':
        //std::cout << "triangles" << std::endl;
        infile >> a;
        infile >> b;
        infile >> c;
        trianglesvertices.push_back({ a - offset, b - offset, c - offset}); // offset 1 
        break;      
      default:
        break;
      }
    }
    infile.close();

    auto surf = std::make_shared<SurfaceMesh>(vertexPositions, trianglesvertices);
    if(write2file)
        write(surf, "");
    return surf;
  }

  bool SURFMeshIO::write(const std::shared_ptr<SurfaceMesh> p_Mesh, const std::string & p_filePath)
  {
    std::ofstream meshfile("sphere.srf");
    if (!meshfile.is_open())
    {
      std::cout << "Unable to create or open file.";
      return false;
    }
    meshfile << "### Asian Dragon surfacemesh" << std::endl;
    meshfile << "### Contains " << p_Mesh->getNumVertices() << " vertices" << std::endl;
    meshfile << "### Contains " << p_Mesh->getTrianglesVertices().size() << " triangles" << std::endl;

    auto& vertices = p_Mesh->getVertexPositions();
    for (auto& v : vertices)
    {
      meshfile << "v " << v.x() << ' ' << v.y() << ' ' << v.z() << std::endl;
    }

    auto& tetras = p_Mesh->getTrianglesVertices();
    for (auto& tri : tetras)
    {
      meshfile << "t " << tri[0] << ' ' << tri[1] << ' ' << tri[2] << std::endl;
    }

    meshfile.close();
    return true;
  }

}