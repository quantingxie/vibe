#ifndef CPDMESHIO_H
#define CPDMESHIO_H

#include <memory>
#include <array>

#include "imstkcpdMesh.h"

namespace cpd
{
  class MeshIO
  {
    enum class MeshFileType
    {
      TET,
      VEG,
      UNDEFINED
    };

  public:

    MeshIO() = default;
    ~MeshIO() = default;

    static std::shared_ptr<Mesh> read(const std::string& p_filePath);

    static bool write(const std::shared_ptr<Mesh> p_Mesh, const std::string& p_filePath);

    static const MeshFileType getFileType(const std::string& p_filePath);

  };
  using MeshIOPtr = std::shared_ptr<MeshIO>;

  class TETMeshIO
  {
  public:
    static std::shared_ptr<TetrahedronMesh> read(const std::string& p_filePath);
    static bool write(const std::shared_ptr<TetrahedronMesh> p_Mesh, const std::string& p_filePath);
  };

  class SURFMeshIO
  {
  public:
    static std::shared_ptr<SurfaceMesh> read(const std::string& p_filePath);
    static bool write(const std::shared_ptr<SurfaceMesh> p_Mesh, const std::string& p_filePath);
  };

}

#endif // !CPDMESHIO_H

