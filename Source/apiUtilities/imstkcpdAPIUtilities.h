#ifndef CPDAPIUTILITIES_H
#define CPDAPIUTILITIES_H

#include "imstkcpdScene.h"
#include "imstkcpdParticleObject.h"
#include "imstkcpdGeometry.h"

namespace cpd
{
  void createForces(const ScenePtr p_scene, const std::vector<std::array<double, 3>>& forces,
    const std::vector<bool>& isDistributed);

  void createTimeIntegrator(const ScenePtr p_scene, const double timestep,
    const BaseTimeIntegrator::Type type = BaseTimeIntegrator::Type::SIE);

  ParticleObjectPtr createPlane(double p_width, const ScenePtr p_scene, const std::string& p_resourceDir);

  ParticleObjectPtr createPlane1(double p_width, const ScenePtr p_scene, const std::string& p_resourceDir);

  ParticleObjectPtr createCylinder(const ScenePtr p_scene, const std::string& p_resourceDir);
  
  ParticleObjectPtr create2DBeam(const ScenePtr p_scene, const std::string& p_resourceDir);

  ParticleObjectPtr create3DBeam(const ScenePtr p_scene, const std::string& p_resourceDir);

  ParticleObjectPtr createChain(const ScenePtr p_scene, const std::string& p_resourceDir);

  ParticleObjectPtr createHarmonicOscillator(const ScenePtr p_scene, const std::string& p_resourceDir);

  ParticleObjectPtr createObjectFromMeshIO(const ScenePtr p_scene, const std::string& p_resourceDir);

  ParticleObjectPtr createObjectFromTetraMesh(const ScenePtr p_scene, const TetrahedronMeshPtr p_tetraMesh);

  ParticleObjectPtr createObjectFromSurfMesh(const ScenePtr p_scene, const std::string& p_resourceDir);

  void writeTetra(TetrahedronMeshPtr p_mesh);
  void writeSurf(SurfaceMeshPtr p_mesh);

}
#endif // !CPDAPIUTILITIES_H