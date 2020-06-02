#ifndef CPDPARTICLEOBJECT_H
#define CPDPARTICLEOBJECT_H

#include <vector>
#include <memory>

//#include "imstkcpdParticle.h"
#include "imstkcpdUtilities.h"
#include "imstkcpdParticleProperty.h"
#include "imstkcpdConstraintBase.h"

#include "imstkcpdMesh.h"
#include "imstkcpdGeometry.h"

#include "imstkcpdUtilities.h"

namespace cpd
{
  class ConstraintSet;
  using ConstraintSetPtr = std::shared_ptr<ConstraintSet>;

  class ParticleObject
  {
  public:
    enum class Mode
    {
      POS,
      DISP,
      Both
    };

    enum class ParticlesInitTpye
    {
      PARTICLES,
      SURFACEMESH,
      VOLUMETRICMESH
    };

  public:

    ParticleObject();
    //ParticleObject(GeometryPtr p_geometry);
    ParticleObject(StdVectorOfVec3d& p_positions);
    ParticleObject(VolumeMeshPtr p_volumeMesh);
    ParticleObject(SurfaceMeshPtr p_surfaceMesh);

    void initialize();
	void initializeDerivativeVariables();

    void setInitPositions(StdVectorOfVec3d p_particles);

    void setParticleCollided(unsigned p_index) { m_collided[p_index] = true; }
    unsigned getParticleCount() const { return m_positions.size(); }

    void setParticleProperty(ParticlePropertyPtr p_property);
    ParticlePropertyPtr getParticleProperty() const { return m_particleProperty; }
    void setParticleInitType(ParticlesInitTpye p_particleInitTpye) { m_particleInitTpye = p_particleInitTpye; }
    ParticlesInitTpye getParticleInitType() const { return m_particleInitTpye; }

    void setMassList(std::vector<double>& p_massList) { m_massList = p_massList; }
    std::vector<double>& getMassList() { return m_massList; }

    void setProperties(double p_property1, double p_property2 = 0.4);
    double getProperties(unsigned p_id) const { return m_properties[p_id]; }

    float getParticleRadius() const { return m_particleProperty->getProperty(PropertyType::RADIUS); }
    float getParticleProximity() const { return m_particleProperty->getProperty(PropertyType::PROXIMITY); }
    //void enableCollision() { if (!m_CollisionEnabled) m_CollisionEnabled = true; }
    //void disableCollision() { if (m_CollisionEnabled) m_CollisionEnabled = false; }
    //void enableSelfCollision() { if (!m_SelfCollisionEnabled) m_SelfCollisionEnabled = true; }
    //void disableSelfCollision() { if (m_SelfCollisionEnabled) m_SelfCollisionEnabled = false; }

    void addExternalForce(const Vec3d& p_force) { m_externalForce += p_force; }
    void addDistributedForce(const Vec3d& p_force) { m_distributedForce += p_force; }
    void updateParticleAcceleration();
    void setFixedPoints(std::vector<unsigned> p_index) { m_fixedParticleIndex = p_index; }

	void addExternalForce(const Vec3d& p_force, unsigned p_idx, double p_mass = 1.0);
	void addExternalForce(const Vec3d& p_force, std::vector<unsigned> p_idx, std::vector<double> p_mass);

    void updateTempPositions();
    void updatePositions(bool p_updateVelocity = true);
    StdVectorOfVec3d& getPositions() { return m_positions; }
    StdVectorOfVec3d& getVelocities() { return  m_velocities; }
    StdVectorOfVec3d& getAccelerations() { return  m_accelerations; }

    void updateVelocity();
    void updateBoundingBox();
    BoundingBox getBoundingBox() const { return m_boundingBox; }

    double getTimeStepSize() { return m_dt; }

    void setSurfaceMesh(SurfaceMeshPtr p_surfaceMesh) { m_surfaceMesh = p_surfaceMesh; }
    void setSurfaceMesh(StdVectorOfVec3d p_vertexPositions, std::vector<std::array<size_t, 3>> p_trianglesVertices);
    void setVolumeMesh(VolumeMeshPtr p_VolumeMesh) { m_volumeMesh = p_VolumeMesh; }
    SurfaceMeshPtr getSurfaceMesh();
    VolumeMeshPtr getVolumeMesh() { return m_volumeMesh; }
    void updateMesh();
    void reset();
    bool checkForReset() { return m_toBeReset; }
    void setToReset() { m_toBeReset = true; }

    Vec3d getDisplacement(unsigned p_id);
    Vec3d getTempDisplacement(unsigned p_id);
    VecXd getTempDisplacement();

    void addConstraintType(ConstraintBase::ConstraintType p_type);
    const std::vector<ConstraintBase::ConstraintType>& getConstraintTypes() const { return m_constriantTypes; }
    bool isPBD() const { return m_isPBD; }

    void addConstraintSet(ConstraintSetPtr p_set) { m_constriantSets.push_back(p_set); }
    const std::vector<ConstraintSetPtr>& getConstraintSets() const { return m_constriantSets; }

    void writePositions(Mode p_mode = Mode::DISP);

    void preSimulation(double p_dt);
    void postSimulation();
    void simulate();

    void clearConstraintForces();
    void addToConstraintForce(unsigned p_idx, Vec3d p_force);
    void addToConstraintForce(unsigned p_idx, unsigned p_cid, Vec3d p_force);

    StdVectorOfVec3d& getTemporaryPositions() { return m_temporaryPositions; }

    double getInvMass(unsigned p_idx) { return m_invMass[p_idx]; }
    Vec3d getTemporaryPosition(unsigned p_idx) { return m_temporaryPositions[p_idx]; }
    Vec3d getPosition(unsigned p_idx) { return m_positions[p_idx]; }
    Vec3d getInitPosition(unsigned p_idx) { return m_initPositions[p_idx]; }

    void setTemporaryPosition(unsigned p_idx, Vec3d p_tempPos) { m_temporaryPositions[p_idx] = p_tempPos; }
    void addToTemporaryPosition(unsigned p_idx, Vec3d p_tempPos) { m_temporaryPositions[p_idx] += p_tempPos; }
    bool isForceUpdated(unsigned p_idx) { return m_forceUpdated[p_idx]; }
    Vec3d getConstraintForce(unsigned p_idx, bool p_updateForce = true) { return m_constraintForces[p_idx]; }
    void clearConstraintForce(unsigned p_idx) { m_constraintForces[p_idx].setZero(); m_constraintForceVectors[p_idx].clear(); }
    void updateTempDisplacement();

    double testConvergence();

    Vec6d getLambda(unsigned p_idx) { return m_lambdas[p_idx]; }
    void setLambda(unsigned p_idx, Vec6d p_lambda) { m_lambdas[p_idx] = p_lambda; }

    void setExculdeCollision(std::vector<unsigned> p_index) { m_excludeIDs = p_index; }

  private:
    StdVectorOfVec3d m_positions;
    StdVectorOfVec3d m_initPositions;
    bool m_initPositionGiven = false;
	unsigned m_numNodes = 0;
	bool m_initialized = false;

    StdVectorOfVec3d m_temporaryPositions;
    StdVectorOfVec3d m_velocities;
    StdVectorOfVec3d m_accelerations;
    StdVectorOfVec3d m_constraintForces;
    std::vector<bool> m_collided;
    std::vector<bool> m_forceUpdated;
    std::vector<double> m_invMass;
    StdVectorOfVec3d m_displacements;

    std::vector<std::map<int, Vec3d>> m_constraintForceVectors;
    std::vector<Vec6d> m_lambdas;

    ParticlePropertyPtr m_particleProperty;

    ParticlesInitTpye m_particleInitTpye = ParticlesInitTpye::PARTICLES;

    SurfaceMeshPtr m_surfaceMesh = nullptr;
    VolumeMeshPtr m_volumeMesh = nullptr;
    BoundingBox m_boundingBox;

    std::vector<unsigned> m_fixedParticleIndex;
    std::vector<double> m_massList;

    std::vector<ConstraintBase::ConstraintType> m_constriantTypes;
    std::vector<ConstraintSetPtr> m_constriantSets;

    double m_dt = TIMESTEP;

    Vec3d m_postion = Vec3d(0.0, 0.0, 0.0);
    Vec3d m_velocity = Vec3d(0.0, 0.0, 0.0);

    Vec3d m_externalForce = Vec3d(0.0, 0.0, 0.0);
    Vec3d m_distributedForce = Vec3d(0.0, 0.0, 0.0);

    bool m_isPBD = false;

    std::array<double, 2> m_properties;
    //bool m_CollisionEnabled;
    //bool m_SelfCollisionEnabled;

    bool m_toBeReset = false;
    std::vector<bool> m_excludeCollision;
    std::vector<unsigned> m_excludeIDs;
  };

  using ParticleObjectPtr = std::shared_ptr<ParticleObject>;

}

#endif //  CPDPARTICLEOBJECT_H
