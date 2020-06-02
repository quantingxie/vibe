#ifndef CPDSCENE_H
#define CPDSCENE_H

#include "imstkcpdParticleObject.h"
#include "imstkcpdConstraintSolver.h"
#include "imstkcpdExternalForceSolver.h"
#include "imstkcpdCollisionPair.h"
#include "imstkcpdConstraintBase.h"
#include "imstkcpdUtilities.h"

namespace cpd
{
  class Scene {

  public:
    Scene();

    void addObject(ParticleObjectPtr p_particleObject) { m_objects.push_back(p_particleObject); }
    unsigned getObjectCount() { return m_objects.size(); }
    unsigned getParticleCount();
    std::vector<ParticleObjectPtr> getObjects() { return m_objects; }

    void addExternalForce(ExternalForcePtr p_force) { m_externalForceSolver->addExternalForce(p_force); }
    ExternalForceSolverPtr getExternalForceSolver() { return m_externalForceSolver; }

    void addDistributedForce(ExternalForcePtr p_force) { m_externalForceSolver->addDistributedForce(p_force); }

    void addConstraint(ConstraintSetPtr p_constraintSet, ConstraintBasePtr p_constraint) { m_constraintSolver->addConstraint(p_constraintSet, std::move(p_constraint)); }
    void addCollisionPair(ParticleObjectPtr p_object1, ParticleObjectPtr p_object2);

    float getTimeStepSize() { return m_externalForceSolver->getTimeStepSize(); }
    void setTimeStepSize(float p_timeStep) { m_dt = p_timeStep;}

    unsigned getSolverIteration() { return m_constraintSolver->getSolverIteration(); }
    void setSolverIteration(int p_nbrSolverIteration) { m_constraintSolver->setSolverIteration(p_nbrSolverIteration); }

    unsigned getPreStablizeIteration() { return m_constraintSolver->getPreStabilizeIteration(); }
    void setPreStabilizeIteration(int p_nbrPreStabilizeIteration) { m_constraintSolver->setPreStabilizeIteration(p_nbrPreStabilizeIteration); }

    float getOverRelaxation() { return m_constraintSolver->getOverRelaxation(); }
    void setOverRelaxation(float p_overRelaxation) { m_constraintSolver->setOverRelaxation(p_overRelaxation); }

    void initializeConstraints();

    void preSimulation();
    void postSimulation();
    bool tempPrint(double time, double res);
    bool simulate();

  private:
    float m_dt = TIMESTEP;
    std::vector<ParticleObjectPtr> m_objects;
    ConstraintSolverPtr m_constraintSolver;
    ExternalForceSolverPtr m_externalForceSolver;
    Eigen::VectorXd m_time;
    Eigen::VectorXd m_pos;
  };

  using ScenePtr = std::shared_ptr<Scene>;
}
#endif // !CPDSCENE_H
