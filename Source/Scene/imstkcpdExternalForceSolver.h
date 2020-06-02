#ifndef CPDEXTERNALFORCESOLVER_H
#define CPDEXTERNALFORCESOLVER_H

#include <memory>
#include <vector>
#include <unordered_set>

#include "imstkcpdExternalForce.h"
#include "imstkcpdBaseTimeIntegrator.h"
#include "imstkcpdSemiImplicitEulerIntegrator.h"
#include "imstkcpdCentralDifferenceIntegrator.h"
#include "imstkcpdParticleObject.h"

namespace cpd
{
  class ExternalForceSolver {

  public:
    ExternalForceSolver() {}

    void addExternalForce(ExternalForcePtr p_force) { m_forces.push_back(p_force); }
    void addDistributedForce(ExternalForcePtr p_force) { m_distributedForces.push_back(p_force); }
    const std::vector<ExternalForcePtr>& getDistributedForces() const { return m_distributedForces; }
    const std::vector<ExternalForcePtr>& getExternalForces() const { return m_forces; }

    //void clearExternalForce() { m_forces.clear(); }
    //void clearDistributedForce() { m_distributedForces.clear(); }

    void updateAffectedObject();

    void setTimeIntegrator(BaseTimeIntegratorPtr p_timeIntegrator) { m_timeIntegrator = p_timeIntegrator; }
    const BaseTimeIntegratorPtr getTimeIntegrator() const { return m_timeIntegrator; }

    void setTimeStepSize(float p_timeStep) { m_timeIntegrator->setTimeStepSize(p_timeStep); }
    float getTimeStepSize() const { return m_timeIntegrator->getTimeStepSize(); }

    std::unordered_set<ParticleObjectPtr> getAffectedObjects() const { return m_objects; }

    void exertForces();

  private:
    std::vector<ExternalForcePtr> m_forces;
    std::vector<ExternalForcePtr> m_distributedForces;
    BaseTimeIntegratorPtr m_timeIntegrator;
    std::unordered_set<ParticleObjectPtr> m_objects;
  };

  using ExternalForceSolverPtr = std::shared_ptr<ExternalForceSolver>;
}

#endif // !CPDEXTERNALFORCESOLVER_H
