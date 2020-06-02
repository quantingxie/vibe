#ifndef CPDCONSTRAINTSOLVER
#define CPDCONSTRAINTSOLVER

#include <map>
#include <random>

#include "imstkcpdConstants.h"
#include "imstkcpdConstraintBase.h"
#include "imstkcpdDistanceConstraint.h"
//#include "imstkcpdDihedralConstraint.h"
#include "imstkcpdCollisionPair.h"
#include "imstkcpdConstraintSet.h"
#include <omp.h>

namespace cpd
{
  class ConstraintSolver {
  public:
    ConstraintSolver(bool p_oneStepNG = true) : m_solverIteration(CONSTRAINT_SOLVER_ITERATION),
      m_preStablizeIteration(PRE_STABLIZE_ITERATION), m_overRelaxation(OVER_RELAXATION) {}

    ConstraintSolver(int p_solverIteration, int p_preStablizeIteration, float p_overRelaxation, bool p_oneStepNG = true)
      : m_solverIteration(p_solverIteration), m_preStablizeIteration(p_preStablizeIteration),
      m_overRelaxation(p_overRelaxation) {}

    void setSolverIteration(int p_solverIteration) { m_solverIteration = p_solverIteration; }
    int getSolverIteration() { return m_solverIteration; }

    void setPreStabilizeIteration(int p_preStablizeIteration) { m_preStablizeIteration = p_preStablizeIteration; }
    int getPreStabilizeIteration() { return m_preStablizeIteration; }
    void setOverRelaxation(float p_overRelaxation) { m_overRelaxation = p_overRelaxation; }
    float getOverRelaxation() { return m_overRelaxation; }

    void addConstraint(ConstraintSetPtr p_constraintSet, ConstraintBasePtr p_constraint);
    void addConstraintSet(ConstraintSetPtr p_constraintSet);
    void addCollisionPair(CollisionPairPtr p_collision);

    //void addAffectedParticles();
    ////const std::vector<ParticlePtr>& getAffectedParticles() { return m_affectedParticles; }
    //void clearAffectedParticles();

    void resetLambda();

    //void setParticleConstraints();

    void solveCollision();
    double solveConstraints(double p_error);
    void preStabilize();
    void endTimeStep();
    unsigned solve();

    const std::map<int, std::vector<ConstraintSetPtr>> getConstraintSets() const { return m_constraintSets; }

    void print() { std::cout << "constraint solver." << std::endl; }
  private:
    std::map<int, std::vector<ConstraintSetPtr>> m_constraintSets;
    std::vector<CollisionPairPtr> m_collisionPairs;
    //std::vector<ParticlePtr> m_affectedParticles;
    unsigned m_solverIteration;
    unsigned m_preStablizeIteration;
    float m_overRelaxation;

    std::random_device m_rd;
  };

  using ConstraintSolverPtr = std::shared_ptr<ConstraintSolver>;
}

#endif // !CPDCONSTRAINTSOLVER
