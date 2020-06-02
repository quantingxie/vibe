#ifndef CPDCOLLISIONPAIR_H
#define CPDCOLLISIONPAIR_H

#include <memory>
#include <vector>

#include "imstkcpdConstants.h"
#include "imstkcpdParticleObject.h"
#include "imstkcpdCollisionConstraint.h"

namespace cpd
{
  class CollisionPair {

  public:
    CollisionPair() :m_solverIteration(COLLISION_SOLVER_ITERATION) {}
    CollisionPair(ParticleObjectPtr p_object1, ParticleObjectPtr p_object2, int p_solverIteration)
      :m_object1(p_object1), m_object2(p_object2), m_solverIteration(p_solverIteration) {}

    void setSolverIteration(int p_solverIteration) { m_solverIteration = p_solverIteration; }
    int getSolverIteration() { return m_solverIteration; }

    void initCollisionConstraint();

    void resetConstraints();
    void resetLambda();

    void resolve();
    void detectAndResolve();

    bool broadPhaseCollisionDetection();
    void narrowPhaseCollisionDetetion();

    void printNumberOfCollision() { std::cout << m_collisionConstraints.size() << std::endl; }

  private:
    std::vector<ConstraintBasePtr> m_collisionConstraints;
    ParticleObjectPtr m_object1 = nullptr;
    ParticleObjectPtr m_object2 = nullptr;
    unsigned m_solverIteration;
    unsigned m_pairID;
  };

  using CollisionPairPtr = std::shared_ptr<CollisionPair>;
}

#endif // !CPDCOLLISIONPAIR_H

