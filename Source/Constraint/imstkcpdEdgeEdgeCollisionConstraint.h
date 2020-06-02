#ifndef CPDEDGEEDGECOLLISIONCONSTRAINT_H
#define CPDEDGEEDGECOLLISIONCONSTRAINT_H

#include "imstkcpdCollisionConstraint.h"

namespace cpd
{
  class EdgeEdgeCollisionConstraint : public CollisionConstraintBase<4, 1, 3>
  {

  public:

    EdgeEdgeCollisionConstraint() :CollisionConstraintBase() {}

    ConstraintType getType() const { return ConstraintType(); }

    void initConstraint(ParticleObjectPtr p_object1, const size_t p_idx1, const size_t p_idx2,
      ParticleObjectPtr p_object2, const size_t p_idx3, const size_t p_idx4, double p_stiffness = 1.0);

    void updateConstraint() override;
    void computeDerivative() override;
    bool updateDenominator() override;
    void computeForce() override;

    void integrate(bool p_all) override;
    void initShapeFunction() override;

  private:
    Vec3d m_delta[3];
    Vec3d m_normal;
    double m_st[2];
  };

  using EdgeEdgeCollisionConstraintPtr = std::unique_ptr<EdgeEdgeCollisionConstraint>;
}

#endif // !CPDEDGEEDGECOLLISIONCONSTRAINT_H