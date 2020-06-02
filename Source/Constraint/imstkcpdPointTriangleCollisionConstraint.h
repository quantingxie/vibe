#ifndef CPDPOINTTRIANGLECOLLISIONCONSTRAINT_H
#define CPDPOINTTRIANGLECOLLISIONCONSTRAINT_H

#include "imstkcpdCollisionConstraint.h"

namespace cpd
{
  class PointTriangleCollisionConstraint : public CollisionConstraintBase<4,1,3> {

  public:

    PointTriangleCollisionConstraint() :CollisionConstraintBase() { m_sizeCol[0] = 1; m_sizeCol[1] = 3; }

    ConstraintType getType() const { return ConstraintType();  }

    void initConstraint(ParticleObjectPtr p_object1, const size_t p_idx1, ParticleObjectPtr p_object2, 
      const size_t p_idx2, const size_t p_idx3, const size_t p_idx4, double p_stiffness = 1.0);

    void updateConstraint() override;
    void computeDerivative() override;
    bool updateDenominator() override;
    void computeForce() override;

    void integrate(bool p_all) override;
    void initShapeFunction() override;

    //void reset() override;
    //bool resolve() override;

  //protected:
  //  size_t m_index1[1];
  //  size_t m_index2[3];
  private:
    Vec3d m_delta[3];
  };

  using PointTriangleCollisionConstraintPtr = std::unique_ptr<PointTriangleCollisionConstraint>;

 /* class PointTriangleCollisionConstraint : public CollisionConstraint {

  public:

    PointTriangleCollisionConstraint() :CollisionConstraint() {}

    Type getType() const { return Type::PointTriangle; }

    void initConstraint(ParticleObjectPtr p_object1, const size_t p_idx1,
      ParticleObjectPtr p_object2, const size_t p_idx2, const size_t p_idx3, const size_t p_idx4);

    void reset() override;
    bool resolve() override;

  protected:
    size_t m_index1[1];
    size_t m_index2[3];
  };

  using PointTriangleCollisionConstraintPtr = std::unique_ptr<PointTriangleCollisionConstraint>;
*/
}

#endif // !CPDPOINTTRIANGLECOLLISIONCONSTRAINT_H

