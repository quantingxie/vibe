#ifndef CPDCOLLISIONCONSTRAINT_H
#define CPDCOLLISIONCONSTRAINT_H

#include <vector>

#include "imstkcpdParticleObject.h"
#include "imstkcpdConstraintTemplate.h"

namespace cpd
{
  template <int SIZE, int SIZEC, int DIM>
  class CollisionConstraintBase : public ConstraintTemplate<SIZE, SIZEC, DIM>
  {

  public:
  //  enum class Type
  //  {
  //    EdgeEdge,
  //    PointTriangle
  //  };

    CollisionConstraintBase() : ConstraintTemplate() {}

    ConstraintType getType() const override { return ConstraintType(); }

    void initDerivatives() override {}

    void updateDeltaDisplacement() override;

    void computeForce() override {}
    void updateConstraint() override {}

    void reset() override {}
    bool updateDenominator() override { return true; }

    void computeDerivative() override {}

    void getStrain(std::array<double, SIZEC>& p_straint) const override {}

    void initVariables() override {}

  protected:
    double m_temp = 1.0;
    ParticleObjectPtr m_objectOther;
    unsigned m_sizeCol[2];
  };

  template <int SIZE, int SIZEC, int DIM>
  using CollisionConstraintBasePtr = std::unique_ptr<CollisionConstraintBase<SIZE, SIZEC, DIM>>;


  class CollisionConstraint
  {

  public:
    enum class Type
    {
      EdgeEdge,
      PointTriangle
    };

    CollisionConstraint() = default;
    virtual void reset() = 0;
    virtual bool resolve() = 0;

  protected:
    ParticleObjectPtr m_object1;
    ParticleObjectPtr m_object2;
  };

  using CollisionConstraintPtr = std::unique_ptr<CollisionConstraint>;

}

#endif // !CPDCOLLISIONCONSTRAINT_H

