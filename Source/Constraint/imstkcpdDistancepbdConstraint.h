#ifndef CPDDISTANCEPBDCONSTRAINT_H
#define CPDDISTANCEPBDCONSTRAINT_H

#include "imstkcpdDistanceConstraint.h"

namespace cpd
{
  class DistancepbdConstraint :public DistanceConstraint
  {
  public:
    DistancepbdConstraint(bool p_isCompliant = true) : DistanceConstraint(), m_isCompliant(p_isCompliant) {}
    DistancepbdConstraint(float p_stiffness, bool p_isCompliant = true) : DistanceConstraint(p_stiffness), m_isCompliant(p_isCompliant) {}

    ConstraintType getType() const override;

    void initConstraint(ParticleObjectPtr p_object, unsigned p_idx1, unsigned p_idx2, const double k = CONSTRAINT_STIFFNESS_DISTANCE) override;

    void resolve() override;
    bool updateDenominator() override;

    void computeForce() override;

    bool isCompliant() { return m_isCompliant; }
    void setPBDstiffness(float p_stiffness);

  public:
    static int m_idGenerator;
  private:
    bool m_isCompliant;
  };

  using DistancepbdConstraintPtr = std::unique_ptr<DistancepbdConstraint>;
}

#endif // !CPDDISTANCECONSTRAINT_H
