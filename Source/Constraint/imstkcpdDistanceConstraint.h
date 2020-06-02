#ifndef CPDDISTANCECONSTRAINT_H
#define CPDDISTANCECONSTRAINT_H

#include "imstkcpdConstraintTemplate.h"

namespace cpd
{
  class DistanceConstraint :public ConstraintTemplate<2, 1, 3>
  {
  public:
    DistanceConstraint() : ConstraintTemplate() {}
    DistanceConstraint(float p_stiffness) : ConstraintTemplate(p_stiffness) {}

    ConstraintType getType() const override { return ConstraintType(ConstraintBase::Type::Distance, ConstraintBase::SubType::XCPD); }

    virtual void initConstraint(ParticleObjectPtr p_object, unsigned p_idx1, unsigned p_idx2, const double k = CONSTRAINT_STIFFNESS_DISTANCE);
    void initVariables() override;
    void initDerivatives() override;

    bool updateDenominator() override;

    void computeDerivative() override;
    void computeForce() override;
    void updateConstraint() override;

    void getStrain(std::array<double, 1>& p_straint) const override {}

    void reset() override;


    double getResidual() override;

    void integrate(bool p_all) override;

    void initShapeFunction() override;

  public:
    static int m_idGenerator;
  protected:
    Vec3d m_delta;
  };

  using DistanceConstraintPtr = std::unique_ptr<DistanceConstraint>;
}

#endif // !CPDDISTANCECONSTRAINT_H
