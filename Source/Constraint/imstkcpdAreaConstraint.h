#ifndef CPDAREACONTRAINT_H
#define CPDAREACONTRAINT_H

#include "imstkcpdConstraintTemplate.h"
#include "imstkcpdParticleObject.h"

namespace cpd
{
  class AreaConstraint :public ConstraintTemplate<3, 1, 3>
  {
  public:
    AreaConstraint() : ConstraintTemplate() {}
    AreaConstraint(float p_stiffness) : ConstraintTemplate(p_stiffness) {}

    ConstraintType getType() const override { return ConstraintType(); }

    void initConstraint(ParticleObjectPtr p_object, unsigned p_idx1, unsigned p_idx2, unsigned p_idx3, const double k = CONSTRAINT_STIFFNESS);
    void initVariables() override;
    void initDerivatives() override;

    bool updateDenominator() override;

    void computeDerivative() override;
    void computeForce() override;
    void updateConstraint() override;

    void getStrain(std::array<double, 1>& p_straint) const override {}

    void reset() override;

    void integrate(bool p_all) override;
    void initShapeFunction() override;

  public:
    static int m_idGenerator;
  private:
    Vec3d m_delta[2];
  };

  using AreaConstraintPtr = std::unique_ptr<AreaConstraint>;
}

#endif // !CPDAREACONTRAINT_H