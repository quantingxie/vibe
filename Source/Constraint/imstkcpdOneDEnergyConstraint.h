#ifndef CPDONEDENERGYCONTRAINT_H
#define CPDCSTENERGYCONTRAINT_H

#include "imstkcpdEnergyConstraint.h"

namespace cpd
{
  class OneDEnergyConstraint :public EnergyConstraint<2, 1, 1>
  {
  public:
    OneDEnergyConstraint() : EnergyConstraint() {}

    ConstraintType getType() const override { return ConstraintType(); }

    void initConstraint(ParticleObjectPtr p_object, unsigned p_idx1, unsigned p_idx2, double p_E, double p_mu);
    void initVariables() override;
    void initDerivatives() override;
    void computeDerivative() override;

    void computeForce() override;
    void updateConstraint() override;

    void reset() override;

    void integrate(bool p_all) override;
    void initShapeFunction() override;

  public:
    static int m_idGenerator;
  };

  using OneDEnergyConstraintPtr = std::unique_ptr<OneDEnergyConstraint>;
}

#endif // !CPDONEDENERGYCONTRAINT_H