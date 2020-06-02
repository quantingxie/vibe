#ifndef CPDCSTENERGYCONTRAINT_H
#define CPDCSTENERGYCONTRAINT_H

#include "imstkcpdEnergyConstraint.h"

namespace cpd
{
  class CSTEnergyConstraint :public EnergyConstraint<3, 3, 2>
  {
  public:
    CSTEnergyConstraint() : EnergyConstraint() {}

    ConstraintType getType() const override { return ConstraintType(ConstraintBase::Type::CST2D, ConstraintBase::SubType::XCPD); }

    virtual void initConstraint(ParticleObjectPtr p_object, unsigned p_idx1, unsigned p_idx2, unsigned p_idx3, double p_E, double p_mu);
    void initVariables() override;
    void initDerivatives() override;
    void computeDerivative() override;

    void computeForce() override;
    void updateConstraint() override;
    double getResidual() override;

    void reset() override;

    void integrate(bool p_all) override;
    void initShapeFunction() override;

  public:
    static int m_idGenerator;
  };

  using CSTEnergyConstraintPtr = std::unique_ptr<CSTEnergyConstraint>;
}

#endif // !CPDCSTENERGYCONTRAINT_H