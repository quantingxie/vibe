#ifndef CPDCSVENERGYCONTRAINT_H
#define CPDCSVENERGYCONTRAINT_H

#include "imstkcpdEnergyConstraint.h"

namespace cpd
{
  class CSVEnergyConstraint :public EnergyConstraint<4, 6, 3>
  {
  public:
    CSVEnergyConstraint() : EnergyConstraint() {}

    ConstraintType getType() const override { return ConstraintType(ConstraintBase::Type::CST3D, ConstraintBase::SubType::XCPD); }

    virtual void initConstraint(ParticleObjectPtr p_object, unsigned p_idx1, unsigned p_idx2, unsigned p_idx3, unsigned p_idx4, double p_E, double p_mu);
    void initVariables() override;
    void initDerivatives() override;
    void computeDerivative() override;

    void computeForce() override;
    void updateConstraint() override;

    void reset() override;

    void updateDeno() override;

    void integrate(bool p_all) override;
    void initShapeFunction() override;
    void updateNN() override;

  public:
    static int m_idGenerator;
  private:
    std::vector<std::vector<double>> m_abcd;
  };

  using CSVEnergyConstraintPtr = std::unique_ptr<CSVEnergyConstraint>;
}

#endif // !CPDCSVENERGYCONTRAINT_H