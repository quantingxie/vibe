#ifndef CPDCSVPBDCONTRAINT_H
#define CPDCSVPBDCONTRAINT_H

#include "imstkcpdCSVEnergyConstraint.h"

namespace cpd
{
  class CSVpbdConstraint :public CSVEnergyConstraint
  {
  public:
    CSVpbdConstraint(bool p_isCompliant = true) : CSVEnergyConstraint(), m_isCompliant(p_isCompliant) {}
    CSVpbdConstraint(float p_stiffness, bool p_isCompliant = true) : CSVEnergyConstraint(), m_isCompliant(p_isCompliant) { m_stiffness = p_stiffness; }

    ConstraintType getType() const override;

    void initConstraint(ParticleObjectPtr p_object, unsigned p_idx1, unsigned p_idx2, unsigned p_idx3, unsigned p_idx4, double p_E, double p_mu) override;

    void resolve() override;

    void updateConstraint() override;

    void computeForce() override;

    bool isCompliant() { return m_isCompliant; }
    void setPBDstiffness(float p_stiffness);

  public:
    static int m_idGenerator;
  private:
    bool m_isCompliant;
  };

  using CSVpbdConstraintPtr = std::unique_ptr<CSVpbdConstraint>;
}

#endif // !CPDCSVPBDCONTRAINT_H