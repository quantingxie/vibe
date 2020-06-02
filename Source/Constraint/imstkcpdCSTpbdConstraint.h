#ifndef CPDCSTPBDCONTRAINT_H
#define CPDCSTPBDCONTRAINT_H

#include "imstkcpdCSTEnergyConstraint.h"

namespace cpd
{
  class CSTpbdConstraint :public CSTEnergyConstraint
  {
  public:
    CSTpbdConstraint(bool p_isCompliant = true) : CSTEnergyConstraint(), m_isCompliant(p_isCompliant) {}
    CSTpbdConstraint(float p_stiffness, bool p_isCompliant = true) : CSTEnergyConstraint(), m_isCompliant(p_isCompliant) { m_stiffness = p_stiffness; }

    ConstraintType getType() const override;

    void initConstraint(ParticleObjectPtr p_object, unsigned p_idx1, unsigned p_idx2, unsigned p_idx3, double p_E, double p_mu) override;

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

  using CSTpbdConstraintPtr = std::unique_ptr<CSTpbdConstraint>;
}

#endif // !CPDCSTPBDCONTRAINT_H