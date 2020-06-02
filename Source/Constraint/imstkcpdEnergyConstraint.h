#ifndef CPDENERGYCONSTRAINT_H
#define CPDENERGYCONSTRAINT_H

#include <memory>
#include <array>

#include "imstkcpdConstraintTemplate.h"
#include "imstkcpdParticleObject.h"

namespace cpd
{
  template <int SIZE, int SIZEC, int DIM>
  class EnergyConstraint : public ConstraintTemplate<SIZE, SIZEC, DIM>
  {
  public:
    EnergyConstraint() : ConstraintTemplate() {}

    void initConstraint(ParticleObjectPtr p_object, std::array<unsigned, SIZE> p_idx, double p_E, double p_mu) override;
    bool updateDenominator() override { return true; }

    //void computeDerivative() override {}

    void getStrain(std::array<double, SIZEC>& p_straint) const override;

    VecXd getBI(size_t i) override;

    void updateUnconstrainedDisplacements() override;
    void updateConstraintNew() override;
    void updateStiffnessMatrix() override;

    virtual void updateNN();
  protected:

    double m_YoungsModulus = 1e5;
    double m_PossionRatio = 0.4;

    Eigen::Matrix<double, SIZEC, SIZE*DIM> m_B;
    Eigen::Matrix<double, SIZE*DIM, SIZE*DIM> m_invM;
    Eigen::Matrix<double, SIZEC, SIZEC> m_CMC;

    Eigen::Matrix<double, SIZEC, SIZEC*SIZE> m_NintPre;
    Eigen::Matrix<double, SIZEC, SIZEC*SIZE> m_NintNew;
    Eigen::Matrix<double, SIZEC, SIZEC> m_elasticityM;

    std::array<double, SIZE> m_N;
    std::array<std::array<double, SIZE>, SIZE> m_NN;
    std::array<std::array<double, DIM>, SIZE> m_dNdx;

  };

  template <int SIZE, int SIZEC, int DIM>
  using EnergyConstraintPtr = std::unique_ptr<EnergyConstraint<SIZE, SIZEC, DIM>>;
}
#endif // ! CPDENERGYCONSTRAINT_H