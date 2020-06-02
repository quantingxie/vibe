#ifndef CPDCONSTRAINTTEMPLATE_H
#define CPDCONSTRAINTTEMPLATE_H

#include <memory>
#include <array>

#include "imstkcpdConstraintBase.h"
#include "imstkcpdParticleObject.h"
#include "imstkcpdQuadrature.h"

namespace cpd
{
  template <int SIZE, int SIZEC, int DIM>
  class ConstraintTemplate : public ConstraintBase
  {
  public:
    ConstraintTemplate() :ConstraintBase() {}
    ConstraintTemplate(float p_stiffness) :ConstraintBase(p_stiffness) {}

    virtual void initConstraint(ParticleObjectPtr p_object, std::array<unsigned, SIZE> p_idx, double p_para1, double p_para2 = 0.0);
    virtual void initVariables() = 0;
    virtual void initDerivatives() = 0;

    void evaluate() override;
    void resolve() override;

    virtual void updateDeltaDisplacement();
    virtual void updateConstraint() = 0;

    void clearLamda() override;
    void clearDeltaDisplacement() override;

    virtual void getStrain(std::array<double, SIZEC>& p_straint) const = 0;
    void writeResults(std::string p_fileName) override;

    VecXd evaluateNew(const VecXd& p_lambda) override;
    size_t getConstraintSize() override { return m_sizeC; }
    void setLambda(VecXd& p_lambda) override;
    VecXd getDeno() override;
    void clearConstraintForces();
    VecXd getBI(size_t i) override;
    VecXi getpParticleIDs() override;

    void updateUnconstrainedDisplacements() {}
    VecXd getElasticityMatrix(size_t i) override;
    virtual void updateConstraintNew() {}

    virtual void updateStiffnessMatrix() {}

    virtual void updateDeno() {}

    virtual void initShapeFunction() = 0;

  protected:
    const int m_size = SIZE;
    const int m_sizeC = SIZEC;
    const int m_dim = DIM;

    std::array<unsigned, SIZE>  m_particleIDs;
    ParticleObjectPtr m_object;

    std::array<double, SIZE> m_invMass;
    std::array<bool, SIZE> m_movable;

    std::array<Vec3d, SIZE> m_unconstrainedDisplacements;
    std::array<std::array<double, DIM>, SIZE> m_derivatives;

    Eigen::Matrix< double, SIZEC, 1> m_lamdaList;
    Eigen::Matrix< double, SIZEC, 1> m_deltaLamdaList; // for xPBD
    Eigen::Matrix< double, SIZEC, 1> m_Constraint;

    Eigen::Matrix<double, SIZEC, SIZEC> m_elasticityMatrix;
    Eigen::Matrix<double, SIZEC, SIZEC> m_deno;

    Eigen::Matrix<double, SIZEC, SIZEC> m_stiffnessMatrix;

    std::array<Vec3d, SIZE> m_deltaDisplacement;

    Quadrature m_quadrature;
    VectorFunction m_shapeFunction;

  };

  template <int SIZE, int SIZEC, int DIM>
  using ConstraintTemplatePtr = std::unique_ptr<ConstraintTemplate<SIZE, SIZEC, DIM>>;
}
#endif // ! CPDCONSTRAINTTEMPLATE_H