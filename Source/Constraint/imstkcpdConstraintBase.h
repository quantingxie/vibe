#ifndef CPDCONSTRAINTBASE_H
#define CPDCONSTRAINTBASE_H

#include <memory>
#include <omp.h>

#include "imstkcpdMath.h"
#include "imstkcpdConstants.h"
//#include "imstkcpdParticle.h"


namespace cpd
{
  class ConstraintBase
  {
  public:
    enum class Type
    {
      Distance,
      Dihedral,
      Area,
      Volume,
      Energy,
      CST3D,
      CST2D,
      Linear1D,
      CST3DPBD,
      CST2DPBD
    };

    enum class SubType
    {
      XCPD,
      PBD,
      XPBD
    };

    class ConstraintType
    {
    public:
      ConstraintType() {}
      ConstraintType(Type p_type, SubType p_subType = SubType::XCPD) : m_type(p_type), m_subType(p_subType) {}
      Type getType() const { return m_type; }
      SubType getSubType() const { return m_subType; }
    private:
      Type m_type;
      SubType m_subType;
    };

  public:
    ConstraintBase() {}
    ConstraintBase(float p_stiffness) :m_stiffness(p_stiffness) {}

    virtual ConstraintType getType() const = 0;
    unsigned getID() const { return m_ID; }

    void setTimeStepSize(double p_dt) { m_dt = p_dt; }

    virtual void evaluate() = 0;
    virtual bool updateDenominator() = 0;
    virtual void resolve() = 0;

    virtual void computeDerivative() = 0;
    virtual void computeForce() = 0;

    virtual void clearLamda() = 0;
    virtual void clearDeltaDisplacement() = 0;

    virtual void reset() = 0;

    virtual void writeResults(std::string p_fileName) = 0;

    virtual void updateConstraintNew() = 0;
    virtual VecXd evaluateNew(const VecXd& p_lambda) = 0;
    virtual size_t getConstraintSize() = 0;
    virtual void setLambda(VecXd& p_lambda) = 0;
    virtual VecXd getDeno() = 0;
    virtual VecXd getBI(size_t i) = 0;
    virtual VecXi getpParticleIDs() = 0;

    virtual void updateUnconstrainedDisplacements() = 0;
    virtual VecXd getElasticityMatrix(size_t i) = 0;

    virtual double getResidual() { return 0; }

    virtual void integrate(bool p_all) = 0;


  protected:

    unsigned m_ID;
    double m_volume = 1.0;
    double m_stiffness = 1.0;
    double m_dt = TIMESTEP;

    bool m_updateDerivative = true;
  };

  using ConstraintBasePtr = std::unique_ptr<ConstraintBase>;

}
#endif // ! CPDCONSTRAINTBASE_H