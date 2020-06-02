#ifndef CPDCONSTRAINTSET_H
#define CPDCONSTRAINTSET_H

#include <random>
#include <omp.h>

#include "imstkcpdConstraintBase.h"
#include "imstkcpdDistanceConstraint.h"
#include "imstkcpdDistancepbdConstraint.h"
//#include "imstkcpdDihedralConstraint.h"
#include "imstkcpdAreaConstraint.h"
#include "imstkcpdCSTEnergyConstraint.h"
#include "imstkcpdOneDEnergyConstraint.h"
#include "imstkcpdCSVEnergyConstraint.h"
#include "imstkcpdCSTpbdConstraint.h"
#include "imstkcpdCSVpbdConstraint.h"


namespace cpd
{
  class ConstraintSet
  {
  public:
    ConstraintSet(ParticleObjectPtr p_object) :m_priority(0) { m_object = p_object; }
    ConstraintSet(ParticleObjectPtr p_object, unsigned p_priority) :m_priority(p_priority) { m_object = p_object; }


    ConstraintBase::ConstraintType getType() { return m_type; }
    void setPriority(int p_priority) { m_priority = p_priority; }
    int getPriority() { return m_priority; }

    void setID(int p_ID) { m_setID = p_ID; }
    int getID() { return m_setID; }
    void setTimeStepSize(double p_dt);

    std::vector<ConstraintBasePtr>& getConstraints() { return m_constraints; }
    ConstraintBasePtr getConstraint(int p_index) { return std::move(m_constraints[p_index]); }

    void addConstraint(ConstraintBasePtr p_constraint);
    void clearConstraints();

    //void addAffectedParticles(ParticlePtr p_affectedParticle) { m_affectedParticles.push_back(p_affectedParticle); }
    //const std::vector<ParticlePtr>& getAffectedParticles() { return m_affectedParticles; }
    //void clearAffectedParticles();

    const ParticleObjectPtr& getObject() { return m_object; }

    void resetLambda();
    void shuffle(std::mt19937& p_g);
    void reverse();
    std::vector<unsigned>& getOrders() { return m_orders; }

    bool initializeDistanceConstraints(ConstraintBase::SubType p_subType = ConstraintBase::SubType::XCPD);
    void initDistanceConstraint(size_t p_idx1, size_t p_idx2, ConstraintBase::SubType p_subType, double p_stiffness, bool p_isCompliant);
    bool initializeDihedralConstraints();
    bool initializeAreaConstraints(ConstraintBase::SubType p_subType = ConstraintBase::SubType::XCPD);
    bool initializeCSTEnergyConstraints(ConstraintBase::SubType p_subType = ConstraintBase::SubType::XCPD);
    bool initializeOneDEnergyConstraints();
    bool initializeCSVEnergyConstraints(ConstraintBase::SubType p_subType = ConstraintBase::SubType::XCPD);

    void initializeSystemVariables();
    void print() { std::cout << "constraint set." << std::endl; }

  private:
    void checkTypeAndID() const;

  private:
    std::vector<ConstraintBasePtr> m_constraints;
    std::vector<unsigned> m_orders;
    ParticleObjectPtr m_object;
    ConstraintBase::ConstraintType m_type;
    unsigned m_setID;
    unsigned m_priority;
    //std::vector<ParticlePtr> m_affectedParticles;

    VecXd m_f;
  };

  using ConstraintSetPtr = std::shared_ptr<ConstraintSet>;
}

#endif // !CPDCONSTRAINTSET_H
