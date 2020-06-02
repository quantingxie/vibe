#include "imstkcpdConstraintTemplate.h"
#include "imstkcpdConstraintTemplate.cpp"
#include "imstkcpdEnergyConstraint.cpp"
#include "imstkcpdCollisionConstraint.cpp"

namespace cpd
{
  template class ConstraintTemplate<2, 1, 1>;
  template class ConstraintTemplate<2, 2, 1>;
  template class ConstraintTemplate<2, 1, 3>;
  template class ConstraintTemplate<3, 1, 3>;
  template class ConstraintTemplate<3, 3, 2>;
  template class ConstraintTemplate<4, 4, 2>;
  template class ConstraintTemplate<4, 6, 3>;
  template class ConstraintTemplate<4, 8, 2>;

  template class EnergyConstraint<2, 1, 1>;
  template class EnergyConstraint<2, 2, 1>;
  template class EnergyConstraint<2, 1, 3>;
  template class EnergyConstraint<3, 1, 3>;
  template class EnergyConstraint<3, 3, 2>;
  template class EnergyConstraint<4, 4, 2>;
  template class EnergyConstraint<4, 6, 3>;
  template class EnergyConstraint<4, 8, 2>;

  template class CollisionConstraintBase<4, 1, 3>;
}