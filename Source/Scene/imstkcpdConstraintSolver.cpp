#include "imstkcpdConstraintSolver.h"
#include <chrono>
using Clock = std::chrono::high_resolution_clock;

namespace cpd
{
	void ConstraintSolver::addConstraint(ConstraintSetPtr p_constraintSet, ConstraintBasePtr p_constraint)
	{
		addConstraintSet(p_constraintSet);
		m_constraintSets[p_constraintSet->getPriority()][p_constraintSet->getID()]->addConstraint(std::move(p_constraint));
	}

	void ConstraintSolver::addConstraintSet(ConstraintSetPtr p_constraintSet)
	{
		auto& pair = m_constraintSets.find(p_constraintSet->getPriority());
		if (pair == m_constraintSets.end()) {
			auto& setVector = std::vector<ConstraintSetPtr>();
			setVector.push_back(p_constraintSet);
			m_constraintSets.insert(std::pair<int, std::vector<ConstraintSetPtr>>(p_constraintSet->getPriority(), setVector));
			p_constraintSet->setID(0);
		}
		else
		{
			auto& iter = std::find(pair->second.begin(), pair->second.end(), p_constraintSet);
			if (iter == pair->second.end())
			{
				p_constraintSet->setID(0);
				m_constraintSets[p_constraintSet->getPriority()].push_back(p_constraintSet);
			}
		}
	}

	void ConstraintSolver::addCollisionPair(CollisionPairPtr p_collision)
	{
		m_collisionPairs.push_back(p_collision);
	}

	//void ConstraintSolver::addAffectedParticles()
	//{ // not used for now

	//	for (auto& pair : m_constraintSets)
	//	{
	//		for (auto& set : pair.second)
	//		{
	//			for (auto& constraint : set->getConstraints())
	//			{
	//				// TODO
	//				/*for (auto& p : constraint->getParticles())
	//				{
	//				auto& iter = std::find(getAffectedParticles().begin(), getAffectedParticles().end(), p);
	//				if (iter == getAffectedParticles().end())
	//				m_affectedParticles.push_back(p);
	//				}*/
	//			}
	//		}
	//	}

	//	for (auto& pair : m_collisionPairs)
	//	{
	//		// not implemented yet
	//	}
	//}

	//void ConstraintSolver::clearAffectedParticles()
	//{
	//	m_affectedParticles.clear();
	//}

	void ConstraintSolver::resetLambda()
	{
		for (auto& pair : m_constraintSets)
		{
			for (auto& set : pair.second)
			{
				set->resetLambda();
			}
		}

		for (auto& col : m_collisionPairs)
		{
			col->resetLambda(); // Not doing anything now
		}
	}

	void ConstraintSolver::solveCollision()
	{
		for (auto& pair : m_collisionPairs)
		{
			pair->detectAndResolve();
		}
	}

	double ConstraintSolver::solveConstraints(double p_error)
	{

		double error = 0.0;

		for (auto& pair : m_constraintSets)
		{
			for (auto& set : pair.second)
			{

				auto& constraints = set->getConstraints();
				int n = constraints.size();

#pragma omp parallel for/* num_threads(32)*//*schedule(guided, 40)*/
				for (int idx = 0; idx < n; idx++)
				{
					//constraints[idx]->evaluate();
					constraints[idx]->resolve();
				}

				for (int idx = 0; idx < n; idx++)
				{
					error += 0;// constraints[idx]->getResidual();
				}
			}
		}
		double ratio = error / p_error;
		//std::cout << "error0 = " << p_error << "error = " << error << "r = " << ratio << std::endl;
		return ratio;
	}

	void ConstraintSolver::preStabilize()
	{
	}

	void ConstraintSolver::endTimeStep()
	{
		// not implemented
	}

	unsigned ConstraintSolver::solve()
	{

		unsigned iter = 1;


		//      auto tStart = Clock::now();

		for (auto& pair : m_constraintSets)
			for (auto& set : pair.second)
			{
				set->getObject()->updateTempDisplacement();

				auto& constraints = set->getConstraints();
				int n = constraints.size();
				int thread = 32;

#pragma omp parallel for /*schedule(static, n/thread)*/ /*num_threads(thread)*/
				for (int i = 0; i < n; i++)
				{
					auto& c = constraints[i];
					//c->clearDeltaDisplacement();
					c->clearLamda();
					c->computeDerivative();
					c->updateUnconstrainedDisplacements(); // only EnergyConstraint for now
				}
				set->getObject()->clearConstraintForces();
			}
		//      auto tEnd = Clock::now();


		double error0 = 1.0;
		for (auto& pair : m_constraintSets)
			for (auto& set : pair.second)
			{
				auto& constraints = set->getConstraints();
				int n = set->getOrders().size();
				auto& od = set->getOrders();
				for (int idx = 0; idx < n; idx++)
				{
					error0 += 0;// constraints[od[idx]]->getResidual();
				}
			}

		for (unsigned i = 0; i < iter/*m_solverIteration*/; i++) {
			/*
					for (auto& pair : m_constraintSets)
					for (auto& set : pair.second)
					{
					set->getObject()->updateTempDisplacement();

					auto& constraints = set->getConstraints();
					int n = constraints.size();
					int thread = 32;
					<<<<<<< HEAD
					//#pragma omp parallel for num_threads(thread)
					//            for (int i = 0; i < n; i++)
					//            {
					//constraints[i]->clearDeltaDisplacement();
					//constraints[i]->updateUnconstrainedDisplacements();
					//            }
					//set->getObject()->clearConstraintForces();
					}
					*/
			double r = solveConstraints(error0);
			if (r < 1e-9)
			{
				iter = i + 1;
				break;
			}


		}


		//      auto tEnd2 = Clock::now();

		//      double time1 = std::chrono::duration_cast<std::chrono::nanoseconds>(tEnd - tStart).count();
		//      double time2 = std::chrono::duration_cast<std::chrono::nanoseconds>(tEnd2 - tEnd).count();

		bool collision = true;

		//// pbd dont need to update temp pos
		for (auto& pair : m_constraintSets)
			for (auto& set : pair.second)
			{
				for (auto& pair : m_constraintSets)
					for (auto& set : pair.second)
					{

						set->getObject()->updateTempPositions();
						if (collision)
							set->getObject()->clearConstraintForces();
					}


				if (collision)
				{
					solveCollision();

					for (auto& pair : m_constraintSets)
						for (auto& set : pair.second)
						{
							set->getObject()->updateTempPositions();
						}
				}
			}

		return iter;
	}

}