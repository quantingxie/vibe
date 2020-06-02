#include "imstkcpdScene.h"
#include <chrono>

using Clock = std::chrono::steady_clock;

const int nit = 1000;

namespace cpd
{
  Scene::Scene()
  {
    m_constraintSolver = std::make_shared<ConstraintSolver>();
    m_externalForceSolver = std::make_shared<ExternalForceSolver>();
    m_time.setZero(nit);
    m_pos.setZero(nit);
    //Eigen::initParallel();
  }

  unsigned Scene::getParticleCount()
  {
    unsigned nbrParticle = 0;
    for (auto& object : m_objects)
    {
      nbrParticle += object->getParticleCount();
    }
    return nbrParticle;
  }

  void Scene::addCollisionPair(ParticleObjectPtr p_object1, ParticleObjectPtr p_object2)
  {
    auto& pair = std::make_shared<CollisionPair>(p_object1, p_object2, COLLISION_SOLVER_ITERATION);
    m_constraintSolver->addCollisionPair(pair);
  }

  void Scene::initializeConstraints()
  {
    unsigned priority = 0;

    for (auto& object : m_objects)
    {
      ConstraintSetPtr constraintSet = std::make_shared<ConstraintSet>(object, priority);
      object->addConstraintSet(constraintSet);
      auto& types = object->getConstraintTypes();
      for (auto& type : types)
      {
        auto subtype = type.getSubType();
        switch (type.getType())
        {
        case ConstraintBase::Type::Distance:
          constraintSet->initializeDistanceConstraints(subtype);
          break;
        case ConstraintBase::Type::Dihedral:
          constraintSet->initializeDihedralConstraints();
          break;
        case ConstraintBase::Type::Area:
          constraintSet->initializeAreaConstraints();
          break;
        case ConstraintBase::Type::CST2D:
          constraintSet->initializeCSTEnergyConstraints(subtype);
          break;
        case ConstraintBase::Type::CST3D:
          constraintSet->initializeCSVEnergyConstraints(subtype);
          break;
        case ConstraintBase::Type::Linear1D:
          constraintSet->initializeOneDEnergyConstraints();
          break;
        default:
          break;
        }
      }

      constraintSet->setTimeStepSize(m_dt);
      std::cout << "dt = " << m_dt << std::endl;
      m_constraintSolver->addConstraintSet(constraintSet);
      std::cout << "Number of constraints = " << constraintSet->getConstraints().size() << std::endl;
    }
  }

  void Scene::preSimulation()
  {
    std::cout << "Initialization ..." << std::endl;

    for (auto& object : m_objects)
    {
      object->preSimulation(m_dt);
    }

    m_externalForceSolver->setTimeStepSize(m_dt);

    for (auto& force : m_externalForceSolver->getExternalForces())
    {
      if (force->isForAll())
      {
        for (auto& obj : m_objects)
        {
          force->addAffectedObject(obj);
          obj->addExternalForce(force->getForce());
        }
      }
      else
      {
        for (auto& obj : force->getAffectedObjects())
        {
          obj->addExternalForce(force->getForce());
        }
      }
    }

    for (auto& force : m_externalForceSolver->getDistributedForces())
    {
      if (force->isForAll())
      {
        for (auto& obj : m_objects)
        {
          force->addAffectedObject(obj);
          obj->addDistributedForce(force->getForce());
        }
      }
      else
      {
        for (auto& obj : force->getAffectedObjects())
        {
          obj->addDistributedForce(force->getForce());
        }
      }
    }

    m_externalForceSolver->updateAffectedObject();

    for (auto& object : m_objects)
    {
      object->updateParticleAcceleration();
    }

    initializeConstraints();

    std::cout << "Initialization Finished." << std::endl;
  }

  bool Scene::simulate()
  {
    // reset objects and constraints if requested
    bool resetConstraint = false;
    for (auto& obj : m_objects)
    {
      if (obj->checkForReset())
      {
        obj->reset();
        resetConstraint = true;
      }
    }

    if (resetConstraint) {
      m_constraintSolver->resetLambda();
     }

    auto tStart = Clock::now();

	unsigned iter;

	// time-stepping for external force
	m_externalForceSolver->exertForces();

	// solve constraints and collisions
	iter = m_constraintSolver->solve();

	auto tSolveEnd = Clock::now();

    double pos;
    double res = 0.0;
    // update position and velocity
    for (auto& obj : m_objects)
    {
      //if (!obj->isPBD())
      //  obj->updateTempPositions();

      obj->updatePositions();

    }
    auto tEnd = Clock::now();

    double time = std::chrono::duration_cast<std::chrono::nanoseconds>(tEnd - tStart).count();
    bool finalstate = tempPrint(time, iter);
    return finalstate;
  }

  void Scene::postSimulation()
  {
    std::cout << "Cleaning up ..." << std::endl;

    bool write = true;
    if (write)
    {
      writeVectorMatlabPlot(m_time, "time.m");
      writeVectorMatlabPlot(m_pos, "res.m");
      for (auto& obj : m_objects)
      {
        obj->writePositions();
        auto& sets = obj->getConstraintSets();

        std::ofstream mfile("test.m");
        mfile << "strain = [";
        mfile.close();
        for (auto& s : sets)
        {
          auto& cs = s->getConstraints();
          for (auto& c : cs)
          {
            c->writeResults("test.m");
          }
        }
        mfile.open("test.m", std::ofstream::out | std::ofstream::app);
        mfile << "];";
        mfile.close();
      }
    }
    std::cout << "Clean-up Finished." << std::endl;
  }

  bool Scene::tempPrint(double time, double res)
  {
    //std::cout << '.';
    const int step = 150;
    const int size = 21;
    static double xx[step][size][3];
    static int cstep = 0;
    static int st = 0;
    static int dt = ceil((1.3 / m_dt) / step);
    static bool chain = false;
    if (chain)
    {
      if (st%dt == 0)
      {
        for (auto& obj : m_objects)
        {
          auto& ps = obj->getPositions();//obj->getPrevPositions();
          for (int k = 0; k < size; k++)
          {
            xx[cstep][k][0] = ps[k][0];
            xx[cstep][k][1] = ps[k][1];
            xx[cstep][k][2] = ps[k][2];
          }
        }
        cstep++;
        std::cout << cstep << '/' << st << std::endl;
      }
      st++;
      if (cstep == step)
      {
        chain = false;
        std::ofstream mfile("chaincpd4.m");

        if (!mfile.is_open())
        {
          std::cout << "Unable to create or open file.";
        }

        std::string tempstr = std::string("chaincpd.m");
        mfile << "u = [\n";
        for (int c = 0; c < step; ++c)
        {
          for (auto i = 0; i < size; ++i)
          {
            mfile << std::setprecision(std::numeric_limits<long double>::digits10 + 1) << xx[c][i][0] << ' ';
          }
          mfile << "\n";
        }
        mfile << "];\n";

        mfile << "v = [\n";
        for (int c = 0; c < step; ++c)
        {
          for (auto i = 0; i < size; ++i)
          {
            mfile << std::setprecision(std::numeric_limits<long double>::digits10 + 1) << xx[c][i][1] << ' ';
          }
          mfile << "\n";
        }
        mfile << "];\n";

        mfile << "w = [\n";
        for (int c = 0; c < step; ++c)
        {
          for (auto i = 0; i < size; ++i)
          {
            mfile << std::setprecision(std::numeric_limits<long double>::digits10 + 1) << xx[c][i][2] << ' ';
          }
          mfile << "\n";
        }
        mfile << "];\n";
        mfile.close();
        std::cout << "Write finished.";
      }
    }

    bool printout = true;
    if (printout)
      for (auto& obj : m_objects)
      {
        static int it = 0;

        static double max = 900.0;
        static bool tryit = true;

        // iteration at which max x displacement is reached
        if (tryit)
        {
          if ((obj->getPosition(obj->getParticleCount() - 1)).y() <= max)
            max = (obj->getPosition(obj->getParticleCount() - 1)).y();
          else
          {
            tryit = false;
            std::cout << "it = " << it << ", max = " << max << std::endl;
          }
        }

        static double totaltime = 0;
        static int totaliter = 0;
        if (it < nit)
        {
          m_time[it] = time;
          m_pos[it] = res;
          std::cout << it << std::endl;
          totaltime += time;
          totaliter += res;
        }
        else if (it == nit)
        {
          double average = totaltime / nit;
          for (unsigned ave = 0; ave < nit; ave++)
          {
            if (m_time[ave] > 1.5*average)
              totaltime += average - m_time[ave];
          }
          std::cout << std::endl << "Average time = " << totaltime / nit << std::endl << std::endl;
          std::cout << std::endl << "Average iter = " << totaliter / nit << std::endl << std::endl;
        }

        it++;

        // static state
        static double dist = 1.0;
        static Vec3d pos = Vec3d(0, 0, 0);
        static Vec3d posp = Vec3d(0, 0, 0);
        static bool onetime = true;

        pos = obj->getDisplacement(obj->getParticleCount() - 1);
        static bool first = true;
        if (first)
        {
          first = false;
          pos = Vec3d(1.0, 0.0, 0.0);
        }
        dist = (pos - posp).norm() / pos.norm();


        int n = int(1 / TIMESTEP);
        //n = 1;
        if ((it%n) == 0)
        {
          std::cout << "it = " << it << ':' << dist * n << std::endl;
          printVector(obj->getPosition(obj->getParticleCount() - 1));
          std::cout << (pos - posp).norm() << ',' << pos.norm() << ',' << dist << ", cr = " << dist / TIMESTEP << std::endl;
        }

        posp = pos;

        static bool converge = false;
        static int cv = 0;
        if (onetime) {
          if (dist / TIMESTEP < 1e-6)
          {
            printVector(obj->getPosition(obj->getParticleCount() - 2));
            printVector(obj->getPosition(obj->getParticleCount() - 1));
            std::cout << std::endl;
            if (cv > 10)
            {
              onetime = false;
            }
            converge = true;
            cv++;
            std::cout << "it = " << it << ':' << dist / TIMESTEP << std::endl;
            printVector(obj->getPosition(obj->getParticleCount() - 1));
          }
          else
          {
            converge = false;
            cv = 0;
          }
        }

        if (it > nit)
          onetime = false;

        return onetime;

      }
  }

}
