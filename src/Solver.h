#ifndef SOLVER_H_
#define SOLVER_H_

#include <vector>
#include <map>
#include "GrowthModel.h"
#include <boost/random/mersenne_twister.hpp>
#include <cinttypes>

typedef boost::mt19937 base_generator_type;
//typedef std::mt19937 base_generator_type;

#include "tinyxml2.h"
using namespace tinyxml2;

//! Solver base class.
/*!
 * \class Solver
 *
 * Typical usage:
 * 		- set up model: \link Solver(pars) \endlink
 * 		- run step: \link Setp \endlink
 *
 */
class Solver {
 public:
  Solver() { };
  Solver(std::vector<reaction>) { };

  virtual ~Solver() { };

  virtual void InitializeReactions(GrowthModel * model, int nclones) { };
  //! Perform deterministic solver step
  /*!
   * Solve one step for a list of clones using \p reactions.
   *
   * \param clones list of clones
   * \param population_size total number of clones
   *
   */
  virtual void Step(std::vector<std::uint64_t> &clones, std::uint64_t &population_size) { };
  //! Perform stochastic solver step
  /*!
   * Solve one step for a list of clones using \p reactions.
   *
   * \param clones list of clones
   * \param population_size total number of clones
   * \param gen random generator used by the solver
   *
   */
  virtual void Step(std::vector<std::uint64_t> &clones, std::uint64_t &population_size,
                    std::vector<base_generator_type> gen) { };
  virtual void Step(std::vector<std::uint64_t> &clones, std::uint64_t &population_size,
                    base_generator_type &gen) { };
  //! Print solver parameters
  virtual void PrintPars() { };
  //! Get map with solver parameters
  virtual std::map<std::string, std::string> GetPars() { };

  //! list of reactions for all clones
  std::vector<reaction> reactions;
  bool stochastic;
  bool parallel;
  double dt;

  std::string type;

};

#endif /* SOLVER_H_ */
