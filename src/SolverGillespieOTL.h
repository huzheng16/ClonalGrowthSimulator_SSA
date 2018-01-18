#ifndef SOLVERGILLESPIEOTL_H_
#define SOLVERGILLESPIEOTL_H_

#include "Solver.h"
#include <string>

//! Ordinary tau leaping Gillespie solver
/*!
 * \class SolverGillespieOTL
 *
 * see
 * <a href="http://scitation.aip.org/content/aip/journal/jcp/115/4/10.1063/1.1378322">
 * Gillespie, Daniel T. "Approximate accelerated stochastic simulation of chemically reacting systems."
 * The Journal of Chemical Physics 115.4 (2001): 1716-1733.</a>
 * and <a href="https://en.wikipedia.org/wiki/Tau-leaping">wikipedia</a>
 *
 */
class SolverGillespieOTL: public Solver {
 public:
  SolverGillespieOTL();
  SolverGillespieOTL(std::map<std::string, std::string> pars);
  virtual ~SolverGillespieOTL();
  virtual void Step(std::vector<std::uint64_t> &clones, std::uint64_t &population_size) { };
  virtual void Step(std::vector<std::uint64_t> &clones, std::uint64_t &population_size,
                    base_generator_type &gen);
  virtual void Step(std::vector<std::uint64_t> &clones, std::uint64_t &population_size,
                    std::vector<base_generator_type> gen);
  virtual void InitializeReactions(GrowthModel * model, int nclones);

  virtual std::map<std::string, std::string> GetPars();
  //! Tau leaping step size
  double tau;
};

#endif /* SOLVERGILLESPIEOTL_H_ */
