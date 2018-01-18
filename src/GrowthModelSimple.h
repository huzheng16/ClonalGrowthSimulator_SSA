#ifndef GROWTHMODELSIMPLE_H_
#define GROWTHMODELSIMPLE_H_

#include "GrowthModel.h"

//! Cell growth model with dividing cells.
/*!
 * \class GrowthModelSimple
 * Cells divide according to:
 * \verbatim Cell -> Cell + Cell\endverbatim
 *
 */
class GrowthModelSimple: public GrowthModel {
 public:
  GrowthModelSimple();
  GrowthModelSimple(std::map<std::string, std::string> pars);
  virtual ~GrowthModelSimple();
  virtual std::vector<reaction> GetReactions();
  virtual void PrintPars() { };
  virtual std::map<std::string, std::string> GetPars();

 private:
  //! (mean) rate of cell division
  double division_rate;
  double division_rate_2;
  double division_rate_2_frac;
  //! rate of cell death
  double death_rate;
  //! division rate variation mode
  enum var_mode { none, normal, lognormal, bimodal} division_rate_variation_mode;
  std::vector<std::string> var_mode_strings = { "none", "normal", "lognormal", "bimodal" };

  //! standard deviation of cell division rate
  double division_rate_sd;
  double division_rate_2_sd;

  //! set division rate to vary
  bool vary_division_rate;
  //! set division rate distribution to lognormal (defaults to normal)
  bool vary_division_rate_lognormal;
  //! normal distribution used for division rate variation
  std::normal_distribution<> norm_distr;
  std::normal_distribution<> norm_distr2;

  //! log-normal distribution used for division rate variation
  std::lognormal_distribution<> lognorm_distr;
};

#endif /* GROWTHMODELSIMPLE_H_ */
