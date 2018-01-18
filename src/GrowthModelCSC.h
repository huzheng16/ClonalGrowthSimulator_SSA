#ifndef GROWTHMODELCSC_H_
#define GROWTHMODELCSC_H_

#include "GrowthModel.h"

//! Cell growth model with cancer stem cells (CSC).
/*!
 * \class GrowthModelCSC
 * Cells divide and die at given rates. Division occurs according to
 * the following scheme:
 * \verbatim
   -> SC + SC
   |
SC -> SC + TC1 -> TC2 + TC2 -> ... -> TCN + TCN -> xxx
   |
   -> TC1 + TC1\endverbatim
 *
 */
class GrowthModelCSC: public GrowthModel {
 public:
  GrowthModelCSC();
  //! Initialize model with parameters extracted from xml
  GrowthModelCSC(std::map<std::string, std::string> pars);
  virtual ~GrowthModelCSC();
  virtual void PrintPars();
  virtual std::map<std::string, std::string> GetPars();

  //! Get reactions for a given maximum TC age
  std::vector<reaction> GetReactionsForAge(unsigned int age);
  virtual std::vector<reaction> GetReactions();

 private:
  //! maximum age of transit cells
  unsigned int TC_max_age;
  //! minimum maximum age of transit cells
  unsigned int TC_min_age;
  //! vary age maximum age of transit cells
  bool vary_TC_age;
  //! standard deviation for varying transit cell age
  double TC_age_sd;
  //! division rate of transit cells
  double TC_division_rate;
  //! division rate of stem cells
  double SC_division_rate;
  //! death rate of transit cells
  double TC_death_rate;
  //! death rate of stem cells
  double SC_death_rate;
  //! probability of a stem cell dividing into two stem cells
  double p_sc_2sc;
  //! probability of a stem cell dividing into two transit cells
  double p_sc_2tc;
  //! initial fraction of stem cells
  double init_stem_cell_fraction;
  //! normal distribution used to set varying rates
  std::normal_distribution<> norm_distr;
  //! standard diviation of the division rate
  double division_rate_sd;
  //! vary the division rate
  bool vary_division_rate;
};

#endif /* GROWTHMODELCSC_H_ */
