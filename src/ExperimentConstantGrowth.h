#ifndef EXPERIMENTCONSTANTGROWTH_H_
#define EXPERIMENTCONSTANTGROWTH_H_

#include "Experiment.h"
#include <boost/random/mersenne_twister.hpp>


//! Constant growth experiment.
/*!
 * \class ExperimentConstantGrowth
 *
 * This class implements a constant growth experiment. The cells grow according
 * to the chosen model until \p max_time is reached.
 *
 */
class ExperimentConstantGrowth: public Experiment {
 public:
  //! Create an constant growth Experiment
  ExperimentConstantGrowth();

  //! Create an constant growth Experiment
  /*!
   *
   * \param pars map with parameters as read from the xml
   *
   */
  ExperimentConstantGrowth(std::map<std::string, std::string> pars);

  virtual ~ExperimentConstantGrowth();

  //! Run iterated growth experiment
  virtual void Run();

  //! Create map of simulation settings that can be used to generate the `Experiment` xml context
  virtual std::map<std::string, std::string> GetPars();

  //! simulation time
  double max_time;

  //! frequency for writing output
  int save_freq;

 private:
  double current_time;
};

#endif /* EXPERIMENTCONSTANTGROWTH_H_ */
