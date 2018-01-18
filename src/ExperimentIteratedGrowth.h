#ifndef EXPERIMENTITERATEDGROWTH_H_
#define EXPERIMENTITERATEDGROWTH_H_

#include "Experiment.h"

//! Constant growth experiment.
/*!
 * \class ExperimentIteratedGrowth
 *
 * This class implements an iterated growth experiment. Cells grow until \p max_pop_size is
 * reached. Then, \p number_of_passing_cells are passed to the next generation. This proces is
 * repeat \p number_of_passages times. To prevent infinitely running simulations, growth may be
 * stopped by setting \p max_pass_time. If the growth fase did not result in more cells than
 * \p number_of_passing_cells passage is skipped, or, if stop_sim_if_n_pass_not_reached is set to
 * `true`, the simulation is aborted.
 */
class ExperimentIteratedGrowth: public Experiment {
 public:

  //! Create an iterated growth Experiment
  ExperimentIteratedGrowth();

  //! Create an iterated growth Experiment
  /*!
   *
   * \param pars map with parameters as read from the xml
   *
   */
  ExperimentIteratedGrowth(std::map<std::string, std::string> pars);

  virtual ~ExperimentIteratedGrowth();

  //! Run iterated growth experiment
  virtual void Run();

  //! Print experiment settings
  virtual void PrintInfo();

  //! Create map of simulation settings that can be used to generate the `Experiment` xml context
  virtual std::map<std::string, std::string> GetPars();

  //! Pass subset of cells to the next generation
  /*!
   * Select a given number of cells that will be used in the next generation.
   * For this a list of list of unique cell ids is generated form a uniform
   * distribution and then the corresponding clones are added to the passage
   * population.
   *
   * \param gen Random generator instance used to randomly select the passing cells
   *
   */
  void Passage(base_generator_type &gen, unsigned int nof_pass);
  void PassageUniform(base_generator_type &gen, unsigned int nof_pass);
  void PassageAlt(base_generator_type &gen, unsigned int nof_pass);
  void PassageBiased(base_generator_type &gen, unsigned int nof_pass);

  //! number of passages performed during a simulation
  unsigned int number_of_passages;
  //! number of cells that is passed to the next generation
  unsigned int number_of_passing_cells;
  //! allow to vary #number_of_passing_cells
  bool number_of_passing_cells_vary;
  //! standard deviation #number_of_passing_cells (used for a normal distribution)
  float number_of_passing_cells_sd;
  //! range for  #number_of_passing_cells (used for a uniform distribution)
  std::vector <unsigned int> number_of_passing_cells_range;

  //! size of the initial population
  unsigned int initial_pop_size;
  //! population size at which growth is stopped
  unsigned int max_pop_size;
  //! allow variation #max_pop_size
  bool max_pop_size_vary;
  //! standard deviation #max_pop_size (used for a normal distribution)
  double max_pop_size_sd;
  //! range for  #max_pop_size (used for a uniform distribution)
  std::vector <unsigned int> max_pop_size_range;
  //! time (in days) at which growth is stopped
  double max_pass_time;
  //! stop simulating if the population size after growth < number of passing cells
  bool stop_sim_if_n_pass_not_reached;
  //! passage is biased by clone size
  bool biased_passage;
  bool test_passage;
  double bias_frac;
};

#endif /* EXPERIMENTITERATEDGROWTH_H_ */
