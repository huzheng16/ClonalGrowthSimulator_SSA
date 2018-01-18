#ifndef EXPERIMENT_H_
#define EXPERIMENT_H_

#include "GrowthModel.h"
#include "Solver.h"
#include <boost/random/mersenne_twister.hpp>
#include <cinttypes>
#include <set>

typedef boost::mt19937 base_generator_type;

//! Base class for a simulation experiment.
/*!
 * \class Experiment
 *
 * Typical usage:
 * 		- set up experiment: \link Experiment(pars) \endlink
 * 		- link model and solver instances
 * 		- initialize clones with a given distribution: \link Initialize \endlink
 * 		- run model: \link Run \endlink
 *
 */
class Experiment {

 public:
  Experiment();

  Experiment(std::map<std::string, std::string> pars);

  virtual ~Experiment();

  virtual void Run();

  //! Get a random number generator for each thread
  /*!
   * To initialize the random number generators, random seeds are
   * produced using \p seed.
   */
  std::vector<base_generator_type> GetRandomGenerators();

  //! Initialize clones
  /*!
   * Initialize clones with the clone distribution extracted from the \p init_clones_file:
   *
   * \param init_clones_file file containing a sample distribution of clones
   * \param ncells number of cells that will be initialized
   * \param init_fractions distribution of species over clones
   *
   */
  void InitializeFromFile(std::string init_clones_file, int ncells, std::vector<double> init_fractions);
  void InitializeFromFileMultiSpecies(std::string init_clones_file, int ncells, std::vector<double> init_fractions);
  void InitializeFromFile1Species(std::string init_clones_file, int ncells, std::vector<double> init_fractions);


  //!
  /*!
   * Initialize clones using a uniform distribition.
   *
   * \param nclones number of clones
   * \param ncells number of cells that will be initialized
   * \param init_fractions distribution of species over clones
   *
   */
  void InitializeUniform(unsigned int nclones, int ncells, std::vector<double> init_fractions);

  //!
  /*!
   * Initialize clones based on the values of class parameters \p init_from_file, \p init_clones_file,
   * \p ncells_init, and nclones.
   *
   * \param init_fractions distribution of species over clones
   *
   */
  void Initialize(std::vector<double> init_fractions);

  //!
  /*!
   * Print some information about the experiment setup
   *
   */
  virtual void PrintInfo();

  //!
  /*!
   * Get map with all experiment parameters.
   *
   */
  virtual std::map<std::string, std::string> GetPars() { };

  //!
  /*!
   * Save number of cells per clone and per species at the current moment. Each row in the output file represents the counts
   * at a given time step. Each column represents the count for given clone and species.
   *
   * \b Example: if there are three clones, (A,B,C) and 2 species (1,2), then the columns are: `A1,B1,C1,A2,B2,C2`.
   */
  void SaveState();

  //!
  /*!
   * Save number of cells per clone and per species at the current moment. Each row in the output file represents the counts
   * at a given time step. Each column represents the count for given clone and species.
   *
   * \b Example: if there are three clones, (A,B,C) and 2 species (1,2), then the columns are: `A1,B1,C1,A2,B2,C2`.
   *
   * \param time time step to be put in the first column
   */
  void SaveState(double time);

  //!
  /*!
   * Save number of cells per clone and per species at the current moment. Each row in the output file represents the counts
   * at a given time step. Each column represents the count for given clone and species.
   *
   * \b Example: if there are three clones, (A,B,C) and 2 species (1,2), then the columns are: `A1,B1,C1,A2,B2,C2`.
   *
   * \param time time step to be put in the first column
   * \param postfix of the output file
   */
  void SaveState(double time, std::string postfix);

  //!
  /*!
   * Write population size to a file with the following columns: `time population_size n_species_1, ..., n_species_n`.
   *
   * \param time time stamp for the first column
   *
   */
  void SavePopulationSize(double time);

  //!
  /*!
   * Write population size to a file with the following columns: `time population_size n_species_1, ..., n_species_n`.
   *
   */
  void SavePopulationSize();


  //!
  /*!
   * Write rate of each transition to file.
   *
   */
  void SaveRates();

  //!
  /*!
   * Write some population stats to file
   *
   */
  void SavePopulationStats(double time);

  //!
  /*!
   * Gzip all output files
   *
   */
  void GzipOutput();

  //!
  /*!
   * Set up folder for simulation results
   *
   */
  void SetupSimulationFolder();

  //! Random generator used to set up the experiment
  boost::random::mt19937 gen;
  //! Seed for simulation algorithm (growth and passage)
  unsigned seed;
  //! Seed for initialization
  unsigned init_seed;
  //! Filename of file with initial clone distribtion
  std::string init_clones_file;
  //! Initialize clones using a stored distribution
  bool init_from_file;
  //! Number of cells to initialize
  unsigned int ncells_init;
  //! Growth model that contains all possible transitions
  GrowthModel *model;
  //! Solver used to run the experiment
  Solver *solver;
  //! Vector containing clone sizes
  std::vector<std::uint64_t> clones;
  //! Number of clones
  unsigned int nclones;
  //! Number of cores used to run the experiment
  int ncores;
  //! Experiment type
  std::string type;
  //! Save simulation settings to xml
  bool save_xml;
  //! Filename of output xml
  std::string xml_out;
  //! Simulation identifier
  std::string simulation_name;
  //! Path to output data
  std::string data_path;
  //! Gzip generated data
  bool gzip;
  //! List of output files
  std::set <std::string> out_files;

 protected:
  //!
  /*!
   * Read initial clone distribution from file
   *
   * \param init_clones_file name of file with input data
   * \return clone fraction (cnt/n_reads) per clone
   *
   */
  std::vector<double> ReadInitialClones(std::string init_clones_file);
  //! Store results in separate folder
  bool use_simdir;
  //! Store simulation results
  bool write_output;
  //! Store population size
  bool write_pop_size;
  //! Sum of all clone frequencies
  std::uint64_t population_size;
  //! Current simulation step
  unsigned int current_step;
  //! Write only population stats
  bool write_only_stats;


 private:
  //! Number of clones at initialization
  int nclones0;
  bool is_number(const std::string& s);
};

#endif /* EXPERIMENT_H_ */
