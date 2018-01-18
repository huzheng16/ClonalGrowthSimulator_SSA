#include <boost/random/uniform_int_distribution.hpp>
#include <boost/random/discrete_distribution.hpp>
#include "ExperimentIteratedGrowth.h"
#include <fstream>      // std::ofstream
#include <set>
#include <chrono>
#include <cmath>
#include <boost/filesystem.hpp>
#include <boost/algorithm/string.hpp>

ExperimentIteratedGrowth::ExperimentIteratedGrowth() : Experiment() {
  max_pop_size = 4148265;
  number_of_passages = 24;
  number_of_passing_cells = 300000;
  ncells_init = 0;
  write_output = false;
  data_path = ".";
  init_clones_file = "/home/mpalm/projects/tumhet/src/growth_simulator/SRR1211210.CC_final.txt";
  init_from_file = true;
  type = "IteratedGrowth";
  gzip = false;
  use_simdir = false;
  max_pass_time = -1;
  bias_frac = 0;
  test_passage = false;
  write_only_stats = false;
}


ExperimentIteratedGrowth::ExperimentIteratedGrowth(std::map<std::string, std::string> pars) :
    ExperimentIteratedGrowth() {
  std::cout << "Set up iterated growth experiment\n";
  for (auto par : pars) {
    if (par.first.compare("Name") == 0) {
      simulation_name = par.second;
      write_output = true;
    } else if (par.first.compare("SaveXML") == 0) {
      save_xml = true;
      xml_out = par.second;
    } else if (par.first.compare("OutPath") == 0)
      data_path = par.second;
    else if (par.first.compare("CriticalPopulationSize") == 0)
      max_pop_size = (unsigned int) std::stoi(par.second);
    else if (par.first.compare("NumberOfCellsToKeep") == 0)
      number_of_passing_cells = (unsigned int) std::stoi(par.second);
    else if (par.first.compare("InitialPopulationSize") == 0)
      ncells_init = (unsigned int) std::stoi(par.second);
    else if (par.first.compare("NumberOfPassages") == 0)
      number_of_passages = (unsigned int) std::stoi(par.second);
    else if (par.first.compare("InitFile") == 0)
      init_clones_file = par.second;
    else if (par.first.compare("InitUniform") == 0) {
      nclones = (unsigned int) std::stoi(par.second);
      init_from_file = false;
    } else if (par.first.compare("Seed") == 0)
      seed = (unsigned int) std::stoi(par.second);
    else if (par.first.compare("InitSeed") == 0)
      init_seed = (unsigned int) std::stoi(par.second);
    else if (par.first.compare("WritePopSize") == 0)
      write_pop_size = true;
    else if (par.first.compare("SaveOnlyStats") == 0)
      write_only_stats = true;
    else if (par.first.compare("gzip") == 0)
      gzip = true;
    else if (par.first.compare("UseSimDir") == 0)
      use_simdir = true;
    else if (par.first.compare("MaxPassTime") == 0)
      max_pass_time = std::stod(par.second);
    else if (par.first.compare("StopSimIfNPassNotReached") == 0)
      stop_sim_if_n_pass_not_reached = true;
    else if (par.first.compare("BiasedPassage") == 0) {
      biased_passage = true;
      bias_frac = std::stod(par.second);
    }
    else if (par.first.compare("TestPassage") == 0)
      test_passage = true;
    else if (par.first.compare("VaryNumberOfCellsToKeepNormal") == 0) {
      number_of_passing_cells_vary = true;
      number_of_passing_cells_range = {0,0};
      number_of_passing_cells_sd = (unsigned int) std::stoi(par.second);
    }
    else if (par.first.compare("VaryNumberOfCellsToKeepUniform") == 0){
      number_of_passing_cells_vary = true;
      number_of_passing_cells_sd = 0;
      std::vector<std::string> strs;
      boost::split(strs, par.second, boost::is_any_of(", "));
      number_of_passing_cells_range = {(unsigned int)std::stoi(strs[0]),(unsigned int)std::stoi(strs[1])};
      std::cout << number_of_passing_cells_range[0] << " - " << number_of_passing_cells_range[1] << std::endl;
    }
    else if (par.first.compare("VaryCriticalPopulationSizeNormal") == 0) {
      max_pop_size_vary = true;
      max_pop_size_sd = (unsigned int) std::stoi(par.second);
      max_pop_size_range = {0,0};
    }
    else if (par.first.compare("VaryCriticalPopulationSizeUniform") == 0){
      max_pop_size_vary = true;
      max_pop_size_sd = 0;
      std::vector<std::string> strs;
      boost::split(strs, par.second, boost::is_any_of(", "));
      max_pop_size_range = {(unsigned int)std::stoi(strs[0]),(unsigned int)std::stoi(strs[1])};
      std::cout << max_pop_size_range[0] << " - " << max_pop_size_range[1] << std::endl;

    }
    else
      std::cout << "Unknown element " << par.first << " in Experiment Iterated Growth" << std::endl;
  }
  if (ncells_init == 0)
    ncells_init = number_of_passing_cells;
  if (init_seed == 0) { init_seed = seed; }
  if (use_simdir){ SetupSimulationFolder(); }
}

std::map<std::string, std::string> ExperimentIteratedGrowth::GetPars() {
  std::map<std::string, std::string> pars = {{"Name", simulation_name},
                                             {"OutPath", data_path},
                                             {"CriticalPopulationSize", std::to_string(max_pop_size)},
                                             {"NumberOfCellsToKeep", std::to_string(number_of_passing_cells)},
                                             {"NumberOfPassages", std::to_string(number_of_passages)},
                                             {"Seed", std::to_string(seed)},
                                             {"InitSeed", std::to_string(init_seed)}};
  if (init_from_file){ pars["InitFile"] = init_clones_file;}
  else{ pars["InitUniform"] = std::to_string(nclones); }
  if (write_only_stats){ pars["SaveOnlyStats"] = "";}
  if (biased_passage) {pars["BiasedPassage"] = ""; }
  if (write_pop_size) { pars["WritePopSize"] = ""; }
  if (test_passage) {pars["TestPassage"] = "";}
  if (gzip) { pars["gzip"] = ""; }
  if (use_simdir) { pars["UseSimDir"] = ""; }
  if (save_xml) { pars["SaveXML"] = xml_out; }
  return pars;
}

ExperimentIteratedGrowth::~ExperimentIteratedGrowth() { }


void ExperimentIteratedGrowth::Passage(base_generator_type &gen, unsigned int nof_pass) {
  if (biased_passage)
    PassageBiased(gen, nof_pass);
  else if (test_passage)
    PassageAlt(gen, nof_pass);
  else
    PassageUniform(gen, nof_pass);
}

void ExperimentIteratedGrowth::PassageAlt(base_generator_type &gen, unsigned int nof_pass){
  std::vector<int> cells;
  for (unsigned int i = 0; i < clones.size(); i++){
    for (unsigned int j = 0; j < clones[i]; j++){ cells.push_back(i); }
  }
  std::shuffle(cells.begin(), cells.end(), gen);
  std::vector<int> keep(cells.begin(), cells.begin()+nof_pass);
  std::vector<std::uint64_t> passed(clones.size(), 0);
  for (auto k : keep){
    passed[k] += 1; }
  std::cout << "passed " << std::accumulate(passed.begin(), passed.end(), 0) << " cells\n";

  clones = passed;
  population_size = (std::uint64_t) std::accumulate(clones.begin(), clones.end(), 0);
}

void ExperimentIteratedGrowth::PassageUniform(base_generator_type &gen, unsigned int nof_pass) {
  std::cout << "Uniform passage\n";
  // Construct list of cells that will be passed to the next generation
  boost::random::uniform_int_distribution<int> d(0, (int) population_size - 1);
  std::set<int> rands;
  while (rands.size() < nof_pass)
    rands.insert(d(gen));
  //LOOP OVER CLONES TO DETERMINE HOW MANY CELLS ARE SELECTED
  std::uint64_t CloneStart = 0;        //number of 1st cell in clone
  unsigned int NBefore = 0;
  unsigned int NAfter = 0;
  std::vector<std::uint64_t> passed(clones.size(), 0);
  for (unsigned int cnt = 0; cnt < clones.size(); cnt++) {
    unsigned int n = 0;    //number of cells in current clone after passaging
    NBefore += clones[cnt];
    std::uint64_t CloneEnd = CloneStart + clones[cnt];
    if (clones[cnt] != 0) {
      for (std::uint64_t i = CloneStart; i < CloneEnd; i++)
        n += rands.count((const int) i);
    }
    if (n > clones[cnt]) { std::cout << "ERROR: increment clone size at passaging\n"; }
    passed[cnt] = n;
    NAfter += clones[cnt];
    CloneStart = CloneEnd;
  }
  clones = passed;
  population_size = (std::uint64_t) std::accumulate(clones.begin(), clones.end(), 0);
}

void ExperimentIteratedGrowth::PassageBiased(base_generator_type &gen, unsigned int nof_pass) {
  std::cout << "Biased passage\n";
  std::vector<double> probs;
  auto num = (double)((clones.size()-1)*population_size);
  auto n = (double)clones.size();
  auto N = (double)population_size;
  std::cout << "Biased passage with bias_frac = " << bias_frac << std::endl;
  for (auto cnt : clones){ probs.push_back((bias_frac*N+cnt*(n*(1-bias_frac)-1))/num);}
  boost::random::discrete_distribution<> dist(probs);
  std::vector<std::uint64_t> passed(clones.size(), 0);
  unsigned int n_passed = 0;
  int clone;
  while (n_passed < nof_pass){
    clone = dist(gen);
    if (passed[clone] < clones[clone]){
      passed[clone] += 1;
      n_passed += 1;
    }
  }
  clones = passed;
  population_size = (std::uint64_t) std::accumulate(clones.begin(), clones.end(), 0);
}


void ExperimentIteratedGrowth::Run() {
  std::vector<base_generator_type> generators = GetRandomGenerators();
  base_generator_type gen = generators[0];
  std::cout << "Run Iterated Growth\n";
  current_step = 0;
  double time = 0;
  // clean old output files
  std::vector<std::string> f_extensions = {"clones_init", "clones_before_passage", "clones_after_passage"};
  for (auto fe : f_extensions) {
    std::string fn = data_path + "/" + simulation_name + "_" + fe + ".txt";
    boost::filesystem::path fpath(fn);
    if (boost::filesystem::exists(fpath)) { boost::filesystem::remove(fpath); }
  }
  // write initial state
  if (write_output) {
    std::cout << "\tSave passage " << current_step << " to " << data_path + "/" + simulation_name + "_clones.txt" <<
        std::endl;
    if (!write_only_stats){SaveState(time, "clones_init");}
    SavePopulationStats(time);
  }
  if (write_pop_size) { SavePopulationSize(); }
  double step_time;
  // Set up distribution for experimentator bias
  unsigned int MaxPop = max_pop_size;
  unsigned int NPass = number_of_passing_cells;
  while ((current_step < number_of_passages) && (population_size > 0)) {
    step_time = 0;
    current_step++;
    if (population_size == 0) {
      std::cout << "All cells are dead :'(\n";
      return;
    }
    if (max_pop_size_vary) {
      if (max_pop_size_sd > 0) {
        std::normal_distribution<> dn(max_pop_size,max_pop_size_sd);
        MaxPop = (unsigned int) dn(gen);
      }
      else {
        std::uniform_int_distribution<int> du(max_pop_size_range[0], max_pop_size_range[1]);
        MaxPop = (unsigned int) round(du(gen));
      }
    }
    if (number_of_passing_cells_vary){
      if (number_of_passing_cells_sd > 0){
        std::normal_distribution<> dn(number_of_passing_cells, number_of_passing_cells_sd);
        NPass = (unsigned int) round(dn(gen));
      }
      else {
        std::uniform_int_distribution<int> du(number_of_passing_cells_range[0],number_of_passing_cells_range[1]);
        NPass = (unsigned int) du(gen);
      }
    }
    std::cout << "Change critical population size from " << max_pop_size << " to " << MaxPop << std::endl;
    std::cout << "Change number of passed cells from " << number_of_passing_cells << " to " << NPass << std::endl;
    std::cout << "Passage " << current_step << " of " << number_of_passages << std::endl;
    std::cout << "\tGrow population from " << population_size << " to " << MaxPop << " cells\n";
    auto start = std::chrono::steady_clock::now();
    while (population_size < MaxPop) {
      if (population_size == 0) {
        std::cout << "All cells are dead :'(\n";
        return;
      }
      if (solver->stochastic) {
        if ((solver->parallel) && (ncores > 1))
          solver->Step(clones, population_size, generators);
        else
          solver->Step(clones, population_size, gen);
      }
      else
        solver->Step(clones, population_size);
      step_time += solver->dt;
      if ((max_pass_time > 0) && (max_pass_time < step_time)){
        std::cout << "Finished growth step after " << step_time << " days with a population of " << population_size << " cells\n";
        break;
      }
    }
    time += step_time;
    auto end = std::chrono::steady_clock::now();
    int msec = (int) std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
    std::cout << "\tsimulating growth took " << msec / 1000. << " seconds\n";
    std::cout << "\ttime = " << time << std::endl;
    std::cout << "\tPass " << NPass << " of " << population_size << " cells\n";
    if (write_output) {
      if (!write_only_stats){SaveState(time, "clones_before_passage");}
      SavePopulationStats(time);
    }
    if (population_size >= NPass)
      Passage(gen,NPass);
    else{
      std::cout << "Skipped passage because the population size is smaller than " << number_of_passing_cells << std::endl;
      if (stop_sim_if_n_pass_not_reached) {
        std::cout << "Stop population because the population size after " << step_time;
        std::cout << " is smaller than " << number_of_passing_cells << std::endl;
        return;
      }
    }
    if (write_output && (!write_only_stats)) {
      std::cout << "\tSave passage " << current_step << " to " <<
          data_path + "/" + simulation_name + "_clones.txt" << std::endl;
      SaveState(time, "clones_after_passage");
    }
    if (write_pop_size) { SavePopulationSize(); }
  }

}

void ExperimentIteratedGrowth::PrintInfo() {
  Experiment::PrintInfo();
}
