#include <boost/random/uniform_int_distribution.hpp>
#include "ExperimentConstantGrowth.h"
#include <set>


ExperimentConstantGrowth::ExperimentConstantGrowth() : Experiment() {
  max_time = 1;
  save_freq = 1;
  ncells_init = 1000;
  write_output = false;
  data_path = ".";
  init_clones_file = "/home/mpalm/projects/tumhet/src/growth_simulator/SRR1211210.CC_final.txt";
  init_from_file = true;
  type = "IteratedGrowth";
  gzip = false;
  use_simdir = false;
  write_only_stats = false;
}


ExperimentConstantGrowth::ExperimentConstantGrowth(std::map<std::string, std::string> pars) :
    ExperimentConstantGrowth() {
  std::cout << "Set up iterated growth experiment\n";
  for (auto par : pars) {
    if (par.first.compare("Name") == 0) {
      simulation_name = par.second;
      write_output = true;
    }
    else if (par.first.compare("SaveXML") == 0) {
      save_xml = true;
      xml_out = par.second;
    }
    else if (par.first.compare("SimulationTime") == 0)
      max_time = std::stod(par.second);
    else if (par.first.compare("SaveFreq") == 0)
      save_freq = std::stoi(par.second);
    else if (par.first.compare("OutPath") == 0)
      data_path = par.second;
    else if (par.first.compare("SaveOnlyStats") == 0)
      write_only_stats = true;
    else if (par.first.compare("InitialPopulationSize") == 0)
      ncells_init = (unsigned int) std::stoi(par.second);
    else if (par.first.compare("InitFile") == 0)
      init_clones_file = par.second;
    else if (par.first.compare("InitUniform") == 0) {
      nclones = (unsigned int) std::stoi(par.second);
      init_from_file = false;
    }
    else if (par.first.compare("Seed") == 0)
      seed = (unsigned int) std::stoi(par.second);
    else if (par.first.compare("InitSeed") == 0)
      init_seed = (unsigned int) std::stoi(par.second);
    else if (par.first.compare("gzip") == 0)
      gzip = true;
    else if (par.first.compare("UseSimDir") == 0)
      use_simdir = true;
    else if (par.first.compare("WritePopSize") == 0)
      write_pop_size = true;
    else
      std::cout << "Unknown element " << par.first << " in Experiment Constant Growth" << std::endl;
  }
  if (init_seed == 0) { init_seed = seed; }
  if (use_simdir){ SetupSimulationFolder(); }
}

std::map<std::string, std::string> ExperimentConstantGrowth::GetPars() {
  std::map<std::string, std::string> pars = {{"SimulationTime", std::to_string(max_time)},
                                             {"SaveFreq", std::to_string(save_freq)},
                                             {"OutPath", data_path},
                                             {"InitialPopulationSize", std::to_string(ncells_init)},
                                             {"Seed", std::to_string(seed)},
                                             {"InitSeed", std::to_string(init_seed)}};
  if (init_from_file){ pars["InitFile"] = init_clones_file;}
  else{ pars["InitUniform"] = std::to_string(nclones); }
  if (write_only_stats){ pars["SaveOnlyStats"] = "";}
  if (save_xml) { pars["SaveXML"] = xml_out; }
  if (gzip) { pars["gzip"] = ""; }
  if (use_simdir) { pars["UseSimDir"] = ""; }
  return pars;
}


ExperimentConstantGrowth::~ExperimentConstantGrowth() { }


void ExperimentConstantGrowth::Run() {
  std::vector<base_generator_type> generators = GetRandomGenerators();
  base_generator_type gen = generators[0];
  std::cout << "Run Constant Growth\n";
  current_step = 0;
  current_time = 0;
  unsigned int next_save_step = 0;
  double next_print_step = 0;
  if (write_output) {
    if (!write_only_stats){SaveState(current_time);}
    SavePopulationStats(current_time);
    if (write_pop_size) { SavePopulationSize(); }
    next_save_step += save_freq;
  }
  while (current_time <= max_time) {
    current_step++;
    if (solver->stochastic) {
      if (solver->parallel)
        solver->Step(clones, population_size, generators);
      else
        solver->Step(clones, population_size, gen);
    }
    else
      solver->Step(clones, population_size);
    current_time += solver->dt;
    if (current_time >= next_print_step) {
      std::cout << (int) (100 * current_time / max_time) << "%: " << population_size << " cells\n";
      next_print_step += 0.1 * max_time;
    }
    if ((write_output) && (current_step >= next_save_step)) {
      if (!write_only_stats){SaveState(current_time);}
      SavePopulationStats(current_time);
      if (write_pop_size) { SavePopulationSize(); }
      next_save_step += save_freq;
    }
  }
  std::cout << "100%: " << population_size << " cells\n";
}
