#include <iostream>
#include <fstream>      // std::ofstream
#include "Experiment.h"
#include <chrono>
#include <algorithm>    // std::min_element, std::max_element
#include <boost/iostreams/filtering_streambuf.hpp>
#include <boost/iostreams/copy.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <boost/filesystem.hpp>
#include <limits>       // std::numeric_limits


Experiment::Experiment() {
  seed = (unsigned int) std::chrono::system_clock::now().time_since_epoch().count();
  init_seed = 0;
  current_step = 0;
  gzip = false;
  use_simdir = false;
}

Experiment::Experiment(std::map<std::string, std::string> pars) :
    Experiment() { }

Experiment::~Experiment() { }

void Experiment::SaveState() { SaveState(current_step); }

void Experiment::SaveRates() {
  std::ofstream of;
  of.open(data_path + "/" + simulation_name + "_reaction_rates.txt", std::ios::trunc);
  out_files.insert(data_path + "/" + simulation_name + "_reaction_rates.txt");
  for (auto r : solver->reactions)
    of << r.str() << "\t" << r.rate << "\n";
  of.close();
}

void Experiment::SaveState(double time) {
  SaveState(time, "clones");
}


void Experiment::SaveState(double time, std::string postfix) {
  std::ofstream of;
  of.open(data_path + "/" + simulation_name + "_" + postfix + ".txt", std::ios::app);
  out_files.insert(data_path + "/" + simulation_name + "_" + postfix + ".txt");
  std::string s = "\n" + std::to_string(time);
  for (std::uint64_t cnt : clones) { s += "\t" + std::to_string(cnt); }
  of << s;
  of.close();
}

void Experiment::SavePopulationStats(double time){
  // Create list of frequencies for the population (lumping all species into one)
  std::vector <uint64_t> freq(nclones,0);
  if (model->nspecies == 1)
    freq = clones;
  else {
    for (unsigned int s = 0; s < nclones; s++) {
      for (unsigned int i = s; i < clones.size(); i += nclones) { freq[s] += clones[i]; }
    }
  }

  // compute the fraction of clones that is left, compared to the inital state
  unsigned int nleft = 0;
  for (auto cnt : freq){ if (cnt > 0){ nleft++; } }
  if (current_step == 0){ nclones0 = nleft;}
  double fleft = nleft/((double)nclones0);

  // compute fraction of clones that makes up the first 50 percent of the population
  std::vector <uint64_t> freq_nz;
  for (auto cnt : freq){ if (cnt > 0){ freq_nz.push_back(cnt); }}
  std::sort(freq_nz.begin(), freq_nz.end(), std::greater<int>());
  std::vector <uint64_t> sfreq(nleft,0);
  sfreq[0] = freq_nz[0];
  for (unsigned int i = 1; i < nleft; i++){ sfreq[i] = sfreq[i-1]+freq_nz[i];}
  double f50 = 100*sfreq[(int)(.5*nleft)]/((double)population_size);

  // compute gini coefficient
  std::vector <double> frac;
  for (auto cnt : freq_nz){frac.push_back(cnt/((double)population_size));}
  double sdiff = 0;
  for (auto f : frac){
    for (auto ff : frac){ sdiff += fabs(f-ff); }
  }
  double gini = sdiff/(2.*nleft);
  // compute population size
  int popsize = std::accumulate(clones.begin() + 0, clones.begin() + nclones, 0);
  // write results to file
  std::ofstream of;
  if (current_step == 0) {
    out_files.insert(data_path + "/" + simulation_name + "_popstats.txt");
    of.open(data_path + "/" + simulation_name + "_popstats.txt", std::ios::trunc);
    of << "#time\tpopsize\tfleft\tf50\tgini";
  }
  else
    of.open(data_path + "/" + simulation_name + "_popstats.txt", std::ios::app);
  of << "\n" << time << "\t" << popsize << "\t" << fleft << "\t" << f50 << "\t" << gini;
  of.close();
}

void Experiment::SavePopulationSize() { SavePopulationSize(current_step); }

void Experiment::SavePopulationSize(double time) {
  std::ofstream of;
  if (current_step == 0) {
    out_files.insert(data_path + "/" + simulation_name + "_popsize.txt");
    of.open(data_path + "/" + simulation_name + "_popsize.txt", std::ios::trunc);
    of << "#time\ttotal";
    for (unsigned int i = 0; i < model->nspecies; i++) { of << "\t" << model->species_names[i]; }
  }
  else
    of.open(data_path + "/" + simulation_name + "_popsize.txt", std::ios::app);
  nclones = (unsigned int) clones.size() / model->nspecies;
  of << "\n" << time << "\t" << population_size;
  for (unsigned int i = 0; i < model->nspecies; i++)
    of << "\t" << std::accumulate(clones.begin() + i * nclones, clones.begin() + (i + 1) * nclones, 0);
  of.close();
}

void Experiment::GzipOutput() {
  for (auto f : out_files){
    std::cout << "gzip " << f << std::endl;
    std::ifstream inStream(f, std::ios_base::in);
    std::ofstream outStream(f+".gz", std::ios_base::out);
    boost::iostreams::filtering_streambuf< boost::iostreams::input> in;
    in.push( boost::iostreams::gzip_compressor());
    in.push( inStream );
    boost::iostreams::copy(in, outStream);
    remove(f.c_str());
   }
}

void Experiment::SetupSimulationFolder(){
  std::string data_path_new = data_path + "/" + simulation_name + "/";
  boost::filesystem::path p (data_path_new);
  // clean up old simulation files
  boost::filesystem::create_directory(p);
  data_path = data_path_new;
}

void Experiment::Run() {
  std::cout << "Run undefined experiment\n";
}

std::vector<base_generator_type> Experiment::GetRandomGenerators() {
  std::vector<base_generator_type> generators;
  std::cout << "Create " << ncores << " random generators\n";
  if (ncores == 1) {
    static thread_local base_generator_type *g = nullptr;
    if (!g) g = new base_generator_type(seed);
    generators.push_back(*g);
  }
  else {
    std::minstd_rand0 master_gen(seed);
    std::uniform_int_distribution<int> distr(0, std::numeric_limits<int>::max());
    for (int i = 0; i < ncores; i++) {
      static thread_local base_generator_type *g = nullptr;
      if (!g) g = new base_generator_type(distr(master_gen));
      generators.push_back(*g);
    }
  }
  return generators;
}

void Experiment::Initialize(std::vector<double> init_fractions){
  if (init_from_file)
    InitializeFromFile(init_clones_file,ncells_init,init_fractions);
  else
    InitializeUniform(nclones,ncells_init,init_fractions);
}

void Experiment::InitializeUniform(unsigned int nclones, int ncells, std::vector<double> init_fractions) {
  std::cout << "Initialize clones from uniformly" << std::endl;
  clones.resize(nclones * model->nspecies);
  gen.seed(init_seed);
  std::vector<int> idx;
  std::uniform_real_distribution<double> d(0.0,1.0);
  std::uint64_t nleft,n_per_bc;
  unsigned int min_n_per_bc = (unsigned int)floor(ncells/((double)nclones));
  double p_min = ncells/((double)nclones)-min_n_per_bc;
  double p;
  for (unsigned int i = 0; i < model->nspecies; i++) { idx.push_back(i); }
  for (unsigned int i = 0; i < nclones; i++) {
    // Compute number of cells per clone
    n_per_bc = min_n_per_bc;
    p = d(gen);
    if (p < p_min)
      n_per_bc += 1;
    nleft = n_per_bc;
    // distribute barcodes over species
    std::random_shuffle(idx.begin(), idx.end());
    unsigned int j = 0;
    while ((nleft > 0) && (j < model->nspecies)) {
      std::uint64_t n = (std::uint64_t) ceil(init_fractions[idx[j]] * n_per_bc);
      if (n <= nleft) {
        clones[i + idx[j] * nclones] = n;
        nleft -= n;
      }
      else {
        clones[i + idx[j] * nclones] = nleft;
        nleft = 0;
      }
      j++;
    }
  }
  population_size = (std::uint64_t) std::accumulate(clones.begin(), clones.end(), 0);
  std::cout << "Initialized " << population_size << " cells divided over " << nclones << " clones and ";
  std::cout << model->nspecies << " species\n";
}


void Experiment::InitializeFromFile(std::string init_clones_file, int ncells, std::vector<double> init_fractions) {
  std::cout << "Initialize clones from " << init_clones_file << std::endl;
//  InitializeFromFileMultiSpecies(init_clones_file, ncells, init_fractions);
  if (model->nspecies > 1){ InitializeFromFileMultiSpecies(init_clones_file, ncells, init_fractions); }
  else { InitializeFromFile1Species(init_clones_file, ncells, init_fractions);}
  population_size = (std::uint64_t) std::accumulate(clones.begin(), clones.end(), 0);
  std::cout << "Initialized " << population_size << " cells divided over " << nclones << " clones and ";
  std::cout << model->nspecies << " species\n";
}

void Experiment::InitializeFromFileMultiSpecies(std::string init_clones_file, int ncells, std::vector<double> init_fractions) {
  std::cout << "Initialize clones from " << init_clones_file << std::endl;
  for (unsigned int i = 0; i < model->nspecies; i++) { std::cout << init_fractions[i] << std::endl; }
  std::vector<double> init_clones = ReadInitialClones(init_clones_file);
  nclones = (unsigned int) init_clones.size();
  clones.resize(nclones * model->nspecies);
  std::cout << "Initialize clones with seed = " << init_seed << std::endl;
  std::srand(init_seed);
  std::vector<int> idx;
  for (unsigned int i = 0; i < model->nspecies; i++) { idx.push_back(i); }
  for (unsigned int i = 0; i < nclones; i++) {
    std::uint64_t nleft = (std::uint64_t) round(ncells * init_clones[i]);
    std::random_shuffle(idx.begin(), idx.end());
    unsigned int j = 0;
    // distribute cells evenly over species
    while ((nleft > 0) && (j < model->nspecies)) {
      std::uint64_t n = (std::uint64_t) round(init_fractions[idx[j]] * ncells * init_clones[i]);
      if (n <= nleft) {
        clones[i + idx[j] * nclones] = n;
        nleft -= n;
      }
      else {
        clones[i + idx[j] * nclones] = nleft;
        nleft = 0;
      }
      j++;
    }
    // if there are any 'leftovers' distribute them in the order of shuffled index vector
    j = 0;
    while (nleft > 0){
      clones[i + idx[j]*nclones] += 1;
      nleft -= 1;
      j = (j < model->nspecies) ? j+1 : 0;
    }
  }
  population_size = (std::uint64_t) std::accumulate(clones.begin(), clones.end(), 0);
  std::cout << "Initialized " << population_size << " cells divided over " << nclones << " clones and ";
  std::cout << model->nspecies << " species\n";
}

void Experiment::InitializeFromFile1Species(std::string init_clones_file, int ncells, std::vector<double> init_fractions) {
  std::cout << "Initialize clones from " << init_clones_file << std::endl;
  std::vector<double> init_clones = ReadInitialClones(init_clones_file);
  nclones = (unsigned int) init_clones.size();
  clones.resize(nclones);
  std::cout << "Initialize clones with seed = " << init_seed << std::endl;
  for (unsigned int i = 0; i < nclones; i++) {
    clones[i] = (std::uint64_t) round(ncells * init_clones[i]);
  }
  population_size = (std::uint64_t) std::accumulate(clones.begin(), clones.end(), 0);
  std::cout << "Initialized " << population_size << " cells divided over " << nclones << " clones and ";
  std::cout << model->nspecies << " species\n";
}

bool Experiment::is_number(const std::string& s) {
  return !s.empty() && std::find_if(s.begin(),
                                    s.end(), [](char c) { return !std::isdigit(c); }) == s.end();
}


std::vector<double> Experiment::ReadInitialClones(std::string fn) {
  std::ifstream init_file(fn);
  std::string Barcode;
  std::string CloneSize;
  std::vector<int> reads;
  std::vector<double> init_clones;
  long double nof_reads = 0;
  while (init_file >> Barcode >> CloneSize) {
    if (is_number(CloneSize)) {
      reads.push_back(atoi(CloneSize.c_str()));
      nof_reads += atoi(CloneSize.c_str());
    }
  }
  for (int cnt : reads) {
    init_clones.push_back(cnt / nof_reads);
  }
  return init_clones;
}

void Experiment::PrintInfo() {
  solver->PrintPars();
  model->PrintPars();
}
