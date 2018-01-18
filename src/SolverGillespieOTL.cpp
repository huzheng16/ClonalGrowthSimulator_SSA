#include "SolverGillespieOTL.h"
#include <boost/random/poisson_distribution.hpp>
#include <fstream>      // std::ofstream
#include <algorithm>    // std::max

#ifdef _OPENMP
#include <omp.h>
#else
#define omp_get_thread_num() 0
#endif
SolverGillespieOTL::SolverGillespieOTL() {
  type = "GillespieOTL";
  tau = 0.0005;
  stochastic = true;
  parallel = false;
  dt = tau;
}

SolverGillespieOTL::SolverGillespieOTL(std::map<std::string, std::string> pars) : SolverGillespieOTL() {
  std::cout << "Set up Gillespie OTL solver\n";
  for (auto par : pars) {
    std::cout << par.first << " " << par.second << std::endl;
    if (par.first.compare("tau") == 0)
      tau = std::stod(par.second);
    else
      std::cout << "Unknown parameter " << par.first << " for solver type GillespieOTL" << std::endl;
  }
  std::cout << "Finished setting up Gillespie OTL solver\n";
  std::cout << "time step = " << tau << std::endl;
  dt = tau;
}

void SolverGillespieOTL::InitializeReactions(GrowthModel * model, int nclones) {
  for (int i = 0; i < nclones; i++) {
    std::vector < reaction> spec_reac = model->GetReactions();
    for (unsigned int j = 0; j < spec_reac.size(); j++) {
      reaction R;
      R.rate = spec_reac[j].rate;
      for (auto s : spec_reac[j].source) { R.source.push_back({s.idx * nclones + i, s.coeff}); }
      for (auto s : spec_reac[j].product) { R.product.push_back({s.idx * nclones + i, s.coeff}); }
      reactions.push_back(R);
    }
  }
  std::cout << "Created " << reactions.size() << " reactions\n";
};


void SolverGillespieOTL::Step(std::vector<std::uint64_t> &clones,
                              std::uint64_t &population_size,
                              std::vector<base_generator_type> gen) {
  dt = tau;
  std::vector<std::vector<int> > change(gen.size(), std::vector<int>(clones.size(), 0));
  std::vector<std::uint64_t> pop_size_change(gen.size(), 0);
  // fix population_size
#pragma omp parallel num_threads(gen.size())
  {
#pragma omp for
    for (unsigned int i = 0; i < reactions.size(); i++) {
      // compute propensity
      double propensity = reactions[i].rate;
      std::uint64_t k_max = clones[reactions[i].source[0].idx] / reactions[i].source[0].coeff;
      for (auto r : reactions[i].source) {
        propensity *= (clones[r.idx] > r.coeff) ? r.coeff * clones[r.idx] : 0;
        k_max = std::max(k_max, (std::uint64_t) r.coeff);
      }
      // stop if there is no source
      if (propensity == 0)
        continue;

      // compute how often the reaction will occur
      boost::random::poisson_distribution<std::uint64_t> d(propensity * tau);

      int id = omp_get_thread_num();
      std::uint64_t k = std::min(k_max, d(gen[id]));
      // stop if the reaction does not occur
      if (k == 0)
        continue;
      // remove sources
      if (k > k_max)
        std::cout << "This can not be >:(\n \t k = " << k << " k_max = " << k_max << "\n";
      for (auto item : reactions[i].source) {
        change[id][item.idx] -= k * item.coeff;
        pop_size_change[id] -= k * item.coeff;
      }
      // add products
      for (auto item : reactions[i].product) {
        change[id][item.idx] += k * item.coeff;
        pop_size_change[id] += k * item.coeff;
      }
    }
  }
  for (unsigned int t = 0; t < gen.size(); t++) {
    population_size += pop_size_change[t];
    for (unsigned int i = 0; i < clones.size(); i++)
      clones[i] = (clones[i] + change[t][i] > 0) ? clones[i] + change[t][i] : 0;
  }
//	for (unsigned int i = 0; i < clones.size(); i++){clones[i] = (clones[i]+change[i] > 0) ? clones[i]+change[i] : 0;}
//		clones[i] = std::max(clones[i]+change[i],0);}
}

void SolverGillespieOTL::Step(std::vector<std::uint64_t> &clones,
                              std::uint64_t &population_size,
                              base_generator_type &gen) {
  double propensity;
  std::uint64_t k_max, k;
  std::vector<int> change(clones.size(), 0);
  for (unsigned int i = 0; i < reactions.size(); i++) {
    // compute propensity
    propensity = reactions[i].rate;
    k_max = clones[reactions[i].source[0].idx] / reactions[i].source[0].coeff;
    for (auto r : reactions[i].source) {
      propensity *= (clones[r.idx] > r.coeff) ? r.coeff * clones[r.idx] : 0;
      k_max = std::max(k_max, (std::uint64_t) r.coeff);
    }
    // stop if there is no source
    if (propensity == 0)
      continue;

    // compute how often the reaction will occur
    boost::random::poisson_distribution<std::uint64_t> d(propensity * tau);
    k = std::min(k_max, d(gen));
    //		k=d(gen);
    // stop if the reaction does not occur
    if (k == 0)
      continue;
    // remove sources
    if (k > k_max)
      std::cout << "This can not be >:(\n \t k = " << k << " k_max = " << k_max << "\n";
    for (auto item : reactions[i].source) {
      change[item.idx] -= k * item.coeff;
      population_size -= k * item.coeff;
    }
    // add products
    for (auto item : reactions[i].product) {
      change[item.idx] += k * item.coeff;
      population_size += k * item.coeff;
    }
  }
  for (unsigned int i = 0; i < clones.size(); i++) {
    clones[i] = (clones[i] + change[i] > 0) ? clones[i] + change[i] : 0;
  }
}

SolverGillespieOTL::~SolverGillespieOTL() { }

std::map<std::string, std::string> SolverGillespieOTL::GetPars() {
  std::map<std::string, std::string> pars = {{"tau", std::to_string(tau)}};
  return pars;
}
