#include "GrowthModelCSC.h"
#include <boost/algorithm/string.hpp>
#include <chrono>

GrowthModelCSC::GrowthModelCSC() {
  type = "CSC";
  TC_min_age = 2;
  TC_max_age = 2;
  vary_TC_age = false;
  TC_age_sd = 2;
  seed = (unsigned int) std::chrono::system_clock::now().time_since_epoch().count();
  nspecies = TC_max_age + 1;
  species_names[0] = "SC";
  for (unsigned int i = 1; i < nspecies; i++) { species_names[i] = "TC" + std::to_string(i); }
  TC_division_rate = 1;
  SC_division_rate = 1;
  TC_death_rate = 0;
  SC_death_rate = 0;
  p_sc_2sc = 0.5;
  p_sc_2tc = 0.1;
  init_stem_cell_fraction = -1;
  for (unsigned int i = 0; i < nspecies; i++) { init_fractions.push_back(1. / nspecies); }
}

GrowthModelCSC::GrowthModelCSC(std::map<std::string, std::string> pars) : GrowthModelCSC() {
  std::cout << "Set up CSC model\n";
  for (auto par : pars) {
    if (par.first.compare("MaxAgeTransitCell") == 0)
      TC_max_age = (unsigned int)std::stoi(par.second);
    if (par.first.compare("VaryMaxAgeTransitCellSD") == 0) {
      vary_TC_age = true;
      TC_age_sd = std::stod(par.second);
    }
    else if (par.first.compare("DivisionRateStemCell") == 0)
      SC_division_rate = std::stod(par.second);
    else if (par.first.compare("DivisionRateTransitCell") == 0)
      TC_division_rate = std::stod(par.second);
    else if (par.first.compare("DeathRateStemCell") == 0)
      SC_death_rate = std::stod(par.second);
    else if (par.first.compare("DeathRateTransitCell") == 0)
      TC_death_rate = std::stod(par.second);
    else if (par.first.compare("ProbSymmStemCellDivision") == 0)
      p_sc_2sc = std::stod(par.second);
    else if (par.first.compare("ProbSymmStemCellDifferentiation") == 0)
      p_sc_2tc = std::stod(par.second);
    else if (par.first.compare("DivisionRateSD") == 0) {
      division_rate_sd = std::stod(par.second);
      vary_division_rate = true;
    }
    else if (par.first.compare("InitialFractions") == 0) {
      std::vector<std::string> strs;
      boost::split(strs, par.second, boost::is_any_of(","));
      init_fractions.clear();
      for (unsigned int i = 0; i < strs.size(); i++) { init_fractions.push_back(std::stod(strs[i])); }
    }
    else if (par.first.compare("InitialStemCellFraction") == 0) {
      init_stem_cell_fraction = std::stod(par.second);
    }
    else if (par.first.compare("Seed") == 0)
      seed = (unsigned int)(std::stoi(par.second));
    else
      std::cout << "Unknown parameter " << par.first << " for model type CSC" << std::endl;
  }
  nspecies = TC_max_age + 1;
  species_names[0] = "SC";
  if (init_stem_cell_fraction >= 0) {
    init_fractions.clear();
    init_fractions.push_back(init_stem_cell_fraction);
    for (unsigned int i = 1; i < nspecies; i++)
      init_fractions.push_back((1 - init_stem_cell_fraction) / (nspecies - 1));
    std::cout << init_fractions.size() << " " << nspecies << std::endl;
  }
  for (unsigned int i = 1; i < nspecies; i++) { species_names[i] = "TC" + std::to_string(i); }
  if (init_fractions.size() != nspecies) {
    std::cout << "InitialFractions contains " << init_fractions.size() << " elements, but there are " << nspecies
        << " speceis\n";
    init_fractions.clear();
    for (unsigned int i = 0; i < nspecies; i++) { init_fractions.push_back(1. / nspecies); }
  }
  gen.seed(seed);

}

std::map<std::string, std::string> GrowthModelCSC::GetPars() {
  std::map<std::string, std::string> pars = {{"MaxAgeTransitCell", std::to_string(TC_max_age)},
                                             {"DivisionRateStemCell", std::to_string(SC_division_rate)},
                                             {"DivisionRateTransitCell", std::to_string(TC_division_rate)},
                                             {"DeathRateStemCell", std::to_string(SC_death_rate)},
                                             {"DeathRateTransitCell", std::to_string(TC_death_rate)},
                                             {"ProbSymmStemCellDivision", std::to_string(p_sc_2sc)},
                                             {"ProbSymmStemCellDifferentiation", std::to_string(p_sc_2tc)},
                                             {"InitialStemCellFraction", std::to_string(init_stem_cell_fraction)}};
  pars["InitialFractions"] = std::to_string(init_fractions[0]);
  for (unsigned int i = 1; i < init_fractions.size(); i++)
    pars["InitialFractions"] += "," + std::to_string(init_fractions[i]);
  return pars;
}

std::vector<reaction> GrowthModelCSC::GetReactions() {
  unsigned int age;
  if (vary_TC_age) {
    age = std::max(TC_min_age,(unsigned int)round(.5*(TC_max_age+TC_min_age)+TC_age_sd*norm_distr(gen)));
    age = std::min(age,TC_max_age);
  }
  else {
    age = TC_max_age;
  }
  return GetReactionsForAge(age);
}

std::vector<reaction> GrowthModelCSC::GetReactionsForAge(unsigned int age) {
//  unsigned int nreac = 2 * (nspecies + 1);
  unsigned int nreac = 2 * (age + 2)-1;
  std::vector<reaction> reac(nreac);
  // reactions 0 - 2 = stem cell division
  // reactions 3 - age+1 = TC division
  // reactions age+2 = cell death at t = tmax
  // reactions age+3 = stem cell death
  // reactions age+4 - 2*age+2 = TC cell death
  // SC -> SC + SC
  double fac = 1.0;
  if (vary_division_rate)
    fac = std::max(1+division_rate_sd*norm_distr(gen), 0.0);
  reac[0] = reaction(fac*SC_division_rate * p_sc_2sc, {{0, 1}}, {{0, 2}});
  // SC -> SC + TC
  reac[1] = reaction(fac*SC_division_rate * (1 - p_sc_2sc - p_sc_2tc), {{0, 1}}, {{0, 1}, {1, 1}});
  // SC -> TC + TC
  reac[2] = reaction(fac*SC_division_rate * p_sc_2tc, {{0, 1}}, {{1, 2}});
  // TC(i) -> TC(i+1) + TC(i+1)
  for (unsigned int j = 1; j < age; j++)
    reac[(j + 2)] = reaction(fac*TC_division_rate, {{j, 1}}, {{(j + 1), 2}});
  // TC(N) -> xxx
  reac[age + 2] = reaction(fac*TC_division_rate, {{age, 1}}, {});
  // cell death
  reac[age + 3] = reaction(SC_death_rate, {{0, 1}}, {});
  for (unsigned int j = 1; j < age; j++)
    reac[j + age + 3] = reaction(TC_death_rate, {{j, 1}}, {});
  return reac;
}


GrowthModelCSC::~GrowthModelCSC() { }

void GrowthModelCSC::PrintPars() {
  GrowthModel::PrintPars();
}
