#include "GrowthModelSimple.h"
#include <chrono>

GrowthModelSimple::GrowthModelSimple() {
  type = "Simple";
  nspecies = 1;
  init_fractions.push_back(1);
  species_names[0] = "SC";
  division_rate = 1;
  division_rate_sd = .1;
  division_rate_2 = 1;
  division_rate_2_sd = .1;
  division_rate_2_frac = 0;
  death_rate = 0;
  division_rate_variation_mode = none;
  seed = (unsigned int) std::chrono::system_clock::now().time_since_epoch().count();
}

GrowthModelSimple::GrowthModelSimple(std::map<std::string, std::string> pars) : GrowthModelSimple() {
  for (auto par : pars) {
    if (par.first.compare("division_rate") == 0)
      division_rate = std::stod(par.second);
    if (par.first.compare("division_rate2") == 0)
      division_rate_2 = std::stod(par.second);
    else if (par.first.compare("DivisionRate") == 0)
      division_rate = std::stod(par.second);
    else if (par.first.compare("DivisionRate2") == 0)
      division_rate_2 = std::stod(par.second);
    else if (par.first.compare("DeathRate") == 0)
      death_rate = std::stod(par.second);
    else if (par.first.compare("DivisionRateSD") == 0) {
      division_rate_sd = std::stod(par.second);
      vary_division_rate = true;
    }
    else if (par.first.compare("DivisionRate2SD") == 0)
      division_rate_2_sd = std::stod(par.second);
    else if (par.first.compare("DivisionRate2Frac") == 0)
      division_rate_2_frac = std::stod(par.second);
    else if (par.first.compare("DivisionRateLogNormal") == 0)
      division_rate_variation_mode = lognormal;
    else if (par.first.compare("DivisionRateVariationMode") == 0){
      if (par.second.compare("none") == 0){ division_rate_variation_mode = none;}
      else if (par.second.compare("normal") == 0){ division_rate_variation_mode = normal;}
      else if (par.second.compare("lognormal") == 0){ division_rate_variation_mode = lognormal;}
      else if (par.second.compare("bimodal") == 0){ division_rate_variation_mode = bimodal;}
    }
    else if (par.first.compare("DivisionRateBimodal") == 0)
      division_rate_variation_mode = bimodal;
    else if (par.first.compare("Seed") == 0)
      seed = (unsigned int)(std::stoi(par.second));
    else
      std::cout << "Unknown parameter " << par.first << " for model type Simple" << std::endl;
  }
  if ((vary_division_rate) && (division_rate_variation_mode == none)){
    division_rate_variation_mode = normal;
  }
  std::cout << "Initialize rng for GrowthModel with seed = " << seed << std::endl;
  gen.seed(seed);
  if ((division_rate_variation_mode == normal) || (division_rate_variation_mode == bimodal))
    norm_distr = std::normal_distribution<>(1,division_rate_sd);
  if (division_rate_variation_mode == bimodal)
    norm_distr2 = std::normal_distribution<>(1,division_rate_2_sd);
  if (division_rate_variation_mode == lognormal) {
    double mu = log(1 / sqrt(1 + division_rate_sd * division_rate_sd));
    double sigma = sqrt(log(1 + division_rate_sd * division_rate_sd));
    lognorm_distr = std::lognormal_distribution<>(mu, sigma);
  }
}

std::map<std::string, std::string> GrowthModelSimple::GetPars() {
  std::map<std::string, std::string> pars = {{"DivisionRate", std::to_string(division_rate)},
                                             {"DeathRate", std::to_string(death_rate)},
                                             {"DivisionRateVariationMode", var_mode_strings[division_rate_variation_mode]}};
  return pars;
}

std::vector<reaction> GrowthModelSimple::GetReactions(){
  double rate = 0;
  if (division_rate_variation_mode == none){
    rate = division_rate;
  }
  else if (division_rate_variation_mode == normal){
    rate = std::max(division_rate * norm_distr(gen), 0.0);
  }
  else if (division_rate_variation_mode == lognormal){
    rate = division_rate * lognorm_distr(gen);
  }
  else if (division_rate_variation_mode == bimodal) {
    std::uniform_real_distribution<double> unif_distr(0.0,1.0);
    if (unif_distr(gen) < division_rate_2_frac){
      rate = std::max(division_rate_2 * norm_distr2(gen), 0.0);
    }
    else{
      rate = std::max(division_rate * norm_distr(gen), 0.0);
    }
  }

  if (death_rate > 0) { reactions = {reaction(rate, {{0, 1}}, {{0, 2}}), reaction(death_rate, {{0, 1}}, {{0, 0}})}; }
  else { reactions = {reaction(rate, {{0, 1}}, {{0, 2}})}; }
  return reactions;
}


GrowthModelSimple::~GrowthModelSimple() { }
