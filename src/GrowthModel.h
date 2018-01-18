#ifndef GROWTHMODEL_H_
#define GROWTHMODEL_H_

#include <string>
#include <random>
#include <vector>
#include <map>
#include <iostream>
#include <boost/random/mersenne_twister.hpp>
typedef boost::mt19937 base_generator_type;
//! Reaction species
/*!
 * \struct reaction_species
 *
 * Object that holds the index and reaction coefficient of a reactant or product.
 *
 */
struct reaction_species { unsigned int idx; unsigned int coeff; };
//struct reaction{double rate; std::vector<reaction_species> source; std::vector<reaction_species> product;};

//! Reaction object.
/*!
 * Container class for a transition.
 *
 * \class reaction
 *
 *
 */
class reaction {
 public:
  /*! reaction rate */
  double rate;
  /*! reaction sources */
  std::vector<reaction_species> source;
  /*! reaction products */
  std::vector<reaction_species> product;

  reaction() {
    rate = 0;
  };


  //! Create a reaction instance
  /*!
   *
   * \param r: reaction rate
   * \param s: list of \link reaction_species \endlink that are used in the reaction
   * \param p: list of \link reaction_species \endlink that are produced by the reaction
   *
   */
  reaction(double r, std::vector<reaction_species> s, std::vector<reaction_species> p) {
    rate = r;
    source = s;
    product = p;
  };

  //! print formatted reaction using species indices
  std::string str() {
    std::map<unsigned int, std::string> species_names;
    for (auto s : source) { species_names[s.idx] = "species[" + std::to_string(s.idx) + "]"; }
    for (auto s : product) { species_names[s.idx] = "species[" + std::to_string(s.idx) + "]"; }
    return str(species_names);
  };

  //! print formatted reaction using species names
  /*!
   *
   * \param species_names: map linking the species indices to species names
   *
   */
  std::string str(std::map<unsigned int, std::string> species_names) {
    std::string s;
    for (unsigned int i = 0; i < source.size(); i++) {
      if (i > 0) { s += " + "; }
      if (source[i].coeff > 1) { s += std::to_string(source[i].coeff); }
      s += species_names[source[i].idx];
    }
    char srate[50];
    sprintf(srate, "%.3f", rate);
    s += " ---(" + std::string(srate) + ")---> ";
    for (unsigned int i = 0; i < product.size(); i++) {
      if (i > 0) { s += " + "; }
      if (product[i].coeff > 1) { s += std::to_string(product[i].coeff); }
      s += species_names[product[i].idx];
    }
    return s;
  };
};

//! Base class for cell growth model.
/*!
 * \class GrowthModel
 *
 * Typical usage:
 * 		- create model with parameters from xml: \link GrowthModel(pars) \endlink
 * 		- create a list of reactions: \link SetReactions \endlink
 *
 */
class GrowthModel {
 public:

  GrowthModel() { };
  virtual ~GrowthModel() { };
  //! Print model paramters
  virtual void PrintPars() {
    std::cout << "Reactions:\n";
    for (auto r : reactions) { std::cout << "\t" << r.str(species_names) << "\n"; };
  };
  //! Get list of transitions in the model
  virtual std::vector<reaction> GetReactions() { };
  //! Get map of model parameters
  virtual std::map<std::string, std::string> GetPars() { };

  //! Initialize model with parameters extracted from xml
  GrowthModel(std::map<std::string, std::string> pars) { };

  //! number of species
  unsigned int nspecies;
  //! list of model reactions
  std::vector<reaction> reactions;
  //! mapping of species names
  std::map<unsigned int, std::string> species_names;
  //! initial distribution of the species
  std::vector<double> init_fractions;
  //! Type of model
  std::string type;
  //! Random generator instance used to vary the division rates
  base_generator_type gen;
  //! Seed for \link gen \endlink
  unsigned int seed;
};

#endif /* GROWTHMODEL_H_ */
