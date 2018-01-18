#include <iostream>		//std::cout
#include <map>
#include <chrono>
#include <time.h>

#include "tinyxml2.h"
using namespace tinyxml2;

#include "Experiment.h"
#include "ExperimentIteratedGrowth.h"
#include "ExperimentConstantGrowth.h"
#include "SolverGillespieOTL.h"
#include "GrowthModelCSC.h"
#include "GrowthModelSimple.h"

//! Simulation builder
/*!
 * \class Builder
 *
 * The builder parses the xml, creates the model, solver and experiment
 * instances, and performs all actions needed before the experiment can run.
 */
class Builder {
 public:
  GrowthModel *model;
  Solver *solver;
  Experiment *experiment;
  int ncores;

  Builder() {
    experiment = NULL;
    model = NULL;
    solver = NULL;
    ncores = 1;
  }


  Builder(std::string fn) {
    experiment = NULL;
    model = NULL;
    solver = NULL;
    ncores = 1;
    ParseFile(fn);
    if (model) {
      experiment->model = model;
      experiment->solver = solver;
      experiment->ncores = ncores;
      std::cout << "Experiment seed = " << experiment->seed << std::endl;
      // Initialize cell population
      experiment->Initialize(model->init_fractions);
      // Initialize reactions
      solver->InitializeReactions(model,experiment->nclones);
      experiment->SaveRates();
    }
    std::string xml_out = experiment->xml_out.empty() ? experiment->data_path+"/"+experiment->simulation_name+".xml" : experiment->xml_out;
    this->ToXML(xml_out);
    experiment->out_files.insert(xml_out);
  }


  void ToXML(std::string fn) {
    tinyxml2::XMLDocument doc;
    tinyxml2::XMLNode *pRoot = doc.NewElement("Simulation");
    tinyxml2::XMLElement *base_el;

    time_t rawtime;
    char timestamp [80];
    time (&rawtime);
    struct tm * timeinfo = localtime (&rawtime);
    strftime (timestamp,80,"%Y-%m-%d %H:%M:%S",timeinfo);
    (pRoot->ToElement())->SetAttribute("run",timestamp);

    // Write model parameters
    base_el = doc.NewElement("Model");
    base_el->SetAttribute("type", model->type.c_str());
    for (auto parset : model->GetPars()) {
      tinyxml2::XMLElement *el = doc.NewElement(parset.first.c_str());
      if (parset.second.size() > 0) { el->SetText(parset.second.c_str()); }
      base_el->InsertEndChild(el);
    }
    pRoot->InsertEndChild(base_el);
    // Write solver parameters
    base_el = doc.NewElement("Solver");
    base_el->SetAttribute("type", solver->type.c_str());
    for (auto parset : solver->GetPars()) {
      tinyxml2::XMLElement *el = doc.NewElement(parset.first.c_str());
      if (parset.second.size() > 0) { el->SetText(parset.second.c_str()); }
      base_el->InsertEndChild(el);
    }
    pRoot->InsertEndChild(base_el);
    // Write experiment parameters
    base_el = doc.NewElement("Experiment");
    base_el->SetAttribute("type", experiment->type.c_str());
    for (auto parset : experiment->GetPars()) {
      tinyxml2::XMLElement *el = doc.NewElement(parset.first.c_str());
      if (parset.second.size() > 0) { el->SetText(parset.second.c_str()); }
      base_el->InsertEndChild(el);
    }
    pRoot->InsertEndChild(base_el);
    doc.InsertFirstChild(pRoot);
    XMLError eResult = doc.SaveFile(fn.c_str());
  }

  std::map<std::string, std::string> ParsePars(const tinyxml2::XMLNode *node) {
    std::map<std::string, std::string> pars;
    const tinyxml2::XMLElement *el = node->FirstChildElement();
    while (el) {
      if (el->FirstChild())
        pars[el->Name()] = el->GetText();
      else
        pars[el->Name()] = "";
      el = el->NextSiblingElement();
    }
    return pars;
  }

  void BuildModel(const tinyxml2::XMLNode *node) {
    std::string node_name;
    if ((node->ToElement())->Attribute("type", "CSC")) {
      model = new GrowthModelCSC(ParsePars(node));
    }
    else if ((node->ToElement())->Attribute("type", "Simple")) {
      model = new GrowthModelSimple(ParsePars(node));
    }
  }

  void BuildSolver(const tinyxml2::XMLNode *node) {
    std::string node_name;
    if ((node->ToElement())->Attribute("type", "GillespieOTL")) {
      solver = new SolverGillespieOTL(ParsePars(node));
    }
  }

  void BuildExperiment(const tinyxml2::XMLNode *node) {
    if ((node->ToElement())->Attribute("type", "IteratedGrowth")) {
      experiment = new ExperimentIteratedGrowth(ParsePars(node));
    }
    else if ((node->ToElement())->Attribute("type", "ConstantGrowth")) {
      experiment = new ExperimentConstantGrowth(ParsePars(node));
    }
  }

  void ParseFile(std::string fn) {
    std::cout << "Load file " << fn << std::endl;
    tinyxml2::XMLDocument xmlDoc;
    xmlDoc.LoadFile(fn.c_str());
    std::cout << "finished loading xml\n";
    tinyxml2::XMLNode *root = xmlDoc.LastChild();
    (root->ToElement())->QueryIntAttribute("ncores", &ncores);
    const tinyxml2::XMLNode *node = root->FirstChild();
    std::string node_name;
    while (node) {
      node_name = std::string(node->Value());
      std::cout << node_name << std::endl;
      if (node_name.compare("Model") == 0)
        BuildModel(node);
      else if (node_name.compare("Solver") == 0)
        BuildSolver(node);
      else if (node_name.compare("Experiment") == 0)
        BuildExperiment(node);
      node = node->NextSibling();
    }
    std::cout << "Finished parsing xml\n";
  }
};

std::map<std::string, std::string> ParsePars(const tinyxml2::XMLNode *node) {
  std::map<std::string, std::string> pars;
  const tinyxml2::XMLElement *el = node->FirstChildElement();
  while (el) {
    pars[el->Name()] = el->GetText();
    el = el->NextSiblingElement();
  }
  return pars;
}

int main(int argc, char *argv[]) {
  auto start = std::chrono::steady_clock::now();
  bool verbose = true;
  // find parameter file
  if (argc < 2) {
    std::cout << "Settings file missing\n";
    exit(EXIT_FAILURE);
  }
  if (argc == 3) {
    if ((std::string(argv[1])).compare("-q")) { verbose = false; }
    else if ((std::string(argv[1])).compare("-v")) { verbose = true; }

  }
  std::string fn = argv[argc - 1];

  Builder builder = Builder(fn);
  if (verbose)
    builder.experiment->PrintInfo();

  // Run simulation
  builder.experiment->Run();
  if (builder.experiment->gzip)
    builder.experiment->GzipOutput();


  // Get simulation time
  if (verbose) {
    auto end = std::chrono::steady_clock::now();
    int tsec = (int) std::chrono::duration_cast<std::chrono::seconds>(end - start).count();
    int tmin = (int) std::chrono::duration_cast<std::chrono::minutes>(end - start).count();
    int thour = (int) std::chrono::duration_cast<std::chrono::hours>(end - start).count();
    std::cout << "Simulation time: ";
    if (thour > 0)
      std::cout << thour << " hours ";
    if ((tmin > 0) || (thour > 0))
      std::cout << tmin - (thour * 60) << " minutes ";
    std::cout << tsec - (tmin * 60) << " seconds\n";
  }
}
