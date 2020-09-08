#include <iostream>
#include <memory>
#include <fstream>
#include <sstream>
#include <string>

#include "UHH2/TstarTstar/include/NeuralNetworkModules.h"

#include "UHH2/core/include/AnalysisModule.h"
#include "UHH2/core/include/Selection.h"
#include "UHH2/core/include/Utils.h"
#include "UHH2/core/include/Event.h"

#include "UHH2/common/include/EventHists.h"
#include "UHH2/common/include/ElectronHists.h"
#include "UHH2/common/include/MuonHists.h"
#include "UHH2/common/include/MuonIds.h"
#include "UHH2/common/include/ElectronIds.h"
#include "UHH2/common/include/ObjectIdUtils.h"
#include "UHH2/common/include/JetIds.h"
#include "UHH2/common/include/JetCorrections.h"
#include "UHH2/common/include/CleaningModules.h"
#include "UHH2/common/include/CommonModules.h"
#include "UHH2/common/include/JetCorrections.h"
#include "UHH2/common/include/MCWeight.h"
#include "UHH2/common/include/EventVariables.h"
#include "UHH2/common/include/CleaningModules.h"
#include "UHH2/common/include/EventVariables.h"
#include "UHH2/common/include/NSelections.h"
#include "UHH2/common/include/AdditionalSelections.h"
#include "UHH2/common/include/LumiSelection.h"
#include "UHH2/common/include/TriggerSelection.h"
#include "UHH2/common/include/Utils.h"
#include "UHH2/common/include/TTbarReconstruction.h"


using namespace std;

float inv_mass(const LorentzVector& p4){ return p4.isTimelike() ? p4.mass() : -sqrt(-p4.mass2()); }

// A quick functions to get the DNN mean and std for normalisation
// a vector containing the values is returned
std::vector<double> get_DNN_vals(const std::string& file) {
  std::vector<double> vals;
  std::ifstream infile(file);
  std::string line;
  while (std::getline(infile, line)) {
    std::istringstream iss(line);
    vals.push_back(std::stod(line));
  }
  return vals;
}

NeuralNetworkModule::NeuralNetworkModule(Context& ctx, const std::string & ModelName, const std::string& ConfigName): NeuralNetworkBase(ctx, ModelName, ConfigName) {}

bool NeuralNetworkModule::process(Event & event, std::vector<double> values){
  CreateInputs(event, values);
  std::vector<std::string> inputs;
  inputs.reserve(LayerInputs.size());
  for (auto const& imap: LayerInputs) inputs.push_back(imap.first);

  tensorflow::run(NNsession, inputs, NNInputs, LayerOutputs, &NNoutputs);
  return true;
}

void NeuralNetworkModule::CreateInputs(Event & event, std::vector<double> values){
  // clear possible previous values
  NNInputs.clear();
  NNoutputs.clear();
  // add tensors to vector
  NNInputs.push_back( tensorflow::Tensor(tensorflow::DT_FLOAT, {1, values.size()}));
  for(uint i = 0; i < values.size(); i++){
    NNInputs.at(0).tensor<float, 2>()(0,i) = values.at(i);
  }

  // sanity check when including
  if (NNInputs.size()!=LayerInputs.size()) throw logic_error("NeuralNetworkIncluder.cxx: Create a number of inputs diffetent wrt. LayerInputs.size()="+to_string(LayerInputs.size()));
}

// ######################################
// ######################################
// ######################################

NeuralNetworkInputCreator::NeuralNetworkInputCreator(Context& ctx) {
  h_primlep = ctx.get_handle<FlavorParticle>("PrimaryLepton");
  h_do_masspoint = ctx.get_handle<bool>("do_masspoint");
  h_is_muevt = ctx.get_handle<bool>("is_muevt");

  // Check if MC or Data
  is_MC = ctx.get("dataset_type") == "MC";
  // Init random for masspoint
  std::srand(std::time(nullptr)); // Initialize random
}

std::vector<double> NeuralNetworkInputCreator::getInputs() {
  return DNNInputs;
}

bool NeuralNetworkInputCreator::setInputs(std::vector<double> vec) {
  DNNInputs = vec;
  return true;
}

bool NeuralNetworkInputCreator::createInputs(Event& event) {
  // push values to vector
  const FlavorParticle& lepton = event.get(h_primlep);
  std::vector<double> values;

  // Lepton
  values.push_back(lepton.pt());
  values.push_back(lepton.eta());
  values.push_back(lepton.phi());
  if(event.get(h_is_muevt)) values.push_back(event.muons->at(0).relIso());
  else values.push_back(event.electrons->at(0).relIso());

  // HOTVR jet 1
  uint i = 0;
  for (const auto & topjet : *event.topjets){
    if(i == 3) break;
    values.push_back(topjet.pt());
    values.push_back(topjet.eta());
    values.push_back(topjet.phi());
    values.push_back(topjet.tau1_groomed());
    values.push_back(topjet.tau2_groomed());
    values.push_back(topjet.tau3_groomed());
    values.push_back(topjet.subjets().size());
    // matching AK4 jets
    double radius = 600/topjet.pt();
    if(radius > 1.5) radius = 1.5;
    else if (radius < 0.1) radius = 0.1;
    //std::vector<jet> matched_jets; //unused atm, but may be useful at some point?
    double maxbtag = 0;
    for (const auto & jet : *event.jets) {
      if(deltaR(jet, topjet) < radius){
        //matched_jets.push_back(jet);
        if(jet.btag_DeepCSV() > maxbtag) maxbtag = jet.btag_DeepCSV();
      }
    }
    values.push_back(maxbtag);
    i++;
  }
  for(;i < 3;i++){
    values.push_back(0); // pt
    values.push_back(3.5); // eta
    values.push_back(4); // phi
    values.push_back(0); // tau1
    values.push_back(0); // tau2
    values.push_back(0); // tau3
    values.push_back(0); // subjets
    values.push_back(0); // max btag
  }

  // Neutrino
  std::vector<LorentzVector> neutrinos = NeutrinoReconstruction(lepton.v4(), event.met->v4());
  values.push_back(neutrinos.at(0).pt());
  values.push_back(neutrinos.at(0).eta());
  values.push_back(neutrinos.at(0).phi());

  // Event variables
  values.push_back(event.jets->size());
  values.push_back(event.topjets->size());

  if(values.size() != 33) {
    std::cout << "You are an idiot that can not count! " << values.size() << endl;
    throw "You are an idiot that can not count!";
  }

  DNNInputs = values;

  return true;
}

// ######################################
// ######################################
// ######################################

NeuralNetworkInputNormalizer::NeuralNetworkInputNormalizer(Context& ctx, const std::string& path) {
  // Define the path for the normalisations
  string DNN_path = path;
  string MeansName = DNN_path+"/means.txt";
  string StdsName = DNN_path+"/stds.txt";

  // get vector vector
  DNNInputs_mean = get_DNN_vals(MeansName);
  DNNInputs_std = get_DNN_vals(StdsName);
}

std::vector<double> NeuralNetworkInputNormalizer::normalizeInputs(uhh2::Event& event, std::vector<double> DNNInputs) {

  assert(DNNInputs.size() == DNNInputs_mean.size());
  assert(DNNInputs.size() == DNNInputs_std.size());


  for(uint i = 0; i < DNNInputs.size(); i++){
    DNNInputs.at(i) = (DNNInputs.at(i)-DNNInputs_mean.at(i))/DNNInputs_std.at(i);
  }

  return DNNInputs;

}

// ######################################
// ######################################
// ######################################

NeuralNetworkIncluder::NeuralNetworkIncluder(Context& ctx, bool parametrized) {
  is_parametrized = parametrized;
  path = "/nfs/dust/cms/user/flabe/CMSSW/CMSSW_10_2_10/src/UHH2/MLCorner/TstarNN/";
  if(parametrized) path += "Parametric";
  else path += "NonParametric";
  NNInputCreator.reset(new NeuralNetworkInputCreator(ctx));
  NNInputNormalizer.reset(new NeuralNetworkInputNormalizer(ctx, path));
  NNModule.reset(new NeuralNetworkModule(ctx, path+"/bestModel/model.pb", path+"/bestModel/model.config.pbtxt"));
  h_masspoint = ctx.get_handle<double>("masspoint");
  h_DNN_output = ctx.declare_event_output<double>("DNN_output");
}

bool NeuralNetworkIncluder::process(Event& event) {
  NNInputCreator->createInputs(event);
  std::vector<double> inputs = NNInputCreator->getInputs();
  if(is_parametrized) inputs.push_back(event.get(h_masspoint));

  // Normalizing
  inputs = NNInputNormalizer->normalizeInputs(event, inputs);

  NNInputCreator->setInputs(inputs);
  NNModule->process(event, NNInputCreator->getInputs());
  std::vector<tensorflow::Tensor> NNOutputs = NNModule->GetOutputs();
  event.set(h_DNN_output, (double)(NNOutputs.at(0).tensor<float, 2>()(0,0)));
  return true;
}

// ######################################
// ######################################
// ######################################

NeuralNetworkInputWriter::NeuralNetworkInputWriter(Context& ctx) {
  NNInputCreator.reset(new NeuralNetworkInputCreator(ctx));

  for(uint i = 0; i < 33; i++){
    handles.push_back(ctx.declare_event_output<double>("DNN_Input_"+std::to_string(i)));
  }

  // Define Handle for DNN
  h_DNN_Inputs = ctx.declare_event_output<std::vector<double>>("DNN_Inputs");
  h_masspoint = ctx.declare_event_output<double>("masspoint");

  is_MC = ctx.get("dataset_type") == "MC";
}

bool NeuralNetworkInputWriter::process(uhh2::Event& event) {
  NNInputCreator->createInputs(event);
  std::vector<double> Inputs = NNInputCreator->getInputs();
  event.set(h_DNN_Inputs, Inputs);

  int i = 0;
  for(auto & handle : handles){
    event.set(handle, Inputs.at(i));
    i++;
  }

  double masspoint = -1;
  if(is_MC){
    assert(event.genparticles);
    for(const GenParticle & gp : *event.genparticles){
      if((gp.pdgId() == 9000005 || gp.pdgId() == -9000005) && (gp.status()==23 || gp.status()==22)){
        masspoint = inv_mass(gp.v4());
      }
    }
  }
  if(masspoint == -1){ // if no TstarFound (eg Data or Background) set random between 700 and 1600
    masspoint = std::rand() % 1801 + 200; // random between 200 and 2000, inclusively!
  }
  event.set(h_masspoint, masspoint);


  return true;
}
