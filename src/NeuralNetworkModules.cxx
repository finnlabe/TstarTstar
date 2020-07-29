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
  values.push_back(lepton.energy());
  if(event.get(h_is_muevt)) values.push_back(event.muons->at(0).relIso());
  else values.push_back(event.electrons->at(0).relIso());
  double min_deltaR = 999;
  double pt_rel = -1;
  for(auto &jet : *event.jets){
    double cur_deltaR = deltaR(jet, lepton);
    if(cur_deltaR < min_deltaR) {
      min_deltaR = cur_deltaR;
      pt_rel = lepton.pt()/jet.pt();
    }
  }
  values.push_back(min_deltaR);
  values.push_back(pt_rel);

  // AK4 jet 1
  auto jet1 = event.jets->at(0);
  values.push_back(jet1.pt());
  values.push_back(jet1.eta());
  values.push_back(jet1.phi());
  values.push_back(inv_mass(jet1.v4()));
  values.push_back(jet1.btag_DeepCSV());
  // AK4 jet 2
  auto jet2 = event.jets->at(1);
  values.push_back(jet2.pt());
  values.push_back(jet2.eta());
  values.push_back(jet2.phi());
  values.push_back(inv_mass(jet2.v4()));
  values.push_back(jet2.btag_DeepCSV());
  // AK4 jet 3
  auto jet3 = event.jets->at(2);
  values.push_back(jet3.pt());
  values.push_back(jet3.eta());
  values.push_back(jet3.phi());
  values.push_back(inv_mass(jet3.v4()));
  values.push_back(jet3.btag_DeepCSV());
  // AK4 jet 4
  auto jet4 = event.jets->at(3);
  values.push_back(jet4.pt());
  values.push_back(jet4.eta());
  values.push_back(jet4.phi());
  values.push_back(inv_mass(jet4.v4()));
  values.push_back(jet4.btag_DeepCSV());

  // HOTVR jet 1
  TopJet topjet1 = event.topjets->at(0);
  values.push_back(topjet1.pt());
  values.push_back(topjet1.eta());
  values.push_back(topjet1.phi());
  values.push_back(inv_mass(topjet1.v4()));
  values.push_back(topjet1.tau1_groomed());
  values.push_back(topjet1.tau2_groomed());
  values.push_back(topjet1.tau3_groomed());
  values.push_back(topjet1.subjets().size());
  // HOTVR jet 1
  if(event.topjets->size() > 1){ // may not be guaranteed
    TopJet topjet2 = event.topjets->at(1);
    values.push_back(topjet2.pt());
    values.push_back(topjet2.eta());
    values.push_back(topjet2.phi());
    values.push_back(inv_mass(topjet2.v4()));
    values.push_back(topjet2.tau1_groomed());
    values.push_back(topjet2.tau2_groomed());
    values.push_back(topjet2.tau3_groomed());
    values.push_back(topjet2.subjets().size());
  }
  else { // fill dummy values
    values.push_back(0); // pt
    values.push_back(3.5); // eta
    values.push_back(4); // phi
    values.push_back(0); // inv mass
    values.push_back(0); // tau1
    values.push_back(0); // tau2
    values.push_back(0); // tau2
    values.push_back(0); // subjets
  }

  // Neutrino
  std::vector<LorentzVector> neutrinos = NeutrinoReconstruction(lepton.v4(), event.met->v4());
  values.push_back(neutrinos.at(0).pt());
  values.push_back(neutrinos.at(0).eta());
  values.push_back(neutrinos.at(0).phi());

  // Event variables
  values.push_back(event.jets->size());
  values.push_back(event.topjets->size());
  // sphericity
  double bottom_sum = 0;
  for (const auto & jet : *event.jets){
    bottom_sum += jet.v4().P2();
  }
  double upper_sum;
  upper_sum = 0;
  for (const auto & jet : *event.jets){ // px px
    upper_sum += jet.v4().px()*jet.v4().px();
  }
  values.push_back(upper_sum/bottom_sum);
  upper_sum = 0;
  for (const auto & jet : *event.jets){ // px py
    upper_sum += jet.v4().px()*jet.v4().py();
  }
  values.push_back(upper_sum/bottom_sum);
  upper_sum = 0;
  for (const auto & jet : *event.jets){ // px pz
    upper_sum += jet.v4().px()*jet.v4().pz();
  }
  values.push_back(upper_sum/bottom_sum);
  upper_sum = 0;
  for (const auto & jet : *event.jets){ // py py
    upper_sum += jet.v4().py()*jet.v4().py();
  }
  values.push_back(upper_sum/bottom_sum);
  upper_sum = 0;
  for (const auto & jet : *event.jets){ // py pz
    upper_sum += jet.v4().py()*jet.v4().pz();
  }
  values.push_back(upper_sum/bottom_sum);
  upper_sum = 0;
  for (const auto & jet : *event.jets){ // pz pz
    upper_sum += jet.v4().pz()*jet.v4().pz();
  }
  values.push_back(upper_sum/bottom_sum);

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
  NNModule.reset(new NeuralNetworkModule(ctx, path+"/model.pb", path+"/model.config.pbtxt"));
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

  for(uint i = 0; i < 54; i++){
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
