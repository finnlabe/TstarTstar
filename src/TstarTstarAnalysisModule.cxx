#include <iostream>
#include <memory>
#include <string>

// UHH2 stuff
#include "UHH2/core/include/AnalysisModule.h"
#include "UHH2/core/include/Event.h"
#include "UHH2/common/include/CommonModules.h"
#include "UHH2/common/include/LuminosityHists.h"
#include "UHH2/common/include/CleaningModules.h"
#include "UHH2/common/include/ElectronHists.h"
#include "UHH2/common/include/NSelections.h"
#include "UHH2/common/include/ElectronIds.h"
#include "UHH2/common/include/PhotonIds.h"
#include <UHH2/common/include/MuonIds.h>
#include <UHH2/common/include/TriggerSelection.h>
#include "UHH2/common/include/TTbarGen.h"
#include "UHH2/common/include/TopJetIds.h"
#include "UHH2/common/include/MCWeight.h"

// TstarTstar stuff
#include "UHH2/TstarTstar/include/TstarTstarSelections.h"
#include "UHH2/TstarTstar/include/TstarTstarHists.h"
#include "UHH2/TstarTstar/include/TstarTstarDNNHists.h"
#include "UHH2/TstarTstar/include/TstarTstarDNNInputHists.h"
#include "UHH2/TstarTstar/include/TstarTstarRecoTstarHists.h"
#include "UHH2/TstarTstar/include/TstarTstarAllGenHists.h"
#include "UHH2/TstarTstar/include/TstarTstarGenHists.h"
#include "UHH2/TstarTstar/include/TstarTstarGenRecoMatchedHists.h"
#include "UHH2/TstarTstar/include/TstarTstarReconstructionModules.h"
#include "UHH2/TstarTstar/include/ReconstructionTstarHypothesis.h"
#include "UHH2/TstarTstar/include/TstarTstarGenMatch.h"
#include "UHH2/TstarTstar/include/NeuralNetworkModules.h"

// other stuff
#include "UHH2/HOTVR/include/HOTVRIds.h"

using namespace std;
using namespace uhh2;

namespace uhh2 {

// quick method to calculate inv_mass
float inv_mass(const LorentzVector& p4){ return p4.isTimelike() ? p4.mass() : -sqrt(-p4.mass2()); }

/** \brief Module for the T*T*->ttbar gg MC based study
 *
 * This is the central class which calls other AnalysisModules, Hists or Selection classes.
 * This AnalysisModule, in turn, is called (via AnalysisModuleRunner) by SFrame.
 */
class TstarTstarAnalysisModule: public AnalysisModule {
public:

  explicit TstarTstarAnalysisModule(Context & ctx);
  virtual bool process(Event & event) override;

private:

  // ###### Modules ######
  unique_ptr<Selection> toptagevt_sel;
  unique_ptr<NeuralNetworkInputWriter> DNN_InputWriter;

  std::unique_ptr<uhh2::AnalysisModule> ttgenprod;
  std::unique_ptr<uhh2::AnalysisModule> reco_primlep;


  // ##### Histograms #####
  std::unique_ptr<Hists> h_main, h_main_ttag, h_main_nottag, h_main_mu, h_main_mu_lowpt, h_main_mu_highpt, h_main_ele, h_main_ele_lowpt, h_main_ele_highpt;
  std::unique_ptr<Hists> h_crosscheck;
  std::unique_ptr<Hists> h_main_gen;

  std::unique_ptr<Hists> h_DNN_Inputs;


  // ###### Handles ######
  uhh2::Event::Handle<double> h_evt_weight;
  uhh2::Event::Handle<FlavorParticle> h_primlep;
  uhh2::Event::Handle<double> h_ST;
  uhh2::Event::Handle<LorentzVector> h_neutrino;

  // for reconstruction
  uhh2::Event::Handle<TTbarGen> h_ttbargen;
  uhh2::Event::Handle<int> h_flag_toptagevent;
  uhh2::Event::Handle<int> h_flag_muonevent;

  // for DNN output
  uhh2::Event::Handle<bool> h_do_masspoint;

  // ###### Control Switches ######
  bool debug = false;
  bool outputDNNvalues = true;
  bool do_masspoint = false;


  // ###### other needed definitions ######
  bool isTrigger;
  bool is_MC;

  bool is_TTbar;
  bool is_Signal;
  bool is_Data;

  TH1D* ST_ratio;
  uhh2::Event::Handle<double> h_ST_weight;
};


TstarTstarAnalysisModule::TstarTstarAnalysisModule(Context & ctx){

  // debug messagt
  if(debug) {
    cout << "Hello World from TstarTstarAnalysisModule!" << endl;
    // If running in SFrame, the keys "dataset_version", "dataset_type", "dataset_lumi",
    // and "target_lumi" are set to the according values in the xml file. For CMSSW, these are
    // not set automatically, but can be set in the python config file.
    for(auto & kv : ctx.get_all()){
        cout << " " << kv.first << " = " << kv.second << endl;
    }
  }

  // ###### 0. Setting Variables ######
  // MC or real data
  is_MC = ctx.get("dataset_type") == "MC";

  is_TTbar = (ctx.get("dataset_version").find("TT") != std::string::npos);
  is_Signal = (ctx.get("dataset_version").find("Tstar") != std::string::npos);


  // ###### 1. Set up modules ######
  // primary lepton
  reco_primlep.reset(new PrimaryLepton(ctx));

  // GEN things
  if(is_MC){
    ttgenprod.reset(new TTbarGenProducer(ctx, "ttbargen", false));
  }

  // top tag definition
  TopJetId topjetID = AndId<TopJet>(HOTVRTopTag(), Tau32Groomed(0.56));
  toptagevt_sel.reset(new TopTagEventSelection(topjetID));;


  // ###### 3. Set up histograms ######
  // before Reconstruction
  h_crosscheck.reset(new TstarTstarHists(ctx, "crosscheck"));
  h_main.reset(new TstarTstarHists(ctx, "main"));
  h_main_ttag.reset(new TstarTstarHists(ctx, "main_ttag"));
  h_main_nottag.reset(new TstarTstarHists(ctx, "main_nottag"));
  h_main_mu.reset(new TstarTstarHists(ctx, "main_mu"));
  h_main_mu_lowpt.reset(new TstarTstarHists(ctx, "main_mu_lowpt"));
  h_main_mu_highpt.reset(new TstarTstarHists(ctx, "main_mu_highpt"));
  h_main_ele.reset(new TstarTstarHists(ctx, "main_ele"));
  h_main_ele_lowpt.reset(new TstarTstarHists(ctx, "main_ele_lowpt"));
  h_main_ele_highpt.reset(new TstarTstarHists(ctx, "main_ele_highpt"));
  h_main_gen.reset(new TstarTstarGenHists(ctx, "main_gen"));

  // DNN hists
  h_DNN_Inputs.reset(new TstarTstarDNNInputHists(ctx, "DNN_Inputs"));


  // ###### 4. Init handles ######
  h_evt_weight = ctx.get_handle<double>("evt_weight");
  h_primlep = ctx.get_handle<FlavorParticle>("PrimaryLepton");
  h_flag_muonevent = ctx.declare_event_output<int>("flag_muonevent");
  h_flag_toptagevent = ctx.declare_event_output<int>("flag_toptagevent");
  h_neutrino = ctx.declare_event_output<LorentzVector>("neutrino");

  if(is_MC) h_ttbargen = ctx.get_handle<TTbarGen>("ttbargen");

  // DNN output
  if(outputDNNvalues){
    DNN_InputWriter.reset(new NeuralNetworkInputWriter(ctx));
    h_do_masspoint = ctx.get_handle<bool>("do_masspoint");
    h_ST = ctx.declare_event_output<double>("ST");
    h_ST_weight = ctx.declare_event_output<double>("ST_weight");
  }

  TFile *f = new TFile("/nfs/dust/cms/user/flabe/MLCorner/TstarNN/reweightingApproach/output/data/ST_weights.root");
  ST_ratio = (TH1D*)f->Get("ST_ratio");

}


bool TstarTstarAnalysisModule::process(Event & event) {

  // debug message
  if(debug){cout << endl << "TstarTstarAnalysisModule: Starting to process event (runid, eventid) = (" << event.run << ", " << event.event << "); weight = " << event.weight << endl;}

  // reapply weights
  event.weight = event.get(h_evt_weight);
  if(debug) cout << "weights applied." << endl;

  h_crosscheck->fill(event);

  // set lepton channel
  const bool muon_evt = (event.muons->size() == 1);
  event.set(h_flag_muonevent, int(muon_evt));

  // set primlep
  reco_primlep->process(event);

  // check ttag
  const bool pass_ttag = toptagevt_sel->passes(event);
  event.set(h_flag_toptagevent, int(pass_ttag));

  // ttgen
  if(is_MC) ttgenprod->process(event);

  // neutrinoreconstruction
  const Particle& lepton = event.get(h_primlep); // Primary Lepton has to be set
  std::vector<LorentzVector> neutrinos = NeutrinoReconstruction(lepton.v4(), event.met->v4());
  if(debug) std::cout << "We have this many neutrino options: " << neutrinos.size() << std::endl;
  for(auto &ntr : neutrinos) {
    if(debug) std::cout << "Neutrino pt: " << ntr.pt() << std::endl;
    double Wmass = inv_mass(ntr+lepton.v4());
    if(debug) std::cout << "W mass: " << Wmass << std::endl;
  }
  event.set(h_neutrino, neutrinos.at(0));
  // TODO find some better way to select best neutrino reconstruction?

  if(debug) std::cout << "Hists before everything" << std::endl;
  // fill hists before things have happened
  h_main->fill(event);
  h_main_gen->fill(event);
  if(pass_ttag) h_main_ttag->fill(event);
  else h_main_nottag->fill(event);
  if(event.get(h_flag_muonevent)){
    h_main_mu->fill(event);
    if(event.get(h_primlep).pt()<60) h_main_mu_lowpt->fill(event);
    else h_main_mu_highpt->fill(event);
  }
  else {
    h_main_ele->fill(event);
    if(event.get(h_primlep).pt()<120) h_main_ele_lowpt->fill(event);
    else h_main_ele_highpt->fill(event);
  }

  // ########################################
  // ########### DNN Preparation ############
  // ########################################

  if(debug) cout << "Start DNN stuff" << endl;

  // Filling output for DNN
  if(outputDNNvalues){
    event.set(h_do_masspoint, do_masspoint);
    if(debug) cout << "Write inputs" << endl;
    DNN_InputWriter->process(event);
    if(debug) cout << "plot inputs" << endl;
    h_DNN_Inputs->fill(event);
    double st_jets = 0.;
    for(const auto & jet : *event.topjets) st_jets += jet.pt();
    for(const auto & lepton : *event.electrons) st_jets += lepton.pt();
    for(const auto & lepton : *event.muons) st_jets += lepton.pt();
    event.set(h_ST, st_jets);
    double ST_weight = 1;
    if(is_TTbar) ST_weight = 1/(ST_ratio->GetBinContent(ST_ratio->GetXaxis()->FindBin(st_jets)));
    event.set(h_ST_weight, ST_weight);
  }

  if(debug){cout << "Done ##################################" << endl;}
  return true;

}

// as we want to run the ExampleCycleNew directly with AnalysisModuleRunner,
// make sure the TstarTstarAnalysisModule is found by class name. This is ensured by this macro:
UHH2_REGISTER_ANALYSIS_MODULE(TstarTstarAnalysisModule)

}
