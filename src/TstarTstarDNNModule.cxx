#include <iostream>
#include <memory>
#include <string>

// UHH2 stuff
#include "UHH2/core/include/AnalysisModule.h"
#include "UHH2/core/include/Event.h"
#include "UHH2/common/include/ElectronHists.h"
#include "UHH2/common/include/NSelections.h"
#include "UHH2/common/include/ElectronIds.h"
#include "UHH2/common/include/PhotonIds.h"
#include <UHH2/common/include/MuonIds.h>
#include "UHH2/common/include/TTbarGen.h"
#include "UHH2/common/include/TopJetIds.h"
#include "UHH2/common/include/MCWeight.h"

// TstarTstar
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

// other
#include "UHH2/HOTVR/include/HOTVRIds.h"

using namespace std;
using namespace uhh2;

namespace uhh2 {

/** \brief Module for the T*T*->ttbar gg MC based study
 *
 * This is the central class which calls other AnalysisModules, Hists or Selection classes.
 * This AnalysisModule, in turn, is called (via AnalysisModuleRunner) by SFrame.
 */
class TstarTstarDNNModule: public AnalysisModule {
public:

  explicit TstarTstarDNNModule(Context & ctx);
  virtual bool process(Event & event) override;

private:

  // ###### Modules ######
  std::unique_ptr<uhh2::AnalysisModule> reco_primlep;
  std::unique_ptr<NeuralNetworkIncluder> DNN_Includer;


  // ##### Histograms #####
  std::unique_ptr<Hists> h_topcheck, h_topcheck_reweighted;
  std::unique_ptr<Hists> h_STreweighted;

  std::unique_ptr<Hists> h_AfterDNNcut_02, h_AfterDNNcut_03, h_AfterDNNcut_04, h_AfterDNNcut_05, h_AfterDNNcut_06, h_AfterDNNcut_07, h_AfterDNNcut_08;
  std::unique_ptr<Hists> h_notDNNcut_02,   h_notDNNcut_03,   h_notDNNcut_04,   h_notDNNcut_05,   h_notDNNcut_06,   h_notDNNcut_07,   h_notDNNcut_08;

  std::unique_ptr<Hists> h_DNN, h_DNN_reweighted, h_DNN_reweighted_2, h_DNN_lowpt, h_DNN_medpt, h_DNN_highpt;
  std::unique_ptr<Hists> h_DNN_lowST, h_DNN_highST, h_DNN_lowDNN, h_DNN_highDNN, h_DNN_highST_lowDNN, h_DNN_highST_highDNN;
  std::unique_ptr<Hists> h_AfterDNN, h_AfterDNN_lowST, h_AfterDNN_highST, h_AfterDNN_lowDNN, h_AfterDNN_highDNN, h_AfterDNN_highST_lowDNN, h_AfterDNN_highST_highDNN;


  // ###### Control switches ######
  bool debug = false;
  bool do_masspoint = false;
  bool doAddInputs = false;


  // ###### Handles ######
  uhh2::Event::Handle<int> h_flag_toptagevent;
  uhh2::Event::Handle<int> h_flag_muonevent;
  uhh2::Event::Handle<double> h_evt_weight;

  uhh2::Event::Handle<double> h_ST_weight;

  uhh2::Event::Handle<double> h_DNN_output;
  uhh2::Event::Handle<bool> h_do_masspoint;
  uhh2::Event::Handle<double> h_ST;
  uhh2::Event::Handle<bool> h_DoAddInputs;


  // ###### other parameters ######
  bool is_MC;
  bool is_TTbar;
  bool is_Signal;
  bool is_Data;
};


TstarTstarDNNModule::TstarTstarDNNModule(Context & ctx){

  // debug message
  if(debug) {
    cout << "Hello World from TstarTstarDNNModule!" << endl;
    // If running in SFrame, the keys "dataset_version", "dataset_type", "dataset_lumi",
    // and "target_lumi" are set to the according values in the xml file. For CMSSW, these are
    // not set automatically, but can be set in the python config file.
    for(auto & kv : ctx.get_all()){
        cout << " " << kv.first << " = " << kv.second << endl;
    }
  }

  // ###### 0. Setting variables ######
  // MC or real data
  is_MC = ctx.get("dataset_type") == "MC";

  is_TTbar = (ctx.get("dataset_version").find("TT") != std::string::npos);
  is_Signal = (ctx.get("dataset_version").find("Tstar") != std::string::npos);

  // ###### 1. set up modules ######
  // primlep
  reco_primlep.reset(new PrimaryLepton(ctx));

  // DNN modules
  DNN_Includer.reset(new NeuralNetworkIncluder(ctx, do_masspoint));


  // 3. Set up Hists classes:
  h_topcheck.reset(new TstarTstarAllGenHists(ctx, "topcheck"));
  h_topcheck_reweighted.reset(new TstarTstarAllGenHists(ctx, "topcheck_reweighted"));

  h_STreweighted.reset(new TstarTstarHists(ctx, "topcheck"));

  h_AfterDNNcut_02.reset(new TstarTstarHists(ctx, "AfterDNNcut_02"));
  h_notDNNcut_02.reset(new TstarTstarHists(ctx, "notDNNcut_02"));
  h_AfterDNNcut_03.reset(new TstarTstarHists(ctx, "AfterDNNcut_03"));
  h_notDNNcut_03.reset(new TstarTstarHists(ctx, "notDNNcut_03"));
  h_AfterDNNcut_04.reset(new TstarTstarHists(ctx, "AfterDNNcut_04"));
  h_notDNNcut_04.reset(new TstarTstarHists(ctx, "notDNNcut_04"));
  h_AfterDNNcut_05.reset(new TstarTstarHists(ctx, "AfterDNNcut_05"));
  h_notDNNcut_05.reset(new TstarTstarHists(ctx, "notDNNcut_05"));
  h_AfterDNNcut_06.reset(new TstarTstarHists(ctx, "AfterDNNcut_06"));
  h_notDNNcut_06.reset(new TstarTstarHists(ctx, "notDNNcut_06"));
  h_AfterDNNcut_07.reset(new TstarTstarHists(ctx, "AfterDNNcut_07"));
  h_notDNNcut_07.reset(new TstarTstarHists(ctx, "notDNNcut_07"));
  h_AfterDNNcut_08.reset(new TstarTstarHists(ctx, "AfterDNNcut_08"));
  h_notDNNcut_08.reset(new TstarTstarHists(ctx, "notDNNcut_08"));

  h_DNN.reset(new TstarTstarDNNHists(ctx, "DNN"));
  h_DNN_reweighted.reset(new TstarTstarDNNHists(ctx, "DNN_reweighted"));
  h_DNN_reweighted_2.reset(new TstarTstarDNNHists(ctx, "DNN_reweighted_2"));
  h_DNN_lowpt.reset(new TstarTstarDNNHists(ctx, "DNN_lowpt"));
  h_DNN_medpt.reset(new TstarTstarDNNHists(ctx, "DNN_medpt"));
  h_DNN_highpt.reset(new TstarTstarDNNHists(ctx, "DNN_highpt"));
  h_DNN_lowST.reset(new TstarTstarDNNHists(ctx, "DNN_lowST"));
  h_DNN_highST.reset(new TstarTstarDNNHists(ctx, "DNN_highST"));
  h_DNN_lowDNN.reset(new TstarTstarDNNHists(ctx, "DNN_lowDNN"));
  h_DNN_highDNN.reset(new TstarTstarDNNHists(ctx, "DNN_highDNN"));
  h_DNN_highST_lowDNN.reset(new TstarTstarDNNHists(ctx, "DNN_highST_lowDNN"));
  h_DNN_highST_highDNN.reset(new TstarTstarDNNHists(ctx, "DNN_highST_highDNN"));

  h_AfterDNN.reset(new TstarTstarHists(ctx, "AfterDNN"));
  h_AfterDNN_lowST.reset(new TstarTstarHists(ctx, "AfterDNN_lowST"));
  h_AfterDNN_highST.reset(new TstarTstarHists(ctx, "AfterDNN_highST"));
  h_AfterDNN_lowDNN.reset(new TstarTstarHists(ctx, "AfterDNN_lowDNN"));
  h_AfterDNN_highDNN.reset(new TstarTstarHists(ctx, "AfterDNN_highDNN"));
  h_AfterDNN_highST_lowDNN.reset(new TstarTstarHists(ctx, "AfterDNN_highST_lowDNN"));
  h_AfterDNN_highST_highDNN.reset(new TstarTstarHists(ctx, "AfterDNN_highST_highDNN"));

  // ###### 4. init handles ######
  h_evt_weight = ctx.get_handle<double>("evt_weight");
  h_flag_toptagevent = ctx.get_handle<int>("flag_toptagevent");
  h_flag_muonevent = ctx.get_handle<int>("flag_muonevent");
  h_DoAddInputs = ctx.declare_event_output<bool>("doAddInputs");

  h_ST_weight = ctx.declare_event_output<double>("ST_weight");

  h_DNN_output = ctx.get_handle<double>("DNN_output");

}


bool TstarTstarDNNModule::process(Event & event) {

  // debug message
  if(debug){cout << endl << "TstarTstarDNNModule: Starting to process event (runid, eventid) = (" << event.run << ", " << event.event << "); weight = " << event.weight << endl;}

  // reapply weights
  event.weight = event.get(h_evt_weight);
  if(debug) cout << "weights applied." << endl;

  // get ST weights
  double ST_weight = event.get(h_ST_weight);

  // set primary lepton
  reco_primlep->process(event);

  // setting addInputs
  event.set(h_DoAddInputs, doAddInputs);

  event.weight = ST_weight * event.get(h_evt_weight);
  h_STreweighted->fill(event);
  event.weight = event.get(h_evt_weight);

  // claculating ST
  double st_jets = 0;
  for(const auto & jet : *event.topjets) st_jets += jet.pt();
  for(const auto & lepton : *event.electrons) st_jets += lepton.pt();
  for(const auto & lepton : *event.muons) st_jets += lepton.pt();

  // ################
  // ### DNN Part ###
  // ################

  // including trained DNN model
  if(debug) cout << "Include DNN model" << endl;
  DNN_Includer->process(event);

  // hists
  if(debug) std::cout << "Plotting" << endl;
  h_DNN->fill(event);
  h_topcheck->fill(event);
  if(is_TTbar) event.weight *= ST_weight;
  h_DNN_reweighted->fill(event);
  h_topcheck_reweighted->fill(event);
  event.weight = event.get(h_evt_weight);
  if(event.topjets->at(0).pt() < 500) h_DNN_lowpt->fill(event);
  else if (event.topjets->at(0).pt() < 1000) h_DNN_medpt->fill(event);
  else h_DNN_highpt->fill(event);

  if(is_MC) { // blinding!!!
    // do DNN cuts
    if(event.get(h_DNN_output) > 0.2) h_AfterDNNcut_02->fill(event);
    else h_notDNNcut_02->fill(event);
    if(event.get(h_DNN_output) > 0.3) h_AfterDNNcut_03->fill(event);
    else h_notDNNcut_03->fill(event);
    if(event.get(h_DNN_output) > 0.4) h_AfterDNNcut_04->fill(event);
    else h_notDNNcut_04->fill(event);
    if(event.get(h_DNN_output) > 0.5) h_AfterDNNcut_05->fill(event);
    else h_notDNNcut_05->fill(event);
    if(event.get(h_DNN_output) > 0.6) h_AfterDNNcut_06->fill(event);
    else h_notDNNcut_06->fill(event);
    if(event.get(h_DNN_output) > 0.7) h_AfterDNNcut_07->fill(event);
    else h_notDNNcut_07->fill(event);
    if(event.get(h_DNN_output) > 0.8) h_AfterDNNcut_08->fill(event);
    else h_notDNNcut_08->fill(event);
  }


  // some more plotting
  double DNNoutput = event.get(h_DNN_output);
  h_AfterDNN->fill(event);

  if(st_jets>500){
    h_DNN_highST->fill(event);
    h_AfterDNN_highST->fill(event);
    if(DNNoutput > 0.6) {
      h_DNN_highST_highDNN->fill(event);
      h_AfterDNN_highST_highDNN->fill(event);
    }
    else {
      h_DNN_highST_lowDNN->fill(event);
      h_AfterDNN_highST_lowDNN->fill(event);
    }
  }
  else {
    h_DNN_lowST->fill(event);
    h_AfterDNN_lowST->fill(event);
  }

  if(DNNoutput > 0.6) {
    h_DNN_highDNN->fill(event);
    h_AfterDNN_highDNN->fill(event);
  }
  else {
    h_DNN_lowDNN->fill(event);
    h_AfterDNN_lowDNN->fill(event);
  }

  if(debug){cout << "Done ##################################" << endl;}
  return true;

}

// as we want to run the ExampleCycleNew directly with AnalysisModuleRunner,
// make sure the TstarTstarDNNModule is found by class name. This is ensured by this macro:
UHH2_REGISTER_ANALYSIS_MODULE(TstarTstarDNNModule)

}
