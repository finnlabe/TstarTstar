#include <iostream>
#include <memory>
#include <string>

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
#include "UHH2/TstarTstar/include/TstarTstarSelections.h"
#include "UHH2/TstarTstar/include/TstarTstarHists.h"
#include "UHH2/TstarTstar/include/TstarTstarDNNHists.h"
#include "UHH2/TstarTstar/include/TstarTstarAfterDNNHists.h"
#include "UHH2/TstarTstar/include/TstarTstarDNNInputHists.h"
#include "UHH2/TstarTstar/include/TstarTstarRecoTstarHists.h"
#include "UHH2/TstarTstar/include/TstarTstarAllGenHists.h"
#include "UHH2/TstarTstar/include/TstarTstarGenHists.h"
#include "UHH2/TstarTstar/include/TstarTstarGenRecoMatchedHists.h"
#include "UHH2/TstarTstar/include/TstarTstarReconstructionModules.h"
#include "UHH2/TstarTstar/include/ReconstructionTstarHypothesis.h"
#include "UHH2/TstarTstar/include/TstarTstarGenMatch.h"
#include "UHH2/HOTVR/include/HOTVRIds.h"

using namespace std;
using namespace uhh2;

namespace uhh2 {

/** \brief Module for the T*T*->ttbar gg MC based study
 *
 * This is the central class which calls other AnalysisModules, Hists or Selection classes.
 * This AnalysisModule, in turn, is called (via AnalysisModuleRunner) by SFrame.
 */
class TstarTstarAfterDNNModule: public AnalysisModule {
public:

  explicit TstarTstarAfterDNNModule(Context & ctx);
  virtual bool process(Event & event) override;

private:

  // ##### Histograms #####
  // Store the Hists collection as member variables. Again, use unique_ptr to avoid memory leaks.
  std::unique_ptr<Hists> h_DNN, h_DNN_lowST, h_DNN_highST, h_DNN_lowDNN, h_DNN_highDNN, h_DNN_highST_lowDNN, h_DNN_highST_highDNN;
  std::unique_ptr<Hists> h_DNNeval, h_DNNeval_lowST, h_DNNeval_highST, h_DNNeval_lowDNN, h_DNNeval_highDNN, h_DNNeval_highST_lowDNN, h_DNNeval_highST_highDNN;
  std::unique_ptr<Hists> h_Hists, h_Hists_lowST, h_Hists_highST, h_Hists_lowDNN, h_Hists_highDNN, h_Hists_highST_lowDNN, h_Hists_highST_highDNN;


  // Bools for Debugging/Options
  bool debug = false;

  // DNN stuff
  uhh2::Event::Handle<double> h_DNN_output;
  uhh2::Event::Handle<double> h_ST;

  std::unique_ptr<uhh2::AnalysisModule> reco_primlep;

  // Handles
  uhh2::Event::Handle<int> h_flag_toptagevent;
  uhh2::Event::Handle<int> h_flag_muonevent;

  uhh2::Event::Handle<double> h_evt_weight;

  // bools for channel and stuff. will be read in later
  bool is_MC;
  bool is_TTbar;
  bool is_Signal;
  bool is_Data;
};


TstarTstarAfterDNNModule::TstarTstarAfterDNNModule(Context & ctx){

  if(debug) {
    cout << "Hello World from TstarTstarAfterDNNModule!" << endl;

    // If running in SFrame, the keys "dataset_version", "dataset_type", "dataset_lumi",
    // and "target_lumi" are set to the according values in the xml file. For CMSSW, these are
    // not set automatically, but can be set in the python config file.
    for(auto & kv : ctx.get_all()){
        cout << " " << kv.first << " = " << kv.second << endl;
    }
   }

  // 0. Reading in whether MC and if so, which channel
  is_MC = ctx.get("dataset_type") == "MC";

  reco_primlep.reset(new PrimaryLepton(ctx));

  //TopTag
  h_flag_toptagevent = ctx.get_handle<int>("flag_toptagevent");

  // Check which lepton is present and save in
  h_flag_muonevent = ctx.get_handle<int>("flag_muonevent");

  // 3. Set up Hists classes:
  h_DNN.reset(new TstarTstarDNNHists(ctx, "DNN"));
  h_DNN_lowST.reset(new TstarTstarDNNHists(ctx, "DNN_lowST"));
  h_DNN_highST.reset(new TstarTstarDNNHists(ctx, "DNN_highST"));
  h_DNN_lowDNN.reset(new TstarTstarDNNHists(ctx, "DNN_lowDNN"));
  h_DNN_highDNN.reset(new TstarTstarDNNHists(ctx, "DNN_highDNN"));
  h_DNN_highST_lowDNN.reset(new TstarTstarDNNHists(ctx, "DNN_highST_lowDNN"));
  h_DNN_highST_highDNN.reset(new TstarTstarDNNHists(ctx, "DNN_highST_highDNN"));

  h_DNNeval.reset(new TstarTstarAfterDNNHists(ctx, "DNNeval"));
  h_DNNeval_lowST.reset(new TstarTstarAfterDNNHists(ctx, "DNNeval_lowST"));
  h_DNNeval_highST.reset(new TstarTstarAfterDNNHists(ctx, "DNNeval_highST"));
  h_DNNeval_lowDNN.reset(new TstarTstarAfterDNNHists(ctx, "DNNeval_lowDNN"));
  h_DNNeval_highDNN.reset(new TstarTstarAfterDNNHists(ctx, "DNNeval_highDNN"));
  h_DNNeval_highST_lowDNN.reset(new TstarTstarAfterDNNHists(ctx, "DNNeval_highST_lowDNN"));
  h_DNNeval_highST_highDNN.reset(new TstarTstarAfterDNNHists(ctx, "DNNeval_highST_highDNN"));

  h_Hists.reset(new TstarTstarHists(ctx, "Hists"));
  h_Hists_lowST.reset(new TstarTstarHists(ctx, "Hists_lowST"));
  h_Hists_highST.reset(new TstarTstarHists(ctx, "Hists_highST"));
  h_Hists_lowDNN.reset(new TstarTstarHists(ctx, "Hists_lowDNN"));
  h_Hists_highDNN.reset(new TstarTstarHists(ctx, "Hists_highDNN"));
  h_Hists_highST_lowDNN.reset(new TstarTstarHists(ctx, "Hists_highST_lowDNN"));
  h_Hists_highST_highDNN.reset(new TstarTstarHists(ctx, "Hists_highST_highDNN"));

  // 5. Handles for DNN
  h_DNN_output = ctx.get_handle<double>("DNN_output");
  h_evt_weight = ctx.get_handle<double>("evt_weight");

  is_TTbar = (ctx.get("dataset_version").find("TT") != std::string::npos);
  is_Signal = (ctx.get("dataset_version").find("Tstar") != std::string::npos);
}


bool TstarTstarAfterDNNModule::process(Event & event) {

  if(debug){cout << endl << "TstarTstarAfterDNNModule: Starting to process event (runid, eventid) = (" << event.run << ", " << event.event << "); weight = " << event.weight << endl;}

  // reapply weights
  event.weight = event.get(h_evt_weight);
  if(debug) cout << "weights applied." << endl;

  reco_primlep->process(event);//set "primary lepton"

  double st_jets = 0;
  for(const auto & jet : *event.topjets) st_jets += jet.pt();
  double DNNoutput = event.get(h_DNN_output);

  h_DNN->fill(event);
  h_DNNeval->fill(event);
  h_Hists->fill(event);

  if(st_jets>500){
    h_DNN_highST->fill(event);
    h_DNNeval_highST->fill(event);
    h_Hists_highST->fill(event);
    if(DNNoutput > 0.6) {
      h_DNN_highST_highDNN->fill(event);
      h_DNNeval_highST_highDNN->fill(event);
      h_Hists_highST_highDNN->fill(event);
    }
    else {
      h_DNN_highST_lowDNN->fill(event);
      h_DNNeval_highST_lowDNN->fill(event);
      h_Hists_highST_lowDNN->fill(event);
    }
  }
  else {
    h_DNN_lowST->fill(event);
    h_DNNeval_lowST->fill(event);
    h_Hists_lowST->fill(event);
  }

  if(DNNoutput > 0.6) {
    h_DNN_highDNN->fill(event);
    h_DNNeval_highDNN->fill(event);
    h_Hists_highDNN->fill(event);
  }
  else {
    h_DNN_lowDNN->fill(event);
    h_DNNeval_lowDNN->fill(event);
    h_Hists_lowDNN->fill(event);
  }

  return true;

}

// as we want to run the ExampleCycleNew directly with AnalysisModuleRunner,
// make sure the TstarTstarAfterDNNModule is found by class name. This is ensured by this macro:
UHH2_REGISTER_ANALYSIS_MODULE(TstarTstarAfterDNNModule)

}
