#include <iostream>
#include <memory>

#include "UHH2/core/include/AnalysisModule.h"
#include "UHH2/core/include/Event.h"
#include "UHH2/common/include/CommonModules.h"
#include "UHH2/common/include/CleaningModules.h"
#include "UHH2/common/include/ElectronHists.h"
#include "UHH2/common/include/NSelections.h"
#include "UHH2/common/include/ElectronIds.h"
#include "UHH2/common/include/PhotonIds.h"
#include "UHH2/common/include/MuonIds.h"
#include "UHH2/common/include/TriggerSelection.h"
#include "UHH2/TstarTstar/include/ModuleBASE.h"
#include "UHH2/TstarTstar/include/TstarTstarSelections.h"
#include "UHH2/TstarTstar/include/TstarTstarHists.h"
#include "UHH2/TstarTstar/include/TstarTstarGenHists.h"
#include "UHH2/TstarTstar/include/TstarTstarGenRecoMatchedHists.h"
#include "UHH2/common/include/TTbarGen.h"
#include "UHH2/TstarTstar/include/TstarTstarRecoTstarHists.h"
#include "UHH2/TstarTstar/include/TstarTstarReconstructionModules.h"
#include "UHH2/HOTVR/include/HOTVRIds.h"

#include "UHH2/common/include/MCWeight.h"

using namespace std;
using namespace uhh2;

//namespace uhh2 {

/** \brief Module for the full selection in T*T*->ttbar gg/gamma search
 *
 * All objects are expected to be corrected in PreSelection stage
 *
 */
class TstarTstarSelectionModule: public AnalysisModule {
public:

    explicit TstarTstarSelectionModule(Context & ctx);
    virtual bool process(Event & event) override;

private:

  // Declare the Selections to use. Use unique_ptr to ensure automatic call of delete in the destructor,
  // to avoid memory leaks.

  unique_ptr<AnalysisModule> LumiWeight_module;
  unique_ptr<uhh2::AnalysisModule> reco_primlep;
  uhh2::Event::Handle<FlavorParticle> h_primlep;

  // Store the Hists collection as member variables. Again, use unique_ptr to avoid memory leaks.
  std::unique_ptr<Hists> h_AK4cleaning,     h_AK8cleaning,     h_METsel,     h_2Dcut,     h_ttagsel,     h_phisel;
  std::unique_ptr<Hists> h_AK4cleaning_gen, h_AK8cleaning_gen, h_METsel_gen, h_2Dcut_gen, h_ttagsel_gen, h_phisel_gen;

  std::unique_ptr<Hists> h_AK8jetsel;
  std::unique_ptr<Hists> h_AK8jetsel_gen;

  unique_ptr<Selection> met_sel;
  unique_ptr<Selection> st_sel;
  unique_ptr<Selection> twodcut_sel;
  unique_ptr<Selection> toptagevt_sel;

  unique_ptr<JetCleaner> AK4cleaner;
  unique_ptr<TopJetCleaner> AK8cleaner;

  // GEN stuff
  std::unique_ptr<uhh2::AnalysisModule> ttgenprod;
  uhh2::Event::Handle<TTbarGen> h_ttbargen;

  bool debug = false;

  // bools for channel and stuff. will be read in later
  bool is_MC;
  bool is_tgtg, is_tgtgamma;

  TopJetId topjetID;

};

TstarTstarSelectionModule::TstarTstarSelectionModule(Context & ctx){

  // 0. Reading in whether MC and if so, which channel
  is_MC = ctx.get("dataset_type") == "MC";

  if(debug) {
    cout << "Hello World from TstarTstarSelectionModule!" << endl;

    // If running in SFrame, the keys "dataset_version", "dataset_type", "dataset_lumi",
    // and "target_lumi" are set to the according values in the xml file. For CMSSW, these are
    // not set automatically, but can be set in the python config file.
    for(auto & kv : ctx.get_all()){
      cout << " " << kv.first << " = " << kv.second << endl;
    }
  }

  // 1. set up lumi rewitghting
  LumiWeight_module.reset(new MCLumiWeight(ctx));

  if(is_MC){
    // Prepare GEN
    ttgenprod.reset(new TTbarGenProducer(ctx, "ttbargen", false));
  }

  // 2. set up selections
  if(debug) cout << "Setting up Selections." << endl;

  // 2D cut
  twodcut_sel.reset(new TwoDCut(0.4, 25.0));  // The same as in Z'->ttbar semileptonic

  // Top Tagging
  /**
  topjetID = AndId<TopJet>(HOTVRTopTag(), Tau32Groomed(0.56));
  toptagevt_sel.reset(new TopTagEventSelection(topjetID));
  **/

  // 4. Set up Hists
  if(debug) cout << "Setting up Hists." << endl;
  h_2Dcut.reset(new TstarTstarHists(ctx, "After2D"));
  //h_ttagsel.reset(new TstarTstarHists(ctx, "AfterTtagsel"));

  h_2Dcut_gen.reset(new TstarTstarGenHists(ctx, "After2D_gen"));
  //h_ttagsel_gen.reset(new TstarTstarGenHists(ctx, "AfterTtagsel_gen"));

}


bool TstarTstarSelectionModule::process(Event & event) {

  if(debug) cout << "TstarTstarSelectionModule: Starting to process event (runid, eventid) = (" << event.run << ", " << event.event << "); weight = " << event.weight << endl;

  LumiWeight_module->process(event); // apply correct weights
  if(debug) cout << "Lumi weights applied." << endl;

  //Fill ttgen object for correct matching check, etc
  if(is_MC){
    ttgenprod->process(event);
    if(debug) cout << "ttgen produced." << endl;
  }

  // Lepton-2Dcut
  for(auto& muo : *event.muons){
    if(debug) cout<<"AFTER Muon (pt,eta): "<<muo.pt()<<", "<<muo.eta()<<endl;
    float    dRmin, pTrel;
    std::tie(dRmin, pTrel) = drmin_pTrel(muo, *event.jets);
    muo.set_tag(Muon::twodcut_dRmin, dRmin);
    muo.set_tag(Muon::twodcut_pTrel, pTrel);
  }
  for(auto& ele : *event.electrons){
    if(debug) cout<<"Electron (pt,eta): "<<ele.pt()<<", "<<ele.eta()<<endl;
    float    dRmin, pTrel;
    std::tie(dRmin, pTrel) = drmin_pTrel(ele, *event.jets);
    ele.set_tag(Electron::twodcut_dRmin, dRmin);
    ele.set_tag(Electron::twodcut_pTrel, pTrel);
  }
  const bool pass_twodcut = twodcut_sel->passes(event);
  if(!pass_twodcut) return false;
  h_2Dcut->fill(event);
  h_2Dcut_gen->fill(event);
  if(debug) cout << "Passed 2D cut." << endl;

  // TopTagEventSelection
  /**
  bool pass_ttag = (toptagevt_sel->passes(event));
  if(!pass_ttag) return false;
  h_ttagsel->fill(event);
  h_ttagsel_gen->fill(event);
  if(debug) cout << "Filled hists after ttagsel" << endl;
  **/

  return true;
}

// as we want to run the ExampleCycleNew directly with AnalysisModuleRunner,
// make sure the TstarTstarSelectionModule is found by class name. This is ensured by this macro:
UHH2_REGISTER_ANALYSIS_MODULE(TstarTstarSelectionModule)
