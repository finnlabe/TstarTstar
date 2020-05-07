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
  std::unique_ptr<Hists> h_2Dcut, h_trigger;
  std::unique_ptr<Hists> h_2Dcut_gen, h_trigger_gen;
  std::unique_ptr<Hists> h_2Dcut_ele, h_trigger_ele;
  std::unique_ptr<Hists> h_2Dcut_mu, h_trigger_mu;
  std::unique_ptr<Hists> h_triggerSingleJet, h_triggerSingleLeptonMu, h_triggerSingleLeptonEle, h_triggerHT, h_triggerPFHT;
  std::unique_ptr<Hists> h_triggerSingleJet_mu, h_triggerHT_mu, h_triggerPFHT_mu;
  std::unique_ptr<Hists> h_triggerSingleJet_ele, h_triggerHT_ele, h_triggerPFHT_ele;

  unique_ptr<Selection> met_sel;
  unique_ptr<Selection> st_sel;
  unique_ptr<Selection> twodcut_sel;
  unique_ptr<Selection> toptagevt_sel;

  unique_ptr<Selection> triggerSingleJet450_sel;
  unique_ptr<Selection> triggerSingleLeptonEle1_sel;
  unique_ptr<Selection> triggerSingleLeptonEle2_sel;
  unique_ptr<Selection> triggerSingleLeptonEle3_sel;
  unique_ptr<Selection> triggerSingleLeptonMu1_sel;
  unique_ptr<Selection> triggerSingleLeptonMu2_sel;
  unique_ptr<Selection> triggerSingleLeptonMu3_sel;
  unique_ptr<Selection> triggerSingleLeptonMu4_sel;
  unique_ptr<Selection> triggerHT1_sel, triggerHT2_sel, triggerHT3_sel, triggerHT4_sel, triggerHT5_sel,  triggerHT6_sel;
  unique_ptr<Selection> triggerPFHT_sel;

  unique_ptr<JetCleaner> AK4cleaner;
  unique_ptr<TopJetCleaner> AK8cleaner;

  // GEN stuff
  std::unique_ptr<uhh2::AnalysisModule> ttgenprod;
  uhh2::Event::Handle<TTbarGen> h_ttbargen;
  uhh2::Event::Handle<bool> h_is_muevt;

  bool debug = false;
  bool isTrigger = true;

  // bools for channel and stuff. will be read in later
  bool is_MC;

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

  //trigger studies
  triggerSingleJet450_sel.reset(new TriggerSelection("HLT_PFJet450_v*"));

  // TODO check 2017 2018 threshholds
  // until 120 GeV
  triggerSingleLeptonEle1_sel.reset(new TriggerSelection("HLT_Ele27_WPTight_Gsf_v*"));
  // above 120 GeV
  triggerSingleLeptonEle2_sel.reset(new TriggerSelection("HLT_Photon175_v*"));
  triggerSingleLeptonEle3_sel.reset(new TriggerSelection("HLT_Ele115_CaloIdVT_GsfTrkIdT_v*"));

  // until 27 GeV
  triggerSingleLeptonMu1_sel.reset(new TriggerSelection("HLT_IsoMu24_v*"));
  triggerSingleLeptonMu2_sel.reset(new TriggerSelection("HLT_IsoTkMu24_v*"));
  // above 60 GeV
  triggerSingleLeptonMu3_sel.reset(new TriggerSelection("HLT_Mu50_v*"));
  triggerSingleLeptonMu4_sel.reset(new TriggerSelection("HLT_Mu55_v*"));

  triggerHT1_sel.reset(new TriggerSelection("HLT_HT430to450_v*"));
  triggerHT2_sel.reset(new TriggerSelection("HLT_HT450to470_v*"));
  triggerHT3_sel.reset(new TriggerSelection("HLT_HT470to500_v*"));
  triggerHT4_sel.reset(new TriggerSelection("HLT_HT500to550_v*"));
  triggerHT5_sel.reset(new TriggerSelection("HLT_HT550to650_v*"));
  triggerHT6_sel.reset(new TriggerSelection("HLT_HT650_v*"));

  triggerPFHT_sel.reset(new TriggerSelection("HLT_PFHT900_v*"));

  // 4. Set up Hists
  if(debug) cout << "Setting up Hists." << endl;
  h_2Dcut.reset(new TstarTstarHists(ctx, "After2D"));
  h_2Dcut_ele.reset(new TstarTstarHists(ctx, "After2D_ele"));
  h_2Dcut_mu.reset(new TstarTstarHists(ctx, "After2D_mu"));
  h_trigger.reset(new TstarTstarHists(ctx, "AfterTrigger"));
  h_trigger_mu.reset(new TstarTstarHists(ctx, "AfterTrigger_mu"));
  h_trigger_ele.reset(new TstarTstarHists(ctx, "AfterTrigger_ele"));
  //h_ttagsel.reset(new TstarTstarHists(ctx, "AfterTtagsel"));

  h_2Dcut_gen.reset(new TstarTstarGenHists(ctx, "After2D_gen"));
  h_trigger_gen.reset(new TstarTstarGenHists(ctx, "AfterTrigger_gen"));
  //h_ttagsel_gen.reset(new TstarTstarGenHists(ctx, "AfterTtagsel_gen"));

  h_triggerSingleJet.reset(new TstarTstarHists(ctx, "triggerSingleJet"));
  h_triggerSingleLeptonMu.reset(new TstarTstarHists(ctx, "triggerSingleLeptonMu"));
  h_triggerSingleLeptonEle.reset(new TstarTstarHists(ctx, "triggerSingleLeptonEle"));
  h_triggerHT.reset(new TstarTstarHists(ctx, "triggerHT"));
  h_triggerPFHT.reset(new TstarTstarHists(ctx, "triggerPFHT"));

  h_triggerSingleJet_mu.reset(new TstarTstarHists(ctx, "triggerSingleJet_mu"));
  h_triggerHT_mu.reset(new TstarTstarHists(ctx, "triggerHT_mu"));
  h_triggerPFHT_mu.reset(new TstarTstarHists(ctx, "triggerPFHT_mu"));

  h_triggerSingleJet_ele.reset(new TstarTstarHists(ctx, "triggerSingleJet_ele"));
  h_triggerHT_ele.reset(new TstarTstarHists(ctx, "triggerHT_ele"));
  h_triggerPFHT_ele.reset(new TstarTstarHists(ctx, "triggerPFHT_ele"));

  h_is_muevt = ctx.get_handle<bool>("is_muevt");

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
  if(event.get(h_is_muevt)) h_2Dcut_mu->fill(event);
  else h_2Dcut_ele->fill(event);
  if(debug) cout << "Passed 2D cut." << endl;

  // TopTagEventSelection
  /**
  bool pass_ttag = (toptagevt_sel->passes(event));
  if(!pass_ttag) return false;
  h_ttagsel->fill(event);
  h_ttagsel_gen->fill(event);
  if(debug) cout << "Filled hists after ttagsel" << endl;
  **/

  // Trigger studies
  if(isTrigger){
    if(debug) cout << "is Trigger" << endl;
    // SingleJet
    bool pass_trigger_SingleJet = (triggerSingleJet450_sel->passes(event) && event.jets->at(0).pt()>450);
    if(pass_trigger_SingleJet){
      h_triggerSingleJet->fill(event);
      if(event.get(h_is_muevt)) h_triggerSingleJet_mu->fill(event);
      else h_triggerSingleJet_ele->fill(event);
    }
    if(debug) cout << "done SingleJet" << endl;
    bool pass_trigger_SingleMu = (triggerSingleLeptonMu1_sel->passes(event) || triggerSingleLeptonMu2_sel->passes(event)
				  || triggerSingleLeptonMu3_sel->passes(event) || triggerSingleLeptonMu4_sel->passes(event));
    if(pass_trigger_SingleMu && (event.muons->size() == 1)){
      h_triggerSingleLeptonMu->fill(event);
    }
    if(debug) cout << "done SingleMu" << endl;
    bool pass_trigger_SingleEle = (triggerSingleLeptonEle1_sel->passes(event) || triggerSingleLeptonEle2_sel->passes(event) || triggerSingleLeptonEle3_sel->passes(event));
    if(pass_trigger_SingleEle && (event.electrons->size() == 1)){
      h_triggerSingleLeptonEle->fill(event);
    }
    if(debug) cout << "done SingleEle" << endl;
    bool pass_trigegr_HT = triggerHT1_sel->passes(event) || triggerHT2_sel->passes(event) || triggerHT3_sel->passes(event)
      || triggerHT4_sel->passes(event) || triggerHT5_sel->passes(event) || triggerHT6_sel->passes(event);
    if(pass_trigegr_HT){
      h_triggerHT->fill(event);
      if(event.get(h_is_muevt)) h_triggerHT_mu->fill(event);
      else h_triggerHT_ele->fill(event);
    }
    if(debug) cout << "done HT" << endl;
    bool pass_trigegr_PFHT = triggerPFHT_sel->passes(event);
    if(pass_trigegr_PFHT){
      h_triggerPFHT->fill(event);
      if(event.get(h_is_muevt)) h_triggerPFHT_mu->fill(event);
      else h_triggerPFHT_ele->fill(event);
    }
    if(debug) cout << "done PFHT" << endl;
  }

  // Trigger
  bool pass_trigger = false;
  bool pass_trigger_SingleMu = (triggerSingleLeptonMu1_sel->passes(event) || triggerSingleLeptonMu2_sel->passes(event) || triggerSingleLeptonMu3_sel->passes(event) || triggerSingleLeptonMu4_sel->passes(event));
  if(pass_trigger_SingleMu && (event.muons->size() == 1)){ pass_trigger = true; }
  bool pass_trigger_SingleEle = (triggerSingleLeptonEle1_sel->passes(event) || triggerSingleLeptonEle2_sel->passes(event) || triggerSingleLeptonEle3_sel->passes(event));
  if(pass_trigger_SingleEle && (event.electrons->size() == 1)){ pass_trigger = true; }
  if(!pass_trigger) return false;
  h_trigger->fill(event);
  h_trigger_gen->fill(event);
  if(event.get(h_is_muevt)) h_trigger_mu->fill(event);
  else h_trigger_ele->fill(event);
  if(debug) cout<<"Filled hists after Trigger"<<endl;

  return true;
}

// as we want to run the ExampleCycleNew directly with AnalysisModuleRunner,
// make sure the TstarTstarSelectionModule is found by class name. This is ensured by this macro:
UHH2_REGISTER_ANALYSIS_MODULE(TstarTstarSelectionModule)
