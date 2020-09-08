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
  unique_ptr<AnalysisModule> sf_muon_ID, sf_muon_trig_lowpt, sf_muon_trig_highpt, sf_muon_iso, sf_ele_ID, sf_ele_trig, sf_ele_iso, sf_toptag;
  unique_ptr<uhh2::AnalysisModule> reco_primlep;
  uhh2::Event::Handle<FlavorParticle> h_primlep;

  // Store the Hists collection as member variables. Again, use unique_ptr to avoid memory leaks.
  std::unique_ptr<Hists> h_2Dcut, h_trigger;
  std::unique_ptr<Hists> h_2Dcut_gen, h_trigger_gen;
  std::unique_ptr<Hists> h_2Dcut_ele, h_trigger_ele;
  std::unique_ptr<Hists> h_2Dcut_ele_lowpt, h_trigger_ele_lowpt;
  std::unique_ptr<Hists> h_2Dcut_ele_highpt, h_trigger_ele_highpt;
  std::unique_ptr<Hists> h_2Dcut_mu, h_trigger_mu;
  std::unique_ptr<Hists> h_2Dcut_mu_lowpt, h_trigger_mu_lowpt;
  std::unique_ptr<Hists> h_2Dcut_mu_highpt, h_trigger_mu_highpt;
  std::unique_ptr<Hists> h_triggerSingleJet, h_triggerSingleLeptonMu, h_triggerSingleLeptonEle, h_triggerHT, h_triggerPFHT;
  std::unique_ptr<Hists> h_triggerSingleJet_mu, h_triggerHT_mu, h_triggerPFHT_mu;
  std::unique_ptr<Hists> h_triggerSingleJet_ele, h_triggerHT_ele, h_triggerPFHT_ele;
  std::unique_ptr<Hists> h_btagcut, h_btagcut_mu, h_btagcut_ele, h_btagcut_mu_lowpt, h_btagcut_ele_lowpt, h_btagcut_mu_highpt, h_btagcut_ele_highpt;
  std::unique_ptr<Hists> h_beginSel, h_beginSel_mu, h_beginSel_mu_lowpt, h_beginSel_mu_highpt, h_beginSel_ele, h_beginSel_ele_lowpt, h_beginSel_ele_highpt;
  std::unique_ptr<Hists> h_beginSel_gen, h_btagcut_gen;

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
  unique_ptr<Selection> triggerHT1_sel, triggerHT2_sel, triggerHT3_sel, triggerHT4_sel, triggerHT5_sel,  triggerHT6_sel;
  unique_ptr<Selection> triggerPFHT_sel;

  unique_ptr<JetCleaner> AK4cleaner;
  unique_ptr<TopJetCleaner> AK8cleaner;

  // handles
  std::unique_ptr<uhh2::AnalysisModule> ttgenprod;
  uhh2::Event::Handle<TTbarGen> h_ttbargen;
  uhh2::Event::Handle<bool> h_is_muevt;

  uhh2::Event::Handle<double> h_evt_weight;
  uhh2::Event::Handle<double> h_ST;

  bool debug = false;
  bool isTrigger = false;

  // bools for channel and stuff. will be read in later
  bool is_MC;
  bool data_isMu = false;
  bool data_is2017B = false;

  TString year;

  TopJetId topjetID;

};

TstarTstarSelectionModule::TstarTstarSelectionModule(Context & ctx){

  // 0. Reading in whether MC and if so, which channel
  is_MC = ctx.get("dataset_type") == "MC";
  TString year = ctx.get("year", "<not set>");
  if(debug) cout << "Year is " << year << "." << endl;

  if(!is_MC) data_isMu = (ctx.get("dataset_version").find("SingleMuon") != std::string::npos);
  if(!is_MC) data_is2017B = (ctx.get("dataset_version").find("SingleElectron2017_RunB") != std::string::npos);

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
  if(is_MC || !data_isMu) {
    // The following exist for both 2016 and 2017
    // until 120 GeV, except for 2017B, there for whole range
    triggerSingleLeptonEle1_sel.reset(new TriggerSelection("HLT_Ele27_WPTight_Gsf_v*"));
    // above 120 GeV
    triggerSingleLeptonEle2_sel.reset(new TriggerSelection("HLT_Photon175_v*"));
    if(!data_is2017B) triggerSingleLeptonEle3_sel.reset(new TriggerSelection("HLT_Ele115_CaloIdVT_GsfTrkIdT_v*"));
  }
  if(is_MC || data_isMu){
    // until 27 GeV
    triggerSingleLeptonMu1_sel.reset(new TriggerSelection("HLT_IsoMu24_v*"));
    if(year == "2016") triggerSingleLeptonMu2_sel.reset(new TriggerSelection("HLT_IsoTkMu24_v*"));
    // above 60 GeV
    triggerSingleLeptonMu3_sel.reset(new TriggerSelection("HLT_Mu50_v*"));
  }

  // 4. Set up Hists
  if(debug) cout << "Setting up Hists." << endl;
  h_2Dcut.reset(new TstarTstarHists(ctx, "After2D"));
  h_2Dcut_ele.reset(new TstarTstarHists(ctx, "After2D_ele"));
  h_2Dcut_ele_lowpt.reset(new TstarTstarHists(ctx, "After2D_ele_lowpt"));
  h_2Dcut_ele_highpt.reset(new TstarTstarHists(ctx, "After2D_ele_highpt"));
  h_2Dcut_mu.reset(new TstarTstarHists(ctx, "After2D_mu"));
  h_2Dcut_mu_lowpt.reset(new TstarTstarHists(ctx, "After2D_mu_lowpt"));
  h_2Dcut_mu_highpt.reset(new TstarTstarHists(ctx, "After2D_mu_highpt"));
  h_btagcut.reset(new TstarTstarHists(ctx, "AfterBtag"));
  h_btagcut_ele.reset(new TstarTstarHists(ctx, "AfterBtag_ele"));
  h_btagcut_ele_lowpt.reset(new TstarTstarHists(ctx, "AfterBtag_ele_lowpt"));
  h_btagcut_ele_highpt.reset(new TstarTstarHists(ctx, "AfterBtag_ele_highpt"));
  h_btagcut_mu.reset(new TstarTstarHists(ctx, "AfterBtag_mu"));
  h_btagcut_mu_lowpt.reset(new TstarTstarHists(ctx, "AfterBtag_mu_lowpt"));
  h_btagcut_mu_highpt.reset(new TstarTstarHists(ctx, "AfterBtag_mu_highpt"));
  h_trigger.reset(new TstarTstarHists(ctx, "AfterTrigger"));
  h_trigger_mu.reset(new TstarTstarHists(ctx, "AfterTrigger_mu"));
  h_trigger_mu_lowpt.reset(new TstarTstarHists(ctx, "AfterTrigger_mu_lowpt"));
  h_trigger_mu_highpt.reset(new TstarTstarHists(ctx, "AfterTrigger_mu_highpt"));
  h_trigger_ele.reset(new TstarTstarHists(ctx, "AfterTrigger_ele"));
  h_trigger_ele_lowpt.reset(new TstarTstarHists(ctx, "AfterTrigger_ele_lowpt"));
  h_trigger_ele_highpt.reset(new TstarTstarHists(ctx, "AfterTrigger_ele_highpt"));

  h_beginSel.reset(new TstarTstarHists(ctx, "beginSel"));
  h_beginSel_mu.reset(new TstarTstarHists(ctx, "beginSel_mu"));
  h_beginSel_mu_lowpt.reset(new TstarTstarHists(ctx, "beginSel_mu_lowpt"));
  h_beginSel_mu_highpt.reset(new TstarTstarHists(ctx, "beginSel_mu_highpt"));
  h_beginSel_ele.reset(new TstarTstarHists(ctx, "beginSel_ele"));
  h_beginSel_ele_lowpt.reset(new TstarTstarHists(ctx, "beginSel_ele_lowpt"));
  h_beginSel_ele_highpt.reset(new TstarTstarHists(ctx, "beginSel_ele_highpt"));
  //h_ttagsel.reset(new TstarTstarHists(ctx, "AfterTtagsel"));

  h_beginSel_gen.reset(new TstarTstarGenHists(ctx, "beginSel_gen"));
  h_2Dcut_gen.reset(new TstarTstarGenHists(ctx, "After2D_gen"));
  h_btagcut_gen.reset(new TstarTstarGenHists(ctx, "AfterBtag_gen"));
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
  h_evt_weight = ctx.get_handle<double>("evt_weight");

  h_ST = ctx.declare_event_output<double>("ST");

  h_primlep = ctx.get_handle<FlavorParticle>("PrimaryLepton");
  reco_primlep.reset(new PrimaryLepton(ctx));

  // Scale factor stuff
  // Muon
  //sf_muon_trig_lowpt.reset(new MCMuonScaleFactor(ctx, "/nfs/dust/cms/user/flabe/CMSSW/CMSSW_10_2_10/src/UHH2/common/data/2016/MuonTrigger_EfficienciesAndSF_average_RunBtoH.root", "IsoMu24_OR_IsoTkMu24_PtEtaBins", 0.5, "muon_trigger_lowpt", true, "nominal"));
  //sf_muon_trig_highpt.reset(new MCMuonScaleFactor(ctx, "/nfs/dust/cms/user/flabe/CMSSW/CMSSW_10_2_10/src/UHH2/common/data/2016/MuonTrigger_EfficienciesAndSF_average_RunBtoH.root", "IsoMu50_OR_IsoTkMu_50_PtEtaBins", 0.5, "muon_trigger_highpt", true, "nominal"));
  //sf_muon_ID.reset(new MCMuonScaleFactor(ctx, "/nfs/dust/cms/user/flabe/CMSSW/CMSSW_10_2_10/src/UHH2/common/data/2016/MuonID_EfficienciesAndSF_average_RunBtoH.root", "NUM_TightID_DEN_genTracks_eta_pt", 1, "muon_tightID", true, "nominal"));
  //sf_muon_iso.reset(new MCMuonScaleFactor(ctx, "/nfs/dust/cms/user/flabe/CMSSW/CMSSW_10_2_10/src/UHH2/common/data/2016/MuonIso_EfficienciesAndSF_average_RunBtoH.root", "NUM_TightRelIso_DEN_TightIDandIPCut_eta_pt", 1, "muon_iso", true, "nominal"));

  // Electron
  //sf_ele_trig.reset();
  //sf_ele_ID.reset();
  //sf_ele_iso.reset();
}


bool TstarTstarSelectionModule::process(Event & event) {

  if(debug) cout << "TstarTstarSelectionModule: Starting to process event (runid, eventid) = (" << event.run << ", " << event.event << "); weight = " << event.weight << endl;

  //LumiWeight_module->process(event); // apply correct weights
  event.weight = event.get(h_evt_weight);
  if(debug) cout << "weights applied." << endl;

  //Fill ttgen object for correct matching check, etc
  if(is_MC){
    ttgenprod->process(event);
    if(debug) cout << "ttgen produced." << endl;
  }

  reco_primlep->process(event);//set "primary lepton"

  if(debug) std::cout << "Fill Crosscheck hists" << endl;
  h_beginSel->fill(event);
  h_beginSel_gen->fill(event);
  if(event.get(h_is_muevt)){
    h_beginSel_mu->fill(event);
    if(event.muons->at(0).pt()<=60) h_beginSel_mu_lowpt->fill(event);
    else h_beginSel_mu_highpt->fill(event);
  }
  else {
    h_beginSel_ele->fill(event);
    if(event.electrons->at(0).pt()<=120) h_beginSel_ele_lowpt->fill(event);
    else h_beginSel_ele_highpt->fill(event);
  }

  // Btag Selection
  bool pass_btagcut = false;
  for (const auto & jet: *event.jets){
    if(jet.btag_DeepCSV() > 0.2219) pass_btagcut = true;
  }
  if(!pass_btagcut) return false;
  h_btagcut->fill(event);
  h_btagcut_gen->fill(event);
  if(event.get(h_is_muevt)){
    h_btagcut_mu->fill(event);
    if(event.muons->at(0).pt()<=60) h_btagcut_mu_lowpt->fill(event);
    else h_btagcut_mu_highpt->fill(event);
  }
  else {
    h_btagcut_ele->fill(event);
    if(event.electrons->at(0).pt()<=120) h_btagcut_ele_lowpt->fill(event);
    else h_btagcut_ele_highpt->fill(event);
  }
  if(debug) cout << "Passed btag selection cut." << endl;

  // Lepton-2Dcut
  if(debug) std::cout << "Start 2D cut" << endl;
  bool pass_2D = true;
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
  if(event.muons->size()==1){if(event.muons->at(0).pt()>60) pass_2D = pass_twodcut;}
  else if(event.electrons->size()==1){if(event.electrons->at(0).pt()>120) pass_2D = pass_twodcut;}
  else std::cout << "How did this happen???" << endl;
  if(!pass_2D) return false;
  h_2Dcut->fill(event);
  h_2Dcut_gen->fill(event);
  if(event.get(h_is_muevt)){
    h_2Dcut_mu->fill(event);
    if(event.muons->at(0).pt()<=60) h_2Dcut_mu_lowpt->fill(event);
    else h_2Dcut_mu_highpt->fill(event);
  }
  else {
    h_2Dcut_ele->fill(event);
    if(event.electrons->at(0).pt()<=120) h_2Dcut_ele_lowpt->fill(event);
    else h_2Dcut_ele_highpt->fill(event);
  }
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
  // if(isTrigger){
  //   if(debug) cout << "is Trigger" << endl;
  //   // SingleJet
  //   bool pass_trigger_SingleJet = (triggerSingleJet450_sel->passes(event) && event.jets->at(0).pt()>450);
  //   if(pass_trigger_SingleJet){
  //     h_triggerSingleJet->fill(event);
  //     if(event.get(h_is_muevt)) h_triggerSingleJet_mu->fill(event);
  //     else h_triggerSingleJet_ele->fill(event);
  //   }
  //   if(debug) cout << "done SingleJet" << endl;
  //   bool pass_trigegr_HT = triggerHT1_sel->passes(event) || triggerHT2_sel->passes(event) || triggerHT3_sel->passes(event)
  //     || triggerHT4_sel->passes(event) || triggerHT5_sel->passes(event) || triggerHT6_sel->passes(event);
  //   if(pass_trigegr_HT){
  //     h_triggerHT->fill(event);
  //     if(event.get(h_is_muevt)) h_triggerHT_mu->fill(event);
  //     else h_triggerHT_ele->fill(event);
  //   }
  //   if(debug) cout << "done HT" << endl;
  //   bool pass_trigegr_PFHT = triggerPFHT_sel->passes(event);
  //   if(pass_trigegr_PFHT){
  //     h_triggerPFHT->fill(event);
  //     if(event.get(h_is_muevt)) h_triggerPFHT_mu->fill(event);
  //     else h_triggerPFHT_ele->fill(event);
  //   }
  //   if(debug) cout << "done PFHT" << endl;
  // }

  // deltaR cut to suppress QCD
  FlavorParticle primary_lepton = event.get(h_primlep);
  double min_deltaR = 999;
  for(auto &jet : *event.jets){
    double cur_deltaR = deltaR(jet, primary_lepton);
    if(cur_deltaR < min_deltaR) min_deltaR = cur_deltaR;
  }
  bool pass_dR = false;
  if(event.get(h_is_muevt)){
    if(event.muons->at(0).pt()<=60) pass_dR = min_deltaR > 0.4;
    else pass_dR = min_deltaR > 0.2;
  }
  else {
    if(event.electrons->at(0).pt()<=120) pass_dR = min_deltaR > 0.4;
    else pass_dR = min_deltaR > 0.2;
  }
  if(!pass_dR) return false;
  if(debug) cout << "Passed dR cut." << endl;

  // Trigger
  bool pass_trigger = false;
  bool pass_trigger_SingleMu_lowpt = false;
  bool pass_trigger_SingleMu_highpt = false;
  if(is_MC || data_isMu){
    if(year == "2016") {
      pass_trigger_SingleMu_lowpt = (triggerSingleLeptonMu1_sel->passes(event) || triggerSingleLeptonMu2_sel->passes(event));
    }
    else pass_trigger_SingleMu_lowpt = triggerSingleLeptonMu1_sel->passes(event);
    pass_trigger_SingleMu_highpt = triggerSingleLeptonMu3_sel->passes(event);
  }
  bool pass_trigger_SingleEle_lowpt = false;
  bool pass_trigger_SingleEle_highpt = false;
  if(is_MC || !data_isMu){
    pass_trigger_SingleEle_lowpt = triggerSingleLeptonEle1_sel->passes(event);
    if(data_is2017B) pass_trigger_SingleEle_highpt = (triggerSingleLeptonEle2_sel->passes(event) || triggerSingleLeptonEle1_sel->passes(event));
    else pass_trigger_SingleEle_highpt = (triggerSingleLeptonEle2_sel->passes(event) || triggerSingleLeptonEle3_sel->passes(event));
  }
  if(pass_trigger_SingleMu_lowpt && (event.muons->size() == 1)){if(event.muons->at(0).pt()<=60) pass_trigger = true; }
  if(pass_trigger_SingleMu_highpt && (event.muons->size() == 1)){if(event.muons->at(0).pt()>60) pass_trigger = true; }
  if(pass_trigger_SingleEle_lowpt && (event.electrons->size() == 1)){if(event.electrons->at(0).pt()<=120)pass_trigger = true; }
  if(pass_trigger_SingleEle_highpt && (event.electrons->size() == 1)){if(event.electrons->at(0).pt()>120)pass_trigger = true; }
  if(!pass_trigger) return false;
  h_trigger->fill(event);
  h_trigger_gen->fill(event);
  if(event.get(h_is_muevt)){
    h_trigger_mu->fill(event);
    if(event.muons->at(0).pt()<=60) h_trigger_mu_lowpt->fill(event);
    else h_trigger_mu_highpt->fill(event);
    h_triggerSingleLeptonMu->fill(event);
  }
  else {
    h_trigger_ele->fill(event);
    if(event.electrons->at(0).pt()<=120) h_trigger_ele_lowpt->fill(event);
    else h_trigger_ele_highpt->fill(event);
    h_triggerSingleLeptonEle->fill(event);
  }
  if(debug) cout<<"Filled hists after Trigger"<<endl;

  // Scale factors.
  // A lot of stuff will happen here. hopefully.


  // Outputting ST
  double st_jets = 0.;
  for(const auto & jet : *event.topjets) st_jets += jet.pt();
  event.set(h_ST, st_jets);

  return true;
}

// as we want to run the ExampleCycleNew directly with AnalysisModuleRunner,
// make sure the TstarTstarSelectionModule is found by class name. This is ensured by this macro:
UHH2_REGISTER_ANALYSIS_MODULE(TstarTstarSelectionModule)
