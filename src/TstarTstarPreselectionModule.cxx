#include <iostream>
#include <memory>

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
#include "UHH2/common/include/MCWeight.h"
#include <UHH2/common/include/TriggerSelection.h>

#include "UHH2/TstarTstar/include/TstarTstarCustomIds.h"
#include "UHH2/TstarTstar/include/TstarTstarSelections.h"
#include "UHH2/TstarTstar/include/TstarTstarHists.h"
#include "UHH2/TstarTstar/include/TstarTstarGenHists.h"
#include "UHH2/TstarTstar/include/TstarTstarGenRecoMatchedHists.h"

#include "UHH2/HOTVR/include/HOTVRJetCorrectionModule.h"
#include "UHH2/HOTVR/include/HadronicTop.h"
#include "UHH2/HOTVR/include/HOTVRScaleFactor.h"


using namespace std;
using namespace uhh2;

namespace uhh2 {

/** \brief Module for the T*T*->ttbar gg/gamma preselection
 *
 * Corrects all objects via CommonModules and applies some loose cuts
 *
 */
class TstarTstarPreselectionModule: public AnalysisModule {
public:

  explicit TstarTstarPreselectionModule(Context & ctx);
  virtual bool process(Event & event) override;

private:

  // Apply common modules: JetPFid, JEC, JER, MET corrections, etc
  std::unique_ptr<CommonModules> common;
  std::unique_ptr<AnalysisModule> HOTVRCorr;

  std::unique_ptr<AnalysisModule> HadronicTopFinder;
  std::unique_ptr<AnalysisModule> HOTVRScale;

  // Declare the Selections to use. Use unique_ptr to ensure automatic call of delete in the destructor, to avoid memory leaks.
  unique_ptr<TopJetCleaner> HOTVRcleaner;
  unique_ptr<MuonCleaner> MuCleaner_lowpt;
  unique_ptr<MuonCleaner> MuCleaner_highpt;
  unique_ptr<ElectronCleaner> EleCleaner_lowpt;
  unique_ptr<ElectronCleaner> EleCleaner_highpt;
  unique_ptr<Selection> met_sel;

  unique_ptr<Selection> triggerSingleJet450_sel;
  unique_ptr<Selection> triggerSingleLeptonEle1_sel;
  unique_ptr<Selection> triggerSingleLeptonEle2_sel;
  unique_ptr<Selection> triggerSingleLeptonEle3_sel;
  unique_ptr<Selection> triggerSingleLeptonMu1_sel;
  unique_ptr<Selection> triggerSingleLeptonMu2_sel;
  unique_ptr<Selection> triggerSingleLeptonMu3_sel;
  unique_ptr<Selection> triggerHT1_sel, triggerHT2_sel, triggerHT3_sel, triggerHT4_sel, triggerHT5_sel,  triggerHT6_sel;
  unique_ptr<Selection> triggerPFHT_sel;

  // ##### Histograms #####
  // Store the Hists collection as member variables. Again, use unique_ptr to avoid memory leaks.
  std::unique_ptr<Hists> h_nocuts,     h_common,     h_lepsel,     h_fatjetsel,     h_METsel;
  std::unique_ptr<Hists> h_nocuts_gen, h_common_gen, h_lepsel_gen, h_fatjetsel_gen, h_METsel_gen;
  std::unique_ptr<Hists> h_lepsel_ele,     h_fatjetsel_ele,     h_METsel_ele;
  std::unique_ptr<Hists> h_lepsel_ele_lowpt,     h_fatjetsel_ele_lowpt,     h_METsel_ele_lowpt;
  std::unique_ptr<Hists> h_lepsel_ele_highpt,     h_fatjetsel_ele_highpt,     h_METsel_ele_highpt;
  std::unique_ptr<Hists> h_lepsel_mu,     h_fatjetsel_mu,     h_METsel_mu;
  std::unique_ptr<Hists> h_lepsel_mu_lowpt,     h_fatjetsel_mu_lowpt,     h_METsel_mu_lowpt;
  std::unique_ptr<Hists> h_lepsel_mu_highpt,     h_fatjetsel_mu_highpt,     h_METsel_mu_highpt;
  std::unique_ptr<Hists> h_triggerSingleLeptonMu, h_triggerSingleLeptonEle;
  std::unique_ptr<Hists> h_trigger;
  std::unique_ptr<Hists> h_trigger_gen;
  std::unique_ptr<Hists> h_trigger_ele;
  std::unique_ptr<Hists> h_trigger_ele_lowpt;
  std::unique_ptr<Hists> h_trigger_ele_highpt;
  std::unique_ptr<Hists> h_trigger_mu;
  std::unique_ptr<Hists> h_trigger_mu_lowpt;
  std::unique_ptr<Hists> h_trigger_mu_highpt;
  std::unique_ptr<uhh2::AnalysisModule> ttgenprod;
  uhh2::Event::Handle<TTbarGen> h_ttbargen;
  uhh2::Event::Handle<bool> h_is_muevt;

  uhh2::Event::Handle<double> h_evt_weight;

  std::unique_ptr<LuminosityHists> lumihist_beforeSel, lumihist_afterTrigger, lumihist_afterLepSel, lumihist_afterNJets, lumihist_afterMET;


  bool debug = false;
  bool doTriggerSel = true;

  // bools for channel and stuff. will be read in later
  bool is_MC;
  bool data_isMu = false;
  bool data_is2017B = false;

  TString year;

};


TstarTstarPreselectionModule::TstarTstarPreselectionModule(Context & ctx){


  if(debug) {
    cout << "Hello World from TstarTstarPreselectionModule!" << endl;

    // If running in SFrame, the keys "dataset_version", "dataset_type", "dataset_lumi",
    // and "target_lumi" are set to the according values in the xml file. For CMSSW, these are
    // not set automatically, but can be set in the python config file.
    for(auto & kv : ctx.get_all()){
      cout << " " << kv.first << " = " << kv.second << endl;
    }

  }

  year = ctx.get("year", "<not set>");
  if(year == "<not set>"){
    if(ctx.get("dataset_version").find("2016") != std::string::npos) year = "2016";
    else if(ctx.get("dataset_version").find("2017") != std::string::npos) year = "2017";
    else if(ctx.get("dataset_version").find("2018") != std::string::npos) year = "2018";
    else throw "No year found in dataset name!";
  }
  if(true) cout << "Year is " << year << "." << endl;

  // 0. Reading in whether MC and if so, which channel
  is_MC = ctx.get("dataset_type") == "MC";

  if(!is_MC) data_isMu = (ctx.get("dataset_version").find("SingleMuon") != std::string::npos);
  if(!is_MC) data_is2017B = (ctx.get("dataset_version").find("SingleElectron2017_RunB") != std::string::npos);

  if(debug) cout << "trying..." << endl;
  ctx.get("HOTVRTopTagSFs");
  if(debug) cout << "done." << endl;

  // 1. setup modules. CommonModules
  common.reset(new CommonModules());
  common->switch_metcorrection();

  HOTVRCorr.reset(new HOTVRJetCorrectionModule(ctx));
  HadronicTopFinder.reset(new HadronicTop(ctx));
  TopJetId topjetID = AndId<TopJet>(HOTVRTopTag(), Tau32Groomed(0.56)); // Top Tag that is used later
  HOTVRScale.reset(new HOTVRScaleFactor(ctx, topjetID));

  // Electron
  ElectronId eleID_lowpt = ElectronID_Summer16_tight;
  ElectronId eleID_highpt = ElectronID_Summer16_tight_noIso;
  double electron_pt_lowpt(30.);
  double electron_pt_highpt(120.);
  common->set_electron_id(OrId<Electron>( AndId<Electron>(PtEtaSCCut(electron_pt_lowpt, 2.4), eleID_lowpt, EleMaxPtCut(120.)),  AndId<Electron>(PtEtaSCCut(electron_pt_highpt, 2.4), eleID_highpt)));

  //Muon
  MuonId muID_lowpt = AndId<Muon>(MuonID(Muon::CutBasedIdTight), MuonID(Muon::PFIsoTight));
  MuonId muID_highpt = MuonID(Muon::CutBasedIdTight);
  double muon_pt_lowpt(27.);
  double muon_pt_highpt(60.);
  common->set_muon_id(OrId<Muon>( AndId<Muon>(PtEtaCut(muon_pt_lowpt, 2.4), muID_lowpt, MuMaxPtCut(60.)),  AndId<Muon>(PtEtaCut(muon_pt_highpt, 2.4), muID_highpt)));

  // Jets
  common->switch_jetlepcleaner();
  common->switch_jetPtSorter();
  double jet_pt(30.);
  common->set_jet_id(AndId<Jet>(PtEtaCut(jet_pt, 2.5), JetPFID(JetPFID::WP_TIGHT_PUPPI)));
  HOTVRcleaner.reset(new TopJetCleaner(ctx, PtEtaCut(150.0, 2.5)));

  // init common.
  common->init(ctx);

  if(is_MC){
    // Prepare GEN
    ttgenprod.reset(new TTbarGenProducer(ctx, "ttbargen", false));
  }

  // 2. set up selections

  // MET selection
  met_sel.reset(new METCut  (50.,1e9));

  // trigger defs
  if(is_MC || !data_isMu) {
    // The following exist for both 2016 and 2017
    // until 120 GeV, except for 2017B, there for whole range
    if(year == "2018") triggerSingleLeptonEle1_sel.reset(new TriggerSelection("HLT_Ele32_WPTight_Gsf_v*"));
    else triggerSingleLeptonEle1_sel.reset(new TriggerSelection("HLT_Ele27_WPTight_Gsf_v*"));
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

  // 3. Set up Hists classes:
  h_nocuts.reset(new TstarTstarHists(ctx, "NoCuts"));
  h_common.reset(new TstarTstarHists(ctx, "AfterCommon"));
  h_lepsel.reset(new TstarTstarHists(ctx, "AfterLepSel"));
  h_fatjetsel.reset(new TstarTstarHists(ctx, "AfterAK8jets"));
  h_METsel.reset(new TstarTstarHists(ctx, "AfterMET"));

  h_lepsel_ele.reset(new TstarTstarHists(ctx, "AfterLepSel_ele"));
  h_fatjetsel_ele.reset(new TstarTstarHists(ctx, "AfterAK8jets_ele"));
  h_METsel_ele.reset(new TstarTstarHists(ctx, "AfterMET_ele"));

  h_lepsel_ele_lowpt.reset(new TstarTstarHists(ctx, "AfterLepSel_ele_lowpt"));
  h_fatjetsel_ele_lowpt.reset(new TstarTstarHists(ctx, "AfterAK8jets_ele_lowpt"));
  h_METsel_ele_lowpt.reset(new TstarTstarHists(ctx, "AfterMET_ele_lowpt"));

  h_lepsel_ele_highpt.reset(new TstarTstarHists(ctx, "AfterLepSel_ele_highpt"));
  h_fatjetsel_ele_highpt.reset(new TstarTstarHists(ctx, "AfterAK8jets_ele_highpt"));
  h_METsel_ele_highpt.reset(new TstarTstarHists(ctx, "AfterMET_ele_highpt"));

  h_lepsel_mu.reset(new TstarTstarHists(ctx, "AfterLepSel_mu"));
  h_fatjetsel_mu.reset(new TstarTstarHists(ctx, "AfterAK8jets_mu"));
  h_METsel_mu.reset(new TstarTstarHists(ctx, "AfterMET_mu"));

  h_lepsel_mu_lowpt.reset(new TstarTstarHists(ctx, "AfterLepSel_mu_lowpt"));
  h_fatjetsel_mu_lowpt.reset(new TstarTstarHists(ctx, "AfterAK8jets_mu_lowpt"));
  h_METsel_mu_lowpt.reset(new TstarTstarHists(ctx, "AfterMET_mu_lowpt"));

  h_lepsel_mu_highpt.reset(new TstarTstarHists(ctx, "AfterLepSel_mu_highpt"));
  h_fatjetsel_mu_highpt.reset(new TstarTstarHists(ctx, "AfterAK8jets_mu_highpt"));
  h_METsel_mu_highpt.reset(new TstarTstarHists(ctx, "AfterMET_mu_highpt"));

  h_trigger.reset(new TstarTstarHists(ctx, "AfterTrigger"));
  h_trigger_mu.reset(new TstarTstarHists(ctx, "AfterTrigger_mu"));
  h_trigger_mu_lowpt.reset(new TstarTstarHists(ctx, "AfterTrigger_mu_lowpt"));
  h_trigger_mu_highpt.reset(new TstarTstarHists(ctx, "AfterTrigger_mu_highpt"));

  h_trigger_ele.reset(new TstarTstarHists(ctx, "AfterTrigger_ele"));
  h_trigger_ele_lowpt.reset(new TstarTstarHists(ctx, "AfterTrigger_ele_lowpt"));
  h_trigger_ele_highpt.reset(new TstarTstarHists(ctx, "AfterTrigger_ele_highpt"));

  h_nocuts_gen.reset(new TstarTstarGenHists(ctx, "NoCuts_gen"));
  h_common_gen.reset(new TstarTstarGenHists(ctx, "AfterCommon_gen"));
  h_trigger_gen.reset(new TstarTstarGenHists(ctx, "AfterTrigger_gen"));
  h_lepsel_gen.reset(new TstarTstarGenHists(ctx, "AfterLepSel_gen"));
  h_fatjetsel_gen.reset(new TstarTstarGenHists(ctx, "AfterAK8jets_gen"));
  h_METsel_gen.reset(new TstarTstarGenHists(ctx, "AfterMET_gen"));

  h_triggerSingleLeptonMu.reset(new TstarTstarHists(ctx, "triggerSingleLeptonMu"));
  h_triggerSingleLeptonEle.reset(new TstarTstarHists(ctx, "triggerSingleLeptonEle"));

  lumihist_beforeSel.reset(new LuminosityHists(ctx, "lumihist_beforeSel"));
  lumihist_afterTrigger.reset(new LuminosityHists(ctx, "lumihist_afterTrigger"));
  lumihist_afterLepSel.reset(new LuminosityHists(ctx, "lumihist_afterLepSel"));
  lumihist_afterNJets.reset(new LuminosityHists(ctx, "lumihist_afterNJets"));
  lumihist_afterMET.reset(new LuminosityHists(ctx, "lumihist_afterMET"));



  h_is_muevt = ctx.declare_event_output<bool>("is_muevt");
  h_evt_weight = ctx.declare_event_output<double>("evt_weight");

}


bool TstarTstarPreselectionModule::process(Event & event) {

  if(debug)
    cout << "TstarTstarPreselectionModule: Starting to process event (runid, eventid) = (" << event.run << ", " << event.event << "); weight = " << event.weight << endl;
  if(debug)
    cout<<"N muons = "<<event.muons->size()<<", N electrons = "<<event.electrons->size()<<", N photons = "<<event.photons->size()<<endl;
  for(auto& muo : *event.muons){
    if(debug) cout<<"BEFORE Muon (pt,eta): "<<muo.pt()<<", "<<muo.eta()<<endl;
  }

  if(is_MC){
    //Fill ttgen object for correct matching check, etc
    ttgenprod->process(event);
  }

  h_nocuts->fill(event);
  h_nocuts_gen->fill(event);
  if(debug) cout<<"Filled hists without any cuts, and GEN with no cuts"<<endl;

  // cleaning & common modules
  if(!(common->process(event))) return false;
  if(!(HOTVRCorr->process(event))) return false;
  if(!(HOTVRcleaner->process(event))) return false;

  // Top Tagging scale factors
  HadronicTopFinder->process(event);
  HOTVRScale->process(event);

  h_common->fill(event);
  h_common_gen->fill(event);
  lumihist_beforeSel->fill(event);
  if(debug) cout<<"Filled hists after cleaning"<<endl;

  //---- Preselection

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
  if(pass_trigger_SingleMu_lowpt && (event.muons->size() >= 1)){if(event.muons->at(0).pt()<=60) pass_trigger = true; }
  if(pass_trigger_SingleMu_highpt && (event.muons->size() >= 1)){if(event.muons->at(0).pt()>60) pass_trigger = true; }
  if(pass_trigger_SingleEle_lowpt && (event.electrons->size() >= 1)){if(event.electrons->at(0).pt()<=120)pass_trigger = true; }
  if(pass_trigger_SingleEle_highpt && (event.electrons->size() >= 1)){if(event.electrons->at(0).pt()>120)pass_trigger = true; }
  if(!pass_trigger && doTriggerSel) return false;
  h_trigger->fill(event);
  h_trigger_gen->fill(event);
  lumihist_afterTrigger->fill(event);

  if(debug) cout<<"Filled hists after Trigger"<<endl;

  // Require exactly one muon or one electron
  const bool pass_lep1 = (((event.muons->size() == 1) || (event.electrons->size() == 1)) && (event.electrons->size()+event.muons->size()) == 1);
  if(!pass_lep1) return false;
  if(event.muons->size() == 1) event.set(h_is_muevt, true);
  else event.set(h_is_muevt, false);
  h_lepsel->fill(event);
  h_lepsel_gen->fill(event);
  lumihist_afterLepSel->fill(event);
  if(event.get(h_is_muevt)){
    h_lepsel_mu->fill(event);
    if(event.muons->at(0).pt()<=60) h_lepsel_mu_lowpt->fill(event);
    else h_lepsel_mu_highpt->fill(event);
  }
  else {
    h_lepsel_ele->fill(event);
    if(event.electrons->at(0).pt()<=120) h_lepsel_ele_lowpt->fill(event);
    else h_lepsel_ele_highpt->fill(event);
  }
  if(debug) cout << "Filled hists after lepsel" << endl;

  // fat jet selection (and jet selection because of crash otherwise)
  bool pass_fat_njet = (event.topjets->size()>0);
  bool pass_njet = (event.jets->size()>3);
  if(!pass_fat_njet) return false;
  if(!pass_njet) return false;
  h_fatjetsel->fill(event);
  h_fatjetsel_gen->fill(event);
  lumihist_afterNJets->fill(event);
  if(event.get(h_is_muevt)){
    h_fatjetsel_mu->fill(event);
    if(event.muons->at(0).pt()<=60) h_fatjetsel_mu_lowpt->fill(event);
    else h_fatjetsel_mu_highpt->fill(event);
  }
  else {
    h_fatjetsel_ele->fill(event);
    if(event.electrons->at(0).pt()<=120) h_fatjetsel_ele_lowpt->fill(event);
    else h_fatjetsel_ele_highpt->fill(event);
  }
  if(debug) cout << "Filled hists after fatjetsel" << endl;

  // MET Selection
  bool pass_MET =  met_sel->passes(event);
  if(!pass_MET) return false;
  h_METsel->fill(event);
  h_METsel_gen->fill(event);
  lumihist_afterMET->fill(event);
  if(event.get(h_is_muevt)){
    h_METsel_mu->fill(event);
    if(event.muons->at(0).pt()<=60) h_METsel_mu_lowpt->fill(event);
    else h_METsel_mu_highpt->fill(event);
  }
  else {
    h_METsel_ele->fill(event);
    if(event.electrons->at(0).pt()<=120) h_METsel_ele_lowpt->fill(event);
    else h_METsel_ele_highpt->fill(event);
  }
  if(debug) cout<<"Filled hists after MET"<<endl;

  if(debug) cout << "########### Done with preselection! ###########" << endl << endl;

  event.set(h_evt_weight, event.weight);

  return true;
}

// as we want to run the ExampleCycleNew directly with AnalysisModuleRunner,
// make sure the TstarTstarPreselectionModule is found by class name. This is ensured by this macro:
UHH2_REGISTER_ANALYSIS_MODULE(TstarTstarPreselectionModule)

}
