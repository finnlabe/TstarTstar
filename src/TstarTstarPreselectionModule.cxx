#include <iostream>
#include <memory>

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
#include "UHH2/common/include/MCWeight.h"
#include <UHH2/common/include/TriggerSelection.h>

// TstarTstar custom stuff
#include "UHH2/TstarTstar/include/TstarTstarCustomIds.h"
#include "UHH2/TstarTstar/include/TstarTstarSelections.h"
#include "UHH2/TstarTstar/include/TstarTstarHists.h"
#include "UHH2/TstarTstar/include/TstarTstarGenHists.h"
#include "UHH2/TstarTstar/include/TstarTstarGenRecoMatchedHists.h"

// HOTVR
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

  // ##### Modules #####
  // common modules, corrections etc.
  // examples: JetPFid, JEC, JER, MET corrections, etc
  std::unique_ptr<CommonModules> common;
  std::unique_ptr<AnalysisModule> HOTVRCorr;
  std::unique_ptr<AnalysisModule> HadronicTopFinder;
  std::unique_ptr<AnalysisModule> HOTVRScale;

  // cleaners
  unique_ptr<TopJetCleaner> HOTVRcleaner;
  unique_ptr<MuonCleaner> MuCleaner_lowpt;
  unique_ptr<MuonCleaner> MuCleaner_highpt;
  unique_ptr<ElectronCleaner> EleCleaner_lowpt;
  unique_ptr<ElectronCleaner> EleCleaner_highpt;

  // trigger selections
  unique_ptr<Selection> trg_ele32;
  unique_ptr<Selection> trg_pho175;
  unique_ptr<Selection> trg_ele115;
  unique_ptr<Selection> trg_mu24_iso;
  unique_ptr<Selection> trg_mu24_iso_tk;
  unique_ptr<Selection> trg_mu50;

  // other external selections
  unique_ptr<Selection> met_sel;

  // other stuff
  std::unique_ptr<uhh2::AnalysisModule> ttgenprod;


  // ##### Histograms #####
  // full hists
  std::unique_ptr<Hists> h_nocuts,     h_common,        h_trigger,        h_lepsel,        h_jetsel,        h_fatjetsel,        h_METsel;
  std::unique_ptr<Hists> h_nocuts_gen, h_common_gen,    h_trigger_gen,    h_lepsel_gen,    h_jetsel_gen,    h_fatjetsel_gen,    h_METsel_gen;
  std::unique_ptr<LuminosityHists>     lumihist_common, lumihist_trigger, lumihist_lepsel, lumihist_jetsel, lumihist_fatjetsel, lumihist_METsel;

  // electron channel
  std::unique_ptr<Hists> h_lepsel_ele,        h_jetsel_ele,        h_fatjetsel_ele,        h_METsel_ele;
  std::unique_ptr<Hists> h_lepsel_ele_lowpt,  h_jetsel_ele_lowpt,  h_fatjetsel_ele_lowpt,  h_METsel_ele_lowpt;
  std::unique_ptr<Hists> h_lepsel_ele_highpt, h_jetsel_ele_highpt, h_fatjetsel_ele_highpt, h_METsel_ele_highpt;

  // muon channel
  std::unique_ptr<Hists> h_lepsel_mu,        h_jetsel_mu,        h_fatjetsel_mu,        h_METsel_mu;
  std::unique_ptr<Hists> h_lepsel_mu_lowpt,  h_jetsel_mu_lowpt,  h_fatjetsel_mu_lowpt,  h_METsel_mu_lowpt;
  std::unique_ptr<Hists> h_lepsel_mu_highpt, h_jetsel_mu_highpt, h_fatjetsel_mu_highpt, h_METsel_mu_highpt;


  // ##### Handles #####
  uhh2::Event::Handle<TTbarGen> h_ttbargen;
  uhh2::Event::Handle<bool> h_is_muevt;

  uhh2::Event::Handle<double> h_evt_weight;


  // ##### Control switches #####
  bool debug = false;
  bool doTriggerSel = true;


  // ##### other needed definitions #####
  TString year;
  bool is_MC;
  bool data_isMu = false;
  bool data_is2017B = false;

};


TstarTstarPreselectionModule::TstarTstarPreselectionModule(Context & ctx){

  // debug status message
  if(debug) {
    cout << "Hello World from TstarTstarPreselectionModule!" << endl;
    // If running in SFrame, the keys "dataset_version", "dataset_type", "dataset_lumi",
    // and "target_lumi" are set to the according values in the xml file. For CMSSW, these are
    // not set automatically, but can be set in the python config file.
    for(auto & kv : ctx.get_all()){
      cout << " " << kv.first << " = " << kv.second << endl;
    }
  }

  // ###### 0. Setting some MISC variables ######
  // year of sample
  year = ctx.get("year", "<not set>");
  if(year == "<not set>"){
    if(ctx.get("dataset_version").find("2016") != std::string::npos) year = "2016";
    else if(ctx.get("dataset_version").find("2017") != std::string::npos) year = "2017";
    else if(ctx.get("dataset_version").find("2018") != std::string::npos) year = "2018";
    else throw "No year found in dataset name!";
  }
  if(debug) cout << "Year is " << year << "." << endl;

  // MC or real data
  is_MC = ctx.get("dataset_type") == "MC";

  // muon channel for DATA
  if(!is_MC) data_isMu = (ctx.get("dataset_version").find("SingleMuon") != std::string::npos);

  // check this specific run as it has error
  if(!is_MC) data_is2017B = (ctx.get("dataset_version").find("SingleElectron2017_RunB") != std::string::npos);

  // ###### 1. set up modules ######
  if(debug) cout << "Setting up modules" << endl;
  // CommonModules
  common.reset(new CommonModules());
  common->switch_metcorrection();
  if(debug) cout << "Common done" << endl;

  // HOTVR jets
  // ctx.get("HOTVRTopTagSFs"); // TODO check do I need this line?
  // HOTVRCorr.reset(new HOTVRJetCorrectionModule(ctx)); // crashes
  // HadronicTopFinder.reset(new HadronicTop(ctx));
  // TopJetId topjetID = AndId<TopJet>(HOTVRTopTag(), Tau32Groomed(0.56)); // Top Tag that is used later
  // HOTVRScale.reset(new HOTVRScaleFactor(ctx, topjetID));
  HOTVRcleaner.reset(new TopJetCleaner(ctx, PtEtaCut(150.0, 2.5)));
  if(debug) cout << "HOTVR done" << endl;

  // Electron
  ElectronId eleID_lowpt = ElectronID_Summer16_tight;
  ElectronId eleID_highpt = ElectronID_Summer16_tight_noIso;
  double electron_pt_lowpt(30.);
  double electron_pt_highpt(120.);
  common->set_electron_id(OrId<Electron>( AndId<Electron>(PtEtaSCCut(electron_pt_lowpt, 2.4), eleID_lowpt, EleMaxPtCut(120.)),  AndId<Electron>(PtEtaSCCut(electron_pt_highpt, 2.4), eleID_highpt)));
  if(debug) cout << "Electrons done" << endl;

  // Muon
  MuonId muID_lowpt = AndId<Muon>(MuonID(Muon::CutBasedIdTight), MuonID(Muon::PFIsoTight));
  MuonId muID_highpt = MuonID(Muon::CutBasedIdTight);
  double muon_pt_lowpt(27.);
  double muon_pt_highpt(60.);
  common->set_muon_id(OrId<Muon>( AndId<Muon>(PtEtaCut(muon_pt_lowpt, 2.4), muID_lowpt, MuMaxPtCut(60.)),  AndId<Muon>(PtEtaCut(muon_pt_highpt, 2.4), muID_highpt)));
  if(debug) cout << "Muons done" << endl;

  // AK4 Jets
  common->switch_jetlepcleaner();
  common->switch_jetPtSorter();
  double jet_pt(30.);
  common->set_jet_id(AndId<Jet>(PtEtaCut(jet_pt, 2.5), JetPFID(JetPFID::WP_TIGHT_PUPPI)));
  if(debug) cout << "Jets done" << endl;

  // init common
  common->init(ctx);
  if(debug) cout << "Common init done" << endl;

  // find ttbar for GEN
  if(is_MC){
    // Prepare GEN
    ttgenprod.reset(new TTbarGenProducer(ctx, "ttbargen", false));
  }

  // ###### 2. set up selections ######
  // MET selection
  met_sel.reset(new METCut  (50.,1e9));

  // trigger selections
  if(is_MC || !data_isMu) {
    // The following exist for both 2016 and 2017
    // until 120 GeV, except for 2017B, there for whole range
    if(year == "2018") trg_ele32.reset(new TriggerSelection("HLT_Ele32_WPTight_Gsf_v*"));
    else trg_ele32.reset(new TriggerSelection("HLT_Ele27_WPTight_Gsf_v*"));
    // above 120 GeV
    trg_pho175.reset(new TriggerSelection("HLT_Photon175_v*"));
    if(!data_is2017B) trg_ele115.reset(new TriggerSelection("HLT_Ele115_CaloIdVT_GsfTrkIdT_v*"));
  }
  if(is_MC || data_isMu){
    // until 27 GeV
    trg_mu24_iso.reset(new TriggerSelection("HLT_IsoMu24_v*"));
    if(year == "2016") trg_mu24_iso_tk.reset(new TriggerSelection("HLT_IsoTkMu24_v*"));
    // above 60 GeV
    trg_mu50.reset(new TriggerSelection("HLT_Mu50_v*"));
  }

  // ###### 3. set up hists ######
  // general
  h_nocuts.reset(new TstarTstarHists(ctx, "NoCuts"));
  h_common.reset(new TstarTstarHists(ctx, "AfterCommon"));
  h_trigger.reset(new TstarTstarHists(ctx, "AfterTrigger"));
  h_lepsel.reset(new TstarTstarHists(ctx, "AfterLep"));
  h_jetsel.reset(new TstarTstarHists(ctx, "AfterJets"));
  h_fatjetsel.reset(new TstarTstarHists(ctx, "AfterFatJets"));
  h_METsel.reset(new TstarTstarHists(ctx, "AfterMET"));

  // GEN hists
  h_nocuts_gen.reset(new TstarTstarGenHists(ctx, "NoCuts_gen"));
  h_common_gen.reset(new TstarTstarGenHists(ctx, "AfterCommon_gen"));
  h_trigger_gen.reset(new TstarTstarGenHists(ctx, "AfterTrigger_gen"));
  h_lepsel_gen.reset(new TstarTstarGenHists(ctx, "AfterLep_gen"));
  h_jetsel_gen.reset(new TstarTstarGenHists(ctx, "AfterJets_gen"));
  h_fatjetsel_gen.reset(new TstarTstarGenHists(ctx, "AfterFatJets_gen"));
  h_METsel_gen.reset(new TstarTstarGenHists(ctx, "AfterMET_gen"));

  // Lumi hists
  lumihist_common.reset(new LuminosityHists(ctx, "lumihist_AfterCommon"));
  lumihist_trigger.reset(new LuminosityHists(ctx, "lumihist_AfterTrigger"));
  lumihist_lepsel.reset(new LuminosityHists(ctx, "lumihist_AfterLep"));
  lumihist_jetsel.reset(new LuminosityHists(ctx, "lumihist_AfterJets"));
  lumihist_fatjetsel.reset(new LuminosityHists(ctx, "lumihist_AfterFatJets"));
  lumihist_METsel.reset(new LuminosityHists(ctx, "lumihist_AfterMET"));

  // electron channel
  h_lepsel_ele.reset(new TstarTstarHists(ctx, "AfterLepSel_ele"));
  h_jetsel_ele.reset(new TstarTstarHists(ctx, "AfterJets_ele"));
  h_fatjetsel_ele.reset(new TstarTstarHists(ctx, "AfterFatJets_ele"));
  h_METsel_ele.reset(new TstarTstarHists(ctx, "AfterMET_ele"));

  h_lepsel_ele_lowpt.reset(new TstarTstarHists(ctx, "AfterLepSel_ele_lowpt"));
  h_jetsel_ele_lowpt.reset(new TstarTstarHists(ctx, "AfterJets_ele_lowpt"));
  h_fatjetsel_ele_lowpt.reset(new TstarTstarHists(ctx, "AfterFatJets_ele_lowpt"));
  h_METsel_ele_lowpt.reset(new TstarTstarHists(ctx, "AfterMET_ele_lowpt"));

  h_lepsel_ele_highpt.reset(new TstarTstarHists(ctx, "AfterLepSel_ele_highpt"));
  h_jetsel_ele_highpt.reset(new TstarTstarHists(ctx, "AfterJets_ele_highpt"));
  h_fatjetsel_ele_highpt.reset(new TstarTstarHists(ctx, "AfterAK8jets_ele_highpt"));
  h_METsel_ele_highpt.reset(new TstarTstarHists(ctx, "AfterMET_ele_highpt"));

  // muon channel
  h_lepsel_mu.reset(new TstarTstarHists(ctx, "AfterLepSel_mu"));
  h_jetsel_mu.reset(new TstarTstarHists(ctx, "AfterJets_mu"));
  h_fatjetsel_mu.reset(new TstarTstarHists(ctx, "AfterFatJets_mu"));
  h_METsel_mu.reset(new TstarTstarHists(ctx, "AfterMET_mu"));

  h_lepsel_mu_lowpt.reset(new TstarTstarHists(ctx, "AfterLepSel_mu_lowpt"));
  h_jetsel_mu_lowpt.reset(new TstarTstarHists(ctx, "AfterJets_mu_lowpt"));
  h_fatjetsel_mu_lowpt.reset(new TstarTstarHists(ctx, "AfterFatJets_mu_lowpt"));
  h_METsel_mu_lowpt.reset(new TstarTstarHists(ctx, "AfterMET_mu_lowpt"));

  h_lepsel_mu_highpt.reset(new TstarTstarHists(ctx, "AfterLepSel_mu_highpt"));
  h_jetsel_mu_highpt.reset(new TstarTstarHists(ctx, "AfterJets_mu_highpt"));
  h_fatjetsel_mu_highpt.reset(new TstarTstarHists(ctx, "AfterAK8jets_mu_highpt"));
  h_METsel_mu_highpt.reset(new TstarTstarHists(ctx, "AfterMET_mu_highpt"));


  // ###### 4. init handles ######
  h_is_muevt = ctx.declare_event_output<bool>("is_muevt");
  h_evt_weight = ctx.declare_event_output<double>("evt_weight");

}


bool TstarTstarPreselectionModule::process(Event & event) {

  // debug status messages
  if(debug){
    cout << "TstarTstarPreselectionModule: Starting to process event (runid, eventid) = (" << event.run << ", " << event.event << "); weight = " << event.weight << endl;
    cout<<"N muons = "<<event.muons->size()<<", N electrons = "<<event.electrons->size()<<", N photons = "<<event.photons->size()<<endl;
    for(auto& muo : *event.muons) cout<<"BEFORE Muon (pt,eta): "<<muo.pt()<<", "<<muo.eta()<<endl;
  }

  // Fill ttbar object for gen
  if(is_MC){
    ttgenprod->process(event);
  }

  // hists before anything happened
  h_nocuts->fill(event);
  h_nocuts_gen->fill(event);
  if(debug) cout<<"Filled hists without any cuts, and GEN with no cuts"<<endl;

  // ###### common modules, corrections & cleaning ######
  if(!(common->process(event))) return false;
  //if(!(HOTVRCorr->process(event))) return false;
  if(!(HOTVRcleaner->process(event))) return false;

  // ###### scale factors ######
  // Top Tagging scale factors
  // HadronicTopFinder->process(event);
  // HOTVRScale->process(event);

  // hists before selection
  h_common->fill(event);
  h_common_gen->fill(event);
  lumihist_common->fill(event);
  if(debug) cout<<"Filled hists after cleaning"<<endl;

  // #################
  // ### Selection ###
  // #################

  // ###### Trigger selection ######
  bool pass_trigger = false;
  bool pass_trigger_SingleMu_lowpt = false;
  bool pass_trigger_SingleMu_highpt = false;
  if(is_MC || data_isMu){
    if(year == "2016") {
      pass_trigger_SingleMu_lowpt = (trg_mu24_iso->passes(event) || trg_mu24_iso_tk->passes(event));
    }
    else pass_trigger_SingleMu_lowpt = trg_mu24_iso->passes(event);
    pass_trigger_SingleMu_highpt = trg_mu50->passes(event);
  }
  bool pass_trigger_SingleEle_lowpt = false;
  bool pass_trigger_SingleEle_highpt = false;
  if(is_MC || !data_isMu){
    pass_trigger_SingleEle_lowpt = trg_ele32->passes(event);
    if(data_is2017B) pass_trigger_SingleEle_highpt = (trg_pho175->passes(event) || trg_ele32->passes(event));
    else pass_trigger_SingleEle_highpt = (trg_pho175->passes(event) || trg_ele115->passes(event));
  }
  if(pass_trigger_SingleMu_lowpt && (event.muons->size() >= 1)){if(event.muons->at(0).pt()<=60) pass_trigger = true; }
  if(pass_trigger_SingleMu_highpt && (event.muons->size() >= 1)){if(event.muons->at(0).pt()>60) pass_trigger = true; }
  if(pass_trigger_SingleEle_lowpt && (event.electrons->size() >= 1)){if(event.electrons->at(0).pt()<=120)pass_trigger = true; }
  if(pass_trigger_SingleEle_highpt && (event.electrons->size() >= 1)){if(event.electrons->at(0).pt()>120)pass_trigger = true; }
  if(!pass_trigger && doTriggerSel) return false;

  // hists
  h_trigger->fill(event);
  h_trigger_gen->fill(event);
  lumihist_trigger->fill(event);
  if(debug) cout<<"Filled hists after Trigger"<<endl;


  // ###### Lepton selection ######
  const bool pass_lep1 = (((event.muons->size() == 1) || (event.electrons->size() == 1)) && (event.electrons->size()+event.muons->size()) == 1);
  if(!pass_lep1) return false;

  // setting muevt handle
  if(event.muons->size() == 1) event.set(h_is_muevt, true);
  else event.set(h_is_muevt, false);

  // hists
  h_lepsel->fill(event);
  h_lepsel_gen->fill(event);
  lumihist_lepsel->fill(event);
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

  // ###### jet selection ######
  bool pass_njet = (event.jets->size()>3);
  if(!pass_njet) return false;

  // hists
  h_jetsel->fill(event);
  h_jetsel_gen->fill(event);
  lumihist_jetsel->fill(event);
  if(event.get(h_is_muevt)){
    h_jetsel_mu->fill(event);
    if(event.muons->at(0).pt()<=60) h_jetsel_mu_lowpt->fill(event);
    else h_jetsel_mu_highpt->fill(event);
  }
  else {
    h_jetsel_ele->fill(event);
    if(event.electrons->at(0).pt()<=120) h_jetsel_ele_lowpt->fill(event);
    else h_jetsel_ele_highpt->fill(event);
  }
  if(debug) cout << "Filled hists after fatjetsel" << endl;

  // ###### fat jet selection ######
  bool pass_fat_njet = (event.topjets->size()>0);
  if(!pass_fat_njet) return false;

  // hists
  h_fatjetsel->fill(event);
  h_fatjetsel_gen->fill(event);
  lumihist_fatjetsel->fill(event);
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


  // ###### MET Selection ######
  bool pass_MET =  met_sel->passes(event);
  if(!pass_MET) return false;

  // hists
  h_METsel->fill(event);
  h_METsel_gen->fill(event);
  lumihist_METsel->fill(event);
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


  // outputting event weight for following modules
  event.set(h_evt_weight, event.weight);


  if(debug) cout << "########### Done with preselection! ###########" << endl << endl;
  return true;
}

// as we want to run the ExampleCycleNew directly with AnalysisModuleRunner,
// make sure the TstarTstarPreselectionModule is found by class name. This is ensured by this macro:
UHH2_REGISTER_ANALYSIS_MODULE(TstarTstarPreselectionModule)

}
