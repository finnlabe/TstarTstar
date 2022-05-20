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
#include "UHH2/common/include/MCWeight.h"

// TstarTstar custom stuff
#include "UHH2/TstarTstar/include/TstarTstarCustomIds.h"
#include "UHH2/TstarTstar/include/TstarTstarSelections.h"
#include "UHH2/TstarTstar/include/TstarTstarHists.h"
#include "UHH2/TstarTstar/include/TstarTstarGenHists.h"
#include "UHH2/TstarTstar/include/TstarTstarGenRecoMatchedHists.h"
#include "UHH2/TstarTstar/include/TstarTstarSFHists.h"

// HOTVR
#include "UHH2/HOTVR/include/HOTVRJetCorrectionModule.h"


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

  // cleaners
  unique_ptr<TopJetCleaner> HOTVRcleaner;

  // other external selections
  unique_ptr<Selection> met_sel;

  // other stuff
  std::unique_ptr<uhh2::AnalysisModule> ttgenprod;


  // ##### Histograms #####
  // full hists
  std::unique_ptr<Hists> h_nocuts,     h_common,         h_lepsel,        h_jetsel,        h_fatjetsel,        h_METsel;
  std::unique_ptr<Hists> h_nocuts_gen, h_common_gen,     h_lepsel_gen,    h_jetsel_gen,    h_fatjetsel_gen,    h_METsel_gen;
  std::unique_ptr<LuminosityHists> lumihist_common, lumihist_lepsel, lumihist_jetsel, lumihist_fatjetsel, lumihist_METsel;

  // electron channel
  std::unique_ptr<Hists> h_lepsel_ele,        h_jetsel_ele,        h_fatjetsel_ele,        h_METsel_ele;
  std::unique_ptr<Hists> h_lepsel_ele_lowpt,  h_jetsel_ele_lowpt,  h_fatjetsel_ele_lowpt,  h_METsel_ele_lowpt;
  std::unique_ptr<Hists> h_lepsel_ele_highpt, h_jetsel_ele_highpt, h_fatjetsel_ele_highpt, h_METsel_ele_highpt;

  // muon channel
  std::unique_ptr<Hists> h_lepsel_mu,        h_jetsel_mu,        h_fatjetsel_mu,        h_METsel_mu;
  std::unique_ptr<Hists> h_lepsel_mu_lowpt,  h_jetsel_mu_lowpt,  h_fatjetsel_mu_lowpt,  h_METsel_mu_lowpt;
  std::unique_ptr<Hists> h_lepsel_mu_highpt, h_jetsel_mu_highpt, h_fatjetsel_mu_highpt, h_METsel_mu_highpt;

  std::unique_ptr<Hists> h_afterSelection_gen, h_afterSelection_genmatch;

  // ##### Handles #####
  uhh2::Event::Handle<TTbarGen> h_ttbargen;
  uhh2::Event::Handle<bool> h_is_muevt;
  uhh2::Event::Handle<bool> h_is_highpt;
  uhh2::Event::Handle<double> h_evt_weight;

  uhh2::Event::Handle<bool> h_MC_isfake2017B;
  uhh2::Event::Handle<bool> h_MC_isfake2016B;

  // ##### Control switches #####
  bool debug = false;

  // ##### other needed definitions #####
  TString year;
  bool is_MC;
  bool data_isMu = false;
  bool data_isPhoton = false;
  bool data_is2017B = false;
  bool data_is2016B = false;
  bool isTriggerSFMeasurement = false;

};


TstarTstarPreselectionModule::TstarTstarPreselectionModule(Context & ctx){

  // setting debug from xml file
  if(ctx.get("debug", "<not set>") == "true") debug = true;

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
    else if(ctx.get("dataset_version").find("UL16preVFP") != std::string::npos) year = "UL16preVFP";
    else if(ctx.get("dataset_version").find("UL16postVFP") != std::string::npos) year = "UL16postVFP";
    else if(ctx.get("dataset_version").find("UL17") != std::string::npos) year = "UL17";
    else if(ctx.get("dataset_version").find("UL18") != std::string::npos) year = "UL18";
    else throw "No year found in dataset name!";
  }
  if(debug) cout << "Year is " << year << "." << endl;

  // MC or real data
  is_MC = ctx.get("dataset_type") == "MC";

  // muon channel for DATA
  if(!is_MC) data_isMu = (ctx.get("dataset_version").find("SingleMuon") != std::string::npos);
  if(data_isMu) std::cout << "this data sample is a muon sample" << std::endl;
  if(!is_MC) data_isPhoton = (ctx.get("dataset_version").find("SinglePhoton") != std::string::npos);
  if(data_isPhoton) std::cout << "this data sample is a photon sample" << std::endl;

  // check this specific run as it has error
  if(!is_MC) data_is2017B = (ctx.get("dataset_version").find("SingleElectron_RunB_UL17") != std::string::npos) || (ctx.get("dataset_version").find("SingleMuon_RunB_UL17") != std::string::npos) || (ctx.get("dataset_version").find("SinglePhoton_RunB_UL17") != std::string::npos);
  if(!is_MC) data_is2016B = (ctx.get("dataset_version").find("SingleElectron_RunB_UL16preVFP") != std::string::npos) || (ctx.get("dataset_version").find("SingleMuon_RunB_UL16preVFP") != std::string::npos) || (ctx.get("dataset_version").find("SinglePhoton_RunB_UL16preVFP") != std::string::npos);

  // ###### 1. set up modules ######
  if(debug) cout << "Setting up modules" << endl;

  // CommonModules
  common.reset(new CommonModules());
  if(debug) cout << "Common ini" << endl;

  // HOTVR jets
  HOTVRcleaner.reset(new TopJetCleaner(ctx, AndId<Jet>(PtEtaCut(200, 2.5), JetPFID(JetPFID::WP_TIGHT_PUPPI)) ));
  HOTVRCorr.reset(new HOTVRJetCorrectionModule(ctx));

  if(debug) cout << "HOTVR done" << endl;

  // Electron
  ElectronId eleID_lowpt = ElectronTagID(Electron::mvaEleID_Fall17_iso_V2_wp90);
  ElectronId eleID_highpt = ElectronTagID(Electron::mvaEleID_Fall17_noIso_V2_wp90);
  double electron_pt_lowpt(40.);
  double electron_pt_highpt(120.);
  common->set_electron_id(OrId<Electron>( AndId<Electron>(PtEtaSCCut(electron_pt_lowpt, 2.4), eleID_lowpt, EleMaxPtCut(120.)),  AndId<Electron>(PtEtaSCCut(electron_pt_highpt, 2.4), eleID_highpt)));
  if(debug) cout << "Electrons done" << endl;

  // Muon
  MuonId muID_lowpt = AndId<Muon>(MuonID(Muon::CutBasedIdTight), MuonID(Muon::PFIsoTight));
  MuonId muID_highpt = MuonID(Muon::CutBasedIdGlobalHighPt);
  double muon_pt_lowpt(30.);
  double muon_pt_highpt(55.);
  common->set_muon_id(OrId<Muon>( AndId<Muon>(PtEtaCut(muon_pt_lowpt, 2.4), muID_lowpt, MuMaxPtCut(muon_pt_highpt)), AndId<Muon>(PtEtaCut(muon_pt_highpt, 2.4), muID_highpt)));
  if(debug) cout << "Muons done" << endl;

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
    if(debug) cout << "TTbarGenProducer done" << endl;
  }

  // ###### 2. set up selections ######
  // MET selection
  met_sel.reset(new METCut  (50.,1e9));

  // ###### 3. set up hists ######
  // general
  h_nocuts.reset(new TstarTstarHists(ctx, "NoCuts"));
  h_common.reset(new TstarTstarHists(ctx, "AfterCommon"));
  h_lepsel.reset(new TstarTstarHists(ctx, "AfterLep"));
  h_jetsel.reset(new TstarTstarHists(ctx, "AfterJets"));
  h_fatjetsel.reset(new TstarTstarHists(ctx, "AfterFatJets"));
  h_METsel.reset(new TstarTstarHists(ctx, "AfterMET"));

  // GEN hists
  h_nocuts_gen.reset(new TstarTstarGenHists(ctx, "NoCuts_gen"));
  h_common_gen.reset(new TstarTstarGenHists(ctx, "AfterCommon_gen"));
  h_lepsel_gen.reset(new TstarTstarGenHists(ctx, "AfterLep_gen"));
  h_jetsel_gen.reset(new TstarTstarGenHists(ctx, "AfterJets_gen"));
  h_fatjetsel_gen.reset(new TstarTstarGenHists(ctx, "AfterFatJets_gen"));
  h_METsel_gen.reset(new TstarTstarGenHists(ctx, "AfterMET_gen"));

  // Lumi hists
  lumihist_common.reset(new LuminosityHists(ctx, "lumihist_AfterCommon"));
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

  h_afterSelection_gen.reset(new TstarTstarGenHists(ctx, "AfterSel_gen"));
  h_afterSelection_genmatch.reset(new TstarTstarGenRecoMatchedHists(ctx, "AfterSel_genmatch"));

  // ###### 4. init handles ######
  h_is_muevt = ctx.declare_event_output<bool>("is_muevt");
  h_is_highpt = ctx.declare_event_output<bool>("is_highpt");
  h_evt_weight = ctx.declare_event_output<double>("evt_weight");

  string IsTriggerSFMeasurement = ctx.get("IsTriggerSFMeasurement", "False");
  if(IsTriggerSFMeasurement == "True") isTriggerSFMeasurement = true;
  if(isTriggerSFMeasurement) std::cout << "Changing the lepton requirement as this is Electron trigger SF measurement run!" << std::endl;

  h_MC_isfake2017B = ctx.declare_event_output<bool>("MC_isfake2017B");
  h_MC_isfake2016B = ctx.declare_event_output<bool>("MC_isfake2016B");

  std::srand(std::time(nullptr));


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
  if(!(HOTVRCorr->process(event))) return false;
  if(!(HOTVRcleaner->process(event))) return false;

  if(debug) cout<<"common modules done"<<endl;

  // hists before selection
  h_common->fill(event);
  h_common_gen->fill(event);
  lumihist_common->fill(event);
  if(debug) cout<<"Filled hists after cleaning"<<endl;

  // #################
  // ### Selection ###
  // #################

  bool MC_isfake2017B = false;
  bool MC_isfake2016B = false;

  if(year == "UL17" && is_MC) {
    // in UL17 we need to fake 11.6% (run B part) to a different trigger

    float r = static_cast <float> (rand()) / static_cast <float> (RAND_MAX);
    if(r < 0.116) MC_isfake2017B = true;

    event.set(h_MC_isfake2017B, MC_isfake2017B);
    event.set(h_MC_isfake2016B, false);
    if(debug && MC_isfake2017B) std::cout << "This MC event was set to be part of a fake 2017B run. r was " << r << "." << std::endl;

  } else if(year == "UL16preVFP" && is_MC) {
    // in UL16preVFP we need to fake 14.29% (run B part) to a different trigger

    float r = static_cast <float> (rand()) / static_cast <float> (RAND_MAX);
    if(r < 0.1429) MC_isfake2016B = true;

    event.set(h_MC_isfake2017B, false);
    event.set(h_MC_isfake2016B, MC_isfake2016B);
    if(debug && MC_isfake2016B) std::cout << "This MC event was set to be part of a fake 2017B run. r was " << r << "." << std::endl;

  } else {
    event.set(h_MC_isfake2017B, false);
    event.set(h_MC_isfake2016B, false);
  }

  // bools
  bool pass_lep1;
  bool pass_njet;
  bool pass_fat_njet;
  bool pass_MET;

  // ###### Lepton selection ######
  pass_lep1 = (((event.muons->size() == 1) || (event.electrons->size() == 1)) && (event.electrons->size()+event.muons->size()) == 1);
  if(isTriggerSFMeasurement) pass_lep1 = ((event.muons->size() == 1) && (event.electrons->size() == 1)); // this ensures orthogonal dataset
  if(!pass_lep1) return false;

  // setting muevt handle
  if(event.muons->size() == 1) event.set(h_is_muevt, true);
  else event.set(h_is_muevt, false);

  // set is_highpt
  if(event.get(h_is_muevt)){
    if(event.muons->at(0).pt()<=55) event.set(h_is_highpt, false);
    else event.set(h_is_highpt, true);
  }
  else {
    if(event.electrons->at(0).pt()<=120) event.set(h_is_highpt, false);
    else event.set(h_is_highpt, true);
  }

  // hists
  h_lepsel->fill(event);
  h_lepsel_gen->fill(event);
  lumihist_lepsel->fill(event);
  if(event.get(h_is_muevt)){
    h_lepsel_mu->fill(event);
    if(event.get(h_is_highpt)) h_lepsel_mu_lowpt->fill(event);
    else h_lepsel_mu_highpt->fill(event);
  }
  else {
    h_lepsel_ele->fill(event);
    if(event.get(h_is_highpt)) h_lepsel_ele_lowpt->fill(event);
    else h_lepsel_ele_highpt->fill(event);
  }
  if(debug) cout << "Filled hists after lepsel" << endl;

  // ###### jet selection ######
  pass_njet = (event.jets->size()>3);
  if(isTriggerSFMeasurement) pass_njet = (event.jets->size()>1); // only need to require 1 jet for trigger SF measurement
  if(!pass_njet) return false;

  // hists
  h_jetsel->fill(event);
  h_jetsel_gen->fill(event);
  lumihist_jetsel->fill(event);
  if(event.get(h_is_muevt)){
    h_jetsel_mu->fill(event);
    if(event.get(h_is_highpt)) h_jetsel_mu_lowpt->fill(event);
    else h_jetsel_mu_highpt->fill(event);
  }
  else {
    h_jetsel_ele->fill(event);
    if(event.get(h_is_highpt)) h_jetsel_ele_lowpt->fill(event);
    else h_jetsel_ele_highpt->fill(event);
  }
  if(debug) cout << "Filled hists after fatjetsel" << endl;

  // ###### fat jet selection ######
  pass_fat_njet = (event.topjets->size()>0);
  if(!pass_fat_njet) return false;

  // hists
  h_fatjetsel->fill(event);
  h_fatjetsel_gen->fill(event);
  lumihist_fatjetsel->fill(event);
  if(event.get(h_is_muevt)){
    h_fatjetsel_mu->fill(event);
    if(event.get(h_is_highpt)) h_fatjetsel_mu_lowpt->fill(event);
    else h_fatjetsel_mu_highpt->fill(event);
  }
  else {
    h_fatjetsel_ele->fill(event);
    if(event.get(h_is_highpt)) h_fatjetsel_ele_lowpt->fill(event);
    else h_fatjetsel_ele_highpt->fill(event);
  }
  if(debug) cout << "Filled hists after fatjetsel" << endl;

  // ###### MET Selection ######
  pass_MET =  met_sel->passes(event);
  if(!pass_MET) return false;

  // hists
  h_METsel->fill(event);
  h_METsel_gen->fill(event);
  lumihist_METsel->fill(event);
  if(event.get(h_is_muevt)){
    h_METsel_mu->fill(event);
    if(event.get(h_is_highpt)) h_METsel_mu_lowpt->fill(event);
    else h_METsel_mu_highpt->fill(event);
  }
  else {
    h_METsel_ele->fill(event);
    if(event.get(h_is_highpt)) h_METsel_ele_lowpt->fill(event);
    else h_METsel_ele_highpt->fill(event);
  }
  if(debug) cout<<"Filled hists after MET"<<endl;

  // some gen check hists
  h_afterSelection_gen->fill(event);
  h_afterSelection_genmatch->fill(event);

  // outputting event weight for following modules
  event.set(h_evt_weight, event.weight);

  if(debug) cout << "########### Done with preselection! ###########" << endl << endl;
  return true;
}

// as we want to run the ExampleCycleNew directly with AnalysisModuleRunner,
// make sure the TstarTstarPreselectionModule is found by class name. This is ensured by this macro:
UHH2_REGISTER_ANALYSIS_MODULE(TstarTstarPreselectionModule)

}
