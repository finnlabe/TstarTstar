#include <iostream>
#include <memory>
#include <chrono>

// TODO clean. Do i need all those?
// UHH2 stuff
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
#include "UHH2/common/include/TTbarGen.h"
#include "UHH2/common/include/MCWeight.h"
#include <UHH2/common/include/DetectorCleaning.h>
#include "UHH2/common/include/JetIds.h"

// TstarTstar stuff
#include "UHH2/TstarTstar/include/ModuleBASE.h"
#include "UHH2/TstarTstar/include/TstarTstarSelections.h"
#include "UHH2/TstarTstar/include/TstarTstarHists.h"
#include "UHH2/TstarTstar/include/TstarTstarGenHists.h"
#include "UHH2/TstarTstar/include/TstarTstarGenRecoMatchedHists.h"
#include "UHH2/TstarTstar/include/TstarTstarRecoTstarHists.h"
#include "UHH2/TstarTstar/include/TstarTstarReconstructionModules.h"

// other stuff
#include "UHH2/HOTVR/include/HOTVRIds.h"
#include "UHH2/HOTVR/include/HadronicTop.h"
#include "UHH2/HOTVR/include/HOTVRScaleFactor.h"

using namespace std;
using namespace uhh2;

namespace uhh2 {

// quick method to calculate inv_mass
float inv_mass(const LorentzVector& p4){ return p4.isTimelike() ? p4.mass() : -sqrt(-p4.mass2()); }

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

  // ###### Modules ######
  // general
  unique_ptr<uhh2::AnalysisModule> reco_primlep;
  unique_ptr<uhh2::AnalysisModule> ttgenprod;

  // scale factors
  std::unique_ptr<AnalysisModule> HadronicTopFinder;
  std::unique_ptr<AnalysisModule> HOTVRScale;
  std::unique_ptr<AnalysisModule> ScaleFactor_muon_ID;
  std::unique_ptr<AnalysisModule> ScaleFactor_muon_iso;
  std::unique_ptr<AnalysisModule> ScaleFactor_muon_trigger;
  std::unique_ptr<AnalysisModule> ScaleFactor_ele_ID;
  std::unique_ptr<AnalysisModule> ScaleFactor_ele_ID_noiso;
  std::unique_ptr<AnalysisModule> ScaleFactor_ele_trigger;
  std::unique_ptr<AnalysisModule> ScaleFactor_btagging;

  // selections
  unique_ptr<Selection> twodcut_sel;
  unique_ptr<Selection> toptagevt_sel;
  unique_ptr<HEMCleanerSelection> HEMCleaner;

  // triggers
  // TODO clean this
  unique_ptr<Selection> triggerSingleJet450_sel;
  unique_ptr<Selection> triggerSingleLeptonEle1_sel;
  unique_ptr<Selection> triggerSingleLeptonEle2_sel;
  unique_ptr<Selection> triggerSingleLeptonEle3_sel;
  unique_ptr<Selection> triggerSingleLeptonMu1_sel;
  unique_ptr<Selection> triggerSingleLeptonMu2_sel;
  unique_ptr<Selection> triggerSingleLeptonMu3_sel;
  unique_ptr<Selection> triggerHT1_sel, triggerHT2_sel, triggerHT3_sel, triggerHT4_sel, triggerHT5_sel,  triggerHT6_sel;
  unique_ptr<Selection> triggerPFHT_sel;

  // ###### Histograms ######
  std::unique_ptr<Hists> h_beginSel,            h_btagcut,             h_2Dcut,               h_dRcut,               h_STcut            , h_corrections;
  std::unique_ptr<Hists> h_beginSel_gen,        h_btagcut_gen,         h_2Dcut_gen,           h_dRcut_gen,           h_STcut_gen        , h_corrections_gen;
  std::unique_ptr<Hists> h_beginSel_ele,        h_btagcut_ele,         h_2Dcut_ele,           h_dRcut_ele,           h_STcut_ele        , h_corrections_ele;
  std::unique_ptr<Hists> h_beginSel_ele_lowpt,  h_btagcut_ele_lowpt,   h_2Dcut_ele_lowpt,     h_dRcut_ele_lowpt,     h_STcut_ele_lowpt  , h_corrections_ele_lowpt;
  std::unique_ptr<Hists> h_beginSel_ele_highpt, h_btagcut_ele_highpt,  h_2Dcut_ele_highpt,    h_dRcut_ele_highpt,    h_STcut_ele_highpt , h_corrections_ele_highpt;
  std::unique_ptr<Hists> h_beginSel_mu,         h_btagcut_mu,          h_2Dcut_mu,            h_dRcut_mu,            h_STcut_mu         , h_corrections_mu;
  std::unique_ptr<Hists> h_beginSel_mu_lowpt,   h_btagcut_mu_lowpt,    h_2Dcut_mu_lowpt,      h_dRcut_mu_lowpt,      h_STcut_mu_lowpt   , h_corrections_mu_lowpt;
  std::unique_ptr<Hists> h_beginSel_mu_highpt,  h_btagcut_mu_highpt,   h_2Dcut_mu_highpt,     h_dRcut_mu_highpt,     h_STcut_mu_highpt  , h_corrections_mu_highpt;
  std::unique_ptr<Hists> h_beginSel_nobtag,     h_btagcut_nobtag,      h_2Dcut_nobtag,        h_dRcut_nobtag,        h_STcut_nobtag     , h_corrections_nobtag;
  std::unique_ptr<Hists> h_STcut_nobtag_ele, h_corrections_nobtag_ele;

  std::unique_ptr<Hists> h_afterSelection_gen, h_afterSelection_genmatch;
  std::unique_ptr<Hists> h_afterSelection, h_afterHEMcleaning, h_afterHEMcleaning_ele, h_afterHEMcleaning_mu;

  // TODO better trigger plots!
  std::unique_ptr<Hists> h_trigger, h_trigger_mu, h_trigger_ele;
  std::unique_ptr<Hists> h_trigger_gen;

  std::unique_ptr<Hists> h_beforeBcorrections, h_afterBcorrections, h_dRcut_nobtag_ele;

  // ###### Handles ######
  uhh2::Event::Handle<FlavorParticle> h_primlep;
  uhh2::Event::Handle<TTbarGen> h_ttbargen;
  uhh2::Event::Handle<bool> h_is_muevt;
  uhh2::Event::Handle<bool> h_is_highpt;
  uhh2::Event::Handle<double> h_evt_weight;
  uhh2::Event::Handle<LorentzVector> h_neutrino;
  uhh2::Event::Handle<double> h_ST;
  uhh2::Event::Handle<bool> h_is_triggered;
  uhh2::Event::Handle<bool> h_is_btagevent;

  // SF handles
  uhh2::Event::Handle<float> h_weight_sfmu_id;
  uhh2::Event::Handle<float> h_weight_sfmu_id_down;
  uhh2::Event::Handle<float> h_weight_sfmu_id_up;
  uhh2::Event::Handle<float> h_weight_sfmu_isolation;
  uhh2::Event::Handle<float> h_weight_sfmu_isolation_down;
  uhh2::Event::Handle<float> h_weight_sfmu_isolation_up;
  uhh2::Event::Handle<float> h_weight_sfele_id;
  uhh2::Event::Handle<float> h_weight_sfele_id_down;
  uhh2::Event::Handle<float> h_weight_sfele_id_up;

  // ###### Control Switches ######
  bool debug = false;

  // ###### other needed definitions ######
  bool is_MC;
  bool data_isMu = false;
  bool data_is2017B = false;
  TString year;
  TopJetId topjetID;

  // b-tagging handles to try to fix memory leak
  uhh2::Event::Handle<float> h_weight_btagdisc_central;
  uhh2::Event::Handle<float> h_weight_btagdisc_jesup;
  uhh2::Event::Handle<float> h_weight_btagdisc_jesdown;
  uhh2::Event::Handle<float> h_weight_btagdisc_lfup;
  uhh2::Event::Handle<float> h_weight_btagdisc_lfdown;
  uhh2::Event::Handle<float> h_weight_btagdisc_hfup;
  uhh2::Event::Handle<float> h_weight_btagdisc_hfdown;
  uhh2::Event::Handle<float> h_weight_btagdisc_hfstats1up;
  uhh2::Event::Handle<float> h_weight_btagdisc_hfstats1down;
  uhh2::Event::Handle<float> h_weight_btagdisc_hfstats2up;
  uhh2::Event::Handle<float> h_weight_btagdisc_hfstats2down;
  uhh2::Event::Handle<float> h_weight_btagdisc_lfstats1up;
  uhh2::Event::Handle<float> h_weight_btagdisc_lfstats1down;
  uhh2::Event::Handle<float> h_weight_btagdisc_lfstats2up;
  uhh2::Event::Handle<float> h_weight_btagdisc_lfstats2down;
  uhh2::Event::Handle<float> h_weight_btagdisc_cferr1up;
  uhh2::Event::Handle<float> h_weight_btagdisc_cferr1down;
  uhh2::Event::Handle<float> h_weight_btagdisc_cferr2up;
  uhh2::Event::Handle<float> h_weight_btagdisc_cferr2down;

};


TstarTstarSelectionModule::TstarTstarSelectionModule(Context & ctx):
  h_weight_btagdisc_central (ctx.declare_event_output<float>("weight_btagdisc_central")),
  h_weight_btagdisc_jesup (ctx.declare_event_output<float>("weight_btagdisc_jesup")),
  h_weight_btagdisc_jesdown (ctx.declare_event_output<float>("weight_btagdisc_jesdown")),
  h_weight_btagdisc_lfup (ctx.declare_event_output<float>("weight_btagdisc_lfup")),
  h_weight_btagdisc_lfdown (ctx.declare_event_output<float>("weight_btagdisc_lfdown")),
  h_weight_btagdisc_hfup (ctx.declare_event_output<float>("weight_btagdisc_hfup")),
  h_weight_btagdisc_hfdown (ctx.declare_event_output<float>("weight_btagdisc_hfdown")),
  h_weight_btagdisc_hfstats1up (ctx.declare_event_output<float>("weight_btagdisc_hfstats1up")),
  h_weight_btagdisc_hfstats1down (ctx.declare_event_output<float>("weight_btagdisc_hfstats1down")),
  h_weight_btagdisc_hfstats2up (ctx.declare_event_output<float>("weight_btagdisc_hfstats2up")),
  h_weight_btagdisc_hfstats2down (ctx.declare_event_output<float>("weight_btagdisc_hfstats2down")),
  h_weight_btagdisc_lfstats1up (ctx.declare_event_output<float>("weight_btagdisc_lfstats1up")),
  h_weight_btagdisc_lfstats1down (ctx.declare_event_output<float>("weight_btagdisc_lfstats1down")),
  h_weight_btagdisc_lfstats2up (ctx.declare_event_output<float>("weight_btagdisc_lfstats2up")),
  h_weight_btagdisc_lfstats2down (ctx.declare_event_output<float>("weight_btagdisc_lfstats2down")),
  h_weight_btagdisc_cferr1up (ctx.declare_event_output<float>("weight_btagdisc_cferr1up")),
  h_weight_btagdisc_cferr1down (ctx.declare_event_output<float>("weight_btagdisc_cferr1down")),
  h_weight_btagdisc_cferr2up (ctx.declare_event_output<float>("weight_btagdisc_cferr2up")),
  h_weight_btagdisc_cferr2down (ctx.declare_event_output<float>("weight_btagdisc_cferr2down"))
{

  // setting debug from xml file
  if(ctx.get("debug", "<not set>") == "true") debug = true;

  // debug message
  if(debug) {
    cout << "Hello World from TstarTstarSelectionModule!" << endl;
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

  // getting SF file path
  string SF_path = ctx.get("SF_path");

  // year of samples
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

  // muon channel for DATA
  if(!is_MC) data_isMu = (ctx.get("dataset_version").find("SingleMuon") != std::string::npos);

  // check this specific run for trigger
  if(!is_MC) data_is2017B = (ctx.get("dataset_version").find("SingleElectron2017_RunB") != std::string::npos);

  // ###### 1. Set up modules ######
  // ttbar on GEN
  if(is_MC) ttgenprod.reset(new TTbarGenProducer(ctx, "ttbargen", false));

  // primary lepton
  reco_primlep.reset(new PrimaryLepton(ctx));


  // ###### 2. set up selections ######
  if(debug) cout << "Setting up Selections." << endl;
  // 2D cut
  twodcut_sel.reset(new TwoDCut(0.4, 25.0));  // The same as in Z'->ttbar semileptonic

  if(debug) cout << "Setting up HOTVR scale." << endl;
  // HOTVR scale
  topjetID = AndId<TopJet>(HOTVRTopTag(), Tau32Groomed(0.56)); // Top Tag that is used later
  HadronicTopFinder.reset(new HadronicTop(ctx));
  if(debug) cout << "HERE?" << endl;
  //HOTVRScale.reset(new HOTVRScaleFactor(ctx, topjetID));

  if(debug) cout << "Setting up electron scale." << endl;
  // lepton scale
  // electron
  if(year == "2016") {
    ScaleFactor_ele_ID.reset(new MCElecScaleFactor(ctx, SF_path+"/electrons/2016LegacyReReco_ElectronMVA90_Fall17V2.root", 0., "id"));
    ScaleFactor_ele_ID_noiso.reset(new MCElecScaleFactor(ctx, SF_path+"/electrons/2016LegacyReReco_ElectronMVA90noiso_Fall17V2.root", 0., "id"));
  }
  else if(year == "2017") {
    ScaleFactor_ele_ID.reset(new MCElecScaleFactor(ctx, SF_path+"/electrons/2017_ElectronMVA90.root", 0., "id"));
    ScaleFactor_ele_ID_noiso.reset(new MCElecScaleFactor(ctx, SF_path+"/electrons/2017_ElectronMVA90noiso.root", 0., "id"));
  }
  else if(year == "2018" || year == "UL18") {
    ScaleFactor_ele_ID.reset(new MCElecScaleFactor(ctx, SF_path+"/electrons/2018_ElectronMVA90.root", 0., "id"));
    ScaleFactor_ele_ID_noiso.reset(new MCElecScaleFactor(ctx, SF_path+"/electrons/2018_ElectronMVA90noiso.root", 0., "id"));
  }
  h_weight_sfele_id = ctx.get_handle<float>("weight_sfelec_id");
  h_weight_sfele_id_down = ctx.get_handle<float>("weight_sfelec_id_down");
  h_weight_sfele_id_up = ctx.get_handle<float>("weight_sfelec_id_up");

  if(debug) cout << "Setting up muon scale." << endl;

  // muons
  if(year == "2016") {
    ScaleFactor_muon_ID.reset(new MCMuonScaleFactor(ctx, SF_path+"/muons/2016/Efficiencies_muon_generalTracks_Z_Run2016_UL_ID.root", "NUM_TightID_DEN_TrackerMuons_abseta_pt", 0., "id", false));
    ScaleFactor_muon_iso.reset(new MCMuonScaleFactor(ctx, SF_path+"/muons/2016/Efficiencies_muon_generalTracks_Z_Run2016_UL_ISO.root", "NUM_TightRelIso_DEN_TightIDandIPCut_abseta_pt", 0., "isolation", false));
  }
  else if(year == "2017") {
    ScaleFactor_muon_ID.reset(new MCMuonScaleFactor(ctx, SF_path+"/muons/2017/Efficiencies_muon_generalTracks_Z_Run2017_UL_ID.root", "NUM_TightID_DEN_TrackerMuons_abseta_pt", 0., "id", false));
    ScaleFactor_muon_iso.reset(new MCMuonScaleFactor(ctx, SF_path+"/muons/2017/Efficiencies_muon_generalTracks_Z_Run2017_UL_ISO.root", "NUM_TightRelIso_DEN_TightIDandIPCut_abseta_pt", 0., "isolation", false));
  }
  else if(year == "2018" || year == "UL18") {
    ScaleFactor_muon_ID.reset(new MCMuonScaleFactor(ctx, SF_path+"/muons/2018/Efficiencies_muon_generalTracks_Z_Run2018_UL_ID.root", "NUM_TightID_DEN_TrackerMuons_abseta_pt", 0., "id", false));
    ScaleFactor_muon_iso.reset(new MCMuonScaleFactor(ctx, SF_path+"/muons/2018/Efficiencies_muon_generalTracks_Z_Run2018_UL_ISO.root", "NUM_TightRelIso_DEN_TightIDandIPCut_abseta_pt", 0., "isolation", false));
  }
  h_weight_sfmu_id = ctx.get_handle<float>("weight_sfmu_id");
  h_weight_sfmu_id_down = ctx.get_handle<float>("weight_sfmu_id_down");
  h_weight_sfmu_id_up = ctx.get_handle<float>("weight_sfmu_id_up");
  h_weight_sfmu_isolation = ctx.get_handle<float>("weight_sfmu_isolation");
  h_weight_sfmu_isolation_down = ctx.get_handle<float>("weight_sfmu_isolation_down");
  h_weight_sfmu_isolation_up = ctx.get_handle<float>("weight_sfmu_isolation_up");

  if(debug) cout << "Setting up btagging scale." << endl;

  // b-tagging SFs
  ScaleFactor_btagging.reset(new MCBTagDiscriminantReweighting(ctx, BTag::algo::DEEPJET)); // should be enough like this

  if(debug) cout << "Setting up HEM fix." << endl;

  // HEM issue
  HEMCleaner.reset(new HEMCleanerSelection(ctx, "jets", "topjets"));

  // Trigger selections
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

  // ###### 3. Set up Hists ######
  if(debug) cout << "Setting up Hists." << endl;
  h_beginSel.reset(new TstarTstarHists(ctx, "beginSel"));
  h_beginSel_gen.reset(new TstarTstarGenHists(ctx, "beginSel_gen"));
  h_beginSel_mu.reset(new TstarTstarHists(ctx, "beginSel_mu"));
  h_beginSel_mu_lowpt.reset(new TstarTstarHists(ctx, "beginSel_mu_lowpt"));
  h_beginSel_mu_highpt.reset(new TstarTstarHists(ctx, "beginSel_mu_highpt"));
  h_beginSel_ele.reset(new TstarTstarHists(ctx, "beginSel_ele"));
  h_beginSel_ele_lowpt.reset(new TstarTstarHists(ctx, "beginSel_ele_lowpt"));
  h_beginSel_ele_highpt.reset(new TstarTstarHists(ctx, "beginSel_ele_highpt"));

  h_btagcut.reset(new TstarTstarHists(ctx, "AfterBtag"));
  h_btagcut_gen.reset(new TstarTstarGenHists(ctx, "AfterBtag_gen"));
  h_btagcut_ele.reset(new TstarTstarHists(ctx, "AfterBtag_ele"));
  h_btagcut_ele_lowpt.reset(new TstarTstarHists(ctx, "AfterBtag_ele_lowpt"));
  h_btagcut_ele_highpt.reset(new TstarTstarHists(ctx, "AfterBtag_ele_highpt"));
  h_btagcut_mu.reset(new TstarTstarHists(ctx, "AfterBtag_mu"));
  h_btagcut_mu_lowpt.reset(new TstarTstarHists(ctx, "AfterBtag_mu_lowpt"));
  h_btagcut_mu_highpt.reset(new TstarTstarHists(ctx, "AfterBtag_mu_highpt"));

  h_2Dcut.reset(new TstarTstarHists(ctx, "After2D"));
  h_2Dcut_gen.reset(new TstarTstarGenHists(ctx, "After2D_gen"));
  h_2Dcut_ele.reset(new TstarTstarHists(ctx, "After2D_ele"));
  h_2Dcut_ele_lowpt.reset(new TstarTstarHists(ctx, "After2D_ele_lowpt"));
  h_2Dcut_ele_highpt.reset(new TstarTstarHists(ctx, "After2D_ele_highpt"));
  h_2Dcut_mu.reset(new TstarTstarHists(ctx, "After2D_mu"));
  h_2Dcut_mu_lowpt.reset(new TstarTstarHists(ctx, "After2D_mu_lowpt"));
  h_2Dcut_mu_highpt.reset(new TstarTstarHists(ctx, "After2D_mu_highpt"));

  h_dRcut.reset(new TstarTstarHists(ctx, "AfterdR"));
  h_dRcut_gen.reset(new TstarTstarGenHists(ctx, "AfterdR_gen"));
  h_dRcut_ele.reset(new TstarTstarHists(ctx, "AfterdR_ele"));
  h_dRcut_ele_lowpt.reset(new TstarTstarHists(ctx, "AfterdR_ele_lowpt"));
  h_dRcut_ele_highpt.reset(new TstarTstarHists(ctx, "AfterdR_ele_highpt"));
  h_dRcut_mu.reset(new TstarTstarHists(ctx, "AfterdR_mu"));
  h_dRcut_mu_lowpt.reset(new TstarTstarHists(ctx, "AfterdR_mu_lowpt"));
  h_dRcut_mu_highpt.reset(new TstarTstarHists(ctx, "AfterdR_mu_highpt"));

  h_STcut.reset(new TstarTstarHists(ctx, "AfterST"));
  h_STcut_gen.reset(new TstarTstarGenHists(ctx, "AfterST_gen"));
  h_STcut_ele.reset(new TstarTstarHists(ctx, "AfterST_ele"));
  h_STcut_ele_lowpt.reset(new TstarTstarHists(ctx, "AfterST_ele_lowpt"));
  h_STcut_ele_highpt.reset(new TstarTstarHists(ctx, "AfterST_ele_highpt"));
  h_STcut_mu.reset(new TstarTstarHists(ctx, "AfterST_mu"));
  h_STcut_mu_lowpt.reset(new TstarTstarHists(ctx, "AfterST_mu_lowpt"));
  h_STcut_mu_highpt.reset(new TstarTstarHists(ctx, "AfterST_mu_highpt"));
  h_STcut_nobtag.reset(new TstarTstarHists(ctx, "AfterST_nobtag"));
  h_STcut_nobtag_ele.reset(new TstarTstarHists(ctx, "AfterST_nobtag_ele"));

  h_corrections.reset(new TstarTstarHists(ctx, "AfterCorrections"));
  h_corrections_gen.reset(new TstarTstarGenHists(ctx, "AfterCorrections_gen"));
  h_corrections_ele.reset(new TstarTstarHists(ctx, "AfterCorrections_ele"));
  h_corrections_ele_lowpt.reset(new TstarTstarHists(ctx, "AfterCorrections_ele_lowpt"));
  h_corrections_ele_highpt.reset(new TstarTstarHists(ctx, "AfterCorrections_ele_highpt"));
  h_corrections_mu.reset(new TstarTstarHists(ctx, "AfterCorrections_mu"));
  h_corrections_mu_lowpt.reset(new TstarTstarHists(ctx, "AfterCorrections_mu_lowpt"));
  h_corrections_mu_highpt.reset(new TstarTstarHists(ctx, "AfterCorrections_mu_highpt"));
  h_corrections_nobtag.reset(new TstarTstarHists(ctx, "AfterCorrections_nobtag"));
  h_corrections_nobtag_ele.reset(new TstarTstarHists(ctx, "AfterCorrections_nobtag_ele"));

  h_beforeBcorrections.reset(new TstarTstarHists(ctx, "BeforeBCorrections"));
  h_afterBcorrections.reset(new TstarTstarHists(ctx, "AfterBCorrections"));

  h_afterSelection.reset(new TstarTstarHists(ctx, "AfterSel"));
  h_afterSelection_gen.reset(new TstarTstarGenHists(ctx, "AfterSel_gen"));
  h_afterSelection_genmatch.reset(new TstarTstarGenRecoMatchedHists(ctx, "AfterSel_genmatch"));

  h_afterHEMcleaning.reset(new TstarTstarHists(ctx, "AfterHEMcleaning"));
  h_afterHEMcleaning_ele.reset(new TstarTstarHists(ctx, "AfterHEMcleaning_ele"));
  h_afterHEMcleaning_mu.reset(new TstarTstarHists(ctx, "AfterHEMcleaning_mu"));

  // TODO
  h_trigger.reset(new TstarTstarHists(ctx, "TriggerXcheck"));
  h_trigger_mu.reset(new TstarTstarHists(ctx, "TriggerXcheck_mu"));
  h_trigger_ele.reset(new TstarTstarHists(ctx, "TriggerXcheck_ele"));
  h_trigger_gen.reset(new TstarTstarGenHists(ctx, "TriggerXcheck_gen"));


  // ###### 4. Init Handles ######
  h_is_muevt = ctx.get_handle<bool>("is_muevt");
  h_is_highpt = ctx.get_handle<bool>("is_highpt");
  h_evt_weight = ctx.get_handle<double>("evt_weight");
  h_primlep = ctx.get_handle<FlavorParticle>("PrimaryLepton");
  h_neutrino = ctx.declare_event_output<LorentzVector>("neutrino");
  h_ST = ctx.declare_event_output<double>("ST");
  h_is_triggered = ctx.get_handle<bool>("is_triggered");
  h_is_btagevent = ctx.declare_event_output<bool>("is_btagevent");

}


bool TstarTstarSelectionModule::process(Event & event) {

  auto time_start = std::chrono::system_clock::now();

  std::cout << "Muon event: " << event.get(h_is_muevt) << std::endl;

  // try to find memory leak
  event.set(h_weight_btagdisc_central,1.);
  event.set(h_weight_btagdisc_jesup,1.);
  event.set(h_weight_btagdisc_jesdown,1.);
  event.set(h_weight_btagdisc_lfup,1.);
  event.set(h_weight_btagdisc_lfdown,1.);
  event.set(h_weight_btagdisc_hfup,1.);
  event.set(h_weight_btagdisc_hfdown,1.);
  event.set(h_weight_btagdisc_hfstats1up,1.);
  event.set(h_weight_btagdisc_hfstats1down,1.);
  event.set(h_weight_btagdisc_hfstats2up,1.);
  event.set(h_weight_btagdisc_hfstats2down,1.);
  event.set(h_weight_btagdisc_lfstats1up,1.);
  event.set(h_weight_btagdisc_lfstats1down,1.);
  event.set(h_weight_btagdisc_lfstats2up,1.);
  event.set(h_weight_btagdisc_lfstats2down,1.);
  event.set(h_weight_btagdisc_cferr1up,1.);
  event.set(h_weight_btagdisc_cferr1down,1.);
  event.set(h_weight_btagdisc_cferr2up,1.);
  event.set(h_weight_btagdisc_cferr2down,1.);

  // debug status message
  if(debug) cout << "TstarTstarSelectionModule: Starting to process event (runid, eventid) = (" << event.run << ", " << event.event << "); weight = " << event.weight << endl;

  // reapply event weights from handle
  event.weight = event.get(h_evt_weight);
  if(debug) cout << "weights applied." << endl;
  bool is_triggered = event.get(h_is_triggered);

  // Fill ttgen object for correct matching check, etc
  if(is_MC) ttgenprod->process(event);
  if(debug) std::cout << "Filled ttgenprod" << endl;

  // set primary lepton
  reco_primlep->process(event);
  if(debug) std::cout << "Got primary lepton" << endl;

  // hists before anything happened
  if(is_triggered) {
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
  }

  auto time_beforebtagSF = std::chrono::system_clock::now();

  if(is_triggered && !event.get(h_is_muevt)) h_beforeBcorrections->fill(event);

  // b-tagging sfs
  ScaleFactor_btagging->process(event);
  auto time_afterbtagSF = std::chrono::system_clock::now();

  if(is_triggered && !event.get(h_is_muevt)) h_afterBcorrections->fill(event);

  // #################
  // ### Selection ###
  // #################
  if(debug) std::cout << "Starting selection" << endl;

  // ###### Btag Selection ######
  BTag bJetID = BTag(BTag::algo::DEEPJET, BTag::wp::WP_LOOSE);
  bool pass_btagcut = false;
  for (const auto & jet: *event.jets){
    //if(jet.btag_DeepJet() > 0.2219) pass_btagcut = true; // TODO replace this by proper method!!!
    if(bJetID(jet, event)) pass_btagcut = true;
  }
  event.set(h_is_btagevent,pass_btagcut);
  if(pass_btagcut && is_triggered) {
    // hists
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
  } // end of pass_btag = true
  else {
    // TODO control region plots maybe here?
  }

  if(debug) std::cout << "Done b-tag" << endl;
  auto time_afterbtagsel = std::chrono::system_clock::now();

    // ###### Lepton-2Dcut ######
    // if(debug) std::cout << "Start 2D cut" << endl;
    // bool pass_2D = true;
    // for(auto& muo : *event.muons){
    //   if(debug) cout<<"AFTER Muon (pt,eta): "<<muo.pt()<<", "<<muo.eta()<<endl;
    //   float    dRmin, pTrel;
    //   std::tie(dRmin, pTrel) = drmin_pTrel(muo, *event.jets);
    //   muo.set_tag(Muon::twodcut_dRmin, dRmin);
    //   muo.set_tag(Muon::twodcut_pTrel, pTrel);
    // }
    // for(auto& ele : *event.electrons){
    //   if(debug) cout<<"Electron (pt,eta): "<<ele.pt()<<", "<<ele.eta()<<endl;
    //   float    dRmin, pTrel;
    //   std::tie(dRmin, pTrel) = drmin_pTrel(ele, *event.jets);
    //   ele.set_tag(Electron::twodcut_dRmin, dRmin);
    //   ele.set_tag(Electron::twodcut_pTrel, pTrel);
    // }
    // const bool pass_twodcut = twodcut_sel->passes(event);
    // if(event.muons->size()==1){if(event.muons->at(0).pt()>60) pass_2D = pass_twodcut;}
    // else if(event.electrons->size()==1){if(event.electrons->at(0).pt()>120) pass_2D = pass_twodcut;}
    // else std::cout << "How did this happen???" << endl;
    // if(!pass_2D) return false;
    //
    // // hists
    // h_2Dcut->fill(event);
    // h_2Dcut_gen->fill(event);
    // if(event.get(h_is_muevt)){
    //   h_2Dcut_mu->fill(event);
    //   if(event.muons->at(0).pt()<=60) h_2Dcut_mu_lowpt->fill(event);
    //   else h_2Dcut_mu_highpt->fill(event);
    // }
    // else {
    //   h_2Dcut_ele->fill(event);
    //   if(event.electrons->at(0).pt()<=120) h_2Dcut_ele_lowpt->fill(event);
    //   else h_2Dcut_ele_highpt->fill(event);
    // }
    // if(debug) cout << "Passed 2D cut." << endl;

    if(debug) std::cout << "We have muons: " << event.muons->size() << endl;
    if(debug) std::cout << "We have electrons: " << event.electrons->size() << endl;
    if(debug) std::cout << "is_muevt is: " << event.get(h_is_muevt) << endl;


    // ###### dR cut to suppress QCD ######
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

    if(pass_btagcut && is_triggered) { // only fill these for btag cut passes
      // hists
      h_dRcut->fill(event);
      h_dRcut_gen->fill(event);
      if(event.get(h_is_muevt)){
        h_dRcut_mu->fill(event);
        if(event.muons->at(0).pt()<=60) h_dRcut_mu_lowpt->fill(event);
        else h_dRcut_mu_highpt->fill(event);
      }
      else {
        h_dRcut_ele->fill(event);
        if(event.electrons->at(0).pt()<=120) h_dRcut_ele_lowpt->fill(event);
        else h_dRcut_ele_highpt->fill(event);
      }
      if(debug) cout << "Passed dR cut." << endl;
    } else {
      // for later
    }

    if(debug) std::cout << "Done dR" << endl;
    auto time_afterdR = std::chrono::system_clock::now();


    // ST cut to reduce computation time (and file sizes)
    // Neutrin reconstruction
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

    // st calculation
    double st = 0.;
    for(const auto & jet : *event.topjets) st += jet.pt();
    for(const auto & lepton : *event.electrons) st += lepton.pt();
    for(const auto & lepton : *event.muons) st += lepton.pt();
    st += event.met->pt();
    event.set(h_ST, st);

    // st cut
    if(st < 500) return false;

    if(pass_btagcut && is_triggered) { // only fill these for btag cut passes
      // hists
      h_STcut->fill(event);
      h_STcut_gen->fill(event);
      if(event.get(h_is_muevt)){
        h_STcut_mu->fill(event);
        if(event.muons->at(0).pt()<=60) h_STcut_mu_lowpt->fill(event);
        else h_STcut_mu_highpt->fill(event);
      }
      else {
        h_STcut_ele->fill(event);
        if(event.electrons->at(0).pt()<=120) h_STcut_ele_lowpt->fill(event);
        else h_STcut_ele_highpt->fill(event);
      }
      if(debug) cout << "Passed ST cut." << endl;
    } else {
      if (is_triggered) {
        h_STcut_nobtag->fill(event);
        if(!event.get(h_is_muevt)) h_STcut_nobtag_ele->fill(event);
      }
    }

    if(debug) std::cout << "Done ST" << endl;

    auto time_afterST = std::chrono::system_clock::now();

    // ######################
    // ### Scale Factors! ###
    // ######################

    // HOTVR
    HadronicTopFinder->process(event);
    //HOTVRScale->process(event);

    if(debug) std::cout << "Done HOTVR scale" << endl;

    // lepton SFs
    if(event.get(h_is_muevt)){
      ScaleFactor_muon_ID->process(event);
      if(!event.get(h_is_highpt)) ScaleFactor_muon_iso->process(event);
      else {
        event.set(h_weight_sfmu_isolation, 1.);
        event.set(h_weight_sfmu_isolation_down, 1.);
        event.set(h_weight_sfmu_isolation_up, 1.);
      }
      event.set(h_weight_sfele_id, 1.);
      event.set(h_weight_sfele_id_down, 1.);
      event.set(h_weight_sfele_id_up, 1.);
    }
    else {
      if(event.get(h_is_highpt)) ScaleFactor_ele_ID_noiso->process(event);
      ScaleFactor_ele_ID->process(event);

      // fixing empty handles
      event.set(h_weight_sfmu_id, 1.);
      event.set(h_weight_sfmu_id_down, 1.);
      event.set(h_weight_sfmu_id_up, 1.);
      event.set(h_weight_sfmu_isolation, 1.);
      event.set(h_weight_sfmu_isolation_down, 1.);
      event.set(h_weight_sfmu_isolation_up, 1.);
    }

    if(debug) std::cout << "Done Lepton ID, ISO SFs" << endl;

    if(pass_btagcut && is_triggered) { // only fill these for btag cut passes
      // hists
      h_corrections->fill(event);
      h_corrections_gen->fill(event);
      if(event.get(h_is_muevt)){
        h_corrections_mu->fill(event);
        if(event.muons->at(0).pt()<=60) h_STcut_mu_lowpt->fill(event);
        else h_corrections_mu_highpt->fill(event);
      }
      else {
        h_corrections_ele->fill(event);
        if(event.electrons->at(0).pt()<=120) h_STcut_ele_lowpt->fill(event);
        else h_corrections_ele_highpt->fill(event);
      }
      if(debug) cout << "Passed ST cut." << endl;
    } else {
      if(is_triggered) {
        h_corrections_nobtag->fill(event);
        if(!event.get(h_is_muevt)) h_corrections_nobtag_ele->fill(event);
      }
    }

    auto time_afterSFs = std::chrono::system_clock::now();

    // #######################
    // ### Trigger studies ###
    // #######################

    if(pass_btagcut) { // fill these also for "non-triggered" events to have a comparison
      h_trigger->fill(event);
      h_trigger_gen->fill(event);
      if(event.get(h_is_muevt)) h_trigger_mu->fill(event);
      else h_trigger_ele->fill(event);
      if(debug) cout<<"Filled hists after Trigger"<<endl;
    }

    // Fixing the HEM issue
    if(!HEMCleaner->passes(event)) return false;
    if(is_MC) event.weight *= 0.913282; // TODO this will need to be changed once full SingleMuon is available

    if(pass_btagcut) { // fill these also for "non-triggered" events to have a comparison
      h_afterHEMcleaning->fill(event);
      if(event.get(h_is_muevt)) h_afterHEMcleaning_mu->fill(event);
      else h_afterHEMcleaning_ele->fill(event);
      if(debug) cout<<"Filled hists HEM fix"<<endl;
    }

    if(!is_triggered) return false; // now finally reject all non-triggered!

    if(pass_btagcut) { // only fill these for btag cut passes
      // some final plot for comparison
      h_afterSelection->fill(event);
      h_afterSelection_gen->fill(event);
      h_afterSelection_genmatch->fill(event);
    } else {
      //h_nobtagcontrolregion->fill();
    }

    // at the moment control region ends here!
    // if(pass_btagcut) return true;
    auto time_afterAll = std::chrono::system_clock::now();

    // timing output
    std::chrono::duration<double> elapsed_seconds_beforebtagSF = time_beforebtagSF-time_start;
    std::chrono::duration<double> elapsed_seconds_afterbtagSF = time_afterbtagSF-time_start;
    std::chrono::duration<double> elapsed_seconds_afterbtagsel = time_afterbtagsel-time_start;
    std::chrono::duration<double> elapsed_seconds_afterdR = time_afterdR-time_start;
    std::chrono::duration<double> elapsed_seconds_afterST = time_afterST-time_start;
    std::chrono::duration<double> elapsed_seconds_afterSFs = time_afterSFs-time_start;
    std::chrono::duration<double> elapsed_seconds_afterAll = time_afterAll-time_start;

    if(debug) {
      std::cout << "Time before btag SF: " << elapsed_seconds_beforebtagSF.count() << "s\n";
      std::cout << "Time after btag SF: " << elapsed_seconds_afterbtagSF.count() << "s\n";
      std::cout << "Time after btag sel: " << elapsed_seconds_afterbtagsel.count() << "s\n";
      std::cout << "Time after dR sel: " << elapsed_seconds_afterdR.count() << "s\n";
      std::cout << "Time after ST sel: " << elapsed_seconds_afterST.count() << "s\n";
      std::cout << "Time after all SFs: " << elapsed_seconds_afterSFs.count() << "s\n";
      std::cout << "Time after everything: " << elapsed_seconds_afterAll.count() << "s\n";
      std::cout << "###########" << std::endl << std::endl;
    }
    return true;

}

// as we want to run the ExampleCycleNew directly with AnalysisModuleRunner,
// make sure the TstarTstarSelectionModule is found by class name. This is ensured by this macro:
UHH2_REGISTER_ANALYSIS_MODULE(TstarTstarSelectionModule)

}
