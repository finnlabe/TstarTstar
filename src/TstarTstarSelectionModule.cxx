#include <iostream>
#include <memory>

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
#include "UHH2/common/include/LeptonScaleFactors.h"
#include "UHH2/common/include/TopPtReweight.h"

// TstarTstar stuff
#include "UHH2/TstarTstar/include/ModuleBASE.h"
#include "UHH2/TstarTstar/include/TstarTstarSelections.h"
#include "UHH2/TstarTstar/include/TstarTstarHists.h"
#include "UHH2/TstarTstar/include/TstarTstarGenHists.h"
#include "UHH2/TstarTstar/include/TstarTstarGenRecoMatchedHists.h"
#include "UHH2/TstarTstar/include/TstarTstarRecoTstarHists.h"
#include "UHH2/TstarTstar/include/TstarTstarReconstructionModules.h"
#include "UHH2/TstarTstar/include/ElecTriggerSF.h"
#include "UHH2/TstarTstar/include/TstarTstarScaleFactors.h"

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
  std::unique_ptr<AnalysisModule> ScaleFactor_btagging;
  TH2D *eventYieldFactors;
  std::unique_ptr<AnalysisModule> ScaleFactor_NNLO;

  // trigger selections
  std::unique_ptr<Selection> trg_ele_low;
  std::unique_ptr<Selection> trg_pho;
  std::unique_ptr<Selection> trg_ele_high;
  std::unique_ptr<Selection> trg_mu_low_1;
  std::unique_ptr<Selection> trg_mu_low_2;
  std::unique_ptr<Selection> trg_mu_high_1;
  std::unique_ptr<Selection> trg_mu_high_2;
  std::unique_ptr<Selection> trg_mu_high_3;

  // lepton SFs
  std::unique_ptr<AnalysisModule> sf_muon_iso;
  std::unique_ptr<AnalysisModule> sf_muon_ID_lowpt;
  std::unique_ptr<AnalysisModule> sf_muon_ID_highpt;
  std::unique_ptr<AnalysisModule> sf_muon_trigger_lowpt;
  std::unique_ptr<AnalysisModule> sf_muon_trigger_highpt;
  std::unique_ptr<AnalysisModule> sf_muon_iso_DUMMY;
  std::unique_ptr<AnalysisModule> sf_muon_ID_DUMMY;
  std::unique_ptr<AnalysisModule> sf_muon_trigger_DUMMY;

  std::unique_ptr<AnalysisModule> sf_ele_ID_lowpt;
  std::unique_ptr<AnalysisModule> sf_ele_ID_highpt;
  std::unique_ptr<AnalysisModule> sf_ele_reco;
  std::unique_ptr<AnalysisModule> sf_ele_ID_DUMMY;
  std::unique_ptr<AnalysisModule> sf_ele_reco_DUMMY;
  // std::unique_ptr<AnalysisModule> sf_ele_trigger; // TODO add

  // selections
  unique_ptr<Selection> twodcut_sel;
  unique_ptr<Selection> toptagevt_sel;
  unique_ptr<HEMCleanerSelection> HEMCleaner;
  unique_ptr<HEMCleanerMCScale> HEMCleanerMCScaler;

  // ###### Histograms ######
  std::unique_ptr<Hists> h_beginSel,            h_btagcut,             h_2Dcut,               h_dRcut,               h_STcut            ,    h_trigger, h_triggerSF,            h_corrections, h_topptreweighting              ;
  std::unique_ptr<Hists> h_beginSel_gen,        h_btagcut_gen,         h_2Dcut_gen,           h_dRcut_gen,           h_STcut_gen        ,    h_trigger_gen,        h_corrections_gen          ;
  std::unique_ptr<Hists> h_beginSel_ele,        h_btagcut_ele,         h_2Dcut_ele,           h_dRcut_ele,           h_STcut_ele        ,    h_trigger_ele, h_triggerSF_ele,       h_corrections_ele, h_topptreweighting_ele          ;
  std::unique_ptr<Hists> h_beginSel_ele_lowpt,  h_btagcut_ele_lowpt,   h_2Dcut_ele_lowpt,     h_dRcut_ele_lowpt,     h_STcut_ele_lowpt  ,    h_trigger_ele_lowpt, h_triggerSF_ele_lowpt,  h_corrections_ele_lowpt, h_topptreweighting_ele_lowpt    ;
  std::unique_ptr<Hists> h_beginSel_ele_highpt, h_btagcut_ele_highpt,  h_2Dcut_ele_highpt,    h_dRcut_ele_highpt,    h_STcut_ele_highpt ,    h_trigger_ele_highpt, h_triggerSF_ele_highpt, h_corrections_ele_highpt, h_topptreweighting_ele_highpt   ;
  std::unique_ptr<Hists> h_beginSel_mu,         h_btagcut_mu,          h_2Dcut_mu,            h_dRcut_mu,            h_STcut_mu         ,    h_trigger_mu,   h_triggerSF_mu,       h_corrections_mu, h_topptreweighting_mu           ;
  std::unique_ptr<Hists> h_beginSel_mu_lowpt,   h_btagcut_mu_lowpt,    h_2Dcut_mu_lowpt,      h_dRcut_mu_lowpt,      h_STcut_mu_lowpt   ,    h_trigger_mu_lowpt, h_triggerSF_mu_lowpt,  h_corrections_mu_lowpt, h_topptreweighting_mu_lowpt     ;
  std::unique_ptr<Hists> h_beginSel_mu_highpt,  h_btagcut_mu_highpt,   h_2Dcut_mu_highpt,     h_dRcut_mu_highpt,     h_STcut_mu_highpt  ,    h_trigger_mu_highpt,  h_triggerSF_mu_highpt, h_corrections_mu_highpt, h_topptreweighting_mu_highpt    ;
  std::unique_ptr<Hists> h_beginSel_nobtag,     h_btagcut_nobtag,      h_2Dcut_nobtag,        h_dRcut_nobtag,        h_STcut_nobtag     ,    h_trigger_nobtag, h_triggerSF_nobtag,     h_corrections_nobtag, h_topptreweighting_nobtag       ;
  std::unique_ptr<Hists> h_STcut_nobtag_ele, h_corrections_nobtag_ele;

  std::unique_ptr<Hists> h_afterHEMcleaning, h_afterHEMcleaning_ele, h_afterHEMcleaning_mu;

  std::unique_ptr<Hists> h_afterNNLO, h_afterNNLO_ele, h_afterNNLO_mu;

  std::unique_ptr<Hists> h_beforeBcorrections, h_afterBcorrections, h_afterBYieldcorrections;

  std::unique_ptr<Hists> h_notTriggered, h_notTriggered_ele, h_notTriggered_mu;

  // ###### Handles ######
  uhh2::Event::Handle<FlavorParticle> h_primlep;
  uhh2::Event::Handle<TTbarGen> h_ttbargen;
  uhh2::Event::Handle<bool> h_is_muevt;
  uhh2::Event::Handle<bool> h_is_highpt;
  uhh2::Event::Handle<double> h_evt_weight;
  uhh2::Event::Handle<LorentzVector> h_neutrino;
  uhh2::Event::Handle<double> h_ST;
  uhh2::Event::Handle<double> h_STHOTVR;
  uhh2::Event::Handle<bool> h_is_btagevent;
  uhh2::Event::Handle<bool> h_MC_isfake2017B;
  uhh2::Event::Handle<bool> h_MC_isfake2016B;

  // ###### Control Switches ######
  bool debug = false;

  // ###### other needed definitions ######
  bool is_MC;
  bool data_isMu = false;
  bool data_is2017B = false;
  bool data_isPhoton = false;
  bool isTriggerSFMeasurement = false;
  TString year;

  std::unique_ptr<AnalysisModule> TopPtReweighting;

};


TstarTstarSelectionModule::TstarTstarSelectionModule(Context & ctx) {

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

  // channel for DATA
  if(!is_MC) data_isMu = (ctx.get("dataset_version").find("SingleMuon") != std::string::npos);
  if(!is_MC) data_isPhoton = (ctx.get("dataset_version").find("SinglePhoton") != std::string::npos);
  if(data_isMu) std::cout << "This is a muon data sample" << std::endl;
  else if(data_isPhoton) std::cout << "This is a photon data sample" << std::endl;
  else if(!is_MC) std::cout << "This is an electron data sample" << std::endl;

  // check this specific run for trigger
  if(!is_MC) data_is2017B = (ctx.get("dataset_version").find("SingleElectron_RunB_UL17") != std::string::npos) || (ctx.get("dataset_version").find("SingleMuon_RunB_UL17") != std::string::npos) || (ctx.get("dataset_version").find("SinglePhoton_RunB_UL17") != std::string::npos);
  if(data_is2017B) std::cout << "This data sample is from 2017 Run B" << std::endl;

  // ###### 1. Set up modules ######
  // ttbar on GEN
  if(is_MC) ttgenprod.reset(new TTbarGenProducer(ctx, "ttbargen", false));

  // primary lepton
  reco_primlep.reset(new PrimaryLepton(ctx));

  // trigger selections
  if(is_MC || data_isMu) {
    std::cout << "Setting up muon triggers" << std::endl;

    // low pt triggers
    if(year == "2016" || year == "UL16preVFP" || year == "UL16postVFP") {
      trg_mu_low_1.reset(new TriggerSelection("HLT_IsoMu24_v*"));
      trg_mu_low_2.reset(new TriggerSelection("HLT_IsoTkMu24_v*"));
      trg_mu_high_1.reset(new TriggerSelection("HLT_Mu50_v*"));
      trg_mu_high_2.reset(new TriggerSelection("HLT_TkMu50_v*"));
    }
    else if(year == "2017" || year == "UL17") {
      trg_mu_low_1.reset(new TriggerSelection("HLT_IsoMu27_v*"));
      trg_mu_high_1.reset(new TriggerSelection("HLT_Mu50_v*"));
      trg_mu_high_2.reset(new TriggerSelection("HLT_TkMu100_v*"));
      trg_mu_high_3.reset(new TriggerSelection("HLT_OldMu100_v*"));
    }
    else if(year == "2018" || year == "UL18") {
      trg_mu_low_1.reset(new TriggerSelection("HLT_IsoMu24_v*"));
      trg_mu_high_1.reset(new TriggerSelection("HLT_Mu50_v*"));
      trg_mu_high_2.reset(new TriggerSelection("HLT_TkMu100_v*"));
      trg_mu_high_3.reset(new TriggerSelection("HLT_OldMu100_v*"));
    }
  }
  if(is_MC || !data_isMu){
    std::cout << "Setting up ele triggers" << std::endl;

    // low pt triggers
    if(year == "2016" || year == "UL16preVFP" || year == "UL16postVFP") {
      trg_ele_low.reset(new TriggerSelection("HLT_Ele27_WPTight_Gsf_v*"));
      trg_ele_high.reset(new TriggerSelection("HLT_Ele115_CaloIdVT_GsfTrkIdT_v*"));
      if(is_MC || data_isPhoton) trg_pho.reset(new TriggerSelection("HLT_Photon175_v*"));
    }
    else if(year == "2017" || year == "UL17") {
      trg_ele_low.reset(new TriggerSelection("HLT_Ele35_WPTight_Gsf_v*"));
      trg_ele_high.reset(new TriggerSelection("HLT_Ele115_CaloIdVT_GsfTrkIdT_v*"));
      if(is_MC || data_isPhoton) trg_pho.reset(new TriggerSelection("HLT_Photon200_v*"));
    }
    else if(year == "2018" || year == "UL18") {
      trg_ele_low.reset(new TriggerSelection("HLT_Ele32_WPTight_Gsf_v*"));
      trg_ele_high.reset(new TriggerSelection("HLT_Ele115_CaloIdVT_GsfTrkIdT_v*"));
      trg_pho.reset(new TriggerSelection("HLT_Photon200_v*"));
    }

  }

  // ###### 2. set up selections ######
  if(debug) cout << "Setting up Selections." << endl;

  // 2D cut
  twodcut_sel.reset(new TwoDCut(0.4, 25.0));  // The same as in Z'->ttbar semileptonic

  if(debug) cout << "Setting up HOTVR scale." << endl;

  // HOTVR scale
  HadronicTopFinder.reset(new HadronicTop(ctx));

  if(debug) cout << "Setting up electron scale." << endl;

  sf_ele_reco.reset( new uhh2::ElectronRecoScaleFactors(ctx) );
  sf_ele_reco_DUMMY.reset(new uhh2::ElectronRecoScaleFactors(ctx, boost::none, boost::none, boost::none, boost::none, true) );

  sf_ele_ID_lowpt.reset(new uhh2::ElectronIdScaleFactors(ctx, Electron::mvaEleID_Fall17_iso_V2_wp90) );
  sf_ele_ID_highpt.reset( new uhh2::ElectronIdScaleFactors(ctx, Electron::mvaEleID_Fall17_noIso_V2_wp90) );
  sf_ele_ID_DUMMY.reset( new uhh2::ElectronIdScaleFactors(ctx, boost::none, boost::none, boost::none, boost::none, true) );

  if(debug) cout << "Setting up muon scale." << endl;

  sf_muon_ID_lowpt.reset( new uhh2::MuonIdScaleFactors(ctx, Muon::CutBasedIdTight) );
  sf_muon_ID_highpt.reset( new uhh2::MuonIdScaleFactors(ctx, Muon::CutBasedIdGlobalHighPt) );
  sf_muon_ID_DUMMY.reset( new uhh2::MuonIdScaleFactors(ctx, boost::none, boost::none, boost::none, boost::none, true) );

  sf_muon_iso.reset( new uhh2::MuonIsoScaleFactors(ctx, Muon::PFIsoTight, Muon::CutBasedIdTight) ); // only for low pt
  sf_muon_iso_DUMMY.reset( new uhh2::MuonIsoScaleFactors(ctx, boost::none, boost::none, boost::none, boost::none, boost::none, true) );

  sf_muon_trigger_lowpt.reset( new uhh2::MuonTriggerScaleFactors(ctx, false) );
  sf_muon_trigger_highpt.reset( new uhh2::MuonTriggerScaleFactors(ctx, true, false) );
  sf_muon_trigger_DUMMY.reset( new uhh2::MuonTriggerScaleFactors(ctx, boost::none, boost::none, boost::none, boost::none, boost::none, true) );

  if(debug) cout << "Setting up btagging scale." << endl;

  // b-tagging SFs
  ScaleFactor_btagging.reset(new MCBTagDiscriminantReweighting(ctx, BTag::algo::DEEPJET)); // should be enough like this

  if(is_MC) {
    TFile *f = new TFile("/nfs/dust/cms/user/flabe/TstarTstar/ULegacy/CMSSW_10_6_28/src/UHH2/TstarTstar/macros/rootmakros/btagYieldSFs.root");
    TString sample_string = "";
    if(ctx.get("dataset_version").find("TT") != std::string::npos) sample_string = "TTbar";
    else if(ctx.get("dataset_version").find("ST") != std::string::npos) sample_string = "ST";
    else if(ctx.get("dataset_version").find("WJets") != std::string::npos) sample_string = "WJets";
    else if(ctx.get("dataset_version").find("QCD") != std::string::npos) sample_string = "QCD";
    else if(ctx.get("dataset_version").find("Diboson") != std::string::npos) sample_string = "VV";
    if(debug) std::cout << "Apply 2D b-taggin yield SFs for " << sample_string << std::endl;

    if(sample_string != "") eventYieldFactors = (TH2D*)f->Get(sample_string);
    else throw std::runtime_error("Error: can not determine sample type for btagging yield SFs.");
  }

  if(debug) cout << "Setting up HEM fix." << endl;

  // HEM issue
  HEMCleaner.reset(new HEMCleanerSelection(ctx, "jets", "topjets"));
  HEMCleanerMCScaler.reset(new HEMCleanerMCScale(ctx, "jets", "topjets"));

  if(debug) cout << "Setting up NNLO correction." << endl;

  // NNLO corrections
  ScaleFactor_NNLO.reset(new NLOCorrections(ctx));

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

  h_trigger.reset(new TstarTstarHists(ctx, "AfterTrigger"));
  h_trigger_gen.reset(new TstarTstarGenHists(ctx, "AfterTrigger_gen"));
  h_trigger_ele.reset(new TstarTstarHists(ctx, "AfterTrigger_ele"));
  h_trigger_ele_lowpt.reset(new TstarTstarHists(ctx, "AfterTrigger_ele_lowpt"));
  h_trigger_ele_highpt.reset(new TstarTstarHists(ctx, "AfterTrigger_ele_highpt"));
  h_trigger_mu.reset(new TstarTstarHists(ctx, "AfterTrigger_mu"));
  h_trigger_mu_lowpt.reset(new TstarTstarHists(ctx, "AfterTrigger_mu_lowpt"));
  h_trigger_mu_highpt.reset(new TstarTstarHists(ctx, "AfterTrigger_mu_highpt"));
  h_trigger_nobtag.reset(new TstarTstarHists(ctx, "AfterTrigger_nobtag"));

  h_triggerSF.reset(new TstarTstarHists(ctx, "AfterTriggerSF"));
  h_triggerSF_ele.reset(new TstarTstarHists(ctx, "AfterTriggerSF_ele"));
  h_triggerSF_ele_lowpt.reset(new TstarTstarHists(ctx, "AfterTriggerSF_ele_lowpt"));
  h_triggerSF_ele_highpt.reset(new TstarTstarHists(ctx, "AfterTriggerSF_ele_highpt"));
  h_triggerSF_mu.reset(new TstarTstarHists(ctx, "AfterTriggerSF_mu"));
  h_triggerSF_mu_lowpt.reset(new TstarTstarHists(ctx, "AfterTriggerSF_mu_lowpt"));
  h_triggerSF_mu_highpt.reset(new TstarTstarHists(ctx, "AfterTriggerSF_mu_highpt"));
  h_triggerSF_nobtag.reset(new TstarTstarHists(ctx, "AfterTriggerSF_nobtag"));

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

  h_topptreweighting.reset(new TstarTstarHists(ctx, "Aftertopptreweighting"));
  h_topptreweighting_ele.reset(new TstarTstarHists(ctx, "Aftertopptreweighting_ele"));
  h_topptreweighting_ele_lowpt.reset(new TstarTstarHists(ctx, "Aftertopptreweighting_ele_lowpt"));
  h_topptreweighting_ele_highpt.reset(new TstarTstarHists(ctx, "Aftertopptreweighting_ele_highpt"));
  h_topptreweighting_mu.reset(new TstarTstarHists(ctx, "Aftertopptreweighting_mu"));
  h_topptreweighting_mu_lowpt.reset(new TstarTstarHists(ctx, "Aftertopptreweighting_mu_lowpt"));
  h_topptreweighting_mu_highpt.reset(new TstarTstarHists(ctx, "Aftertopptreweighting_mu_highpt"));
  h_topptreweighting_nobtag.reset(new TstarTstarHists(ctx, "Aftertopptreweighting_nobtag"));

  h_afterNNLO.reset(new TstarTstarHists(ctx, "AfterNNLO"));
  h_afterNNLO_ele.reset(new TstarTstarHists(ctx, "AfterNNLO_ele"));
  h_afterNNLO_mu.reset(new TstarTstarHists(ctx, "AfterNNLO_mu"));

  h_beforeBcorrections.reset(new TstarTstarHists(ctx, "BeforeBCorrections"));
  h_afterBcorrections.reset(new TstarTstarHists(ctx, "AfterBCorrections"));
  h_afterBYieldcorrections.reset(new TstarTstarHists(ctx, "AfterBYieldCorrections"));

  h_afterHEMcleaning.reset(new TstarTstarHists(ctx, "AfterHEMcleaning"));
  h_afterHEMcleaning_ele.reset(new TstarTstarHists(ctx, "AfterHEMcleaning_ele"));
  h_afterHEMcleaning_mu.reset(new TstarTstarHists(ctx, "AfterHEMcleaning_mu"));

  h_notTriggered.reset(new TstarTstarHists(ctx, "notTriggered"));
  h_notTriggered_ele.reset(new TstarTstarHists(ctx, "notTriggered_ele"));
  h_notTriggered_mu.reset(new TstarTstarHists(ctx, "notTriggered_mu"));

  // ###### 4. Init Handles ######
  h_is_muevt = ctx.get_handle<bool>("is_muevt");
  h_is_highpt = ctx.get_handle<bool>("is_highpt");
  h_evt_weight = ctx.get_handle<double>("evt_weight");
  h_primlep = ctx.get_handle<FlavorParticle>("PrimaryLepton");
  h_neutrino = ctx.declare_event_output<LorentzVector>("neutrino");
  h_ST = ctx.declare_event_output<double>("ST");
  h_STHOTVR = ctx.declare_event_output<double>("STHOTVR");
  h_is_btagevent = ctx.declare_event_output<bool>("is_btagevent");

  h_MC_isfake2017B = ctx.get_handle<bool>("MC_isfake2017B");

  string IsTriggerSFMeasurement = ctx.get("IsTriggerSFMeasurement", "False");
  if(IsTriggerSFMeasurement == "True") isTriggerSFMeasurement = true;

  TopPtReweighting.reset( new TopPtReweight(ctx, 0.0615, -0.0005, "ttbargen", "weight_ttbar", true) );

}


bool TstarTstarSelectionModule::process(Event & event) {

  // debug status message
  if(debug) cout << "TstarTstarSelectionModule: Starting to process event (runid, eventid) = (" << event.run << ", " << event.event << "); weight = " << event.weight << endl;

  if(debug && event.get(h_MC_isfake2017B)) cout << "This MC event was randomly assigned to 2017 Run B handling" << endl;

  if(debug && event.get(h_is_highpt)) cout << "The event is classified as having high pt" << std::endl;

  // reapply event weights from handle
  event.weight = event.get(h_evt_weight);
  if(debug) cout << "weights applied." << endl;

  // Fill ttgen object for correct matching check, etc
  if(is_MC) ttgenprod->process(event);
  if(debug) std::cout << "Filled ttgenprod" << endl;

  // set primary lepton
  reco_primlep->process(event);
  if(debug) std::cout << "Got primary lepton" << endl;

  // hists before anything happened
  if(debug) std::cout << "Fill Crosscheck hists" << endl;
  h_beginSel->fill(event);
  h_beginSel_gen->fill(event);
  if(event.get(h_is_muevt)){
    h_beginSel_mu->fill(event);
    if(event.get(h_is_highpt)) h_beginSel_mu_highpt->fill(event);
    else h_beginSel_mu_lowpt->fill(event);
  }
  else {
    h_beginSel_ele->fill(event);
    if(event.get(h_is_highpt)) h_beginSel_ele_highpt->fill(event);
    else h_beginSel_ele_lowpt->fill(event);
  }

  // ###### Trigger selection ######
  bool pass_trigger = false;
  bool pass_trigger_SingleMu_lowpt = false;
  bool pass_trigger_SingleMu_highpt = false;
  bool pass_trigger_SingleEle_lowpt = false;
  bool pass_trigger_SingleEle_highpt = false;

  if( (is_MC && event.get(h_is_muevt)) || data_isMu ) {
    if(debug) std::cout << "Entered muon trigger logic" << std::endl;

    if(year == "2016" || year == "UL16preVFP" || year == "UL16postVFP") {
      pass_trigger_SingleMu_lowpt = (trg_mu_low_1->passes(event) || trg_mu_low_2->passes(event));
      pass_trigger_SingleMu_highpt = (trg_mu_high_1->passes(event) || trg_mu_high_2->passes(event));
    }
    else if(year == "2017" || year == "UL17") {
      pass_trigger_SingleMu_lowpt = trg_mu_low_1->passes(event);
      if(data_is2017B || event.get(h_MC_isfake2017B)) pass_trigger_SingleMu_highpt = trg_mu_high_1->passes(event);
      else pass_trigger_SingleMu_highpt = (trg_mu_high_1->passes(event) || trg_mu_high_2->passes(event) || trg_mu_high_3->passes(event));
    }
    else if(year == "2018" || year == "UL18") {
      pass_trigger_SingleMu_lowpt = trg_mu_low_1->passes(event);
      pass_trigger_SingleMu_highpt = (trg_mu_high_1->passes(event) || trg_mu_high_2->passes(event) || trg_mu_high_3->passes(event));
    }

    if(debug) std::cout << "Passed muon trigger logic" << std::endl;
  } else if ( (is_MC && !event.get(h_is_muevt)) || !data_isMu ){
    if(debug) std::cout << "Entered electron trigger logic" << std::endl;

    if(year == "2016" || year == "UL16preVFP" || year == "UL16postVFP") {

      if(!data_isPhoton) pass_trigger_SingleEle_lowpt = trg_ele_low->passes(event);

      if(is_MC) {
        pass_trigger_SingleEle_highpt = (trg_ele_high->passes(event) || trg_pho->passes(event));
      } else {
        if (data_isPhoton) pass_trigger_SingleEle_highpt = (!trg_ele_high->passes(event) && trg_pho->passes(event));
        else pass_trigger_SingleEle_highpt = trg_ele_high->passes(event);
      }

    }
    else if(year == "2017" || year == "UL17") {
      if(!data_isPhoton) pass_trigger_SingleEle_lowpt = trg_ele_low->passes(event);

      if(is_MC) {
        if (event.get(h_MC_isfake2017B)) pass_trigger_SingleEle_highpt = (trg_ele_low->passes(event) || trg_pho->passes(event));
        else pass_trigger_SingleEle_highpt = (trg_ele_high->passes(event) || trg_pho->passes(event));
      } else {
        if (data_is2017B) {
          if(data_isPhoton) pass_trigger_SingleEle_highpt = (!trg_ele_low->passes(event) && trg_pho->passes(event));
          else pass_trigger_SingleEle_highpt = trg_ele_low->passes(event);
        } else {
          if(data_isPhoton) pass_trigger_SingleEle_highpt = (!trg_ele_high->passes(event) && trg_pho->passes(event));
          else pass_trigger_SingleEle_highpt = trg_ele_high->passes(event);
        }
      }

    }
    else if(year == "2018" || year == "UL18") {
      pass_trigger_SingleEle_lowpt = trg_ele_low->passes(event);
      pass_trigger_SingleEle_highpt = (trg_ele_high->passes(event) || trg_pho->passes(event));
    }

    if(debug) std::cout << "Passed electron trigger logic" << std::endl;
  } else {
    throw std::runtime_error("Error: event not fitting any trigger group.");
  }

  // main logic
  if( (is_MC && event.get(h_is_muevt)) || data_isMu ) {
    if(event.get(h_is_highpt)) pass_trigger = pass_trigger_SingleMu_highpt;
    else pass_trigger = pass_trigger_SingleMu_lowpt;
  } else {
    if(event.get(h_is_highpt)) pass_trigger = pass_trigger_SingleEle_highpt;
    else pass_trigger = pass_trigger_SingleEle_lowpt;
  }
  //if(isTriggerSFMeasurement) pass_trigger = pass_trigger_SingleMu_lowpt || pass_trigger_SingleMu_highpt;

  if(!pass_trigger) {
    h_notTriggered->fill(event);
    if(event.get(h_is_muevt)) h_notTriggered_mu->fill(event);
    else h_notTriggered_ele->fill(event);
    if(debug) cout<<"Filled hists for not triggered"<<endl;


    return false;
  }

  {
    // hists
    h_trigger->fill(event);
    h_trigger_gen->fill(event);
    if(event.get(h_is_muevt)){
      h_trigger_mu->fill(event);
      if(event.get(h_is_highpt)) h_trigger_mu_highpt->fill(event);
      else h_trigger_mu_lowpt->fill(event);
    }
    else {
      h_trigger_ele->fill(event);
      if(event.get(h_is_highpt)) h_trigger_ele_highpt->fill(event);
      else h_trigger_ele_lowpt->fill(event);
    }
  }

  // trigger SFs
  if(is_MC && event.get(h_is_muevt) && !isTriggerSFMeasurement){

    if(event.get(h_is_highpt)) {
      sf_muon_trigger_highpt->process(event);
    } else {
      sf_muon_trigger_lowpt->process(event);
    }

  } else {
    sf_muon_trigger_DUMMY->process(event);
  }

  {
    // hists
    h_triggerSF->fill(event);
    if(event.get(h_is_muevt)){
      h_triggerSF_mu->fill(event);
      if(event.get(h_is_highpt)) h_triggerSF_mu_highpt->fill(event);
      else h_triggerSF_mu_lowpt->fill(event);
    }
    else {
      h_triggerSF_ele->fill(event);
      if(event.get(h_is_highpt)) h_triggerSF_ele_highpt->fill(event);
      else h_triggerSF_ele_lowpt->fill(event);
    }
  }

  h_beforeBcorrections->fill(event);

  // b-tagging sfs
  ScaleFactor_btagging->process(event);

  h_afterBcorrections->fill(event);

  if(is_MC) {
    double ht = 0.;
    for(const auto & jet : *event.jets) ht += jet.pt();
    if(ht >= 4000.) ht = 3999.9;

    double btaggingYieldWeight = eventYieldFactors->GetBinContent( eventYieldFactors->GetXaxis()->FindBin(ht),  eventYieldFactors->GetYaxis()->FindBin(event.jets->size()) );
    event.weight *= btaggingYieldWeight;
  }

  h_afterBYieldcorrections->fill(event);

  // #################
  // ### Selection ###
  // #################
  if(debug) std::cout << "Starting selection" << endl;

  // ###### Btag Selection ######
  BTag bJetID = BTag(BTag::algo::DEEPJET, BTag::wp::WP_LOOSE);
  bool pass_btagcut = false;
  for (const auto & jet: *event.jets){
    if(bJetID(jet, event)) pass_btagcut = true;
  }
  event.set(h_is_btagevent,pass_btagcut);
  if(pass_btagcut) {
    // hists
    h_btagcut->fill(event);
    h_btagcut_gen->fill(event);
    if(event.get(h_is_muevt)){
      h_btagcut_mu->fill(event);
      if(event.get(h_is_highpt)) h_btagcut_mu_highpt->fill(event);
      else h_btagcut_mu_lowpt->fill(event);
    }
    else {
      h_btagcut_ele->fill(event);
      if(event.get(h_is_highpt)) h_btagcut_ele_highpt->fill(event);
      else h_btagcut_ele_lowpt->fill(event);
    }
    if(debug) cout << "Passed btag selection cut." << endl;
  } // end of pass_btag = true
  else {
    // TODO control region plots maybe here?
  }

  if(debug) std::cout << "Done b-tag" << endl;

  // ###### Lepton-2Dcut ######
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
  if(event.get(h_is_highpt)) pass_2D = pass_twodcut;
  if(!pass_2D) return false;

  // hists
  h_2Dcut->fill(event);
  h_2Dcut_gen->fill(event);
  if(event.get(h_is_muevt)){
    h_2Dcut_mu->fill(event);
    if(event.get(h_is_highpt)) h_2Dcut_mu_highpt->fill(event);
    else h_2Dcut_mu_lowpt->fill(event);
  }
  else {
    h_2Dcut_ele->fill(event);
    if(event.get(h_is_highpt)) h_2Dcut_ele_highpt->fill(event);
    else h_2Dcut_ele_lowpt->fill(event);
  }
  if(debug) cout << "Passed 2D cut." << endl;

  if(debug) std::cout << "We have muons: " << event.muons->size() << endl;
  if(debug) std::cout << "We have electrons: " << event.electrons->size() << endl;
  if(debug) std::cout << "is_muevt is: " << event.get(h_is_muevt) << endl;

  if(debug) std::cout << "Done dR" << endl;

  // ST cut to reduce computation time (and file sizes)
  // Neutrino reconstruction
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
  for(const auto & lepton : *event.electrons) st += lepton.pt();
  for(const auto & lepton : *event.muons) st += lepton.pt();
  st += event.met->pt();
  double stHOTVR = st;

  for(const auto & jet : *event.jets) st += jet.pt();
  for(const auto & jet : *event.topjets) stHOTVR += jet.pt();

  event.set(h_ST, st);
  event.set(h_STHOTVR, stHOTVR);

  // st cut
  if(st < 500) return false;

  if(pass_btagcut) { // only fill these for btag cut passes
    // hists
    h_STcut->fill(event);
    h_STcut_gen->fill(event);
    if(event.get(h_is_muevt)){
      h_STcut_mu->fill(event);
      if(event.get(h_is_highpt)) h_STcut_mu_highpt->fill(event);
      else h_STcut_mu_lowpt->fill(event);
    }
    else {
      h_STcut_ele->fill(event);
      if(event.get(h_is_highpt)) h_STcut_ele_highpt->fill(event);
      else h_STcut_ele_lowpt->fill(event);
    }
    if(debug) cout << "Passed ST cut." << endl;
  } else {
    h_STcut_nobtag->fill(event);
    if(!event.get(h_is_muevt)) h_STcut_nobtag_ele->fill(event);
  }

  if(debug) std::cout << "Done ST" << endl;

  // Fixing the HEM issue
  if(!HEMCleaner->passes(event)) return false;
  if(is_MC) HEMCleanerMCScaler->process(event);

  if(pass_btagcut) {
    h_afterHEMcleaning->fill(event);
    if(event.get(h_is_muevt)) h_afterHEMcleaning_mu->fill(event);
    else h_afterHEMcleaning_ele->fill(event);
    if(debug) cout<<"Filled hists HEM fix"<<endl;
  }

  // lepton SFs
  if(event.get(h_is_muevt)){

    if(event.get(h_is_highpt)) {
      sf_muon_ID_highpt->process(event);
      sf_muon_iso_DUMMY->process(event);
    } else {
      sf_muon_ID_lowpt->process(event);
      sf_muon_iso->process(event);
    }

    sf_ele_reco_DUMMY->process(event);
    sf_ele_ID_DUMMY->process(event);

  } else {

    sf_ele_reco->process(event);

    if(event.get(h_is_highpt)) {
      sf_ele_ID_highpt->process(event);
    } else {
      sf_ele_ID_lowpt->process(event);
    }

    sf_muon_ID_DUMMY->process(event);
    sf_muon_iso_DUMMY->process(event);

  }

  if(debug) std::cout << "Done Lepton ID, ISO SFs" << endl;

  if(pass_btagcut) { // only fill these for btag cut passes
    // hists
    h_corrections->fill(event);
    h_corrections_gen->fill(event);
    if(event.get(h_is_muevt)){
      h_corrections_mu->fill(event);
      if(event.get(h_is_highpt)) h_corrections_mu_highpt->fill(event);
      else h_corrections_mu_lowpt->fill(event);
    }
    else {
      h_corrections_ele->fill(event);
      if(event.get(h_is_highpt)) h_corrections_ele_highpt->fill(event);
      else h_corrections_ele_lowpt->fill(event);
    }
  } else {
    h_corrections_nobtag->fill(event);
    if(!event.get(h_is_muevt)) h_corrections_nobtag_ele->fill(event);
  }

  /**
  // NNLO corrections
  ScaleFactor_NNLO->process(event);

  if(pass_btagcut) {
    h_afterNNLO->fill(event);
    if(event.get(h_is_muevt)) h_afterNNLO_mu->fill(event);
    else h_afterNNLO_ele->fill(event);
    if(debug) cout<<"Filled hists HEM fix"<<endl;
  }

  if(debug) cout << "NNLO corrections done" << std::endl;

  TopPtReweighting->process(event);

  if(pass_btagcut) { // only fill these for btag cut passes
    // hists
    h_topptreweighting->fill(event);
    if(event.get(h_is_muevt)){
      h_topptreweighting_mu->fill(event);
      if(event.get(h_is_highpt)) h_topptreweighting_mu_highpt->fill(event);
      else h_topptreweighting_mu_lowpt->fill(event);
    }
    else {
      h_topptreweighting_ele->fill(event);
      if(event.get(h_is_highpt)) h_topptreweighting_ele_highpt->fill(event);
      else h_topptreweighting_ele_lowpt->fill(event);
    }
  } else {
    h_topptreweighting_nobtag->fill(event);
  }
  **/

  return true;

}

// as we want to run the ExampleCycleNew directly with AnalysisModuleRunner,
// make sure the TstarTstarSelectionModule is found by class name. This is ensured by this macro:
UHH2_REGISTER_ANALYSIS_MODULE(TstarTstarSelectionModule)

}
