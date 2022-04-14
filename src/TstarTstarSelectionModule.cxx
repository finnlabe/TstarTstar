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

// TstarTstar stuff
#include "UHH2/TstarTstar/include/ModuleBASE.h"
#include "UHH2/TstarTstar/include/TstarTstarSelections.h"
#include "UHH2/TstarTstar/include/TstarTstarHists.h"
#include "UHH2/TstarTstar/include/TstarTstarGenHists.h"
#include "UHH2/TstarTstar/include/TstarTstarGenRecoMatchedHists.h"
#include "UHH2/TstarTstar/include/TstarTstarRecoTstarHists.h"
#include "UHH2/TstarTstar/include/TstarTstarReconstructionModules.h"
#include "UHH2/TstarTstar/include/ElecTriggerSF.h"

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

  // I'll add an "dummy" version here to be used on the "wrong" channel
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
  std::unique_ptr<Hists> h_beginSel,            h_btagcut,             h_2Dcut,               h_dRcut,               h_STcut            , h_corrections,              h_triggercorrections;
  std::unique_ptr<Hists> h_beginSel_gen,        h_btagcut_gen,         h_2Dcut_gen,           h_dRcut_gen,           h_STcut_gen        , h_corrections_gen,          h_triggercorrections_gen;
  std::unique_ptr<Hists> h_beginSel_ele,        h_btagcut_ele,         h_2Dcut_ele,           h_dRcut_ele,           h_STcut_ele        , h_corrections_ele,          h_triggercorrections_ele;
  std::unique_ptr<Hists> h_beginSel_ele_lowpt,  h_btagcut_ele_lowpt,   h_2Dcut_ele_lowpt,     h_dRcut_ele_lowpt,     h_STcut_ele_lowpt  , h_corrections_ele_lowpt,    h_triggercorrections_ele_lowpt;
  std::unique_ptr<Hists> h_beginSel_ele_highpt, h_btagcut_ele_highpt,  h_2Dcut_ele_highpt,    h_dRcut_ele_highpt,    h_STcut_ele_highpt , h_corrections_ele_highpt,   h_triggercorrections_ele_highpt;
  std::unique_ptr<Hists> h_beginSel_mu,         h_btagcut_mu,          h_2Dcut_mu,            h_dRcut_mu,            h_STcut_mu         , h_corrections_mu,           h_triggercorrections_mu;
  std::unique_ptr<Hists> h_beginSel_mu_lowpt,   h_btagcut_mu_lowpt,    h_2Dcut_mu_lowpt,      h_dRcut_mu_lowpt,      h_STcut_mu_lowpt   , h_corrections_mu_lowpt,     h_triggercorrections_mu_lowpt;
  std::unique_ptr<Hists> h_beginSel_mu_highpt,  h_btagcut_mu_highpt,   h_2Dcut_mu_highpt,     h_dRcut_mu_highpt,     h_STcut_mu_highpt  , h_corrections_mu_highpt,    h_triggercorrections_mu_highpt;
  std::unique_ptr<Hists> h_beginSel_nobtag,     h_btagcut_nobtag,      h_2Dcut_nobtag,        h_dRcut_nobtag,        h_STcut_nobtag     , h_corrections_nobtag,       h_triggercorrections_nobtag;
  std::unique_ptr<Hists> h_STcut_nobtag_ele, h_corrections_nobtag_ele,  h_triggercorrections_nobtag_ele, h_triggercorrections_nobtag_mu;

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
  sf_muon_trigger_highpt.reset( new uhh2::MuonTriggerScaleFactors(ctx, true) );
  sf_muon_trigger_DUMMY.reset( new uhh2::MuonTriggerScaleFactors(ctx, boost::none, boost::none, boost::none, boost::none, boost::none, true) );

  if(debug) cout << "Setting up btagging scale." << endl;

  // b-tagging SFs
  ScaleFactor_btagging.reset(new MCBTagDiscriminantReweighting(ctx, BTag::algo::DEEPJET)); // should be enough like this

  if(debug) cout << "Setting up HEM fix." << endl;

  // HEM issue
  HEMCleaner.reset(new HEMCleanerSelection(ctx, "jets", "topjets"));


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

  h_triggercorrections.reset(new TstarTstarHists(ctx, "AfterTriggerSFs"));
  h_triggercorrections_gen.reset(new TstarTstarGenHists(ctx, "AfterTriggerSFs_gen"));
  h_triggercorrections_ele.reset(new TstarTstarHists(ctx, "AfterTriggerSFs_ele"));
  h_triggercorrections_ele_lowpt.reset(new TstarTstarHists(ctx, "AfterTriggerSFs_ele_lowpt"));
  h_triggercorrections_ele_highpt.reset(new TstarTstarHists(ctx, "AfterTriggerSFs_ele_highpt"));
  h_triggercorrections_mu.reset(new TstarTstarHists(ctx, "AfterTriggerSFs_mu"));
  h_triggercorrections_mu_lowpt.reset(new TstarTstarHists(ctx, "AfterTriggerSFs_mu_lowpt"));
  h_triggercorrections_mu_highpt.reset(new TstarTstarHists(ctx, "AfterTriggerSFs_mu_highpt"));
  h_triggercorrections_nobtag.reset(new TstarTstarHists(ctx, "AfterTriggerSFs_nobtag"));
  h_triggercorrections_nobtag_ele.reset(new TstarTstarHists(ctx, "AfterTriggerSFs_nobtag_ele"));
  h_triggercorrections_nobtag_mu.reset(new TstarTstarHists(ctx, "AfterTriggerSFs_nobtag_mu"));


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

  if(is_triggered && !event.get(h_is_muevt)) h_beforeBcorrections->fill(event);

  // b-tagging sfs
  ScaleFactor_btagging->process(event);

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
  if(event.muons->size()==1){if(event.muons->at(0).pt()>60) pass_2D = pass_twodcut;}
  else if(event.electrons->size()==1){if(event.electrons->at(0).pt()>120) pass_2D = pass_twodcut;}
  else std::cout << "How did this happen???" << endl;
  if(!pass_2D) return false;

  // hists
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

  // ######################
  // ### Scale Factors! ###
  // ######################

  // HOTVR
  HadronicTopFinder->process(event);
  //HOTVRScale->process(event);

  if(debug) std::cout << "Done HOTVR scale" << endl;

  // lepton SFs
  if(event.get(h_is_muevt)){

    if(event.get(h_is_highpt)) {
      sf_muon_ID_highpt->process(event);
      sf_muon_iso_DUMMY->process(event);
      sf_muon_trigger_lowpt->process(event);
    } else {
      sf_muon_iso->process(event);
      sf_muon_ID_lowpt->process(event);
      sf_muon_trigger_highpt->process(event);
    }

    sf_ele_reco_DUMMY->process(event);
    sf_ele_ID_DUMMY->process(event);

  }
  else {

    sf_ele_reco->process(event);

    if(event.get(h_is_highpt)) {
      sf_ele_ID_highpt->process(event);
    } else {
      sf_ele_ID_highpt->process(event);
    }

    sf_muon_ID_DUMMY->process(event);
    sf_muon_iso_DUMMY->process(event);
    sf_muon_trigger_DUMMY->process(event);

  }

  if(debug) std::cout << "Done Lepton ID, ISO SFs" << endl;

  if(pass_btagcut && is_triggered) { // only fill these for btag cut passes
    // hists
    h_corrections->fill(event);
    h_corrections_gen->fill(event);
    if(event.get(h_is_muevt)){
      h_corrections_mu->fill(event);
      if(event.muons->at(0).pt()<=60) h_corrections_mu_lowpt->fill(event);
      else h_corrections_mu_highpt->fill(event);
    }
    else {
      h_corrections_ele->fill(event);
      if(event.electrons->at(0).pt()<=120) h_corrections_ele_lowpt->fill(event);
      else h_corrections_ele_highpt->fill(event);
    }
    if(debug) cout << "Passed ST cut." << endl;
  } else {
    if(is_triggered) {
      h_corrections_nobtag->fill(event);
      if(!event.get(h_is_muevt)) h_corrections_nobtag_ele->fill(event);
    }
  }

  if(pass_btagcut && is_triggered) { // only fill these for btag cut passes
    // hists
    h_triggercorrections->fill(event);
    h_triggercorrections_gen->fill(event);
    if(event.get(h_is_muevt)){
      h_triggercorrections_mu->fill(event);
      if(event.muons->at(0).pt()<=60) h_triggercorrections_mu_lowpt->fill(event);
      else h_triggercorrections_mu_highpt->fill(event);
    }
    else {
      h_triggercorrections_ele->fill(event);
      if(event.electrons->at(0).pt()<=120) h_triggercorrections_ele_lowpt->fill(event);
      else h_triggercorrections_ele_highpt->fill(event);
    }
    if(debug) cout << "Passed ST cut." << endl;
  } else {
    if(is_triggered) {
      h_triggercorrections_nobtag->fill(event);
      if(event.get(h_is_muevt)) h_triggercorrections_ele->fill(event);
      else h_triggercorrections_ele->fill(event);
    }
  }

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

  return true;

}

// as we want to run the ExampleCycleNew directly with AnalysisModuleRunner,
// make sure the TstarTstarSelectionModule is found by class name. This is ensured by this macro:
UHH2_REGISTER_ANALYSIS_MODULE(TstarTstarSelectionModule)

}
