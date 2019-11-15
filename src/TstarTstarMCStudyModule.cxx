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
#include <UHH2/common/include/MuonIds.h>
#include <UHH2/common/include/TriggerSelection.h>
#include "UHH2/TstarTstar/include/TstarTstarSelections.h"
#include "UHH2/TstarTstar/include/TstarTstarHists.h"
#include "UHH2/TstarTstar/include/TstarTstarRecoTstarHists.h"
#include "UHH2/TstarTstar/include/TstarTstarGenHists.h"
#include "UHH2/TstarTstar/include/TstarTstarGenRecoMatchedHists.h"
#include "UHH2/TstarTstar/include/TstarTstarReconstructionModules.h"

/**
#include "UHH2/common/include/TTbarReconstruction.h"
#include "UHH2/common/include/ReconstructionHypothesis.h"
#include "UHH2/common/include/ReconstructionHypothesisDiscriminators.h"
**/

using namespace std;
using namespace uhh2;

namespace uhh2 {

/** \brief Module for the T*T*->ttbar gg MC based study  
 *  
 * This is the central class which calls other AnalysisModules, Hists or Selection classes.
 * This AnalysisModule, in turn, is called (via AnalysisModuleRunner) by SFrame.
 */
class TstarTstarMCStudyModule: public AnalysisModule {
public:
    
  explicit TstarTstarMCStudyModule(Context & ctx);
  virtual bool process(Event & event) override;
  
private:
    
  // Apply common modules: JetPFid, JEC, JER, MET corrections, etc
  std::unique_ptr<CommonModules> common;
  
  // Declare the Selections to use. Use unique_ptr to ensure automatic call of delete in the destructor,
  // to avoid memory leaks.
  unique_ptr<Selection> twodcut_sel;   
  //  unique_ptr<Selection> TTbarSemiLepMatchable_selection;
  unique_ptr<TTbarSemiLepMatchableSelection> TTbarSemiLepMatchable_selection;

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
  // Store the Hists collection as member variables. Again, use unique_ptr to avoid memory leaks.
  std::unique_ptr<Hists> h_nocuts, h_common, h_lepsel, h_njetsel, h_nphosel, h_2dcut, h_semilepttbarmatch, h_nosemilepttbarmatch;
  std::unique_ptr<Hists> h_semilepttbarmatch_triggerSingleJet, h_semilepttbarmatch_triggerSingleLeptonMu, h_semilepttbarmatch_triggerSingleLeptonEle, h_semilepttbarmatch_triggerHT, h_semilepttbarmatch_triggerPFHT;
  std::unique_ptr<Hists> h_semilepttbarmatch_gen;
  std::unique_ptr<Hists> h_semilepttbarmatch_genreco,h_semilepttbarmatch_genreco_mu,h_semilepttbarmatch_genreco_ele;
  std::unique_ptr<Hists> h_semilepttbarmatch_triggerSingleJet_genreco, h_semilepttbarmatch_triggerSingleLeptonMu_genreco, h_semilepttbarmatch_triggerSingleLeptonEle_genreco, h_semilepttbarmatch_triggerHT_genreco, h_semilepttbarmatch_triggerPFHT_genreco;
  std::unique_ptr<Hists> h_semilepttbarmatch_triggerSingleJet_genreco_mu, h_semilepttbarmatch_triggerHT_genreco_mu, h_semilepttbarmatch_triggerPFHT_genreco_mu;
  std::unique_ptr<Hists> h_semilepttbarmatch_triggerSingleJet_genreco_ele, h_semilepttbarmatch_triggerHT_genreco_ele, h_semilepttbarmatch_triggerPFHT_genreco_ele;
  std::unique_ptr<Hists> h_After_TstarTstar_Reco, h_After_TstarTstar_Reco_match;
  std::unique_ptr<TstarTstarRecoTstarHists> h_RecoPlots_After_ttbar,h_RecoPlots_After_ttbar_correct_ttbar, h_RecoPlots_After_TstarTstar, h_RecoPlots_After_TstarTstar_match;
  std::unique_ptr<TstarTstarRecoTstarHists>  h_RecoPlots_After_TstarTstar_tgtg,h_RecoPlots_After_TstarTstar_tgtg_correct_ttbar;
  bool isTrigger = false;
  //  bool debug = true;
  bool debug = false;
  std::unique_ptr<uhh2::AnalysisModule> reco_primlep;
  std::unique_ptr<uhh2::AnalysisModule> ttbar_reco;
  std::unique_ptr<ttbarChi2Discriminator> ttbar_discriminator;
  std::unique_ptr<TstarTstar_Reconstruction> TstarTstar_reco;
  std::unique_ptr<TstarTstar_tgluon_tgluon_Reconstruction> TstarTstar_tgluon_tgluon_reco;
  uhh2::Event::Handle<std::vector<ReconstructionHypothesis>> h_ttbar_hyps;
  uhh2::Event::Handle<bool> h_is_ttbar_reconstructed;
  uhh2::Event::Handle<ReconstructionHypothesis> h_recohyp;

  uhh2::Event::Handle<std::vector<ReconstructionTstarHypothesis>> h_tstartstar_hyps;
  //  uhh2::Event::Handle<bool> h_is_tstartstar_reconstructed;
  uhh2::Event::Handle<ReconstructionTstarHypothesis> h_recohyp_tstartstar;


  // uhh2::Event::Handle<float> h_M_Tstar_gluon;
  // uhh2::Event::Handle<float> h_M_Tstar_gamma;
  // uhh2::Event::Handle<float> h_M_Tstar_lep;
  // uhh2::Event::Handle<float> h_M_Tstar_had;
  bool is_tgtg, is_tgtgamma;
  bool check_ttbar_reco;

};

TstarTstarMCStudyModule::TstarTstarMCStudyModule(Context & ctx){
    
 
  if(debug) {
    cout << "Hello World from TstarTstarMCStudyModule!" << endl;
    
    
    // If running in SFrame, the keys "dataset_version", "dataset_type", "dataset_lumi",
    // and "target_lumi" are set to the according values in the xml file. For CMSSW, these are
    // not set automatically, but can be set in the python config file.
    for(auto & kv : ctx.get_all()){
        cout << " " << kv.first << " = " << kv.second << endl;
    }
    
   }
  check_ttbar_reco = false;
  is_tgtg = false; is_tgtgamma = false;
  if(ctx.get("channel") == "tgtg") is_tgtg = true;
  if(ctx.get("channel") == "tgtgamma") is_tgtgamma = true;

    // 1. setup other modules. CommonModules
    common.reset(new CommonModules());
    common->disable_mcpileupreweight();//FixME: PU re-weighting crushes
    common->switch_metcorrection();
    //    common->disable_metfilters();//FixME: met filters missing in some Z' samples, remove after tests with them are done
    common->switch_jetlepcleaner();
    common->switch_jetPtSorter();
    common->set_jet_id(AndId<Jet>(JetPFID(JetPFID::WP_TIGHT_PUPPI), PtEtaCut(30.0,5.2)));
    ElectronId eleID; 
    double electron_pt(20.);
    eleID = ElectronID_Summer16_tight_noIso;
    common->set_electron_id(AndId<Electron>(PtEtaSCCut(electron_pt, 2.5), eleID));

    
    PhotonId phoID; 
    double photon_pt(20.);
    // read in desired photonID from config file
    if (ctx.get("PhotonID") == "cutBasedPhotonIDlooseFall17"){phoID = PhotonTagID(Photon::cutBasedPhotonID_Fall17_94X_V2_loose);}
    else if (ctx.get("PhotonID") == "cutBasedPhotonIDmediumFall17"){phoID = PhotonTagID(Photon::cutBasedPhotonID_Fall17_94X_V2_medium);}
    else if (ctx.get("PhotonID") == "cutBasedPhotonIDtightFall17"){phoID = PhotonTagID(Photon::cutBasedPhotonID_Fall17_94X_V2_tight);}
    else if (ctx.get("PhotonID") == "mvaPhoIDwp90Fall17"){phoID = PhotonTagID(Photon::mvaPhoID_Fall17_iso_V2_wp90);}
    else if (ctx.get("PhotonID") == "mvaPhoIDwp80Fall17"){phoID = PhotonTagID(Photon::mvaPhoID_Fall17_iso_V2_wp80);}
    else {cout << "error: unrecognized photonID "  << ctx.get("PhotonID") << endl;}

    if (ctx.get("PhotonID") != "noPhotonID") {common->set_photon_id(AndId<Photon>(PtEtaCut(photon_pt, 5.2), phoID));}
    else {common->set_photon_id(PtEtaCut(photon_pt, 5.2));}
    

    MuonId muID;
    double muon_pt(20.);
    muID = MuonID(Muon::Highpt);
    //    muID = MuonID(Muon::CutBasedIdTight);
    common->set_muon_id(AndId<Muon>(PtEtaCut(muon_pt, 2.4), muID));
    common->init(ctx);
    
    // 2. set up selections
    ///2D Cut Lepton-Jets
    twodcut_sel.reset(new TwoDCut(0.4, 25.0));  // The same as in Z'->ttbar semileptonic
    TTbarSemiLepMatchable_selection.reset(new TTbarSemiLepMatchableSelection());

    isTrigger = (ctx.get("UseTrigger") == "true");
    //isTrigger =  false; //TEST
    triggerSingleJet450_sel.reset(new TriggerSelection("HLT_PFJet450_v*"));
    //    triggerSingleLeptonEle_sel.reset(new TriggerSelection("HLT_Ele32_eta2p1_WPTight_Gsf_v*"));
    // triggerSingleLeptonMu1_sel.reset(new TriggerSelection("HLT_IsoMu24_v*"));
    // triggerSingleLeptonMu2_sel.reset(new TriggerSelection("HLT_IsoTkMu24_v*"));

    //    triggerSingleLeptonEle_sel.reset(new TriggerSelection("HLT_Ele115_CaloIdVT_GsfTrkIdT_v*"));
    triggerSingleLeptonEle1_sel.reset(new TriggerSelection("HLT_Ele115_CaloIdVT_GsfTrkIdT_v*"));
    triggerSingleLeptonEle2_sel.reset(new TriggerSelection("HLT_Ele25_eta2p1_WPTight_Gsf_v*"));
    triggerSingleLeptonEle3_sel.reset(new TriggerSelection("HLT_Ele32_eta2p1_WPTight_Gsf_v*"));

    triggerSingleLeptonMu1_sel.reset(new TriggerSelection("HLT_Mu50_v*"));
    triggerSingleLeptonMu2_sel.reset(new TriggerSelection("HLT_Mu55_v*"));
    triggerSingleLeptonMu3_sel.reset(new TriggerSelection("HLT_IsoMu24_v*"));
    triggerSingleLeptonMu4_sel.reset(new TriggerSelection("HLT_IsoTkMu24_v*"));

    triggerHT1_sel.reset(new TriggerSelection("HLT_HT430to450_v*"));
    triggerHT2_sel.reset(new TriggerSelection("HLT_HT450to470_v*"));
    triggerHT3_sel.reset(new TriggerSelection("HLT_HT470to500_v*"));
    triggerHT4_sel.reset(new TriggerSelection("HLT_HT500to550_v*"));
    triggerHT5_sel.reset(new TriggerSelection("HLT_HT550to650_v*"));
    triggerHT6_sel.reset(new TriggerSelection("HLT_HT650_v*"));

    triggerPFHT_sel.reset(new TriggerSelection("HLT_PFHT900_v*"));



    // 3. Set up Hists classes:
    h_nocuts.reset(new TstarTstarHists(ctx, "NoCuts"));
    h_common.reset(new TstarTstarHists(ctx, "AfterCommon"));
    h_njetsel.reset(new TstarTstarHists(ctx, "AfterNjets"));
    h_lepsel.reset(new TstarTstarHists(ctx, "AfterLepSel"));
    h_nphosel.reset(new TstarTstarHists(ctx, "AfterNpho"));
    h_2dcut.reset(new TstarTstarHists(ctx, "After2D"));
    h_semilepttbarmatch.reset(new TstarTstarHists(ctx, "SemiLepTTBarMatch"));

    h_semilepttbarmatch_triggerSingleJet.reset(new TstarTstarHists(ctx, "SemiLepTTBarMatch_triggerSingleJet"));
    h_semilepttbarmatch_triggerSingleLeptonMu.reset(new TstarTstarHists(ctx, "SemiLepTTBarMatch_triggerSingleLeptonMu"));
    h_semilepttbarmatch_triggerSingleLeptonEle.reset(new TstarTstarHists(ctx, "SemiLepTTBarMatch_triggerSingleLeptonEle"));
    h_semilepttbarmatch_triggerHT.reset(new TstarTstarHists(ctx, "SemiLepTTBarMatch_triggerHT"));
    h_semilepttbarmatch_triggerPFHT.reset(new TstarTstarHists(ctx, "SemiLepTTBarMatch_triggerPFHT"));

    h_nosemilepttbarmatch.reset(new TstarTstarHists(ctx, "NotSemiLepTTBarMatch"));
    h_semilepttbarmatch_gen.reset(new TstarTstarGenHists(ctx, "SemiLepTTBarMatchGEN"));
    h_semilepttbarmatch_genreco.reset(new TstarTstarGenRecoMatchedHists(ctx, "SemiLepTTBarMatchGENRECO"));
    h_semilepttbarmatch_genreco_mu.reset(new TstarTstarGenRecoMatchedHists(ctx, "SemiLepTTBarMatchGENRECO_mu"));
    h_semilepttbarmatch_genreco_ele.reset(new TstarTstarGenRecoMatchedHists(ctx, "SemiLepTTBarMatchGENRECO_ele"));

    h_semilepttbarmatch_triggerSingleJet_genreco.reset(new TstarTstarGenRecoMatchedHists(ctx, "SemiLepTTBarMatchGENRECO_triggerSingleJet"));
    h_semilepttbarmatch_triggerHT_genreco.reset(new TstarTstarGenRecoMatchedHists(ctx, "SemiLepTTBarMatchGENRECO_triggerHT"));
    h_semilepttbarmatch_triggerSingleJet_genreco_mu.reset(new TstarTstarGenRecoMatchedHists(ctx, "SemiLepTTBarMatchGENRECO_triggerSingleJet_mu"));
    h_semilepttbarmatch_triggerHT_genreco_mu.reset(new TstarTstarGenRecoMatchedHists(ctx, "SemiLepTTBarMatchGENRECO_triggerHT_mu"));
    h_semilepttbarmatch_triggerPFHT_genreco_mu.reset(new TstarTstarGenRecoMatchedHists(ctx, "SemiLepTTBarMatchGENRECO_triggerPFHT_mu"));
    h_semilepttbarmatch_triggerSingleJet_genreco_ele.reset(new TstarTstarGenRecoMatchedHists(ctx, "SemiLepTTBarMatchGENRECO_triggerSingleJet_ele"));
    h_semilepttbarmatch_triggerHT_genreco_ele.reset(new TstarTstarGenRecoMatchedHists(ctx, "SemiLepTTBarMatchGENRECO_triggerHT_ele"));
    h_semilepttbarmatch_triggerPFHT_genreco_ele.reset(new TstarTstarGenRecoMatchedHists(ctx, "SemiLepTTBarMatchGENRECO_triggerPFHT_ele"));
    h_semilepttbarmatch_triggerPFHT_genreco.reset(new TstarTstarGenRecoMatchedHists(ctx, "SemiLepTTBarMatchGENRECO_triggerPFHT"));

    h_semilepttbarmatch_triggerSingleLeptonMu_genreco.reset(new TstarTstarGenRecoMatchedHists(ctx, "SemiLepTTBarMatchGENRECO_triggerSingleLeptonMu"));
    h_semilepttbarmatch_triggerSingleLeptonEle_genreco.reset(new TstarTstarGenRecoMatchedHists(ctx, "SemiLepTTBarMatchGENRECO_triggerSingleLeptonEle"));

    h_After_TstarTstar_Reco.reset(new TstarTstarHists(ctx, "After_TstarTstar_Reco"));
    h_After_TstarTstar_Reco_match.reset(new TstarTstarHists(ctx, "After_TstarTstar_Reco_match"));

    h_RecoPlots_After_ttbar.reset(new TstarTstarRecoTstarHists(ctx, "RecoPlots_After_ttbar"));
    h_RecoPlots_After_ttbar_correct_ttbar.reset(new TstarTstarRecoTstarHists(ctx, "RecoPlots_After_ttbar_correct_ttbar"));
    h_RecoPlots_After_TstarTstar.reset(new TstarTstarRecoTstarHists(ctx, "RecoPlots_After_TstarTstar"));
    h_RecoPlots_After_TstarTstar_match.reset(new TstarTstarRecoTstarHists(ctx, "RecoPlots_After_TstarTstar_match"));

    h_RecoPlots_After_TstarTstar_tgtg.reset(new TstarTstarRecoTstarHists(ctx, "RecoPlots_After_TstarTstar_tgtg"));
    h_RecoPlots_After_TstarTstar_tgtg_correct_ttbar.reset(new TstarTstarRecoTstarHists(ctx, "RecoPlots_After_TstarTstar_tgtg_correct_ttbar"));


    //4. Set up ttbar reconstruction
    const std::string ttbar_hyps_label("TTbarReconstruction");
    const std::string ttbar_chi2_label("Chi2");
    reco_primlep.reset(new PrimaryLepton(ctx));
    ttbar_reco.reset(new HighMassTTbarReconstruction(ctx, NeutrinoReconstruction,ttbar_hyps_label));
    h_ttbar_hyps = ctx.get_handle<std::vector<ReconstructionHypothesis>>(ttbar_hyps_label);
    h_is_ttbar_reconstructed = ctx.get_handle< bool >("is_ttbar_reconstructed_chi2");
    h_recohyp = ctx.declare_event_output<ReconstructionHypothesis>(ttbar_hyps_label+"_best");

    if(is_tgtg){
      const std::string tstartstar_hyps_label("TstarTstar_tgtg");
      h_tstartstar_hyps = ctx.get_handle<std::vector<ReconstructionTstarHypothesis>>(tstartstar_hyps_label);
      h_recohyp_tstartstar = ctx.declare_event_output<ReconstructionTstarHypothesis>(tstartstar_hyps_label+"_best");
    }
    if(is_tgtgamma){
      const std::string tstartstar_hyps_label("TstarTstar_tgtgamma");
      h_tstartstar_hyps = ctx.get_handle<std::vector<ReconstructionTstarHypothesis>>(tstartstar_hyps_label);
      h_recohyp_tstartstar = ctx.declare_event_output<ReconstructionTstarHypothesis>(tstartstar_hyps_label+"_best");
    }

    ttbar_discriminator.reset(new ttbarChi2Discriminator(ctx));

    // h_M_Tstar_gluon = ctx.get_handle< float >("M_Tstar_gluon");
    // h_M_Tstar_gamma = ctx.get_handle< float >("M_Tstar_gamma");
    // h_M_Tstar_lep = ctx.get_handle< float >("M_Tstar_lep");
    // h_M_Tstar_had = ctx.get_handle< float >("M_Tstar_had");

    TstarTstar_reco.reset(new TstarTstar_Reconstruction(ctx));
    TstarTstar_tgluon_tgluon_reco.reset(new TstarTstar_tgluon_tgluon_Reconstruction(ctx));
}


bool TstarTstarMCStudyModule::process(Event & event) {
   
  if(debug)   
    cout << "TstarTstarMCStudyModule: Starting to process event (runid, eventid) = (" << event.run << ", " << event.event << "); weight = " << event.weight << endl;
  // if(debug){
  //   cout<<"BEFORE N photons = "<<event.photons->size()<<endl;
  //   for (const Photon & thisgamma : *event.photons){
  //     cout<<" thisgamma.pt() = "<<thisgamma.pt()<<" thisgamma.eta() = "<<thisgamma.eta()<<endl;
  //   }
  // }

  // event.set(h_M_Tstar_gluon, 0.);
  // event.set(h_M_Tstar_gamma, 0.);
  // event.set(h_M_Tstar_lep, 0.);
  // event.set(h_M_Tstar_had, 0.);

  event.set(h_is_ttbar_reconstructed, false);
  event.set(h_recohyp, ReconstructionHypothesis());
  event.set(h_recohyp_tstartstar,ReconstructionTstarHypothesis());

  if(debug){cout << "Finished initialization of Handle Variables" << endl;}

  h_nocuts->fill(event);
  common->process(event);
  h_common->fill(event);

  // if(debug){
  //   cout<<"AFTER N photons = "<<event.photons->size()<<endl;
  //   for (const Photon & thisgamma : *event.photons){
  //     cout<<" thisgamma.pt() = "<<thisgamma.pt()<<" thisgamma.eta() = "<<thisgamma.eta()<<endl;
  //   }
  // }

  //---- Loose selection
  // Require at least 1 jet
  const bool pass_njet = (event.jets->size()>0);
  if(!pass_njet) return false;
  h_njetsel->fill(event);
  // Require at least one Muon or one Electron
  //  const bool pass_lep1 = ((event.muons->size() >= 1) || (event.electrons->size() >= 1));

  // Require exactly one muon or one electron
  const bool pass_lep1 = (((event.muons->size() == 1) || (event.electrons->size() == 1)) && (event.electrons->size()+event.muons->size()) == 1);

  if(!pass_lep1) return false;
  h_lepsel->fill(event);

  // Require more than one photon
  const bool pass_npho = (event.photons->size()>0);

  if(!pass_npho) return false;
  h_nphosel->fill(event);


  // Lepton-2Dcut variables
  for(auto& muo : *event.muons){
    float    dRmin, pTrel;
    std::tie(dRmin, pTrel) = drmin_pTrel(muo, *event.jets);
    muo.set_tag(Muon::twodcut_dRmin, dRmin);
    muo.set_tag(Muon::twodcut_pTrel, pTrel);
    if(debug) cout<<"Muon.Pt = "<<muo.pt()<<endl;
  }
  for(auto& ele : *event.electrons){
    float    dRmin, pTrel;
    std::tie(dRmin, pTrel) = drmin_pTrel(ele, *event.jets);

    ele.set_tag(Electron::twodcut_dRmin, dRmin);
    ele.set_tag(Electron::twodcut_pTrel, pTrel);
    if(debug) cout<<"Electron.Pt = "<<ele.pt()<<endl;
  }
  const bool pass_twodcut = twodcut_sel->passes(event);
  if(!pass_twodcut) return false;
  h_2dcut->fill(event);
  if(debug) cout<<"passed 2D cut"<<endl;

  //---- Matching to GEN
  const bool pass_ttbarsemilep = TTbarSemiLepMatchable_selection->passes(event);
  if(pass_ttbarsemilep){
    h_semilepttbarmatch->fill(event);
    h_semilepttbarmatch_gen->fill(event);
    h_semilepttbarmatch_genreco->fill(event);
    if(event.muons->size() == 1) h_semilepttbarmatch_genreco_mu->fill(event);
    if(event.electrons->size() == 1) h_semilepttbarmatch_genreco_ele->fill(event);
  }
  else h_nosemilepttbarmatch->fill(event);

  if(isTrigger){
    //    cout<<"AAAAAAAAAa why am I here???"<<endl;
    if(debug) cout<<"N AK4 jets = "<<event.jets->size()<<" N AK8 jets = "<<event.topjets->size()<<endl;
    bool pass_trigger_SingleJet = (triggerSingleJet450_sel->passes(event) && event.jets->at(0).pt()>450);
    if(pass_trigger_SingleJet && pass_ttbarsemilep){
      h_semilepttbarmatch_triggerSingleJet->fill(event);
      h_semilepttbarmatch_triggerSingleJet_genreco->fill(event);
      if(event.muons->size() == 1) h_semilepttbarmatch_triggerSingleJet_genreco_mu->fill(event);
      if(event.electrons->size() == 1) h_semilepttbarmatch_triggerSingleJet_genreco_ele->fill(event);
    }
    bool pass_trigger_SingleMu = (triggerSingleLeptonMu1_sel->passes(event) || triggerSingleLeptonMu2_sel->passes(event)
				  || triggerSingleLeptonMu3_sel->passes(event) || triggerSingleLeptonMu4_sel->passes(event));
    if(pass_trigger_SingleMu && pass_ttbarsemilep && (event.muons->size() == 1)){
      h_semilepttbarmatch_triggerSingleLeptonMu->fill(event);
      h_semilepttbarmatch_triggerSingleLeptonMu_genreco->fill(event);
    }
    bool pass_trigger_SingleEle = (triggerSingleLeptonEle1_sel->passes(event) || triggerSingleLeptonEle2_sel->passes(event) || triggerSingleLeptonEle3_sel->passes(event));
    //    bool pass_trigger_SingleEle = triggerSingleLeptonEle_sel->passes(event);
    if(pass_trigger_SingleEle && pass_ttbarsemilep && (event.electrons->size() == 1)){
      h_semilepttbarmatch_triggerSingleLeptonEle->fill(event);
      h_semilepttbarmatch_triggerSingleLeptonEle_genreco->fill(event);
    }
    bool pass_trigegr_HT = triggerHT1_sel->passes(event) || triggerHT2_sel->passes(event) || triggerHT3_sel->passes(event) 
      || triggerHT4_sel->passes(event) || triggerHT5_sel->passes(event) || triggerHT6_sel->passes(event);
    if(pass_trigegr_HT &&  pass_ttbarsemilep){
      h_semilepttbarmatch_triggerHT->fill(event);
      h_semilepttbarmatch_triggerHT_genreco->fill(event);
      if(debug) cout<<"All HT"<<endl;
      if(event.muons->size() == 1) h_semilepttbarmatch_triggerHT_genreco_mu->fill(event);
      if(debug) cout<<"MU HT"<<endl;
      if(event.electrons->size() == 1) h_semilepttbarmatch_triggerHT_genreco_ele->fill(event);
      if(debug) cout<<"ELE HT"<<endl;
    }
    bool pass_trigegr_PFHT = triggerPFHT_sel->passes(event);
    if(pass_trigegr_PFHT &&  pass_ttbarsemilep){
      h_semilepttbarmatch_triggerPFHT->fill(event);
      if(debug) cout<<"All PFHT RECO"<<endl;
      h_semilepttbarmatch_triggerPFHT_genreco->fill(event);
      if(debug) cout<<"All PFHT GENRECO"<<endl;
      if(event.muons->size() == 1) h_semilepttbarmatch_triggerPFHT_genreco_mu->fill(event);
      if(debug) cout<<"MU PFHT"<<endl;
      if(event.electrons->size() == 1) h_semilepttbarmatch_triggerPFHT_genreco_ele->fill(event);
      if(debug) cout<<"ELE PFHT"<<endl;
    }
    if(debug) cout<<"pass_trigger_SingleJet "<<pass_trigger_SingleJet<<" pass_trigger_SingleMu "<<pass_trigger_SingleMu<<" pass_trigger_SingleEle "<<pass_trigger_SingleEle<<endl;
  }

  reco_primlep->process(event);//set "primary lepton"
  if(debug) {cout << "Starting ttbar reconstruction... ";}\
  ttbar_reco->process(event);//reconstruct ttbar
  
  /** old stuff
  // goal: save only the chi2-best ttbar hypothesis in output sub-ntuple
  std::vector<ReconstructionHypothesis>& hyps = event.get(h_ttbar_hyps); //hyps contains all hypothesises
  if(debug) cout<<"Number of ttbar hyps = "<<hyps.size()<<endl;
  //  const ReconstructionHypothesis* hyp = get_best_hypothesis(hyps, "Chi2"); TODO
  ReconstructionHypothesis hyp;
  if(hyps.size()>0) hyp = hyps.at(0);//FixMe: change placeholder to "best hypothesis" candidate
  double W_mass = hyp.wlep_v4().M();
  if(abs(W_mass-80.399)>10) return false; //store only hypotheses with 2 solutions, which leads to correct W mass
  if(debug) cout<<" Best hypothesis ...  "<<hyp.neutrino_v4()<<" W mass = "<<W_mass<<endl;
  hyps.clear();
  if(debug) cout<<" Best hypothesis copied ...  "<<endl;
  event.set(h_recohyp, hyp); //save "best" hypothesis in event
  **/

  if(debug) {cout << "Finished. Finding best Hypothesis..."<< endl;}   
  ttbar_discriminator->process(event);

  ReconstructionHypothesis hyp = event.get(h_recohyp);  
  if(debug) {cout << "Start TstarTstar reconstruction ..."<< endl;}
  if(event.get(h_is_ttbar_reconstructed)){
    if(is_tgtgamma){ // Tstar+Tstar -> t+g + t+gamma //FixME: the code is broken and most probably won't work. Sorry!
      h_RecoPlots_After_ttbar->fill(event);
      if(TstarTstar_reco->process(event)){
	h_After_TstarTstar_Reco->fill(event);
	h_RecoPlots_After_TstarTstar->fill(event);
	if(pass_ttbarsemilep){   //Check same plots for matching
	  h_After_TstarTstar_Reco_match->fill(event);
	  h_RecoPlots_After_TstarTstar_match->fill(event);
	}
      }
    }
    if(is_tgtg){ // Tstar+Tstar -> t+g + t+g
      bool pass_check_reco_ttbar = false;
      if(pass_ttbarsemilep && check_ttbar_reco){
	std::vector<ReconstructionHypothesis> ttbar_all_hyps = event.get(h_ttbar_hyps);
	double mindRsum = 1e6; int best_match_hyp = -1;
	bool pass_check_reco_ttbar = false;
	  for(unsigned int i=0; i<ttbar_all_hyps.size(); i++){
	    std::pair<bool,double> pass_check_reco_ttbar_pair = TTbarSemiLepMatchable_selection->check_reco(ttbar_all_hyps.at(i));
	    if(pass_check_reco_ttbar_pair.first && pass_check_reco_ttbar_pair.second<mindRsum){
	      mindRsum = pass_check_reco_ttbar_pair.second; 
	      best_match_hyp = i;
	      pass_check_reco_ttbar = true;
	    }
	  
	  if(pass_check_reco_ttbar){
	    h_RecoPlots_After_ttbar_correct_ttbar->fill(event);//fill once per event if at least one good ttbar hypothesis found
	    if(debug) cout<<"Best hyp #"<<best_match_hyp<<" dR_sum = "<<mindRsum<<endl;
	    h_RecoPlots_After_ttbar_correct_ttbar->fill_ttbarhyps(event, ttbar_all_hyps.at(best_match_hyp));
	  }
	}
      }
      //      if(debug) cout<<"pass_ttbarsemilep "<<pass_ttbarsemilep<<" pass_check_reco_ttbar "<<pass_check_reco_ttbar<<endl;
      bool pass_tgluon_tgluon_reco = TstarTstar_tgluon_tgluon_reco->process(event);//modify here to reconstruct with best macth ttbar?
      if(pass_tgluon_tgluon_reco) 
	h_RecoPlots_After_TstarTstar_tgtg->fill(event);
      if(pass_check_reco_ttbar && pass_tgluon_tgluon_reco)
	h_RecoPlots_After_TstarTstar_tgtg_correct_ttbar->fill(event);
    }
  }
  else{
    if(debug){cout << "Event has no best hypothesis!" << endl;}
  }
  if(debug){cout << "Done ##################################" << endl << endl;}
  return true;
}

// as we want to run the ExampleCycleNew directly with AnalysisModuleRunner,
// make sure the TstarTstarMCStudyModule is found by class name. This is ensured by this macro:
UHH2_REGISTER_ANALYSIS_MODULE(TstarTstarMCStudyModule)

}
