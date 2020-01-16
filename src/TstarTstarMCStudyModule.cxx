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
#include "UHH2/common/include/TTbarGen.h"
#include "UHH2/common/include/TopJetIds.h"
#include "UHH2/TstarTstar/include/TstarTstarSelections.h"
#include "UHH2/TstarTstar/include/TstarTstarHists.h"
#include "UHH2/TstarTstar/include/TstarTstarRecoTstarHists.h"
#include "UHH2/TstarTstar/include/TstarTstarGenHists.h"
#include "UHH2/TstarTstar/include/TstarTstarGenRecoMatchedHists.h"
#include "UHH2/TstarTstar/include/TstarTstarReconstructionModules.h"
#include "UHH2/TstarTstar/include/ReconstructionTstarHypothesis.h"
#include "UHH2/TstarTstar/include/TstarTstarGenMatch.h"

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
  
  // Declare the Selections to use. Use unique_ptr to ensure automatic call of delete in the destructor, to avoid memory leaks.
  unique_ptr<Selection> twodcut_sel;   
  unique_ptr<TTbarSemiLepMatchableSelection> TTbarSemiLepMatchable_selection;
  unique_ptr<Selection> met_sel, st_sel; 
  unique_ptr<Selection> topjet_selection;
  unique_ptr<Selection> toptagevt_sel;

  //Trigger Selections
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


  // ##### Histograms #####
  // Store the Hists collection as member variables. Again, use unique_ptr to avoid memory leaks.
  std::unique_ptr<Hists> h_nocuts_gen, h_allcuts_gen;
  std::unique_ptr<Hists> h_nocuts, h_common, h_lepsel, h_njetsel, h_nphosel, h_2dcut, h_allcuts, h_semilepttbarmatch, h_nosemilepttbarmatch;
  std::unique_ptr<Hists> h_semilepttbarmatch_triggerSingleJet, h_semilepttbarmatch_triggerSingleLeptonMu, h_semilepttbarmatch_triggerSingleLeptonEle, h_semilepttbarmatch_triggerHT, h_semilepttbarmatch_triggerPFHT;
  std::unique_ptr<Hists> h_semilepttbarmatch_gen;
  std::unique_ptr<Hists> h_semilepttbarmatch_genreco,h_semilepttbarmatch_genreco_mu,h_semilepttbarmatch_genreco_ele;
  std::unique_ptr<Hists> h_semilepttbarmatch_triggerSingleJet_genreco, h_semilepttbarmatch_triggerSingleLeptonMu_genreco, h_semilepttbarmatch_triggerSingleLeptonEle_genreco, h_semilepttbarmatch_triggerHT_genreco, h_semilepttbarmatch_triggerPFHT_genreco;
  std::unique_ptr<Hists> h_semilepttbarmatch_triggerSingleJet_genreco_mu, h_semilepttbarmatch_triggerHT_genreco_mu, h_semilepttbarmatch_triggerPFHT_genreco_mu;
  std::unique_ptr<Hists> h_semilepttbarmatch_triggerSingleJet_genreco_ele, h_semilepttbarmatch_triggerHT_genreco_ele, h_semilepttbarmatch_triggerPFHT_genreco_ele;

  std::unique_ptr<Hists> h_After_TstarTstar_Reco, h_After_TstarTstar_Reco_match;
  std::unique_ptr<TstarTstarRecoTstarHists> h_RecoPlots_Full, h_RecoPlots_ttag, h_RecoPlots_nottag; 
  std::unique_ptr<TstarTstarRecoTstarHists> h_RecoPlots_After_ttbar_correct_ttbar;
  std::unique_ptr<TstarTstarRecoTstarHists> h_RecoPlots_GEN;
  std::unique_ptr<TstarTstarMergedHists> h_merged;
  
  // Bools for Debugging/Options
  bool debug = false;

  bool check_ttbar_reco = true;
  bool doTopTagged = true;
  bool doNotTopTagged = false;
  

  // Modules
  std::unique_ptr<uhh2::AnalysisModule> ttgenprod;
  std::unique_ptr<uhh2::AnalysisModule> reco_primlep;
  std::unique_ptr<uhh2::AnalysisModule> ttbar_reco;
  std::unique_ptr<ttbarChi2Discriminator> ttbar_discriminator;
  std::unique_ptr<TstarTstarGenMatcher> genmatcher;
  //  std::unique_ptr<ttbarCorrectMatchDiscriminator> ttbar_CorrectMatchDiscriminator;
  std::unique_ptr<CorrectMatchDiscriminator> ttbar_CorrectMatchDiscriminator;
  std::unique_ptr<TstarTstar_tgluon_tgamma_Reconstruction> TstarTstar_tgluon_tgamma_reco;
  std::unique_ptr<TstarTstar_tgluon_tgluon_Reconstruction2> TstarTstar_tgluon_tgluon_reco;
  std::unique_ptr<uhh2::AnalysisModule> ttbar_reco_toptag;

  // Handles
  uhh2::Event::Handle<std::vector<ReconstructionHypothesis>> h_ttbar_hyps;
  uhh2::Event::Handle<bool> h_is_ttbar_reconstructed;
  uhh2::Event::Handle<ReconstructionHypothesis> h_recohyp;
  uhh2::Event::Handle<std::vector<ReconstructionTstarHypothesis>> h_tstartstar_hyps;
  uhh2::Event::Handle<ReconstructionTstarHypothesis> h_recohyp_tstartstar;
  uhh2::Event::Handle<TTbarGen> h_ttbargen;
  uhh2::Event::Handle<int> h_flag_toptagevent;
  uhh2::Event::Handle<int> h_ttag_jet_pos;

  // bools for channel. will be read in later
  bool is_tgtg, is_tgtgamma, isTrigger;
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
  
  // 0. Reading in which channel
  is_tgtg = false; is_tgtgamma = false;
  if(ctx.get("channel") == "tgtg") is_tgtg = true;
  if(ctx.get("channel") == "tgtgamma") is_tgtgamma = true;
  
  // 1. setup other modules. CommonModules
  common.reset(new CommonModules());
  common->disable_mcpileupreweight();//FixME: PU re-weighting crushes
  common->switch_metcorrection();

  // Jets
  common->switch_jetlepcleaner();
  common->switch_jetPtSorter();
  common->set_jet_id(AndId<Jet>(JetPFID(JetPFID::WP_TIGHT_PUPPI), PtEtaCut(30.0,5.2)));
  
  // Electron
  ElectronId eleID; 
  double electron_pt(20.);
  eleID = ElectronID_Summer16_tight_noIso;
  common->set_electron_id(AndId<Electron>(PtEtaSCCut(electron_pt, 2.5), eleID));

  // photon
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
  
  // Muon
  MuonId muID;
  double muon_pt(20.);
  muID = MuonID(Muon::Highpt);
  common->set_muon_id(AndId<Muon>(PtEtaCut(muon_pt, 2.4), muID));
  
  // init common.
  common->init(ctx);



  // Prepare GEN
  ttgenprod.reset(new TTbarGenProducer(ctx, "ttbargen", false));
  h_ttbargen = ctx.get_handle<TTbarGen>("ttbargen");
    


  // 2. set up selections
  // 2D Cut Lepton-Jets
  twodcut_sel.reset(new TwoDCut(0.4, 25.0));  // The same as in Z'->ttbar semileptonic
  TTbarSemiLepMatchable_selection.reset(new TTbarSemiLepMatchableSelection());
  
  // Trigger Stuff
  isTrigger = (ctx.get("UseTrigger") == "true");
  triggerSingleJet450_sel.reset(new TriggerSelection("HLT_PFJet450_v*"));
  // triggerSingleLeptonEle_sel.reset(new TriggerSelection("HLT_Ele32_eta2p1_WPTight_Gsf_v*"));
  // triggerSingleLeptonMu1_sel.reset(new TriggerSelection("HLT_IsoMu24_v*"));
  // triggerSingleLeptonMu2_sel.reset(new TriggerSelection("HLT_IsoTkMu24_v*"));
  
  // triggerSingleLeptonEle_sel.reset(new TriggerSelection("HLT_Ele115_CaloIdVT_GsfTrkIdT_v*"));
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

  //MET selection
  met_sel.reset(new METCut  (50.,1e6));
  
  //ST selection TODO what is this?
  st_sel.reset(new STCut  (500.,1e6));
  
  //Ak8jet selection TODO check whether I really need to cut on min 2 ak8... although it seems alright
  topjet_selection.reset(new NTopJetSelection(3, -1, TopJetId(PtEtaCut(100, 2.1))));

  //TopTag
  const TopJetId topjetID = AndId<TopJet>(TopTagMassWindow(), TopTagSubbtag(DeepCSVBTag::WP_LOOSE),  Tau32(0.65));
  const float minDR_topjet_jet(1.2);
  toptagevt_sel.reset(new TopTagEventSelection(topjetID, minDR_topjet_jet));
  h_flag_toptagevent = ctx.declare_event_output<int>("flag_toptagevent");
 


  // 3. Set up Hists classes:
  h_nocuts_gen.reset(new TstarTstarGenHists(ctx, "NoCuts_GEN"));
  h_allcuts_gen.reset(new TstarTstarGenHists(ctx, "AfterAllcuts_GEN"));

  h_nocuts.reset(new TstarTstarHists(ctx, "NoCuts"));
  h_common.reset(new TstarTstarHists(ctx, "AfterCommon"));
  h_njetsel.reset(new TstarTstarHists(ctx, "AfterNjets"));
  h_lepsel.reset(new TstarTstarHists(ctx, "AfterLepSel"));
  h_nphosel.reset(new TstarTstarHists(ctx, "AfterNpho"));
  h_2dcut.reset(new TstarTstarHists(ctx, "After2D"));
  h_allcuts.reset(new TstarTstarHists(ctx, "AfterAllcuts"));
  h_semilepttbarmatch.reset(new TstarTstarHists(ctx, "SemiLepTTBarMatch"));
  
  h_nosemilepttbarmatch.reset(new TstarTstarHists(ctx, "NotSemiLepTTBarMatch"));
  h_semilepttbarmatch_gen.reset(new TstarTstarGenHists(ctx, "SemiLepTTBarMatchGEN"));
  h_semilepttbarmatch_genreco.reset(new TstarTstarGenRecoMatchedHists(ctx, "SemiLepTTBarMatchGENRECO"));
  h_semilepttbarmatch_genreco_mu.reset(new TstarTstarGenRecoMatchedHists(ctx, "SemiLepTTBarMatchGENRECO_mu"));
  h_semilepttbarmatch_genreco_ele.reset(new TstarTstarGenRecoMatchedHists(ctx, "SemiLepTTBarMatchGENRECO_ele"));

  if(isTrigger){
    h_semilepttbarmatch_triggerSingleJet.reset(new TstarTstarHists(ctx, "SemiLepTTBarMatch_triggerSingleJet"));
    h_semilepttbarmatch_triggerSingleLeptonMu.reset(new TstarTstarHists(ctx, "SemiLepTTBarMatch_triggerSingleLeptonMu"));
    h_semilepttbarmatch_triggerSingleLeptonEle.reset(new TstarTstarHists(ctx, "SemiLepTTBarMatch_triggerSingleLeptonEle"));
    h_semilepttbarmatch_triggerHT.reset(new TstarTstarHists(ctx, "SemiLepTTBarMatch_triggerHT"));
    h_semilepttbarmatch_triggerPFHT.reset(new TstarTstarHists(ctx, "SemiLepTTBarMatch_triggerPFHT"));
    
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
  }

  h_After_TstarTstar_Reco.reset(new TstarTstarHists(ctx, "After_TstarTstar_Reco"));
  h_After_TstarTstar_Reco_match.reset(new TstarTstarHists(ctx, "After_TstarTstar_Reco_match"));

  h_RecoPlots_Full.reset(new TstarTstarRecoTstarHists(ctx, "RecoPlots_Full"));
  h_RecoPlots_ttag.reset(new TstarTstarRecoTstarHists(ctx, "RecoPlots_ttag"));
  h_RecoPlots_nottag.reset(new TstarTstarRecoTstarHists(ctx, "RecoPlots_nottag"));

  h_RecoPlots_GEN.reset(new TstarTstarRecoTstarHists(ctx, "RecoPlots_GEN"));
  h_merged.reset(new TstarTstarMergedHists(ctx, "Merge Check"));



  //4. Set up ttbar reconstruction
  const std::string ttbar_hyps_label("TTbarReconstruction");
  const std::string ttbar_chi2_label("Chi2");
  reco_primlep.reset(new PrimaryLepton(ctx));

  if(is_tgtg){
    ttbar_reco.reset(new HighMassSkipJetsTTbarReconstruction(ctx, NeutrinoReconstruction,ttbar_hyps_label,0));
    h_ttbar_hyps = ctx.get_handle<std::vector<ReconstructionHypothesis>>(ttbar_hyps_label);
    h_is_ttbar_reconstructed = ctx.get_handle< bool >("is_ttbar_reconstructed_chi2");
    h_recohyp = ctx.declare_event_output<ReconstructionHypothesis>(ttbar_hyps_label+"_best");
    
    const std::string tstartstar_hyps_label("TstarTstar_tgtg");
    h_tstartstar_hyps = ctx.get_handle<std::vector<ReconstructionTstarHypothesis>>(tstartstar_hyps_label);
    h_recohyp_tstartstar = ctx.declare_event_output<ReconstructionTstarHypothesis>(tstartstar_hyps_label+"_best");
  }
  if(is_tgtgamma){
    ttbar_reco.reset(new HighMassSkipJetsTTbarReconstruction(ctx, NeutrinoReconstruction,ttbar_hyps_label,0));
    h_ttbar_hyps = ctx.get_handle<std::vector<ReconstructionHypothesis>>(ttbar_hyps_label);
    h_is_ttbar_reconstructed = ctx.get_handle< bool >("is_ttbar_reconstructed_chi2");
    h_recohyp = ctx.declare_event_output<ReconstructionHypothesis>(ttbar_hyps_label+"_best");
    
    const std::string tstartstar_hyps_label("TstarTstar_tgtgamma");
    h_tstartstar_hyps = ctx.get_handle<std::vector<ReconstructionTstarHypothesis>>(tstartstar_hyps_label);
    h_recohyp_tstartstar = ctx.declare_event_output<ReconstructionTstarHypothesis>(tstartstar_hyps_label+"_best");
  }
  
  ttbar_discriminator.reset(new ttbarChi2Discriminator(ctx));
  ttbar_CorrectMatchDiscriminator.reset(new CorrectMatchDiscriminator(ctx,ttbar_hyps_label));
  
  TstarTstar_tgluon_tgamma_reco.reset(new TstarTstar_tgluon_tgamma_Reconstruction(ctx));
  TstarTstar_tgluon_tgluon_reco.reset(new TstarTstar_tgluon_tgluon_Reconstruction2(ctx));

  ttbar_reco_toptag.reset(new TstarTstarTopTagReconstruction(ctx, NeutrinoReconstruction, ttbar_hyps_label, topjetID, minDR_topjet_jet));

  genmatcher.reset(new TstarTstarGenMatcher(ctx));
  
  h_ttag_jet_pos = ctx.get_handle<int>("ttag_jet_pos");

}


bool TstarTstarMCStudyModule::process(Event & event) {
   
  if(debug){cout << endl << "TstarTstarMCStudyModule: Starting to process event (runid, eventid) = (" << event.run << ", " << event.event << "); weight = " << event.weight << endl;}

  // Fix empty handle problem
  // (Dummy values are filled into handles so that they are not empty)
  event.set(h_is_ttbar_reconstructed, false);
  event.set(h_recohyp, ReconstructionHypothesis());
  event.set(h_recohyp_tstartstar, ReconstructionTstarHypothesis());
  std::vector<ReconstructionHypothesis> fixvec;
  event.set(h_ttbar_hyps, fixvec);
  event.set(h_ttag_jet_pos, -1);
  if(debug){cout << "Finished initialization of Handle Variables" << endl;}

  // ##### Selection #####
  h_nocuts->fill(event);
  h_nocuts_gen->fill(event);
  if(!(common->process(event))){return false;}
  h_common->fill(event);

  //Fill ttgen object for correct matching check, etc
  ttgenprod->process(event);

  //---- Loose selection
  // Require at least 1 jet
  const bool pass_njet = (event.jets->size()>0);
  if(!pass_njet) return false;
  h_njetsel->fill(event);

  // Require exactly one muon or one electron
  const bool pass_lep1 = (((event.muons->size() == 1) || (event.electrons->size() == 1)) && (event.electrons->size()+event.muons->size()) == 1);
  if(!pass_lep1) return false;
  h_lepsel->fill(event);
  
  if(is_tgtgamma){
    // Require at least 1 photon
    const bool pass_npho = (event.photons->size()>0);
    if(!pass_npho) return false;
    h_nphosel->fill(event);
  }

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

  // MET Selection
  bool pass_MET =  met_sel->passes(event);
  if(!pass_MET) return false;

  // ST Selection
  bool pass_ST =  st_sel->passes(event);
  if(!pass_ST) return false;

  // AK8 Selection
  bool pass_ak8 = topjet_selection->passes(event);
  if(!pass_ak8) return false;


  if(debug) cout<<"passed all cuts"<<endl;
  h_allcuts->fill(event);
  h_nocuts_gen->fill(event);


  /* TOPTAG-EVENT boolean */
  const bool pass_ttagevt = toptagevt_sel->passes(event);
  /* add flag_toptagevent to output ntuple */
  if(debug){cout << "TTag:" << pass_ttagevt << endl;}
  event.set(h_flag_toptagevent, int(pass_ttagevt));


  // ##### Matching to GEN ######
  const bool pass_ttbarsemilep = TTbarSemiLepMatchable_selection->passes(event); // check whether event is matchable to GEN ttbar
  if(pass_ttbarsemilep){
    if(debug) cout<<"event is matchable"<<endl;
    h_semilepttbarmatch->fill(event);
    h_semilepttbarmatch_gen->fill(event);
    h_semilepttbarmatch_genreco->fill(event);
    if(event.muons->size() == 1) h_semilepttbarmatch_genreco_mu->fill(event);
    if(event.electrons->size() == 1) h_semilepttbarmatch_genreco_ele->fill(event);
  }
  else h_nosemilepttbarmatch->fill(event);


  // ##### Trigger analysis #####
  if(isTrigger){
    if(debug) cout<<"Starting arigger analysis"<<endl;
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



  // ########### TstarTstarReco! ############
  if(debug) cout<<"Starting TstarTstar Reconstruction part"<<endl;
  reco_primlep->process(event);//set "primary lepton"

  if(is_tgtg){

    if(event.get(h_flag_toptagevent) && doTopTagged){ // if we have a TopTag!
      if(debug){cout << "We have an event with TopTag!" << endl;}
      ttbar_reco_toptag->process(event);
    }
    else{
      if(!doNotTopTagged){return false;}
      ttbar_reco->process(event);
    }

    // TODO ttbar_CorrectMatchDiscriminator->process(event);//find matched to ttbar gen hypothesis
      
    if(debug) {cout << "Finished. Finding best Hypothesis..."<< endl;}  
    ttbar_discriminator->process(event); // Chose best ttbar hypothesis
    
    if(debug) {cout << "Start TstarTstar reconstruction ..."<< endl;}

    if(event.get(h_is_ttbar_reconstructed)){
      if(debug){cout << "best ttbar hyp chosen" << endl;}
            
      if(TstarTstar_tgluon_tgluon_reco->process(event)){
	if(debug){cout << "Tstar Reconstructed: filling histogram" << endl;}
	h_After_TstarTstar_Reco->fill(event);
	h_RecoPlots_Full->fill(event);
	h_merged->fill(event);
	if(event.get(h_flag_toptagevent)){h_RecoPlots_ttag->fill(event);}
	else{h_RecoPlots_nottag->fill(event);}

	// TODO Stuff mit pass_ttbarsemilep
      }
      
    }
    else{
      if(debug){cout << "Event has no best hypothesis!" << endl;}
    }
  }
  else if (is_tgtgamma){
    // TODO
  }  
  
  if(false){
    // ##### GEN Matching
    if(check_ttbar_reco){ // Check matching
      ReconstructionTstarHypothesis hyp_tmp = event.get(h_recohyp_tstartstar);
      if(genmatcher->process(event)){
	//cout << "GEN MATCHED!" << endl;
	h_RecoPlots_GEN->fill(event);
	event.set(h_recohyp_tstartstar, hyp_tmp);
      }
    }
  }


  if(debug){cout << "Done ##################################" << endl;}
  return true;
}

// as we want to run the ExampleCycleNew directly with AnalysisModuleRunner,
// make sure the TstarTstarMCStudyModule is found by class name. This is ensured by this macro:
UHH2_REGISTER_ANALYSIS_MODULE(TstarTstarMCStudyModule)

}
