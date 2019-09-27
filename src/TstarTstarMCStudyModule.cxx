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
#include "UHH2/TstarTstar/include/TstarTstarGenHists.h"
#include "UHH2/TstarTstar/include/TstarTstarGenRecoMatchedHists.h"

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
    unique_ptr<Selection> TTbarSemiLepMatchable_selection;
    unique_ptr<Selection> triggerSingleJet450_sel;
    unique_ptr<Selection> triggerSingleLeptonEle_sel;
    unique_ptr<Selection> triggerSingleLeptonMu1_sel;
    unique_ptr<Selection> triggerSingleLeptonMu2_sel;
    unique_ptr<Selection> triggerHT1_sel, triggerHT2_sel, triggerHT3_sel, triggerHT4_sel, triggerHT5_sel,  triggerHT6_sel;
    // Store the Hists collection as member variables. Again, use unique_ptr to avoid memory leaks.
  std::unique_ptr<Hists> h_nocuts, h_common, h_lepsel, h_njetsel, h_2dcut, h_semilepttbarmatch, h_nosemilepttbarmatch;
  std::unique_ptr<Hists> h_semilepttbarmatch_triggerSingleJet, h_semilepttbarmatch_triggerSingleLeptonMu, h_semilepttbarmatch_triggerSingleLeptonEle, h_semilepttbarmatch_triggerHT;
  std::unique_ptr<Hists> h_semilepttbarmatch_gen;
  std::unique_ptr<Hists> h_semilepttbarmatch_genreco,h_semilepttbarmatch_genreco_mu,h_semilepttbarmatch_genreco_ele;
  std::unique_ptr<Hists> h_semilepttbarmatch_triggerSingleJet_genreco, h_semilepttbarmatch_triggerSingleLeptonMu_genreco, h_semilepttbarmatch_triggerSingleLeptonEle_genreco, h_semilepttbarmatch_triggerHT_genreco;
  std::unique_ptr<Hists> h_semilepttbarmatch_triggerSingleJet_genreco_mu, h_semilepttbarmatch_triggerHT_genreco_mu;
  std::unique_ptr<Hists> h_semilepttbarmatch_triggerSingleJet_genreco_ele, h_semilepttbarmatch_triggerHT_genreco_ele;
    bool debug = false;
    bool isTrigger = false;
  //    bool debug = true;
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

    // 1. setup other modules. CommonModules
    common.reset(new CommonModules());
    common->switch_metcorrection();
    common->switch_jetlepcleaner();
    common->switch_jetPtSorter();
    common->set_jet_id(AndId<Jet>(JetPFID(JetPFID::WP_TIGHT_PUPPI), PtEtaCut(30.0,5.2)));
    ElectronId eleID; 
    double electron_pt(20.);
    eleID = ElectronID_Summer16_tight_noIso;
    common->set_electron_id(AndId<Electron>(PtEtaSCCut(electron_pt, 2.5), eleID));

    PhotonId phoID; 
    double photon_pt(20.);
    //    phoID = PhotonTagID(Photon::cutBasedPhotonID_Spring16_V2p2_tight);
    //    phoID = PhotonTagID(Photon::cutBasedPhotonID_Spring16_V2p2_loose);
    phoID = PhotonTagID(Photon::cutBasedPhotonID_Fall17_94X_V2_loose);
    common->set_photon_id(AndId<Photon>(PtEtaCut(photon_pt, 5.2), phoID));

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
    triggerSingleJet450_sel.reset(new TriggerSelection("HLT_PFJet450_v*"));
    //    triggerSingleLeptonEle_sel.reset(new TriggerSelection("HLT_Ele32_eta2p1_WPTight_Gsf_v*"));
    // triggerSingleLeptonMu1_sel.reset(new TriggerSelection("HLT_IsoMu24_v*"));
    // triggerSingleLeptonMu2_sel.reset(new TriggerSelection("HLT_IsoTkMu24_v*"));

    triggerSingleLeptonEle_sel.reset(new TriggerSelection("HLT_Ele115_CaloIdVT_GsfTrkIdT_v*"));
    triggerSingleLeptonMu1_sel.reset(new TriggerSelection("HLT_Mu50_v*"));
    triggerSingleLeptonMu2_sel.reset(new TriggerSelection("HLT_Mu55_v*"));


    triggerHT1_sel.reset(new TriggerSelection("HLT_HT430to450_v*"));
    triggerHT2_sel.reset(new TriggerSelection("HLT_HT450to470_v*"));
    triggerHT3_sel.reset(new TriggerSelection("HLT_HT470to500_v*"));
    triggerHT4_sel.reset(new TriggerSelection("HLT_HT500to550_v*"));
    triggerHT5_sel.reset(new TriggerSelection("HLT_HT550to650_v*"));
    triggerHT6_sel.reset(new TriggerSelection("HLT_HT650_v*"));



    // 3. Set up Hists classes:
    h_nocuts.reset(new TstarTstarHists(ctx, "NoCuts"));
    h_common.reset(new TstarTstarHists(ctx, "AfterCommon"));
    h_njetsel.reset(new TstarTstarHists(ctx, "AfterNjets"));
    h_lepsel.reset(new TstarTstarHists(ctx, "AfterLepSel"));
    h_2dcut.reset(new TstarTstarHists(ctx, "After2D"));
    h_semilepttbarmatch.reset(new TstarTstarHists(ctx, "SemiLepTTBarMatch"));

    h_semilepttbarmatch_triggerSingleJet.reset(new TstarTstarHists(ctx, "SemiLepTTBarMatch_triggerSingleJet"));
    h_semilepttbarmatch_triggerSingleLeptonMu.reset(new TstarTstarHists(ctx, "SemiLepTTBarMatch_triggerSingleLeptonMu"));
    h_semilepttbarmatch_triggerSingleLeptonEle.reset(new TstarTstarHists(ctx, "SemiLepTTBarMatch_triggerSingleLeptonEle"));
    h_semilepttbarmatch_triggerHT.reset(new TstarTstarHists(ctx, "SemiLepTTBarMatch_triggerHT"));

    h_nosemilepttbarmatch.reset(new TstarTstarHists(ctx, "NotSemiLepTTBarMatch"));
    h_semilepttbarmatch_gen.reset(new TstarTstarGenHists(ctx, "SemiLepTTBarMatchGEN"));
    h_semilepttbarmatch_genreco.reset(new TstarTstarGenRecoMatchedHists(ctx, "SemiLepTTBarMatchGENRECO"));
    h_semilepttbarmatch_genreco_mu.reset(new TstarTstarGenRecoMatchedHists(ctx, "SemiLepTTBarMatchGENRECO_mu"));
    h_semilepttbarmatch_genreco_ele.reset(new TstarTstarGenRecoMatchedHists(ctx, "SemiLepTTBarMatchGENRECO_ele"));

    h_semilepttbarmatch_triggerSingleJet_genreco.reset(new TstarTstarGenRecoMatchedHists(ctx, "SemiLepTTBarMatchGENRECO_triggerSingleJet"));
    h_semilepttbarmatch_triggerHT_genreco.reset(new TstarTstarGenRecoMatchedHists(ctx, "SemiLepTTBarMatchGENRECO_triggerHT"));
    h_semilepttbarmatch_triggerSingleJet_genreco_mu.reset(new TstarTstarGenRecoMatchedHists(ctx, "SemiLepTTBarMatchGENRECO_triggerSingleJet_mu"));
    h_semilepttbarmatch_triggerHT_genreco_mu.reset(new TstarTstarGenRecoMatchedHists(ctx, "SemiLepTTBarMatchGENRECO_triggerHT_mu"));
    h_semilepttbarmatch_triggerSingleJet_genreco_ele.reset(new TstarTstarGenRecoMatchedHists(ctx, "SemiLepTTBarMatchGENRECO_triggerSingleJet_ele"));
    h_semilepttbarmatch_triggerHT_genreco_ele.reset(new TstarTstarGenRecoMatchedHists(ctx, "SemiLepTTBarMatchGENRECO_triggerHT_ele"));

    h_semilepttbarmatch_triggerSingleLeptonMu_genreco.reset(new TstarTstarGenRecoMatchedHists(ctx, "SemiLepTTBarMatchGENRECO_triggerSingleLeptonMu"));
    h_semilepttbarmatch_triggerSingleLeptonEle_genreco.reset(new TstarTstarGenRecoMatchedHists(ctx, "SemiLepTTBarMatchGENRECO_triggerSingleLeptonEle"));


}


bool TstarTstarMCStudyModule::process(Event & event) {
   
  if(debug)   
    cout << "TstarTstarMCStudyModule: Starting to process event (runid, eventid) = (" << event.run << ", " << event.event << "); weight = " << event.weight << endl;
  if(debug){
    cout<<"BEFORE N photons = "<<event.photons->size()<<endl;
    for (const Photon & thisgamma : *event.photons){
      cout<<" thisgamma.pt() = "<<thisgamma.pt()<<" thisgamma.eta() = "<<thisgamma.eta()<<endl;
    }
  }

  h_nocuts->fill(event);
  common->process(event);
  h_common->fill(event);

  if(debug){
    cout<<"AFTER N photons = "<<event.photons->size()<<endl;
    for (const Photon & thisgamma : *event.photons){
      cout<<" thisgamma.pt() = "<<thisgamma.pt()<<" thisgamma.eta() = "<<thisgamma.eta()<<endl;
    }
  }

  //---- Loose selection
  // Require at least 1 jet
  const bool pass_njet = (event.jets->size()>0);
  if(!pass_njet) return false;
  h_njetsel->fill(event);
  // Require at least one Muon or one Electron
  const bool pass_lep1 = ((event.muons->size() >= 1) || (event.electrons->size() >= 1));
  if(!pass_lep1) return false;
  h_lepsel->fill(event);
  // Lepton-2Dcut variables
  for(auto& muo : *event.muons){
    float    dRmin, pTrel;
    std::tie(dRmin, pTrel) = drmin_pTrel(muo, *event.jets);
    muo.set_tag(Muon::twodcut_dRmin, dRmin);
    muo.set_tag(Muon::twodcut_pTrel, pTrel);
  }
  for(auto& ele : *event.electrons){
    float    dRmin, pTrel;
    std::tie(dRmin, pTrel) = drmin_pTrel(ele, *event.jets);

    ele.set_tag(Electron::twodcut_dRmin, dRmin);
    ele.set_tag(Electron::twodcut_pTrel, pTrel);
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
    if(debug) cout<<"N jets = "<<event.jets->size()<<endl;
    bool pass_trigger_SingleJet = (triggerSingleJet450_sel->passes(event) && event.jets->at(0).pt()>450);
    if(pass_trigger_SingleJet && pass_ttbarsemilep){
      h_semilepttbarmatch_triggerSingleJet->fill(event);
      h_semilepttbarmatch_triggerSingleJet_genreco->fill(event);
      if(event.muons->size() == 1) h_semilepttbarmatch_triggerSingleJet_genreco_mu->fill(event);
      if(event.electrons->size() == 1) h_semilepttbarmatch_triggerSingleJet_genreco_ele->fill(event);
    }
    bool pass_trigger_SingleMu = (triggerSingleLeptonMu1_sel->passes(event) || triggerSingleLeptonMu2_sel->passes(event));
    if(pass_trigger_SingleMu && pass_ttbarsemilep && (event.muons->size() == 1)){
      h_semilepttbarmatch_triggerSingleLeptonMu->fill(event);
      h_semilepttbarmatch_triggerSingleLeptonMu_genreco->fill(event);
    }
    bool pass_trigger_SingleEle = triggerSingleLeptonEle_sel->passes(event);
    if(pass_trigger_SingleEle && pass_ttbarsemilep && (event.electrons->size() == 1)){
      h_semilepttbarmatch_triggerSingleLeptonEle->fill(event);
      h_semilepttbarmatch_triggerSingleLeptonEle_genreco->fill(event);
    }
    bool pass_trigegr_HT = triggerHT1_sel->passes(event) || triggerHT2_sel->passes(event) || triggerHT3_sel->passes(event) 
      || triggerHT4_sel->passes(event) || triggerHT5_sel->passes(event) || triggerHT6_sel->passes(event);
    if(pass_trigegr_HT &&  pass_ttbarsemilep){
      h_semilepttbarmatch_triggerHT->fill(event);
      h_semilepttbarmatch_triggerHT_genreco->fill(event);
      if(event.muons->size() == 1) h_semilepttbarmatch_triggerHT_genreco_mu->fill(event);
      if(event.electrons->size() == 1) h_semilepttbarmatch_triggerHT_genreco_ele->fill(event);
    }
    if(debug) cout<<"pass_trigger_SingleJet "<<pass_trigger_SingleJet<<" pass_trigger_SingleMu "<<pass_trigger_SingleMu<<" pass_trigger_SingleEle "<<pass_trigger_SingleEle<<endl;
  }
  
  return true;
}

// as we want to run the ExampleCycleNew directly with AnalysisModuleRunner,
// make sure the TstarTstarMCStudyModule is found by class name. This is ensured by this macro:
UHH2_REGISTER_ANALYSIS_MODULE(TstarTstarMCStudyModule)

}
