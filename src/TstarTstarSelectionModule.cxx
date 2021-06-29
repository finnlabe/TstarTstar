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

  // selections
  unique_ptr<Selection> twodcut_sel;
  unique_ptr<Selection> toptagevt_sel;

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
  std::unique_ptr<Hists> h_beginSel,            h_btagcut,             h_2Dcut,               h_dRcut,               h_STcut            ;
  std::unique_ptr<Hists> h_beginSel_gen,        h_btagcut_gen,         h_2Dcut_gen,           h_dRcut_gen,           h_STcut_gen        ;
  std::unique_ptr<Hists> h_beginSel_ele,        h_btagcut_ele,         h_2Dcut_ele,           h_dRcut_ele,           h_STcut_ele        ;
  std::unique_ptr<Hists> h_beginSel_ele_lowpt,  h_btagcut_ele_lowpt,   h_2Dcut_ele_lowpt,     h_dRcut_ele_lowpt,     h_STcut_ele_lowpt  ;
  std::unique_ptr<Hists> h_beginSel_ele_highpt, h_btagcut_ele_highpt,  h_2Dcut_ele_highpt,    h_dRcut_ele_highpt,    h_STcut_ele_highpt ;
  std::unique_ptr<Hists> h_beginSel_mu,         h_btagcut_mu,          h_2Dcut_mu,            h_dRcut_mu,            h_STcut_mu         ;
  std::unique_ptr<Hists> h_beginSel_mu_lowpt,   h_btagcut_mu_lowpt,    h_2Dcut_mu_lowpt,      h_dRcut_mu_lowpt,      h_STcut_mu_lowpt   ;
  std::unique_ptr<Hists> h_beginSel_mu_highpt,  h_btagcut_mu_highpt,   h_2Dcut_mu_highpt,     h_dRcut_mu_highpt,     h_STcut_mu_highpt  ;

  std::unique_ptr<Hists> h_afterSelection_gen, h_afterSelection_genmatch;
  std::unique_ptr<Hists> h_afterSelection;

  // TODO better trigger plots!
  std::unique_ptr<Hists> h_trigger, h_trigger_mu, h_trigger_ele;
  std::unique_ptr<Hists> h_trigger_gen;

  // ###### Handles ######
  uhh2::Event::Handle<FlavorParticle> h_primlep;
  uhh2::Event::Handle<TTbarGen> h_ttbargen;
  uhh2::Event::Handle<bool> h_is_muevt;
  uhh2::Event::Handle<double> h_evt_weight;
  uhh2::Event::Handle<LorentzVector> h_neutrino;
  uhh2::Event::Handle<double> h_ST;

  // ###### Control Switches ######
  bool debug = false;
  bool isTrigger = false;

  // ###### other needed definitions ######
  bool is_MC;
  bool data_isMu = false;
  bool data_is2017B = false;
  TString year;
  TopJetId topjetID;

};


TstarTstarSelectionModule::TstarTstarSelectionModule(Context & ctx){

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

  // year of samples
  year = ctx.get("year", "<not set>");
  if(year == "<not set>"){
    if(ctx.get("dataset_version").find("2016") != std::string::npos) year = "2016";
    else if(ctx.get("dataset_version").find("2017") != std::string::npos) year = "2017";
    else if(ctx.get("dataset_version").find("2018") != std::string::npos) year = "2018";
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

  // Top Tagging
  /**
  topjetID = AndId<TopJet>(HOTVRTopTag(), Tau32Groomed(0.56));
  toptagevt_sel.reset(new TopTagEventSelection(topjetID));
  **/

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

  h_afterSelection.reset(new TstarTstarHists(ctx, "AfterSel"));
  h_afterSelection_gen.reset(new TstarTstarGenHists(ctx, "AfterSel_gen"));
  h_afterSelection_genmatch.reset(new TstarTstarGenRecoMatchedHists(ctx, "AfterSel_genmatch"));

  // TODO
  h_trigger.reset(new TstarTstarHists(ctx, "AfterTrigger"));
  h_trigger_mu.reset(new TstarTstarHists(ctx, "AfterTrigger_mu"));
  h_trigger_ele.reset(new TstarTstarHists(ctx, "AfterTrigger_ele"));
  h_trigger_gen.reset(new TstarTstarGenHists(ctx, "AfterTrigger_gen"));


  // ###### 4. Init Handles ######
  h_is_muevt = ctx.get_handle<bool>("is_muevt");
  h_evt_weight = ctx.get_handle<double>("evt_weight");
  h_primlep = ctx.get_handle<FlavorParticle>("PrimaryLepton");
  h_neutrino = ctx.declare_event_output<LorentzVector>("neutrino");
  h_ST = ctx.declare_event_output<double>("ST");

}


bool TstarTstarSelectionModule::process(Event & event) {

  // debug status message
  if(debug) cout << "TstarTstarSelectionModule: Starting to process event (runid, eventid) = (" << event.run << ", " << event.event << "); weight = " << event.weight << endl;

  // reapply event weights from handle
  event.weight = event.get(h_evt_weight);
  if(debug) cout << "weights applied." << endl;

  // Fill ttgen object for correct matching check, etc
  if(is_MC) ttgenprod->process(event);

  // set primary lepton
  reco_primlep->process(event);

  // hists before anything happened
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


  // #################
  // ### Selection ###
  // #################

  // ###### Btag Selection ######
  bool pass_btagcut = false;
  for (const auto & jet: *event.jets){
    if(jet.btag_DeepCSV() > 0.2219) pass_btagcut = true;
  }
  if(pass_btagcut) {
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

    if(pass_btagcut) { // only fill these for btag cut passes
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
      // btagging reverse region plots
    }

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
    st += event.MET.pt();
    event.set(h_ST, st);

    // st cut
    if(st < 500) return false;

    if(pass_btagcut) { // only fill these for btag cut passes
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
      // TODO control region plots here!
    }


    // #######################
    // ### Trigger studies ###
    // #######################

    if(pass_btagcut) { // only fill these for btag cut passes
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
      if(pass_trigger) {
        h_trigger->fill(event);
        h_trigger_gen->fill(event);
        if(event.get(h_is_muevt)) h_trigger_mu->fill(event);
        else h_trigger_ele->fill(event);
        if(debug) cout<<"Filled hists after Trigger"<<endl;
      }
    }

    if(pass_btagcut) { // only fill these for btag cut passes
      // some final plot for comparison
      h_afterSelection->fill(event);
      h_afterSelection_gen->fill(event);
      h_afterSelection_genmatch->fill(event);
    } else {
      //h_nobtagcontrolregion->fill();
    }

    // at the moment control region ends here!
    if(pass_btagcut) return true;
    else return false;

}

// as we want to run the ExampleCycleNew directly with AnalysisModuleRunner,
// make sure the TstarTstarSelectionModule is found by class name. This is ensured by this macro:
UHH2_REGISTER_ANALYSIS_MODULE(TstarTstarSelectionModule)

}
