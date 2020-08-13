#include <iostream>
#include <memory>
#include <string>

// For DNN
#include "UHH2/TstarTstar/include/NeuralNetworkModules.h"

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
#include "UHH2/common/include/MCWeight.h"
#include "UHH2/TstarTstar/include/TstarTstarSelections.h"
#include "UHH2/TstarTstar/include/TstarTstarHists.h"
#include "UHH2/TstarTstar/include/TstarTstarDNNHists.h"
#include "UHH2/TstarTstar/include/TstarTstarDNNInputHists.h"
#include "UHH2/TstarTstar/include/TstarTstarRecoTstarHists.h"
#include "UHH2/TstarTstar/include/TstarTstarAllGenHists.h"
#include "UHH2/TstarTstar/include/TstarTstarGenHists.h"
#include "UHH2/TstarTstar/include/TstarTstarGenRecoMatchedHists.h"
#include "UHH2/TstarTstar/include/TstarTstarReconstructionModules.h"
#include "UHH2/TstarTstar/include/ReconstructionTstarHypothesis.h"
#include "UHH2/TstarTstar/include/TstarTstarGenMatch.h"
#include "UHH2/HOTVR/include/HOTVRIds.h"

using namespace std;
using namespace uhh2;

namespace uhh2 {

  float inv_mass(const LorentzVector& p4){ return p4.isTimelike() ? p4.mass() : -sqrt(-p4.mass2()); }

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

  unique_ptr<TTbarSemiLepMatchableSelection> TTbarSemiLepMatchable_selection;
  unique_ptr<Selection> toptagevt_sel;

  // ##### Histograms #####
  // Store the Hists collection as member variables. Again, use unique_ptr to avoid memory leaks.

  std::unique_ptr<Hists> h_beforeReco, h_beforeReco_ttag, h_beforeReco_nottag, h_afterPrimlep, h_afterHypCreation, h_afterReco_Full, h_afterReco_ttag, h_afterReco_nottag, h_lowchi2, h_highchi2, h_afterGEN, h_afterGEN_onlyttbar, h_afterGEN_onlyttbar_ttag, h_afterGEN_onlyttbar_nottag, h_notReconstructible, h_notReconstructible_ttag, h_notReconstructible_nottag, h_passFatJetSel;
  std::unique_ptr<TstarTstarRecoTstarHists> h_RecoPlots_Full, h_RecoPlots_ttag, h_RecoPlots_nottag, h_RecoPlots_lowchi2, h_RecoPlots_highchi2;
  std::unique_ptr<TstarTstarRecoTstarHists> h_RecoPlots_GEN, h_RecoPlots_GEN_onlyttbar, h_RecoPlots_GEN_onlyttbar_ttag, h_RecoPlots_GEN_onlyttbar_nottag;

  std::unique_ptr<Hists> h_afterPrimlep_mu, h_afterPrimlep_mu_lowpt, h_afterPrimlep_mu_highpt, h_afterPrimlep_ele, h_afterPrimlep_ele_lowpt, h_afterPrimlep_ele_highpt;

  std::unique_ptr<Hists> h_AfterDNNcut, h_notDNNcut, h_ST_reweighted, h_ST_reweighted_2;
  std::unique_ptr<TstarTstarRecoTstarHists> h_RecoPlots_AfterDNNcut, h_RecoPlots_notDNNcut;

  std::unique_ptr<Hists> h_GEN_Hists, h_GEN_Hists_pre;
  std::unique_ptr<Hists> h_DNN_Inputs, h_DNN_Inputs_reweighted, h_DNN_Inputs_reweighted_2;

  std::unique_ptr<Hists> h_top_gluon_checks, h_top_gluon_checks_reweighted, h_top_gluon_checks_reweighted_2;
  std::unique_ptr<TstarTstarDNNHists> h_DNN_Hists, h_DNN_Hists_reweighted, h_DNN_Hists_reweighted_2, h_DNN_Hists_AfterDNNCut, h_DNN_Hists_lowpt, h_DNN_Hists_medpt, h_DNN_Hists_highpt;

  // Bools for Debugging/Options
  bool debug = false;

  // DNN stuff
  bool outputDNNvalues = true;
  bool includeDNNmodel = true;
  bool do_masspoint = false;
  std::unique_ptr<NeuralNetworkInputWriter> DNN_InputWriter;
  std::unique_ptr<NeuralNetworkIncluder> DNN_Includer;
  uhh2::Event::Handle<double> h_DNN_output;
  uhh2::Event::Handle<bool> h_do_masspoint;
  uhh2::Event::Handle<double> h_ST;

  // primlep
  uhh2::Event::Handle<FlavorParticle> h_primlep;

  // Modules
  std::unique_ptr<uhh2::AnalysisModule> ttgenprod;
  std::unique_ptr<uhh2::AnalysisModule> reco_primlep;
  std::unique_ptr<TstarTstarGenDiscriminator> genmatcher;
  std::unique_ptr<TstarTstarGenDiscriminator> genmatcher_onlyttbar;
  std::unique_ptr<TstarTstar_tgtg_TopTag_Reconstruction> TstarTstarHypCreator;
  std::unique_ptr<TstarTstar_Discrimination> TstarTstarHypSelector;

  // Handles
  uhh2::Event::Handle<TTbarGen> h_ttbargen;
  uhh2::Event::Handle<int> h_flag_toptagevent;
  uhh2::Event::Handle<int> h_flag_muonevent;
  uhh2::Event::Handle<std::vector<ReconstructionTstarHypothesis>> h_tstartstar_hyp_vector;
  uhh2::Event::Handle<ReconstructionTstarHypothesis> h_tstartstar_hyp;

  uhh2::Event::Handle<double> h_evt_weight;
  uhh2::Event::Handle<double> h_ST_weight;
  uhh2::Event::Handle<double> h_ST_weight_2;

  // bools for channel and stuff. will be read in later
  bool isTrigger;
  bool is_MC;

  TH1D* ST_ratio;
  TH1D* ST_sig;
  TH1D* ST_bkg;
  bool is_TTbar;
  bool is_Signal;
  bool is_Data;
};


TstarTstarMCStudyModule::TstarTstarMCStudyModule(Context & ctx){

  reco_primlep.reset(new PrimaryLepton(ctx));

  if(debug) {
    cout << "Hello World from TstarTstarMCStudyModule!" << endl;

    // If running in SFrame, the keys "dataset_version", "dataset_type", "dataset_lumi",
    // and "target_lumi" are set to the according values in the xml file. For CMSSW, these are
    // not set automatically, but can be set in the python config file.
    for(auto & kv : ctx.get_all()){
        cout << " " << kv.first << " = " << kv.second << endl;
    }
   }

  // 0. Reading in whether MC and if so, which channel
  is_MC = ctx.get("dataset_type") == "MC";

  // Prepare GEN
  if(is_MC){
    ttgenprod.reset(new TTbarGenProducer(ctx, "ttbargen", false));
    h_ttbargen = ctx.get_handle<TTbarGen>("ttbargen");
    genmatcher.reset(new TstarTstarGenDiscriminator(ctx));
    genmatcher_onlyttbar.reset(new TstarTstarGenDiscriminator(ctx, true));
  }

  //TopTag
  TopJetId topjetID = AndId<TopJet>(HOTVRTopTag(), Tau32Groomed(0.56));
  toptagevt_sel.reset(new TopTagEventSelection(topjetID));
  h_flag_toptagevent = ctx.declare_event_output<int>("flag_toptagevent");

  // Check which lepton is present and save in
  h_flag_muonevent = ctx.declare_event_output<int>("flag_muonevent");

  // 3. Set up Hists classes:
  h_beforeReco.reset(new TstarTstarHists(ctx, "beforeReco"));
  h_beforeReco_ttag.reset(new TstarTstarHists(ctx, "beforeReco_ttag"));
  h_beforeReco_nottag.reset(new TstarTstarHists(ctx, "beforeReco_nottag"));
  h_afterPrimlep.reset(new TstarTstarHists(ctx, "afterPrimlep"));
  h_passFatJetSel.reset(new TstarTstarHists(ctx, "passFatJetSel"));
  h_afterHypCreation.reset(new TstarTstarHists(ctx, "afterHypCreation"));
  h_afterReco_Full.reset(new TstarTstarHists(ctx, "AfterReco_Full"));
  h_afterReco_ttag.reset(new TstarTstarHists(ctx, "AfterReco_ttag"));
  h_afterReco_nottag.reset(new TstarTstarHists(ctx, "AfterReco_nottag"));
  h_lowchi2.reset(new TstarTstarHists(ctx, "Lowchi2"));
  h_highchi2.reset(new TstarTstarHists(ctx, "Highchi2"));
  h_afterGEN.reset(new TstarTstarHists(ctx, "AfterGEN"));
  h_afterGEN_onlyttbar.reset(new TstarTstarHists(ctx, "AfterGEN_onlyttbar"));
  h_afterGEN_onlyttbar_ttag.reset(new TstarTstarHists(ctx, "AfterGEN_onlyttbar_ttag"));
  h_afterGEN_onlyttbar_nottag.reset(new TstarTstarHists(ctx, "AfterGEN_onlyttbar_nottag"));
  h_notReconstructible.reset(new TstarTstarHists(ctx, "notReconstructible"));
  h_notReconstructible_ttag.reset(new TstarTstarHists(ctx, "notReconstructible_ttag"));
  h_notReconstructible_nottag.reset(new TstarTstarHists(ctx, "notReconstructible_nottag"));

  h_ST_reweighted.reset(new TstarTstarHists(ctx, "STreweighted"));
  h_ST_reweighted_2.reset(new TstarTstarHists(ctx, "STreweighted_2"));

  h_afterPrimlep_mu.reset(new TstarTstarHists(ctx, "afterPrimlep_mu"));
  h_afterPrimlep_mu_lowpt.reset(new TstarTstarHists(ctx, "afterPrimlep_mu_lowpt"));
  h_afterPrimlep_mu_highpt.reset(new TstarTstarHists(ctx, "afterPrimlep_mu_highpt"));
  h_afterPrimlep_ele.reset(new TstarTstarHists(ctx, "afterPrimlep_ele"));
  h_afterPrimlep_ele_lowpt.reset(new TstarTstarHists(ctx, "afterPrimlep_ele_lowpt"));
  h_afterPrimlep_ele_highpt.reset(new TstarTstarHists(ctx, "afterPrimlep_ele_highpt"));

  h_RecoPlots_Full.reset(new TstarTstarRecoTstarHists(ctx, "RecoPlots_Full"));
  h_RecoPlots_ttag.reset(new TstarTstarRecoTstarHists(ctx, "RecoPlots_ttag"));
  h_RecoPlots_nottag.reset(new TstarTstarRecoTstarHists(ctx, "RecoPlots_nottag"));
  h_RecoPlots_lowchi2.reset(new TstarTstarRecoTstarHists(ctx, "RecoPlots_lowchi2"));
  h_RecoPlots_highchi2.reset(new TstarTstarRecoTstarHists(ctx, "RecoPlots_highchi2"));
  h_RecoPlots_GEN.reset(new TstarTstarRecoTstarHists(ctx, "RecoPlots_GEN"));
  h_RecoPlots_GEN_onlyttbar.reset(new TstarTstarRecoTstarHists(ctx, "RecoPlots_GEN_onlyttbar"));
  h_RecoPlots_GEN_onlyttbar_ttag.reset(new TstarTstarRecoTstarHists(ctx, "RecoPlots_GEN_onlyttbar_ttag"));
  h_RecoPlots_GEN_onlyttbar_nottag.reset(new TstarTstarRecoTstarHists(ctx, "RecoPlots_GEN_onlyttbar_nottag"));

  h_AfterDNNcut.reset(new TstarTstarHists(ctx, "AfterDNNcut"));
  h_notDNNcut.reset(new TstarTstarHists(ctx, "notDNNcut"));
  h_RecoPlots_AfterDNNcut.reset(new TstarTstarRecoTstarHists(ctx, "RecoPlots_AfterDNNcut"));
  h_RecoPlots_notDNNcut.reset(new TstarTstarRecoTstarHists(ctx, "RecoPlots_notDNNcut"));

  h_GEN_Hists_pre.reset(new TstarTstarGenHists(ctx, "GEN_Hists_beforeReco"));
  h_GEN_Hists.reset(new TstarTstarGenHists(ctx, "GEN_Hists_AfterReco"));

  h_DNN_Hists.reset(new TstarTstarDNNHists(ctx, "DNN_Hists"));
  h_DNN_Hists_reweighted.reset(new TstarTstarDNNHists(ctx, "DNN_Hists_reweighted"));
  h_DNN_Hists_reweighted_2.reset(new TstarTstarDNNHists(ctx, "DNN_Hists_reweighted_2"));
  h_DNN_Hists_AfterDNNCut.reset(new TstarTstarDNNHists(ctx, "DNN_Hists_AfterDNNCut"));
  h_DNN_Hists_lowpt.reset(new TstarTstarDNNHists(ctx, "DNN_Hists_lowpt"));
  h_DNN_Hists_medpt.reset(new TstarTstarDNNHists(ctx, "DNN_Hists_medpt"));
  h_DNN_Hists_highpt.reset(new TstarTstarDNNHists(ctx, "DNN_Hists_highpt"));

  h_DNN_Inputs.reset(new TstarTstarDNNInputHists(ctx, "DNN_Inputs"));
  h_DNN_Inputs_reweighted.reset(new TstarTstarDNNInputHists(ctx, "DNN_Inputs_reweighted"));
  h_DNN_Inputs_reweighted_2.reset(new TstarTstarDNNInputHists(ctx, "DNN_Inputs_reweighted_2"));

  h_top_gluon_checks.reset(new TstarTstarAllGenHists(ctx, "Top_check"));
  h_top_gluon_checks_reweighted.reset(new TstarTstarAllGenHists(ctx, "Top_check_reweighted"));
  h_top_gluon_checks_reweighted_2.reset(new TstarTstarAllGenHists(ctx, "Top_check_reweighted_2"));

  //4. Set up ttbar reconstruction
  h_tstartstar_hyp_vector = ctx.get_handle<std::vector<ReconstructionTstarHypothesis>>("TstarTstar_Hyp_Vector");
  h_tstartstar_hyp = ctx.get_handle<ReconstructionTstarHypothesis>("TstarTstar_Hyp");

  TstarTstarHypCreator.reset(new TstarTstar_tgtg_TopTag_Reconstruction(ctx, NeutrinoReconstruction, topjetID));
  TstarTstarHypSelector.reset(new TstarTstar_Discrimination(ctx));

  // 5. Handles for DNN
  if(outputDNNvalues){
    DNN_InputWriter.reset(new NeuralNetworkInputWriter(ctx));
    h_do_masspoint = ctx.get_handle<bool>("do_masspoint");
    h_ST = ctx.declare_event_output<double>("ST");
  }
  if(includeDNNmodel){
    DNN_Includer.reset(new NeuralNetworkIncluder(ctx, do_masspoint));
    h_DNN_output = ctx.get_handle<double>("DNN_output");
  }

  h_primlep = ctx.get_handle<FlavorParticle>("PrimaryLepton");

  h_evt_weight = ctx.get_handle<double>("evt_weight");
  h_ST_weight = ctx.declare_event_output<double>("ST_weight");
  h_ST_weight_2 = ctx.declare_event_output<double>("ST_weight_flat");

  TFile *f = new TFile("/nfs/dust/cms/user/flabe/CMSSW/CMSSW_10_2_10/src/UHH2/MLCorner/TstarNN/ST_weights.root");
  ST_ratio = (TH1D*)f->Get("ST_ratio");
  ST_sig = (TH1D*)f->Get("ST_sig");
  ST_bkg = (TH1D*)f->Get("ST_bkg");

  is_TTbar = (ctx.get("dataset_version").find("TT") != std::string::npos);
  is_Signal = (ctx.get("dataset_version").find("Tstar") != std::string::npos);
}


bool TstarTstarMCStudyModule::process(Event & event) {

  if(debug){cout << endl << "TstarTstarMCStudyModule: Starting to process event (runid, eventid) = (" << event.run << ", " << event.event << "); weight = " << event.weight << endl;}

  // reapply weights
  event.weight = event.get(h_evt_weight);
  if(debug) cout << "weights applied." << endl;

  h_top_gluon_checks->fill(event);

  // ST reweighting
  double st_jets = 0;
  double ST_weight = 0;
  double ST_weight_2 = 0;
  for(const auto & jet : *event.topjets) st_jets += jet.pt();
  if(is_TTbar){
    if(st_jets < 3000) ST_weight = 1/(ST_ratio->GetBinContent(ST_ratio->GetXaxis()->FindBin(st_jets)));
    if(st_jets < 3000) ST_weight_2 = 1/(ST_bkg->GetBinContent(ST_bkg->GetXaxis()->FindBin(st_jets)));
    event.set(h_ST_weight, ST_weight*event.get(h_evt_weight));
    event.set(h_ST_weight_2, ST_weight_2);
    event.weight *= ST_weight;
    h_ST_reweighted->fill(event);
    h_top_gluon_checks_reweighted->fill(event);
    event.weight = event.get(h_evt_weight);
    event.weight = ST_weight_2;
    h_ST_reweighted_2->fill(event);
    h_top_gluon_checks_reweighted_2->fill(event);
    event.weight = event.get(h_evt_weight);
  }
  else if(is_Signal){
    double ST_weight_2 = 0;
    if(st_jets < 3000) ST_weight_2 = 1/(ST_sig->GetBinContent(ST_sig->GetXaxis()->FindBin(st_jets)));
    event.set(h_ST_weight, event.get(h_evt_weight));
    event.set(h_ST_weight_2, ST_weight_2);
    h_ST_reweighted->fill(event);
    h_top_gluon_checks_reweighted->fill(event);
    event.weight = ST_weight_2;
    h_ST_reweighted_2->fill(event);
    h_top_gluon_checks_reweighted_2->fill(event);
    event.weight = event.get(h_evt_weight);
  }
  else {
    event.set(h_ST_weight, event.get(h_evt_weight));
    event.set(h_ST_weight_2, event.get(h_evt_weight));
    h_ST_reweighted->fill(event);
    h_top_gluon_checks_reweighted->fill(event);
    h_ST_reweighted_2->fill(event);
    h_top_gluon_checks_reweighted_2->fill(event);
  }

  if(is_MC) ttgenprod->process(event);

  // check lepton channel
  const bool muon_evt = (event.muons->size() == 1);
  event.set(h_flag_muonevent, int(muon_evt));

  h_beforeReco->fill(event);
  h_GEN_Hists_pre->fill(event);

  // check wether ttag is present
  const bool pass_ttag = toptagevt_sel->passes(event);
  event.set(h_flag_toptagevent, int(pass_ttag));
  if(debug) cout << "Done checking ttag: " << pass_ttag << endl;

  if(pass_ttag){h_beforeReco_ttag->fill(event);}
  else {h_beforeReco_nottag->fill(event);}

  // Fix empty handle problem
  // (Dummy values are filled into handles so that they are not empty)
  event.set(h_tstartstar_hyp, ReconstructionTstarHypothesis());
  if(debug){cout << "Finished initialization of Handle Variables" << endl;}

  // ########### TstarTstarReco! ############
  if(debug) cout<<"Starting TstarTstar Reconstruction part"<<endl;

  reco_primlep->process(event);//set "primary lepton"

  bool TstarHypsCreated = false;
  bool bestHypFound = false;

  h_afterPrimlep->fill(event);
  if(event.get(h_flag_muonevent)){
    h_afterPrimlep_mu->fill(event);
    if(event.get(h_primlep).pt()<60) h_afterPrimlep_mu_lowpt->fill(event);
    else h_afterPrimlep_mu_highpt->fill(event);
  }
  else {
    h_afterPrimlep_ele->fill(event);
    if(event.get(h_primlep).pt()<120) h_afterPrimlep_ele_lowpt->fill(event);
    else h_afterPrimlep_ele_highpt->fill(event);
  }

  // fat jet selection
  bool pass_fat_njet = (event.topjets->size()>2);

  if(debug){cout << "Starting to construct all TstarTstar Hypothesiseseses" << endl;}
  if(pass_fat_njet){
    h_passFatJetSel->fill(event);
    TstarHypsCreated = TstarTstarHypCreator->process(event);
  }
  if(TstarHypsCreated){
    h_afterHypCreation->fill(event);
    bestHypFound = TstarTstarHypSelector->process(event);
    if(bestHypFound){
      h_RecoPlots_Full->fill(event);
      h_afterReco_Full->fill(event);
      h_GEN_Hists->fill(event);
      if(pass_ttag){
	       h_RecoPlots_ttag->fill(event);
	       h_afterReco_ttag->fill(event);
      }
      else{
	       h_RecoPlots_nottag->fill(event);
	       h_afterReco_nottag->fill(event);
      }
      if(event.get(h_tstartstar_hyp).chi2()<50){
        h_RecoPlots_lowchi2->fill(event);
        h_lowchi2->fill(event);
      }
      else {
        h_RecoPlots_highchi2->fill(event);
        h_highchi2->fill(event);
      }
    }
  }
  if(!TstarHypsCreated || !bestHypFound){
    h_notReconstructible->fill(event);
    if(pass_ttag) h_notReconstructible_ttag->fill(event);
    else h_notReconstructible_nottag->fill(event);
  }

  // ####################
  // ####################
  // ####################

  if(is_MC && TstarHypsCreated && false){
    if(debug){ cout << "Doing GEN matching check" << endl;}
    // ##### GEN Matching
    {
      ReconstructionTstarHypothesis hyp_tmp = event.get(h_tstartstar_hyp);
      if(genmatcher->process(event)){
      	if(debug)cout << "GEN MATCHED!" << endl;
      	h_RecoPlots_GEN->fill(event);
      	h_afterGEN->fill(event);
      	event.set(h_tstartstar_hyp, hyp_tmp);
      }
    }
    {
      ReconstructionTstarHypothesis hyp_tmp = event.get(h_tstartstar_hyp);
      if(genmatcher_onlyttbar->process(event)){
      	if(debug)cout << "GEN MATCHED!" << endl;
      	h_RecoPlots_GEN_onlyttbar->fill(event);
      	h_afterGEN_onlyttbar->fill(event);
      	if(pass_ttag){
      	  h_RecoPlots_GEN_onlyttbar_ttag->fill(event);
      	  h_afterGEN_onlyttbar_ttag->fill(event);
      	}
      	if(pass_ttag){
      	  h_RecoPlots_GEN_onlyttbar_nottag->fill(event);
      	  h_afterGEN_onlyttbar_nottag->fill(event);
      	}
      	event.set(h_tstartstar_hyp, hyp_tmp);
        }
      }
    }

    // ####################
    // ####################
    // ####################

    if(debug) cout << "Start DNN stuff" << endl;
    // Filling output for DNN
    if(outputDNNvalues){
      event.set(h_do_masspoint, do_masspoint);
      double st_jets = 0.;
      for(const auto & jet : *event.topjets) st_jets += jet.pt();
      event.set(h_ST, st_jets);
      DNN_InputWriter->process(event);
      h_DNN_Inputs->fill(event);
      if(is_TTbar) event.weight *= ST_weight;
      h_DNN_Inputs_reweighted->fill(event);
      event.weight = event.get(h_evt_weight);
      if(is_TTbar || is_Signal) event.weight = ST_weight_2;
      h_DNN_Inputs_reweighted_2->fill(event);
      event.weight = event.get(h_evt_weight);
    }
    if(includeDNNmodel){
      DNN_Includer->process(event);
      h_DNN_Hists->fill(event);
      if(is_TTbar) event.weight *= ST_weight;
      h_DNN_Hists_reweighted->fill(event);
      event.weight = event.get(h_evt_weight);
      if(is_TTbar || is_Signal) event.weight = ST_weight_2;
      h_DNN_Hists_reweighted_2->fill(event);
      event.weight = event.get(h_evt_weight);

      if(event.topjets->at(0).pt() < 500) h_DNN_Hists_lowpt->fill(event);
      else if (event.topjets->at(0).pt() < 1000) h_DNN_Hists_medpt->fill(event);
      else h_DNN_Hists_highpt->fill(event);

      if(is_MC) {
        if(event.get(h_DNN_output) > 0.4) {
          h_AfterDNNcut->fill(event);
          h_RecoPlots_AfterDNNcut->fill(event);
          h_DNN_Hists_AfterDNNCut->fill(event);
        }
        else {
          h_notDNNcut->fill(event);
          h_RecoPlots_notDNNcut->fill(event);
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
