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
#include <UHH2/common/include/TriggerSelection.h>
#include "UHH2/common/include/TTbarGen.h"
#include "UHH2/common/include/TopJetIds.h"
#include "UHH2/common/include/MCWeight.h"

// TstarTstar stuff
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
#include "UHH2/TstarTstar/include/NeuralNetworkModules.h"

// other stuff
#include "UHH2/HOTVR/include/HOTVRIds.h"

using namespace std;
using namespace uhh2;

namespace uhh2 {

// quick method to calculate inv_mass
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

  // ###### Modules ######
  unique_ptr<Selection> toptagevt_sel;
  unique_ptr<TTbarSemiLepMatchableSelection> TTbarSemiLepMatchable_selection;
  unique_ptr<NeuralNetworkInputWriter> DNN_InputWriter;

  std::unique_ptr<uhh2::AnalysisModule> ttgenprod;
  std::unique_ptr<uhh2::AnalysisModule> reco_primlep;
  std::unique_ptr<TstarTstarGenDiscriminator> genmatcher;
  std::unique_ptr<TstarTstarGenDiscriminator> genmatcher_onlyttbar;
  std::unique_ptr<TstarTstar_tgtg_TopTag_Reconstruction> TstarTstarHypCreator;
  std::unique_ptr<TstarTstar_tgtg_AK4_Reconstruction> TstarTstarHypCreatorAK4;
  std::unique_ptr<TstarTstar_Discrimination> TstarTstarHypSelector;
  std::unique_ptr<TstarTstar_Discrimination> TstarTstarHypSelectorAK4;


  // ##### Histograms #####
  std::unique_ptr<Hists> h_beforeReco, h_beforeReco_ttag, h_beforeReco_nottag, h_beforeReco_mu, h_beforeReco_mu_lowpt, h_beforeReco_mu_highpt, h_beforeReco_ele, h_beforeReco_ele_lowpt, h_beforeReco_ele_highpt;
  std::unique_ptr<Hists> h_passFatJetSel;
  std::unique_ptr<Hists> h_beforeReco_gen;

  std::unique_ptr<Hists> h_afterHypCreation, h_afterReco_Full, h_afterReco_ttag, h_afterReco_nottag, h_afterReco_lowchi2, h_afterReco_highchi2;
  std::unique_ptr<Hists> h_afterReco_gen;
  std::unique_ptr<Hists> h_RecoPlots_Full, h_RecoPlots_ttag, h_RecoPlots_nottag, h_RecoPlots_lowchi2, h_RecoPlots_highchi2;
  std::unique_ptr<Hists> h_notReconstructible, h_notReconstructible_ttag, h_notReconstructible_nottag;

  std::unique_ptr<Hists> h_afterHypCreation_AK4, h_afterReco_Full_AK4, h_afterReco_ttag_AK4, h_afterReco_nottag_AK4, h_afterReco_lowchi2_AK4, h_afterReco_highchi2_AK4;
  std::unique_ptr<Hists> h_RecoPlots_Full_AK4, h_RecoPlots_ttag_AK4, h_RecoPlots_nottag_AK4, h_RecoPlots_lowchi2_AK4, h_RecoPlots_highchi2_AK4;
  std::unique_ptr<Hists> h_notReconstructible_AK4, h_notReconstructible_ttag_AK4, h_notReconstructible_nottag_AK4;

  std::unique_ptr<Hists> h_afterGEN, h_afterGEN_onlyttbar, h_afterGEN_onlyttbar_ttag, h_afterGEN_onlyttbar_nottag;
  std::unique_ptr<Hists> h_RecoPlots_GEN, h_RecoPlots_GEN_onlyttbar, h_RecoPlots_GEN_onlyttbar_ttag, h_RecoPlots_GEN_onlyttbar_nottag;

  std::unique_ptr<Hists> h_ST_reweighted, h_ST_reweighted_2;
  std::unique_ptr<Hists> h_top_gluon_checks, h_top_gluon_checks_reweighted, h_top_gluon_checks_reweighted_2;
  std::unique_ptr<Hists> h_DNN_Inputs, h_DNN_Inputs_reweighted, h_DNN_Inputs_reweighted_2;


  // ###### Handles ######
  uhh2::Event::Handle<double> h_evt_weight;
  uhh2::Event::Handle<FlavorParticle> h_primlep;
  uhh2::Event::Handle<double> h_ST;

  // for reconstruction
  uhh2::Event::Handle<TTbarGen> h_ttbargen;
  uhh2::Event::Handle<int> h_flag_toptagevent;
  uhh2::Event::Handle<int> h_flag_muonevent;
  uhh2::Event::Handle<std::vector<ReconstructionTstarHypothesis>> h_tstartstar_hyp_vector;
  uhh2::Event::Handle<ReconstructionTstarHypothesis> h_tstartstar_hyp;

  // for DNN output
  uhh2::Event::Handle<bool> h_do_masspoint;
  uhh2::Event::Handle<double> h_ST_weight;
  uhh2::Event::Handle<double> h_ST_weight_2;


  // ###### Control Switches ######
  bool debug = false;
  bool outputDNNvalues = true;
  bool do_masspoint = false;


  // ###### other needed definitions ######
  bool isTrigger;
  bool is_MC;

  TH1D* ST_ratio;
  TH1D* ST_sig;
  TH1D* ST_bkg;
  TH1D* ST_sig_2;
  TH1D* ST_bkg_2;
  bool is_TTbar;
  bool is_Signal;
  bool is_Data;
};


TstarTstarMCStudyModule::TstarTstarMCStudyModule(Context & ctx){

  // debug messagt
  if(debug) {
    cout << "Hello World from TstarTstarMCStudyModule!" << endl;
    // If running in SFrame, the keys "dataset_version", "dataset_type", "dataset_lumi",
    // and "target_lumi" are set to the according values in the xml file. For CMSSW, these are
    // not set automatically, but can be set in the python config file.
    for(auto & kv : ctx.get_all()){
        cout << " " << kv.first << " = " << kv.second << endl;
    }
  }

  // ###### 0. Setting Variables ######
  // MC or real data
  is_MC = ctx.get("dataset_type") == "MC";


  // ###### 1. Set up modules ######
  // primary lepton
  reco_primlep.reset(new PrimaryLepton(ctx));

  // GEN things
  if(is_MC){
    ttgenprod.reset(new TTbarGenProducer(ctx, "ttbargen", false));
    genmatcher.reset(new TstarTstarGenDiscriminator(ctx));
    genmatcher_onlyttbar.reset(new TstarTstarGenDiscriminator(ctx, true));
  }

  // top tag definition
  TopJetId topjetID = AndId<TopJet>(HOTVRTopTag(), Tau32Groomed(0.56));
  toptagevt_sel.reset(new TopTagEventSelection(topjetID));

  // Reconstruction
  TstarTstarHypCreator.reset(new TstarTstar_tgtg_TopTag_Reconstruction(ctx, NeutrinoReconstruction, topjetID));
  TstarTstarHypCreatorAK4.reset(new TstarTstar_tgtg_AK4_Reconstruction(ctx, NeutrinoReconstruction, topjetID));
  TstarTstarHypSelector.reset(new TstarTstar_Discrimination(ctx));
  TstarTstarHypSelectorAK4.reset(new TstarTstar_Discrimination(ctx));


  // ###### 3. Set up histograms ######
  // before Reconstruction
  h_beforeReco.reset(new TstarTstarHists(ctx, "beforeReco"));
  h_beforeReco_ttag.reset(new TstarTstarHists(ctx, "beforeReco_ttag"));
  h_beforeReco_nottag.reset(new TstarTstarHists(ctx, "beforeReco_nottag"));
  h_beforeReco_mu.reset(new TstarTstarHists(ctx, "beforeReco_mu"));
  h_beforeReco_mu_lowpt.reset(new TstarTstarHists(ctx, "beforeReco_mu_lowpt"));
  h_beforeReco_mu_highpt.reset(new TstarTstarHists(ctx, "beforeReco_mu_highpt"));
  h_beforeReco_ele.reset(new TstarTstarHists(ctx, "beforeReco_ele"));
  h_beforeReco_ele_lowpt.reset(new TstarTstarHists(ctx, "beforeReco_ele_lowpt"));
  h_beforeReco_ele_highpt.reset(new TstarTstarHists(ctx, "beforeReco_ele_highpt"));
  h_beforeReco_gen.reset(new TstarTstarGenHists(ctx, "beforeReco_gen"));

  // After additional fat jet selection
  h_passFatJetSel.reset(new TstarTstarHists(ctx, "AfterFatJets"));

  // Reco regular
  h_afterHypCreation.reset(new TstarTstarHists(ctx, "AfterHypCreation"));

  h_afterReco_Full.reset(new TstarTstarHists(ctx, "AfterReco_Full"));
  h_afterReco_ttag.reset(new TstarTstarHists(ctx, "AfterReco_ttag"));
  h_afterReco_nottag.reset(new TstarTstarHists(ctx, "AfterReco_nottag"));
  h_afterReco_lowchi2.reset(new TstarTstarHists(ctx, "AfterReco_lowchi2"));
  h_afterReco_highchi2.reset(new TstarTstarHists(ctx, "AfterReco_highchi2"));
  h_afterReco_gen.reset(new TstarTstarGenHists(ctx, "AfterReco_gen"));

  h_RecoPlots_Full.reset(new TstarTstarRecoTstarHists(ctx, "RecoPlots_Full"));
  h_RecoPlots_ttag.reset(new TstarTstarRecoTstarHists(ctx, "RecoPlots_ttag"));
  h_RecoPlots_nottag.reset(new TstarTstarRecoTstarHists(ctx, "RecoPlots_nottag"));
  h_RecoPlots_lowchi2.reset(new TstarTstarRecoTstarHists(ctx, "RecoPlots_lowchi2"));
  h_RecoPlots_highchi2.reset(new TstarTstarRecoTstarHists(ctx, "RecoPlots_highchi2"));

  h_notReconstructible.reset(new TstarTstarHists(ctx, "notReconstructible"));
  h_notReconstructible_ttag.reset(new TstarTstarHists(ctx, "notReconstructible_ttag"));
  h_notReconstructible_nottag.reset(new TstarTstarHists(ctx, "notReconstructible_nottag"));

  // Reco AK4
  h_afterHypCreation_AK4.reset(new TstarTstarHists(ctx, "afterHypCreation_AK4"));

  h_afterReco_Full_AK4.reset(new TstarTstarHists(ctx, "AfterReco_Full_AK4"));
  h_afterReco_ttag_AK4.reset(new TstarTstarHists(ctx, "AfterReco_ttag_AK4"));
  h_afterReco_nottag_AK4.reset(new TstarTstarHists(ctx, "AfterReco_nottag_AK4"));
  h_afterReco_lowchi2_AK4.reset(new TstarTstarHists(ctx, "AfterReco_lowchi2_AK4"));
  h_afterReco_highchi2_AK4.reset(new TstarTstarHists(ctx, "AfterReco_highchi2_AK4"));

  h_RecoPlots_Full_AK4.reset(new TstarTstarRecoTstarHists(ctx, "RecoPlots_Full_AK4"));
  h_RecoPlots_ttag_AK4.reset(new TstarTstarRecoTstarHists(ctx, "RecoPlots_ttag_AK4"));
  h_RecoPlots_nottag_AK4.reset(new TstarTstarRecoTstarHists(ctx, "RecoPlots_nottag_AK4"));
  h_RecoPlots_lowchi2_AK4.reset(new TstarTstarRecoTstarHists(ctx, "RecoPlots_lowchi2_AK4"));
  h_RecoPlots_highchi2_AK4.reset(new TstarTstarRecoTstarHists(ctx, "RecoPlots_highchi2_AK4"));

  h_notReconstructible_AK4.reset(new TstarTstarHists(ctx, "notReconstructible_AK4"));
  h_notReconstructible_ttag_AK4.reset(new TstarTstarHists(ctx, "notReconstructible_ttag_AK4"));
  h_notReconstructible_nottag_AK4.reset(new TstarTstarHists(ctx, "notReconstructible_nottag_AK4"));

  // GEN matching
  h_RecoPlots_GEN.reset(new TstarTstarRecoTstarHists(ctx, "RecoPlots_GEN"));
  h_afterGEN.reset(new TstarTstarHists(ctx, "AfterGEN"));

  h_afterGEN_onlyttbar.reset(new TstarTstarHists(ctx, "AfterGEN_onlyttbar"));
  h_afterGEN_onlyttbar_ttag.reset(new TstarTstarHists(ctx, "AfterGEN_onlyttbar_ttag"));
  h_afterGEN_onlyttbar_nottag.reset(new TstarTstarHists(ctx, "AfterGEN_onlyttbar_nottag"));
  h_RecoPlots_GEN_onlyttbar.reset(new TstarTstarRecoTstarHists(ctx, "RecoPlots_GEN_onlyttbar"));
  h_RecoPlots_GEN_onlyttbar_ttag.reset(new TstarTstarRecoTstarHists(ctx, "RecoPlots_GEN_onlyttbar_ttag"));
  h_RecoPlots_GEN_onlyttbar_nottag.reset(new TstarTstarRecoTstarHists(ctx, "RecoPlots_GEN_onlyttbar_nottag"));


  // ST reweighting
  h_ST_reweighted.reset(new TstarTstarHists(ctx, "STreweighted"));
  h_ST_reweighted_2.reset(new TstarTstarHists(ctx, "STreweighted_2"));


  // GEN hists that work also for BG
  h_top_gluon_checks.reset(new TstarTstarAllGenHists(ctx, "Top_check"));
  h_top_gluon_checks_reweighted.reset(new TstarTstarAllGenHists(ctx, "Top_check_reweighted"));
  h_top_gluon_checks_reweighted_2.reset(new TstarTstarAllGenHists(ctx, "Top_check_reweighted_2"));


  // DNN hists
  h_DNN_Inputs.reset(new TstarTstarDNNInputHists(ctx, "DNN_Inputs"));
  h_DNN_Inputs_reweighted.reset(new TstarTstarDNNInputHists(ctx, "DNN_Inputs_reweighted"));
  h_DNN_Inputs_reweighted_2.reset(new TstarTstarDNNInputHists(ctx, "DNN_Inputs_reweighted_2"));



  // ###### 4. Init handles ######
  h_evt_weight = ctx.get_handle<double>("evt_weight");
  h_primlep = ctx.get_handle<FlavorParticle>("PrimaryLepton");
  h_flag_toptagevent = ctx.declare_event_output<int>("flag_toptagevent");
  h_flag_muonevent = ctx.declare_event_output<int>("flag_muonevent");

  if(is_MC) h_ttbargen = ctx.get_handle<TTbarGen>("ttbargen");

  // Reconstruction
  h_tstartstar_hyp_vector = ctx.get_handle<std::vector<ReconstructionTstarHypothesis>>("TstarTstar_Hyp_Vector");
  h_tstartstar_hyp = ctx.get_handle<ReconstructionTstarHypothesis>("TstarTstar_Hyp");

  // DNN output
  if(outputDNNvalues){
    DNN_InputWriter.reset(new NeuralNetworkInputWriter(ctx));
    h_do_masspoint = ctx.get_handle<bool>("do_masspoint");
    h_ST = ctx.declare_event_output<double>("ST");
  }

  h_ST_weight = ctx.declare_event_output<double>("ST_weight");
  h_ST_weight_2 = ctx.declare_event_output<double>("ST_weight_flat");



  // ###### 5. other definitions ######
  TFile *f = new TFile("MLCorner/TstarNN/ST_weights.root");
  ST_ratio = (TH1D*)f->Get("ST_ratio");
  ST_sig = (TH1D*)f->Get("ST_sig_split");
  ST_bkg = (TH1D*)f->Get("ST_bkg_split");
  ST_sig_2 = (TH1D*)f->Get("ST_sig_rebinned");
  ST_bkg_2 = (TH1D*)f->Get("ST_bkg_rebinned");

  is_TTbar = (ctx.get("dataset_version").find("TT") != std::string::npos);
  is_Signal = (ctx.get("dataset_version").find("Tstar") != std::string::npos);


}


bool TstarTstarMCStudyModule::process(Event & event) {

  // debug message
  if(debug){cout << endl << "TstarTstarMCStudyModule: Starting to process event (runid, eventid) = (" << event.run << ", " << event.event << "); weight = " << event.weight << endl;}

  // reapply weights
  event.weight = event.get(h_evt_weight);
  if(debug) cout << "weights applied." << endl;

  // set lepton channel
  const bool muon_evt = (event.muons->size() == 1);
  event.set(h_flag_muonevent, int(muon_evt));

  // set primlep
  reco_primlep->process(event);

  // check ttag
  const bool pass_ttag = toptagevt_sel->passes(event);
  event.set(h_flag_toptagevent, int(pass_ttag));

  // ttgen
  if(is_MC) ttgenprod->process(event);

  // prevent empty handles
  event.set(h_tstartstar_hyp, ReconstructionTstarHypothesis());

  // fill hists before things have happened
  event.weight = event.get(h_evt_weight);
  h_beforeReco->fill(event);
  h_beforeReco_gen->fill(event);
  if(pass_ttag) h_beforeReco_ttag->fill(event);
  else h_beforeReco_nottag->fill(event);
  if(event.get(h_flag_muonevent)){
    h_beforeReco_mu->fill(event);
    if(event.get(h_primlep).pt()<60) h_beforeReco_mu_lowpt->fill(event);
    else h_beforeReco_mu_highpt->fill(event);
  }
  else {
    h_beforeReco_ele->fill(event);
    if(event.get(h_primlep).pt()<120) h_beforeReco_ele_lowpt->fill(event);
    else h_beforeReco_ele_highpt->fill(event);
  }


  // ##################################################
  // ########### TstarTstar Reconstruction ############
  // ##################################################

  if(debug) cout<<"Starting TstarTstar Reconstruction part"<<endl;

  // ########### with gluons as HOTVR jets #############
  // needed additional fat jet selection
  bool pass_fat_njet = (event.topjets->size()>2);

  bool TstarHypsCreated = false;
  bool bestHypFound = false;

  if(debug) cout << "Starting to construct all TstarTstar Hypothesiseseses" << endl;
  if(pass_fat_njet){
    h_passFatJetSel->fill(event);
    TstarHypsCreated = TstarTstarHypCreator->process(event);
  }
  if(TstarHypsCreated) h_afterHypCreation->fill(event);

  if(debug) cout << "Starting to find best TstarTstar Hypothesis" << endl;
  if(TstarHypsCreated){
    bestHypFound = TstarTstarHypSelector->process(event);

    if(bestHypFound){
      h_RecoPlots_Full->fill(event);
      h_afterReco_Full->fill(event);
      h_afterReco_gen->fill(event);
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
        h_afterReco_lowchi2->fill(event);
      }
      else {
        h_RecoPlots_highchi2->fill(event);
        h_afterReco_highchi2->fill(event);
      }
    }
  }
  if(!TstarHypsCreated || !bestHypFound){
    h_notReconstructible->fill(event);
    if(pass_ttag) h_notReconstructible_ttag->fill(event);
    else h_notReconstructible_nottag->fill(event);
  }

  // ########### with gluons as AK4 jets #############

  // bool bestHypFoundAK4 = false;
  // bool TstarHypsCreatedAK4;
  //
  // if(debug){cout << "Starting to construct all TstarTstar Hypothesiseseses for AK4 mode" << endl;}
  // TstarHypsCreatedAK4 = TstarTstarHypCreatorAK4->process(event);
  // if(TstarHypsCreatedAK4) h_afterHypCreation_AK4->fill(event);
  //
  // if(debug){cout << "Starting to find best TstarTstar Hypothesis for AK4 mode" << endl;}
  // if(TstarHypsCreatedAK4){
  //
  //   bestHypFoundAK4 = TstarTstarHypSelector->process(event);
  //   if(bestHypFoundAK4){
  //     h_RecoPlots_Full_AK4->fill(event);
  //     h_afterReco_Full_AK4->fill(event);
  //     if(pass_ttag){
  //        h_RecoPlots_ttag_AK4->fill(event);
  //        h_afterReco_ttag_AK4->fill(event);
  //     }
  //     else{
  //        h_RecoPlots_nottag_AK4->fill(event);
  //        h_afterReco_nottag_AK4->fill(event);
  //     }
  //     if(event.get(h_tstartstar_hyp).chi2()<50){
  //       h_RecoPlots_lowchi2_AK4->fill(event);
  //       h_afterReco_lowchi2_AK4->fill(event);
  //     }
  //     else {
  //       h_RecoPlots_highchi2_AK4->fill(event);
  //       h_afterReco_highchi2_AK4->fill(event);
  //     }
  //   }
  // }
  // if(!TstarHypsCreatedAK4 || !bestHypFoundAK4){
  //   h_notReconstructible_AK4->fill(event);
  //   if(pass_ttag) h_notReconstructible_ttag_AK4->fill(event);
  //   else h_notReconstructible_nottag_AK4->fill(event);
  // }

  // ########### GEN matching #############

  if(is_MC && TstarHypsCreated){
    if(debug) cout << "Doing GEN matching check" << endl;
    // For full hypotheses
    {
      ReconstructionTstarHypothesis hyp_tmp = event.get(h_tstartstar_hyp);
      if(genmatcher->process(event)){
      	if(debug) cout << "GEN MATCHED!" << endl;
      	h_RecoPlots_GEN->fill(event);
      	h_afterGEN->fill(event);
      	event.set(h_tstartstar_hyp, hyp_tmp);
      }
    }

    // using only ttbar information
    {
      ReconstructionTstarHypothesis hyp_tmp = event.get(h_tstartstar_hyp);
      if(genmatcher_onlyttbar->process(event)){
      	if(debug) cout << "GEN MATCHED!" << endl;
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

    // ########################################
    // ########### DNN Preparation ############
    // ########################################

    if(debug) cout << "Start DNN stuff" << endl;

    // ###### ST reweighting ######
    double st_jets = 0;
    double ST_weight = 0;
    double ST_weight_2 = 0;
    h_top_gluon_checks->fill(event);
    for(const auto & jet : *event.topjets) st_jets += jet.pt();
    if(is_TTbar){
      if(st_jets < 3000){
        ST_weight = 1/(ST_ratio->GetBinContent(ST_ratio->GetXaxis()->FindBin(st_jets)));
        ST_weight_2 = 1/(500*ST_bkg->GetBinContent(ST_bkg->GetXaxis()->FindBin(st_jets)));
        //ST_weight_2 *= 1/(ST_bkg_2->GetBinContent(ST_bkg_2->GetXaxis()->FindBin(st_jets)));
      }

      event.set(h_ST_weight, ST_weight*event.get(h_evt_weight));
      event.set(h_ST_weight_2, ST_weight_2*event.get(h_evt_weight));
      event.weight *= ST_weight;
      h_ST_reweighted->fill(event);
      h_top_gluon_checks_reweighted->fill(event);
      event.weight = event.get(h_evt_weight);
      event.weight = ST_weight_2*event.get(h_evt_weight);
      h_ST_reweighted_2->fill(event);
      h_top_gluon_checks_reweighted_2->fill(event);
      event.weight = event.get(h_evt_weight);
    }
    else if(is_Signal){
      double ST_weight_2 = 0;
      if(st_jets < 3000){
        ST_weight_2 = 1/(ST_sig->GetBinContent(ST_sig->GetXaxis()->FindBin(st_jets)));
        //ST_weight_2 *= 1/(ST_sig_2->GetBinContent(ST_sig_2->GetXaxis()->FindBin(st_jets)));
      }
      event.set(h_ST_weight, event.get(h_evt_weight));
      event.set(h_ST_weight_2, ST_weight_2*event.get(h_evt_weight));
      h_ST_reweighted->fill(event);
      h_top_gluon_checks_reweighted->fill(event);
      event.weight = ST_weight_2*event.get(h_evt_weight);
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
      event.weight = event.get(h_evt_weight);
    }

    // Filling output for DNN
    if(outputDNNvalues){
      event.set(h_do_masspoint, do_masspoint);
      double st_jets = 0.;
      for(const auto & jet : *event.topjets) st_jets += jet.pt();
      event.set(h_ST, st_jets);
      if(debug) cout << "Write inputs" << endl;
      DNN_InputWriter->process(event);
      if(debug) cout << "plot inputs" << endl;
      h_DNN_Inputs->fill(event);
      if(debug) cout << "plot inputs for reweighted" << endl;
      if(is_TTbar) event.weight *= ST_weight;
      h_DNN_Inputs_reweighted->fill(event);
      event.weight = event.get(h_evt_weight);
      if(is_TTbar || is_Signal) event.weight = ST_weight_2*event.get(h_evt_weight);
      h_DNN_Inputs_reweighted_2->fill(event);
      event.weight = event.get(h_evt_weight);
    }

    // safety
    event.weight = event.get(h_evt_weight);

    if(debug){cout << "Done ##################################" << endl;}
    return true;

}

// as we want to run the ExampleCycleNew directly with AnalysisModuleRunner,
// make sure the TstarTstarMCStudyModule is found by class name. This is ensured by this macro:
UHH2_REGISTER_ANALYSIS_MODULE(TstarTstarMCStudyModule)

}
