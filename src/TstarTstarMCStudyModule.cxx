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
#include "UHH2/common/include/MCWeight.h"
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

  std::unique_ptr<Hists> h_beforeReco, h_afterPrimlep, h_afterHypCreation, h_afterReco_Full, h_afterReco_ttag, h_afterReco_nottag, h_afterMtcut, h_afterGEN;
  std::unique_ptr<TstarTstarRecoTstarHists> h_RecoPlots_Full, h_RecoPlots_ttag, h_RecoPlots_nottag, h_RecoPlots_Mtcut; 
  std::unique_ptr<TstarTstarRecoTstarHists> h_RecoPlots_GEN;
  
  std::unique_ptr<Hists> h_GEN_Hists;

  // Bools for Debugging/Options
  bool debug = false;

  bool doTopTagged = true;
  bool doNotTopTagged = false;
  bool doGen = true;

  // save pt, eta, phi for ttagged jet, leptop, MET, leptopjet, gluon candidates
  bool forDNN = true;
  uhh2::Event::Handle<double> h_DNN_ttaggedjet_pt;
  uhh2::Event::Handle<double> h_DNN_ttaggedjet_eta;
  uhh2::Event::Handle<double> h_DNN_ttaggedjet_phi;
  uhh2::Event::Handle<double> h_DNN_lepton_pt;
  uhh2::Event::Handle<double> h_DNN_lepton_eta;
  uhh2::Event::Handle<double> h_DNN_lepton_phi;
  uhh2::Event::Handle<double> h_DNN_MET_pt;
  uhh2::Event::Handle<double> h_DNN_MET_eta;
  uhh2::Event::Handle<double> h_DNN_MET_phi;
  uhh2::Event::Handle<double> h_DNN_leptopjet_pt;
  uhh2::Event::Handle<double> h_DNN_leptopjet_eta;
  uhh2::Event::Handle<double> h_DNN_leptopjet_phi;
  uhh2::Event::Handle<double> h_DNN_gluon1_pt;
  uhh2::Event::Handle<double> h_DNN_gluon1_eta;
  uhh2::Event::Handle<double> h_DNN_gluon1_phi;
  uhh2::Event::Handle<double> h_DNN_gluon2_pt;
  uhh2::Event::Handle<double> h_DNN_gluon2_eta;
  uhh2::Event::Handle<double> h_DNN_gluon2_phi;

  uhh2::Event::Handle<FlavorParticle> h_primlep;

  // Modules
  std::unique_ptr<uhh2::AnalysisModule> ttgenprod;
  std::unique_ptr<uhh2::AnalysisModule> reco_primlep;
  std::unique_ptr<TstarTstarGenMatcher> genmatcher;
  std::unique_ptr<CorrectMatchDiscriminator> ttbar_CorrectMatchDiscriminator;
  std::unique_ptr<TstarTstar_tgtg_TopTag_Reconstruction> TstarTstarHypCreator;
  std::unique_ptr<TstarTstar_Discrimination> TstarTstarHypSelector;

  std::unique_ptr<uhh2::AnalysisModule> MCWeight;

  // Handles
  uhh2::Event::Handle<TTbarGen> h_ttbargen;
  uhh2::Event::Handle<int> h_flag_toptagevent;
  uhh2::Event::Handle<std::vector<ReconstructionTstarHypothesis>> h_tstartstar_hyp_vector;
  uhh2::Event::Handle<ReconstructionTstarHypothesis> h_tstartstar_hyp;

  // bools for channel and stuff. will be read in later
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

  // Dont save everything!
  //ctx.undeclare_all_event_output();
  
  // 0. Reading in which channel
  is_tgtg = false; is_tgtgamma = false;
  if(ctx.get("channel") == "tgtg") is_tgtg = true;
  if(ctx.get("channel") == "tgtgamma") is_tgtgamma = true;
  
  // Prepare GEN
  if(doGen){
    ttgenprod.reset(new TTbarGenProducer(ctx, "ttbargen", false));
    h_ttbargen = ctx.get_handle<TTbarGen>("ttbargen");
    TTbarSemiLepMatchable_selection.reset(new TTbarSemiLepMatchableSelection());
    genmatcher.reset(new TstarTstarGenMatcher(ctx));
  }
  
  //TopTag
  const TopJetId topjetID = AndId<TopJet>(TopTagMassWindow(), TopTagSubbtag(DeepCSVBTag::WP_LOOSE),  Tau32(0.65));
  const float minDR_topjet_jet(1.2);
  toptagevt_sel.reset(new TopTagEventSelection(topjetID, minDR_topjet_jet));
  h_flag_toptagevent = ctx.get_handle<int>("flag_toptagevent");
 

  // 3. Set up Hists classes:

  h_beforeReco.reset(new TstarTstarHists(ctx, "beforeReco"));
  h_afterPrimlep.reset(new TstarTstarHists(ctx, "afterPrimlep"));
  h_afterHypCreation.reset(new TstarTstarHists(ctx, "afterHypCreation"));
  h_afterReco_Full.reset(new TstarTstarHists(ctx, "AfterReco_Full"));
  h_afterReco_ttag.reset(new TstarTstarHists(ctx, "AfterReco_ttag"));
  h_afterReco_nottag.reset(new TstarTstarHists(ctx, "AfterReco_nottag"));
  h_afterMtcut.reset(new TstarTstarHists(ctx, "AfterMtcut"));
  h_afterGEN.reset(new TstarTstarHists(ctx, "AfterGEN"));

  h_RecoPlots_Full.reset(new TstarTstarRecoTstarHists(ctx, "RecoPlots_Full"));
  h_RecoPlots_ttag.reset(new TstarTstarRecoTstarHists(ctx, "RecoPlots_ttag"));
  h_RecoPlots_nottag.reset(new TstarTstarRecoTstarHists(ctx, "RecoPlots_nottag"));
  h_RecoPlots_Mtcut.reset(new TstarTstarRecoTstarHists(ctx, "RecoPlots_Mtcut"));
 
  h_RecoPlots_GEN.reset(new TstarTstarRecoTstarHists(ctx, "RecoPlots_GEN"));

  h_GEN_Hists.reset(new TstarTstarGenHists(ctx, "GEN_Hists_AfterReco"));
  

  //4. Set up ttbar reconstruction
  reco_primlep.reset(new PrimaryLepton(ctx));

  h_tstartstar_hyp_vector = ctx.get_handle<std::vector<ReconstructionTstarHypothesis>>("TstarTstar_Hyp_Vector");
  h_tstartstar_hyp = ctx.get_handle<ReconstructionTstarHypothesis>("TstarTstar_Hyp");

  TstarTstarHypCreator.reset(new TstarTstar_tgtg_TopTag_Reconstruction(ctx, NeutrinoReconstruction, topjetID, minDR_topjet_jet));
  TstarTstarHypSelector.reset(new TstarTstar_Discrimination(ctx));

  MCWeight.reset(new MCLumiWeight(ctx));

  // 5. Handles for DNN
  if(forDNN){
    h_DNN_ttaggedjet_pt = ctx.declare_event_output<double>("DNN_ttaggedjet_pt");
    h_DNN_ttaggedjet_eta = ctx.declare_event_output<double>("DNN_ttaggedjet_eta");
    h_DNN_ttaggedjet_phi = ctx.declare_event_output<double>("DNN_ttaggedjet_phi");
    h_DNN_lepton_pt = ctx.declare_event_output<double>("DNN_lepton_pt");
    h_DNN_lepton_eta = ctx.declare_event_output<double>("DNN_lepton_eta");
    h_DNN_lepton_phi = ctx.declare_event_output<double>("DNN_lepton_phi");
    h_DNN_MET_pt = ctx.declare_event_output<double>("DNN_MET_pt");
    h_DNN_MET_eta = ctx.declare_event_output<double>("DNN_MET_eta");
    h_DNN_MET_phi = ctx.declare_event_output<double>("DNN_MET_phi");
    h_DNN_leptopjet_pt = ctx.declare_event_output<double>("DNN_leptopjet_pt");
    h_DNN_leptopjet_eta = ctx.declare_event_output<double>("DNN_leptopjet_eta");
    h_DNN_leptopjet_phi = ctx.declare_event_output<double>("DNN_leptopjet_phi");
    h_DNN_gluon1_pt = ctx.declare_event_output<double>("DNN_gluon1_pt");
    h_DNN_gluon1_eta = ctx.declare_event_output<double>("DNN_gluon1_eta");
    h_DNN_gluon1_phi = ctx.declare_event_output<double>("DNN_gluon1_phi");
    h_DNN_gluon2_pt = ctx.declare_event_output<double>("DNN_gluon2_pt");
    h_DNN_gluon2_eta = ctx.declare_event_output<double>("DNN_gluon2_eta");
    h_DNN_gluon2_phi = ctx.declare_event_output<double>("DNN_gluon2_phi");

    reco_primlep.reset(new PrimaryLepton(ctx));
    h_primlep = ctx.get_handle<FlavorParticle>("PrimaryLepton");
  }

}


bool TstarTstarMCStudyModule::process(Event & event) {
   
  if(debug){cout << endl << "TstarTstarMCStudyModule: Starting to process event (runid, eventid) = (" << event.run << ", " << event.event << "); weight = " << event.weight << endl;}

  MCWeight->process(event);
  ttgenprod->process(event);

  h_beforeReco->fill(event);

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

  if(is_tgtg){
    if(debug){cout << "Starting to construct all TstarTstar Hypothesiseseses" << endl;}
    if(event.get(h_flag_toptagevent) && doTopTagged){ // if we have a TopTag!
      TstarHypsCreated = TstarTstarHypCreator->process(event);
      if(TstarHypsCreated){
	h_afterHypCreation->fill(event);
	bestHypFound = TstarTstarHypSelector->process(event);
	if(bestHypFound){
	  h_RecoPlots_Full->fill(event);
	  h_afterReco_Full->fill(event);
	  h_GEN_Hists->fill(event);

	  // Fill a new hist with some cut(s):
	  if(inv_mass(event.get(h_tstartstar_hyp).ttbar_hyp().tophad_v4()) > 165){
	    h_RecoPlots_Mtcut->fill(event);
	    h_afterMtcut->fill(event);
	  }
	}
      }
    }
  }
  else if(is_tgtgamma){
    // TODO
  }

  // Filling output for DNN
  // TODO put this in extra file that you pass an hypothesis
  if(forDNN && bestHypFound){
    if(debug) cout << "Start filling of DNN stuff" << endl;
    ReconstructionTstarHypothesis best_hyp = event.get(h_tstartstar_hyp);
    ReconstructionHypothesis best_hyp_ttbar = best_hyp.ttbar_hyp();

    event.set(h_DNN_ttaggedjet_pt, best_hyp_ttbar.tophad_v4().pt());
    event.set(h_DNN_ttaggedjet_eta, best_hyp_ttbar.tophad_v4().eta());
    event.set(h_DNN_ttaggedjet_phi, best_hyp_ttbar.tophad_v4().phi());
    if(debug) cout << "Done with ttagjet." << endl;

    event.set(h_DNN_lepton_pt, best_hyp_ttbar.lepton().pt());
    event.set(h_DNN_lepton_eta, best_hyp_ttbar.lepton().eta());
    event.set(h_DNN_lepton_phi, best_hyp_ttbar.lepton().phi());
    if(debug) cout << "Done with lepton." << endl;

    event.set(h_DNN_MET_pt, best_hyp_ttbar.neutrino_v4().pt());
    event.set(h_DNN_MET_eta, best_hyp_ttbar.neutrino_v4().eta());
    event.set(h_DNN_MET_phi, best_hyp_ttbar.neutrino_v4().phi());
    if(debug) cout << "Done with neutrino." << endl;

    event.set(h_DNN_leptopjet_pt, best_hyp_ttbar.blep_v4().pt());
    event.set(h_DNN_leptopjet_eta, best_hyp_ttbar.blep_v4().eta());
    event.set(h_DNN_leptopjet_phi, best_hyp_ttbar.blep_v4().phi());
    if(debug) cout << "Done with leptopjet." << endl;

    event.set(h_DNN_gluon1_pt, best_hyp.gluon1_v4().pt());
    event.set(h_DNN_gluon1_eta, best_hyp.gluon1_v4().eta());
    event.set(h_DNN_gluon1_phi, best_hyp.gluon1_v4().phi());
    if(debug) cout << "Done with gluon1." << endl;
    
    event.set(h_DNN_gluon2_pt, best_hyp.gluon2_v4().pt());
    event.set(h_DNN_gluon2_eta, best_hyp.gluon2_v4().eta());
    event.set(h_DNN_gluon2_phi, best_hyp.gluon2_v4().phi());
    if(debug) cout << "Done with gluon2." << endl;
  }
  
  
  if(doGen && TstarHypsCreated){
    if(debug){ cout << "Doing GEN matching check" << endl;}
    // ##### GEN Matching
    ReconstructionTstarHypothesis hyp_tmp = event.get(h_tstartstar_hyp);
    if(genmatcher->process(event)){
	if(debug)cout << "GEN MATCHED!" << endl;
	h_RecoPlots_GEN->fill(event);
	h_afterGEN->fill(event);
	event.set(h_tstartstar_hyp, hyp_tmp);
    }
  }
  
  if(debug){cout << "Done ##################################" << endl;}
  return bestHypFound;
}

// as we want to run the ExampleCycleNew directly with AnalysisModuleRunner,
// make sure the TstarTstarMCStudyModule is found by class name. This is ensured by this macro:
UHH2_REGISTER_ANALYSIS_MODULE(TstarTstarMCStudyModule)

}
