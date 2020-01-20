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

  unique_ptr<TTbarSemiLepMatchableSelection> TTbarSemiLepMatchable_selection;
  unique_ptr<Selection> toptagevt_sel;
    
  // ##### Histograms #####
  // Store the Hists collection as member variables. Again, use unique_ptr to avoid memory leaks.

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
  bool doGen = false;
  

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
  
  // 0. Reading in which channel
  is_tgtg = false; is_tgtgamma = false;
  if(ctx.get("channel") == "tgtg") is_tgtg = true;
  if(ctx.get("channel") == "tgtgamma") is_tgtgamma = true;
  
  // Prepare GEN
  if(doGen){
    ttgenprod.reset(new TTbarGenProducer(ctx, "ttbargen", false));
    h_ttbargen = ctx.get_handle<TTbarGen>("ttbargen");
    TTbarSemiLepMatchable_selection.reset(new TTbarSemiLepMatchableSelection()); // TODO what does this doe???
  }
  
  //TopTag
  const TopJetId topjetID = AndId<TopJet>(TopTagMassWindow(), TopTagSubbtag(DeepCSVBTag::WP_LOOSE),  Tau32(0.65));
  const float minDR_topjet_jet(1.2);
  toptagevt_sel.reset(new TopTagEventSelection(topjetID, minDR_topjet_jet));
  h_flag_toptagevent = ctx.declare_event_output<int>("flag_toptagevent");
 

  // 3. Set up Hists classes:

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

  if(doGen){
    //Fill ttgen object for correct matching check, etc
    ttgenprod->process(event);
  }

  /* TOPTAG-EVENT boolean */
  const bool pass_ttagevt = toptagevt_sel->passes(event);
  /* add flag_toptagevent to output ntuple */
  if(debug){cout << "TTag:" << pass_ttagevt << endl;}
  event.set(h_flag_toptagevent, int(pass_ttagevt));


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
	if(doGen){h_merged->fill(event);}
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
  
  if(doGen){
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
