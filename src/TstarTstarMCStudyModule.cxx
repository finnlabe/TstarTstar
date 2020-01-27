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

  std::unique_ptr<TstarTstarRecoTstarHists> h_RecoPlots_Full, h_RecoPlots_ttag, h_RecoPlots_nottag; 
  std::unique_ptr<TstarTstarRecoTstarHists> h_RecoPlots_GEN;
  
  // Bools for Debugging/Options
  bool debug = false;

  bool check_ttbar_reco = true;
  bool doTopTagged = true;
  bool doNotTopTagged = false;
  bool doGen = false;
  

  // Modules
  std::unique_ptr<uhh2::AnalysisModule> ttgenprod;
  std::unique_ptr<uhh2::AnalysisModule> reco_primlep;
  std::unique_ptr<TstarTstarGenMatcher> genmatcher;
  std::unique_ptr<CorrectMatchDiscriminator> ttbar_CorrectMatchDiscriminator;
  std::unique_ptr<TstarTstar_tgtg_TopTag_Reconstruction> TstarTstarHypCreator;
  std::unique_ptr<TstarTstar_Discrimination> TstarTstarHypSelector;

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
  ctx.undeclare_all_event_output();
  
  // 0. Reading in which channel
  is_tgtg = false; is_tgtgamma = false;
  if(ctx.get("channel") == "tgtg") is_tgtg = true;
  if(ctx.get("channel") == "tgtgamma") is_tgtgamma = true;
  
  // Prepare GEN
  if(doGen){
    ttgenprod.reset(new TTbarGenProducer(ctx, "ttbargen", false));
    h_ttbargen = ctx.get_handle<TTbarGen>("ttbargen");
    TTbarSemiLepMatchable_selection.reset(new TTbarSemiLepMatchableSelection()); // TODO what does this do???
  }
  
  //TopTag
  const TopJetId topjetID = AndId<TopJet>(TopTagMassWindow(), TopTagSubbtag(DeepCSVBTag::WP_LOOSE),  Tau32(0.65));
  const float minDR_topjet_jet(1.2);
  toptagevt_sel.reset(new TopTagEventSelection(topjetID, minDR_topjet_jet));
  h_flag_toptagevent = ctx.declare_event_output<int>("flag_toptagevent");
 

  // 3. Set up Hists classes:
  h_RecoPlots_Full.reset(new TstarTstarRecoTstarHists(ctx, "RecoPlots_Full"));
  h_RecoPlots_ttag.reset(new TstarTstarRecoTstarHists(ctx, "RecoPlots_ttag"));
  h_RecoPlots_nottag.reset(new TstarTstarRecoTstarHists(ctx, "RecoPlots_nottag"));

  h_RecoPlots_GEN.reset(new TstarTstarRecoTstarHists(ctx, "RecoPlots_GEN"));


  //4. Set up ttbar reconstruction
  reco_primlep.reset(new PrimaryLepton(ctx));

  h_tstartstar_hyp_vector = ctx.get_handle<std::vector<ReconstructionTstarHypothesis>>("TstarTstar_Hyp_Vector");
  h_tstartstar_hyp = ctx.get_handle<ReconstructionTstarHypothesis>("TstarTstar_Hyp");
  TstarTstarHypCreator.reset(new TstarTstar_tgtg_TopTag_Reconstruction(ctx, NeutrinoReconstruction, topjetID, minDR_topjet_jet));
  TstarTstarHypSelector.reset(new TstarTstar_Discrimination(ctx));
  
  // GEN stuff
  //genmatcher.reset(new TstarTstarGenMatcher(ctx));
  

}


bool TstarTstarMCStudyModule::process(Event & event) {
   
  if(debug){cout << endl << "TstarTstarMCStudyModule: Starting to process event (runid, eventid) = (" << event.run << ", " << event.event << "); weight = " << event.weight << endl;}

  // Fix empty handle problem
  // (Dummy values are filled into handles so that they are not empty)
  event.set(h_tstartstar_hyp, ReconstructionTstarHypothesis());
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
    if(debug){ cout << "Starting to construct all TstarTstar Hypothesiseseses" << endl;}
    if(event.get(h_flag_toptagevent) && doTopTagged){ // if we have a TopTag!
      if(TstarTstarHypCreator->process(event)){
	if(TstarTstarHypSelector->process(event)){
	  h_RecoPlots_Full->fill(event);
	}
      }
    }
  }
  else if(is_tgtgamma){
    // TODO
  }

  /**
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
  **/

  if(debug){cout << "Done ##################################" << endl;}
  return true;
}

// as we want to run the ExampleCycleNew directly with AnalysisModuleRunner,
// make sure the TstarTstarMCStudyModule is found by class name. This is ensured by this macro:
UHH2_REGISTER_ANALYSIS_MODULE(TstarTstarMCStudyModule)

}
