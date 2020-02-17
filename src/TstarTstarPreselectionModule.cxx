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
#include "UHH2/common/include/MCWeight.h"

#include "UHH2/TstarTstar/include/TstarTstarSelections.h"
#include "UHH2/TstarTstar/include/TstarTstarHists.h"
#include "UHH2/TstarTstar/include/TstarTstarGenHists.h"
#include "UHH2/TstarTstar/include/TstarTstarGenRecoMatchedHists.h"


using namespace std;
using namespace uhh2;

namespace uhh2 {

/** \brief Module for the T*T*->ttbar gg/gamma preselection
 *  
 * Corrects all objects via CommonModules and applies some loose cuts
 *
 */
class TstarTstarPreselectionModule: public AnalysisModule {
public:
    
  explicit TstarTstarPreselectionModule(Context & ctx);
  virtual bool process(Event & event) override;

private:
    
  // Apply common modules: JetPFid, JEC, JER, MET corrections, etc
  std::unique_ptr<CommonModules> common;
  
  // Declare the Selections to use. Use unique_ptr to ensure automatic call of delete in the destructor, to avoid memory leaks.
  unique_ptr<Selection> twodcut_sel;   
  unique_ptr<TTbarSemiLepMatchableSelection> TTbarSemiLepMatchable_selection;
  unique_ptr<Selection> met_sel, st_sel; 
  unique_ptr<Selection> topjet_selection;
  
  // ##### Histograms #####
  // Store the Hists collection as member variables. Again, use unique_ptr to avoid memory leaks.
  std::unique_ptr<Hists> h_nocuts_gen, h_allcuts_gen;
  std::unique_ptr<Hists> h_nocuts, h_common, h_lepsel, h_njetsel, h_nphosel, h_2dcut, h_metcut, h_STcut, h_allcuts;
  std::unique_ptr<Hists> h_semilepttbarmatch, h_nosemilepttbarmatch, h_semilepttbarmatch_genreco,h_semilepttbarmatch_genreco_mu,h_semilepttbarmatch_genreco_ele,h_semilepttbarmatch_gen;

  std::unique_ptr<uhh2::AnalysisModule> ttgenprod;

  bool debug = false;

  // bools for channel and stuff. will be read in later
  bool is_tgtg, is_tgtgamma, isTrigger;

  bool is_MC;

};


TstarTstarPreselectionModule::TstarTstarPreselectionModule(Context & ctx){
  
  
  if(debug) {
    cout << "Hello World from TstarTstarPreselectionModule!" << endl;
    
    
    // If running in SFrame, the keys "dataset_version", "dataset_type", "dataset_lumi",
    // and "target_lumi" are set to the according values in the xml file. For CMSSW, these are
    // not set automatically, but can be set in the python config file.
    for(auto & kv : ctx.get_all()){
      cout << " " << kv.first << " = " << kv.second << endl;
    }
    
  }
  
  // 0. Reading in which channel
  is_MC = ctx.get("dataset_type") == "MC";

  is_tgtg = false; is_tgtgamma = false;
  if(ctx.get("channel") == "tgtg") is_tgtg = true;
  if(ctx.get("channel") == "tgtgamma") is_tgtgamma = true;
  
  // 1. setup other modules. CommonModules
  common.reset(new CommonModules());
  common->disable_mcpileupreweight(); //FixME: PU re-weighting crushes TODO
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

  if(is_tgtgamma){
    // photon
    PhotonId phoID; 
    double photon_pt(20.);
    phoID = PhotonTagID(Photon::mvaPhoID_Fall17_iso_V2_wp90);
    common->set_photon_id(AndId<Photon>(PtEtaCut(photon_pt, 5.2), phoID));
  }

  //Muon
  MuonId muID;
  double muon_pt(20.);
  muID = MuonID(Muon::Highpt);
  common->set_muon_id(AndId<Muon>(PtEtaCut(muon_pt, 2.4), muID));

  // init common.
  common->init(ctx);


  if(is_MC){
    // Prepare GEN
    ttgenprod.reset(new TTbarGenProducer(ctx, "ttbargen", false));
  }


  // 2. set up selections
  // 2D Cut Lepton-Jets
  twodcut_sel.reset(new TwoDCut(0.4, 25.0));  // The same as in Z'->ttbar semileptonic

  // MET selection
  met_sel.reset(new METCut  (50.,1e6));
  
  //ST selection
  st_sel.reset(new STCut  (500.,1e6));

  //Ak8jet selection  TODO maximum cut should be 2 as i need two for photon case later. atm three for code reasons, does not matter for tgtg
  topjet_selection.reset(new NTopJetSelection(3, -1, TopJetId(PtEtaCut(100, 2.1))));


  // 3. Set up Hists classes:
  
  h_nocuts_gen.reset(new TstarTstarGenHists(ctx, "NoCuts_GEN"));
  h_allcuts_gen.reset(new TstarTstarGenHists(ctx, "AfterAllcuts_GEN"));

  h_nocuts.reset(new TstarTstarHists(ctx, "NoCuts"));
  h_common.reset(new TstarTstarHists(ctx, "AfterCommon"));
  h_njetsel.reset(new TstarTstarHists(ctx, "AfterNjets"));
  h_lepsel.reset(new TstarTstarHists(ctx, "AfterLepSel"));
  h_nphosel.reset(new TstarTstarHists(ctx, "AfterNpho"));
  h_2dcut.reset(new TstarTstarHists(ctx, "After2D"));
  h_metcut.reset(new TstarTstarHists(ctx, "AfterMET"));
  h_STcut.reset(new TstarTstarHists(ctx, "AfterST"));
  h_allcuts.reset(new TstarTstarHists(ctx, "AfterAllCuts"));



  if(is_MC){
    // For GEN matching
    h_semilepttbarmatch.reset(new TstarTstarHists(ctx, "SemiLepTTBarMatch"));
    h_nosemilepttbarmatch.reset(new TstarTstarHists(ctx, "NotSemiLepTTBarMatch"));
    h_semilepttbarmatch_gen.reset(new TstarTstarGenHists(ctx, "SemiLepTTBarMatchGEN"));
    h_semilepttbarmatch_genreco.reset(new TstarTstarGenRecoMatchedHists(ctx, "SemiLepTTBarMatchGENRECO"));
    h_semilepttbarmatch_genreco_mu.reset(new TstarTstarGenRecoMatchedHists(ctx, "SemiLepTTBarMatchGENRECO_mu"));
    h_semilepttbarmatch_genreco_ele.reset(new TstarTstarGenRecoMatchedHists(ctx, "SemiLepTTBarMatchGENRECO_ele"));
    
    TTbarSemiLepMatchable_selection.reset(new TTbarSemiLepMatchableSelection());
  }

}


bool TstarTstarPreselectionModule::process(Event & event) {
   
  if(debug)   
    cout << "TstarTstarPreselectionModule: Starting to process event (runid, eventid) = (" << event.run << ", " << event.event << "); weight = " << event.weight << endl;
  if(debug)
    cout<<"N muons = "<<event.muons->size()<<", N electrons = "<<event.electrons->size()<<", N photons = "<<event.photons->size()<<endl;
  for(auto& muo : *event.muons){
    if(debug) cout<<"BEFORE Muon (pt,eta): "<<muo.pt()<<", "<<muo.eta()<<endl;
  }

  if(is_MC){
    //Fill ttgen object for correct matching check, etc
    ttgenprod->process(event);
  }

  h_nocuts->fill(event);
  h_nocuts_gen->fill(event);
  if(debug) cout<<"Filled hists without any cuts, and GEN with no cuts"<<endl;

  if(!(common->process(event))) return false;
  h_common->fill(event);
  if(debug) cout<<"Filled hists after common modules"<<endl;

  //---- Loose selection

  // Require at least 1 jet
  const bool pass_njet = (event.jets->size()>0);
  if(!pass_njet) return false;
  h_njetsel->fill(event);
  if(debug) cout << "Filled hists after njetsel" << endl;

  // Require exactly one muon or one electron
  const bool pass_lep1 = (((event.muons->size() == 1) || (event.electrons->size() == 1)) && (event.electrons->size()+event.muons->size()) == 1);
  if(!pass_lep1) return false;
  h_lepsel->fill(event);
  if(debug) cout << "Filled hists after lepsel" << endl;

   if(is_tgtgamma){
    // Require at least 1 photon
    const bool pass_npho = (event.photons->size()>0);
    if(!pass_npho) return false;
    h_nphosel->fill(event);
    if(debug) cout << "Filled hists after nphosel" << endl;
  }

  // Lepton-2Dcut variables
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
  if(!pass_twodcut) return false;
  h_2dcut->fill(event);
  if(debug) cout<<"Filled hists after 2D"<<endl;

  // MET Selection
  bool pass_MET =  met_sel->passes(event);
  if(!pass_MET) return false;
  h_metcut->fill(event);
  if(debug) cout<<"Filled hists after MET"<<endl;

  // ST Selection
  bool pass_ST =  st_sel->passes(event);
  if(!pass_ST) return false;
  h_STcut->fill(event);
  if(debug) cout<<"Filled hists after ST"<<endl;
  
  // AK8 Selection
  bool pass_ak8 = topjet_selection->passes(event);
  if(!pass_ak8) return false;
  
  if(debug) cout<<"passed all cuts"<<endl;
  h_allcuts->fill(event);
  if(debug) cout<<"error here?"<<endl;
  h_allcuts_gen->fill(event);

  if(is_MC){
    // ##### Matching to GEN ######
    if(debug) cout<<"Matching to GEN..."<<endl;
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
  }

  if(debug) cout << "########### Done with preselection! ###########" << endl << endl;

  return true;
}

// as we want to run the ExampleCycleNew directly with AnalysisModuleRunner,
// make sure the TstarTstarPreselectionModule is found by class name. This is ensured by this macro:
UHH2_REGISTER_ANALYSIS_MODULE(TstarTstarPreselectionModule)

}
