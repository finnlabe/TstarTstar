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
  unique_ptr<Selection> met_sel;
  
  // ##### Histograms #####
  // Store the Hists collection as member variables. Again, use unique_ptr to avoid memory leaks.
  std::unique_ptr<Hists> h_nocuts,     h_common,     h_lepsel,     h_AK8jetsel,     h_AK4jetsel,     h_METsel;
  std::unique_ptr<Hists> h_nocuts_gen, h_common_gen, h_lepsel_gen, h_AK8jetsel_gen, h_AK4jetsel_gen, h_METsel_gen;

  std::unique_ptr<uhh2::AnalysisModule> ttgenprod;
  uhh2::Event::Handle<TTbarGen> h_ttbargen;
  
  bool debug = false;

  // bools for channel and stuff. will be read in later
  bool is_MC;
  bool is_tgtg, is_tgtgamma;

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
  
  // 0. Reading in whether MC and if so, which channel
  is_MC = ctx.get("dataset_type") == "MC";

  if(is_MC){
    is_tgtg = false; is_tgtgamma = false;
    if(ctx.get("channel") == "tgtg") is_tgtg = true;
    if(ctx.get("channel") == "tgtgamma") is_tgtgamma = true;
  }
  
  // 1. setup modules. CommonModules
  common.reset(new CommonModules());
  common->disable_mcpileupreweight(); //FixME: PU re-weighting crushes TODO
  common->switch_metcorrection();

  // Electron
  ElectronId eleID; 
  double electron_pt(20.);
  eleID = ElectronID_Summer16_tight_noIso;
  common->set_electron_id(AndId<Electron>(PtEtaSCCut(electron_pt, 2.5), eleID));

  //Muon
  MuonId muID;
  double muon_pt(20.);
  muID = MuonID(Muon::Highpt);
  common->set_muon_id(AndId<Muon>(PtEtaCut(muon_pt, 2.4), muID));

  // Jets
  common->switch_jetlepcleaner();
  common->switch_jetPtSorter();
  common->set_jet_id(JetPFID(JetPFID::WP_TIGHT_PUPPI));

  // init common.
  common->init(ctx);

  if(is_MC){
    // Prepare GEN
    ttgenprod.reset(new TTbarGenProducer(ctx, "ttbargen", false));
  }

  // 2. set up selections

  // MET selection
  met_sel.reset(new METCut  (30.,1e6));

  // 3. Set up Hists classes:
  h_nocuts.reset(new TstarTstarHists(ctx, "NoCuts"));
  h_common.reset(new TstarTstarHists(ctx, "AfterCommon"));
  h_lepsel.reset(new TstarTstarHists(ctx, "AfterLepSel"));
  h_AK8jetsel.reset(new TstarTstarHists(ctx, "AfterAK8jets"));
  h_AK4jetsel.reset(new TstarTstarHists(ctx, "AfterAK4jets"));
  h_METsel.reset(new TstarTstarHists(ctx, "AfterMET"));

  h_nocuts_gen.reset(new TstarTstarHists(ctx, "NoCuts_gen"));
  h_common_gen.reset(new TstarTstarHists(ctx, "AfterCommon_gen"));
  h_lepsel_gen.reset(new TstarTstarHists(ctx, "AfterLepSel_gen"));
  h_AK8jetsel_gen.reset(new TstarTstarHists(ctx, "AfterAK8jets_gen"));
  h_AK4jetsel_gen.reset(new TstarTstarHists(ctx, "AfterAK4jets_gen"));
  h_METsel_gen.reset(new TstarTstarHists(ctx, "AfterMET_gen"));

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
  h_common_gen->fill(event);
  if(debug) cout<<"Filled hists after cleaning"<<endl;

  
  //---- Preselection

  // Require exactly one muon or one electron
  const bool pass_lep1 = (((event.muons->size() == 1) || (event.electrons->size() == 1)) && (event.electrons->size()+event.muons->size()) == 1);
  if(!pass_lep1) return false;
  h_lepsel->fill(event);
  h_lepsel_gen->fill(event);
  if(debug) cout << "Filled hists after lepsel" << endl;

  // AK8 jet selection
  bool pass_ak8_njet = (event.topjets->size()>1);
  if(!pass_ak8_njet) return false;
  h_AK8jetsel->fill(event);
  h_AK8jetsel_gen->fill(event);
  if(debug) cout << "Filled hists after AK8jetsel" << endl;

  // AK4 jet selection
  const bool pass_ak4_njet = (event.jets->size()>0);
  if(!pass_ak4_njet) return false;
  h_AK4jetsel->fill(event);
  h_AK4jetsel_gen->fill(event);
  if(debug) cout << "Filled hists after AK4jetsel" << endl;

  // MET Selection
  bool pass_MET =  met_sel->passes(event);
  if(!pass_MET) return false;
  h_METsel->fill(event);
  h_METsel_gen->fill(event);
  if(debug) cout<<"Filled hists after MET"<<endl;
  if(debug) cout<<"passed all cuts"<<endl;

  if(debug) cout << "########### Done with preselection! ###########" << endl << endl;

  return true;
}

// as we want to run the ExampleCycleNew directly with AnalysisModuleRunner,
// make sure the TstarTstarPreselectionModule is found by class name. This is ensured by this macro:
UHH2_REGISTER_ANALYSIS_MODULE(TstarTstarPreselectionModule)

}
