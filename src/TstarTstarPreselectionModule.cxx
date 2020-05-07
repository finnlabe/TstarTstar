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
#include <UHH2/common/include/TriggerSelection.h>

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
  unique_ptr<TopJetCleaner> AK8cleaner;
  unique_ptr<TopJetCleaner> HOTVRcleaner;
  unique_ptr<Selection> met_sel;


  // ##### Histograms #####
  // Store the Hists collection as member variables. Again, use unique_ptr to avoid memory leaks.
  std::unique_ptr<Hists> h_nocuts,     h_common,     h_lepsel,     h_fatjetsel,     h_METsel;
  std::unique_ptr<Hists> h_nocuts_gen, h_common_gen, h_lepsel_gen, h_fatjetsel_gen, h_METsel_gen;
  std::unique_ptr<Hists> h_lepsel_ele,     h_fatjetsel_ele,     h_METsel_ele;
  std::unique_ptr<Hists> h_lepsel_mu,     h_fatjetsel_mu,     h_METsel_mu;
  std::unique_ptr<uhh2::AnalysisModule> ttgenprod;
  uhh2::Event::Handle<TTbarGen> h_ttbargen;
  uhh2::Event::Handle<bool> h_is_muevt;

  bool debug = false;

  // bools for channel and stuff. will be read in later
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

  // 0. Reading in whether MC and if so, which channel
  is_MC = ctx.get("dataset_type") == "MC";

  // 1. setup modules. CommonModules
  common.reset(new CommonModules());
  common->switch_metcorrection();

  // Electron
  ElectronId eleID;
  double electron_pt(30.);
  eleID = ElectronID_Summer16_tight_noIso;
  common->set_electron_id(AndId<Electron>(PtEtaSCCut(electron_pt, 2.4), eleID));

  //Muon
  MuonId muID;
  double muon_pt(30.);
  muID = MuonID(Muon::Highpt);
  common->set_muon_id(AndId<Muon>(PtEtaCut(muon_pt, 2.4), muID));

  // Jets
  common->switch_jetlepcleaner();
  common->switch_jetPtSorter();
  double jet_pt(30.);
  common->set_jet_id(AndId<Jet>(PtEtaCut(jet_pt, 2.5), JetPFID(JetPFID::WP_TIGHT_PUPPI)));
  HOTVRcleaner.reset(new TopJetCleaner(ctx, PtEtaCut(150.0, 2.5)));

  // init common.
  common->init(ctx);

  if(is_MC){
    // Prepare GEN
    ttgenprod.reset(new TTbarGenProducer(ctx, "ttbargen", false));
  }

  // 2. set up selections

  // MET selection
  met_sel.reset(new METCut  (50.,1e9));

  // 3. Set up Hists classes:
  h_nocuts.reset(new TstarTstarHists(ctx, "NoCuts"));
  h_common.reset(new TstarTstarHists(ctx, "AfterCommon"));
  h_lepsel.reset(new TstarTstarHists(ctx, "AfterLepSel"));
  h_fatjetsel.reset(new TstarTstarHists(ctx, "AfterAK8jets"));
  h_METsel.reset(new TstarTstarHists(ctx, "AfterMET"));

  h_lepsel_ele.reset(new TstarTstarHists(ctx, "AfterLepSel_ele"));
  h_fatjetsel_ele.reset(new TstarTstarHists(ctx, "AfterAK8jets_ele"));
  h_METsel_ele.reset(new TstarTstarHists(ctx, "AfterMET_ele"));

  h_lepsel_mu.reset(new TstarTstarHists(ctx, "AfterLepSel_mu"));
  h_fatjetsel_mu.reset(new TstarTstarHists(ctx, "AfterAK8jets_mu"));
  h_METsel_mu.reset(new TstarTstarHists(ctx, "AfterMET_mu"));

  h_nocuts_gen.reset(new TstarTstarHists(ctx, "NoCuts_gen"));
  h_common_gen.reset(new TstarTstarHists(ctx, "AfterCommon_gen"));
  h_lepsel_gen.reset(new TstarTstarHists(ctx, "AfterLepSel_gen"));
  h_fatjetsel_gen.reset(new TstarTstarHists(ctx, "AfterAK8jets_gen"));
  h_METsel_gen.reset(new TstarTstarHists(ctx, "AfterMET_gen"));

  h_is_muevt = ctx.declare_event_output<bool>("is_muevt");

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

  // cleaning & common modules
  if(!(common->process(event))) return false;
  if(!(HOTVRcleaner->process(event))) return false;
  h_common->fill(event);
  h_common_gen->fill(event);
  if(debug) cout<<"Filled hists after cleaning"<<endl;


  //---- Preselection

  // Require exactly one muon or one electron
  const bool pass_lep1 = (((event.muons->size() == 1) || (event.electrons->size() == 1)) && (event.electrons->size()+event.muons->size()) == 1);
  if(!pass_lep1) return false;
  if(event.muons->size() == 1) event.set(h_is_muevt, true);
  else event.set(h_is_muevt, false);
  h_lepsel->fill(event);
  h_lepsel_gen->fill(event);
  if(event.get(h_is_muevt)) h_lepsel_mu->fill(event);
  else h_lepsel_ele->fill(event);
  if(debug) cout << "Filled hists after lepsel" << endl;

  // fat jet selection
  bool pass_fat_njet = (event.topjets->size()>0);
  if(!pass_fat_njet) return false;
  h_fatjetsel->fill(event);
  h_fatjetsel_gen->fill(event);
  if(event.get(h_is_muevt)) h_fatjetsel_mu->fill(event);
  else h_fatjetsel_ele->fill(event);
  if(debug) cout << "Filled hists after fatjetsel" << endl;

  // MET Selection
  bool pass_MET =  met_sel->passes(event);
  if(!pass_MET) return false;
  h_METsel->fill(event);
  h_METsel_gen->fill(event);
  if(event.get(h_is_muevt)) h_METsel_mu->fill(event);
  else h_METsel_ele->fill(event);
  if(debug) cout<<"Filled hists after MET"<<endl;

  if(debug) cout << "########### Done with preselection! ###########" << endl << endl;

  return true;
}

// as we want to run the ExampleCycleNew directly with AnalysisModuleRunner,
// make sure the TstarTstarPreselectionModule is found by class name. This is ensured by this macro:
UHH2_REGISTER_ANALYSIS_MODULE(TstarTstarPreselectionModule)

}
