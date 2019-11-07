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
#include "UHH2/common/include/MuonIds.h"
#include "UHH2/common/include/TriggerSelection.h"
#include "UHH2/TstarTstar/include/ModuleBASE.h"
#include "UHH2/TstarTstar/include/TstarTstarSelections.h"
#include "UHH2/TstarTstar/include/TstarTstarHists.h"
#include "UHH2/TstarTstar/include/TstarTstarGenHists.h"
#include "UHH2/TstarTstar/include/TstarTstarGenRecoMatchedHists.h"

#include "UHH2/common/include/MCWeight.h"

using namespace std;
using namespace uhh2;

//namespace uhh2 {

/** \brief Module for the full selection in T*T*->ttbar gg/gamma search 
 *
 * All objects are expected to be corrected in PreSelection stage
 *
 */
// class TstarTstarSelectionModule: public AnalysisModule {
class TstarTstarSelectionModule: public ModuleBASE {
public:
    
    explicit TstarTstarSelectionModule(Context & ctx);
    virtual bool process(Event & event) override;
    void book_histograms(uhh2::Context&, vector<string>);
    void fill_histograms(uhh2::Event&, string, bool);

private:
    

    // Declare the Selections to use. Use unique_ptr to ensure automatic call of delete in the destructor,
    // to avoid memory leaks.

    unique_ptr<Selection> TTbarSemiLepMatchable_selection;
  
    unique_ptr<Selection> triggerPFHT_sel;
    unique_ptr<Selection> triggerSingleJet450_sel;
    unique_ptr<Selection> triggerSingleLeptonEle1_sel;
    unique_ptr<Selection> triggerSingleLeptonEle2_sel;
    unique_ptr<Selection> triggerSingleLeptonEle3_sel;
    unique_ptr<Selection> triggerSingleLeptonMu1_sel;
    unique_ptr<Selection> triggerSingleLeptonMu2_sel;
    unique_ptr<Selection> triggerSingleLeptonMu3_sel;
    unique_ptr<Selection> triggerSingleLeptonMu4_sel;
    unique_ptr<Selection> triggerHT1_sel, triggerHT2_sel, triggerHT3_sel, triggerHT4_sel, triggerHT5_sel,  triggerHT6_sel;
    unique_ptr<Selection>  met_sel, st_sel; //selections defined in UHH2/TstarTstar/include/TstarTstarSelections.h
    unique_ptr<Selection> topjet_selection;

    unique_ptr<AnalysisModule> LumiWeight_module;
    // Store the Hists collection as member variables. Again, use unique_ptr to avoid memory leaks.
  //  std::unique_ptr<Hists> h_nocuts, h_common, h_lepsel, h_2dcut, h_semilepttbarmatch, h_nosemilepttbarmatch;// h_trigger;
  //    std::unique_ptr<Hists> h_semilepttbarmatch_gen;
  //    std::unique_ptr<Hists> h_semilepttbarmatch_genreco;
  bool debug = false;
  //  bool debug = true;
};

void TstarTstarSelectionModule::book_histograms(uhh2::Context& ctx, vector<string> tags){
  for(const auto & tag : tags){
    string  mytag = tag + "_RECO";
    book_HFolder(mytag, new TstarTstarHists(ctx,mytag));
    mytag = tag + "_GEN";
    book_HFolder(mytag, new TstarTstarGenHists(ctx,mytag));
    mytag = tag + "_GENRECO";
    book_HFolder(mytag, new TstarTstarGenRecoMatchedHists(ctx,mytag));
  }
}

  void TstarTstarSelectionModule::fill_histograms(uhh2::Event& event, string tag, bool pass_ttbarsemilep){
    //  for(const auto & tag : tags){
    string mytag = tag + "_RECO";
    HFolder(mytag)->fill(event);
    if(pass_ttbarsemilep){
      mytag = tag + "_GEN";
      HFolder(mytag)->fill(event);
      mytag = tag + "_GENRECO";
      HFolder(mytag)->fill(event);
    }
    //  }
}
  


TstarTstarSelectionModule::TstarTstarSelectionModule(Context & ctx){
    
 
  if(debug) {
    cout << "Hello World from TstarTstarSelectionModule!" << endl;
    
    
    // If running in SFrame, the keys "dataset_version", "dataset_type", "dataset_lumi",
    // and "target_lumi" are set to the according values in the xml file. For CMSSW, these are
    // not set automatically, but can be set in the python config file.
    for(auto & kv : ctx.get_all()){
        cout << " " << kv.first << " = " << kv.second << endl;
    }
    
   }

    // 1. setup common selection  module -> skipped here, because done at preselection
    // 1b. set up lumi rewitghting
    LumiWeight_module.reset(new MCLumiWeight(ctx));
    
    // 2. set up selections
    //Trigger selection
    //

    //HT triggers
    triggerHT1_sel.reset(new TriggerSelection("HLT_HT430to450_v*"));
    triggerHT2_sel.reset(new TriggerSelection("HLT_HT450to470_v*"));
    triggerHT3_sel.reset(new TriggerSelection("HLT_HT470to500_v*"));
    triggerHT4_sel.reset(new TriggerSelection("HLT_HT500to550_v*"));
    triggerHT5_sel.reset(new TriggerSelection("HLT_HT550to650_v*"));
    triggerHT6_sel.reset(new TriggerSelection("HLT_HT650_v*"));

    //PF HT trigger
    triggerPFHT_sel.reset(new TriggerSelection("HLT_PFHT900_v*"));

    //SingleJet trigger
    triggerSingleJet450_sel.reset(new TriggerSelection("HLT_PFJet450_v*"));

    //Lepton trigger
    triggerSingleLeptonEle1_sel.reset(new TriggerSelection("HLT_Ele115_CaloIdVT_GsfTrkIdT_v*"));
    triggerSingleLeptonEle2_sel.reset(new TriggerSelection("HLT_Ele25_eta2p1_WPTight_Gsf_v*"));
    triggerSingleLeptonEle3_sel.reset(new TriggerSelection("HLT_Ele32_eta2p1_WPTight_Gsf_v*"));


    triggerSingleLeptonMu1_sel.reset(new TriggerSelection("HLT_Mu50_v*"));
    triggerSingleLeptonMu2_sel.reset(new TriggerSelection("HLT_Mu55_v*"));
    triggerSingleLeptonMu3_sel.reset(new TriggerSelection("HLT_IsoMu24_v*"));
    triggerSingleLeptonMu4_sel.reset(new TriggerSelection("HLT_IsoTkMu24_v*"));

    //MET selection
    met_sel.reset(new METCut  (50.,1e6));

    //ST selection
    st_sel.reset(new STCut  (500.,1e6));

    //Ak8jet selection
    topjet_selection.reset(new NTopJetSelection(1, -1, TopJetId(PtEtaCut(100, 2.1))));

    //Match to semileptonic ttbar
    TTbarSemiLepMatchable_selection.reset(new TTbarSemiLepMatchableSelection());// for x-checks


    // 3. Set up Hists
    vector<string> histogram_tags = {"PreSelection","PreSelection_mu","PreSelection_ele","AK8sel","AK8sel_mu","AK8sel_ele",
				     "MET","MET_mu","MET_ele","ST","ST_mu","ST_ele","triggerSingleLeptonMu",
				     "triggerSingleLeptonEle","triggerSingleJet_mu","triggerSingleJet_ele",
				     "triggerHT_mu","triggerHT_ele","triggerPFHT_mu","triggerPFHT_ele"};
    book_histograms(ctx, histogram_tags);
   
}


bool TstarTstarSelectionModule::process(Event & event) {
  if(debug)   
    cout << "TstarTstarSelectionModule: Starting to process event (runid, eventid) = (" << event.run << ", " << event.event << "); weight = " << event.weight << endl;
  LumiWeight_module->process(event);
  const bool pass_ttbarsemilep = TTbarSemiLepMatchable_selection->passes(event);
  fill_histograms(event, "PreSelection", pass_ttbarsemilep);
  if(event.muons->size() == 1) fill_histograms(event, "PreSelection_mu", pass_ttbarsemilep);
  if(event.electrons->size() == 1) fill_histograms(event, "PreSelection_ele", pass_ttbarsemilep);
 if(debug)   
   cout << "pass_ttbarsemilep" <<endl;
 if(event.jets->size()<1) return false;//FixMe: why this is happening?


 bool pass_ak8 = topjet_selection->passes(event);
 if(!pass_ak8) return false;

 fill_histograms(event, "AK8sel", pass_ttbarsemilep);
 if(event.muons->size() == 1) fill_histograms(event, "AK8sel_mu", pass_ttbarsemilep);
 if(event.electrons->size() == 1) fill_histograms(event, "AK8sel_ele", pass_ttbarsemilep);


 bool pass_MET =  met_sel->passes(event);
 if(!pass_MET) return false;

 fill_histograms(event, "MET", pass_ttbarsemilep);
 if(event.muons->size() == 1) fill_histograms(event, "MET_mu", pass_ttbarsemilep);
 if(event.electrons->size() == 1) fill_histograms(event, "MET_ele", pass_ttbarsemilep);

 bool pass_ST =  st_sel->passes(event);
 if(!pass_ST) return false;
 fill_histograms(event, "ST", pass_ttbarsemilep);
 if(event.muons->size() == 1) fill_histograms(event, "ST_mu", pass_ttbarsemilep);
 if(event.electrons->size() == 1) fill_histograms(event, "ST_ele", pass_ttbarsemilep);



  bool pass_trigger_SingleJet = (triggerSingleJet450_sel->passes(event) && event.jets->at(0).pt()>500);
  if(pass_trigger_SingleJet){
    if(event.muons->size() == 1) fill_histograms(event, "triggerSingleJet_mu", pass_ttbarsemilep);
    if(event.electrons->size() == 1) fill_histograms(event, "triggerSingleJet_ele", pass_ttbarsemilep);
  }
 if(debug)   
   cout << "pass_trigger_SingleJet" <<endl;

  bool pass_trigger_SingleMu = (triggerSingleLeptonMu1_sel->passes(event) || triggerSingleLeptonMu2_sel->passes(event) 
				|| triggerSingleLeptonMu3_sel->passes(event) || triggerSingleLeptonMu4_sel->passes(event));

  if(pass_trigger_SingleMu){
    if((event.muons->size() == 1) && (event.muons->at(0).pt()>60))
      fill_histograms(event, "triggerSingleLeptonMu", pass_ttbarsemilep); //FixMe: each Muon trigger should have its own threshold
  }
 if(debug)   
   cout << "pass_trigger_SingleMu" <<endl;

 bool pass_trigger_SingleEle = (triggerSingleLeptonEle1_sel->passes(event) || triggerSingleLeptonEle2_sel->passes(event) || triggerSingleLeptonEle3_sel->passes(event));
 if(pass_trigger_SingleEle){
   if((event.electrons->size() == 1) && (event.electrons->at(0).pt()>120)) 
     fill_histograms(event, "triggerSingleLeptonEle", pass_ttbarsemilep); //FixMe: each Electron trigger should have its own threshold
 }
 if(debug)   
   cout << "pass_trigger_SingleEle" <<endl;

 double st_jets = 0.;
 std::vector<Jet>* jets = event.jets;
 for(const auto & jet : *jets) st_jets += jet.pt();

 bool pass_trigegr_HT = (triggerHT1_sel->passes(event) || triggerHT2_sel->passes(event) || triggerHT3_sel->passes(event) 
    || triggerHT4_sel->passes(event) || triggerHT5_sel->passes(event) || triggerHT6_sel->passes(event));
  if(pass_trigegr_HT){
    if((event.muons->size() == 1) && st_jets>650) fill_histograms(event, "triggerHT_mu", pass_ttbarsemilep); //FixME: each HT trigger should have different st_jets threshold
    if((event.electrons->size() == 1) && st_jets>650) fill_histograms(event, "triggerHT_ele", pass_ttbarsemilep);
  }
 if(debug)   
   cout << "pass_trigegr_HT" <<endl;

  bool pass_trigegr_PFHT = triggerPFHT_sel->passes(event);
  if(pass_trigegr_PFHT){
    if((event.muons->size() == 1) && st_jets>900) fill_histograms(event, "triggerPFHT_mu", pass_ttbarsemilep); //FixME: PFHT trigger might have different st_jets threshold
    if((event.electrons->size() == 1) && st_jets>900) fill_histograms(event, "triggerPFHT_ele", pass_ttbarsemilep);
  }
 if(debug)   
   cout << "pass_trigegr_PFHT"<<endl;

  return true;
}

// as we want to run the ExampleCycleNew directly with AnalysisModuleRunner,
// make sure the TstarTstarSelectionModule is found by class name. This is ensured by this macro:
UHH2_REGISTER_ANALYSIS_MODULE(TstarTstarSelectionModule)

//}
