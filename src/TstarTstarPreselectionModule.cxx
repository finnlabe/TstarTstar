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
//#include "UHH2/common/include/TriggerSelection.h"

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
    std::unique_ptr<MuonCleaner>                     muon_cleaner;//DEBUG

    // Declare the Selections to use. Use unique_ptr to ensure automatic call of delete in the destructor,
    // to avoid memory leaks.
    unique_ptr<Selection> twodcut_sel;   
    unique_ptr<Selection> TTbarSemiLepMatchable_selection;
    // std::unique_ptr<uhh2::Selection> trigger40_sel;
    // std::unique_ptr<uhh2::Selection> trigger60_sel;
    // std::unique_ptr<uhh2::Selection> trigger80_sel;
    // std::unique_ptr<uhh2::Selection> trigger140_sel;
    // std::unique_ptr<uhh2::Selection> trigger200_sel;
    // std::unique_ptr<uhh2::Selection> trigger260_sel;
    // std::unique_ptr<uhh2::Selection> trigger320_sel;
    // std::unique_ptr<uhh2::Selection> trigger400_sel;
    // std::unique_ptr<uhh2::Selection> trigger450_sel;
    // std::unique_ptr<uhh2::Selection> trigger500_sel;

    // unique_ptr<Selection> triggerHT1_sel, triggerHT2_sel, triggerHT3_sel, triggerHT4_sel, triggerHT5_sel,  triggerHT6_sel;

    // Store the Hists collection as member variables. Again, use unique_ptr to avoid memory leaks.
  std::unique_ptr<Hists> h_nocuts, h_common, h_lepsel, h_2dcut, h_semilepttbarmatch, h_nosemilepttbarmatch;// h_trigger;
  //    std::unique_ptr<Hists> h_semilepttbarmatch_gen;
  //    std::unique_ptr<Hists> h_semilepttbarmatch_genreco;
  bool debug = false;
  //  bool debug = true;
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

    // 1. setup other modules. CommonModules
    common.reset(new CommonModules());
    common->switch_metcorrection();
    common->switch_jetlepcleaner();
    common->switch_jetPtSorter();
    common->set_jet_id(AndId<Jet>(JetPFID(JetPFID::WP_TIGHT_PUPPI), PtEtaCut(30.0,2.1)));
    ElectronId eleID; 
    double electron_pt(20.);
    eleID = ElectronID_Summer16_tight_noIso;
    common->set_electron_id(AndId<Electron>(PtEtaSCCut(electron_pt, 1.9), eleID));
    PhotonId phoID; 
    double photon_pt(20.);
    //    phoID = PhotonTagID(Photon::cutBasedPhotonID_Spring16_V2p2_tight);
    phoID = PhotonTagID(Photon::cutBasedPhotonID_Fall17_94X_V2_loose);
    common->set_photon_id(AndId<Photon>(PtEtaCut(photon_pt, 2.1), phoID));
    MuonId muID;
    double muon_pt(20.);
    muID = MuonID(Muon::Highpt);
    const MuonId muonID(AndId<Muon>(PtEtaCut(muon_pt, 1.9), muID));
    //    const MuonId muonID(PtEtaCut(muon_pt, 1.9));
    //    muon_cleaner.reset(new MuonCleaner(muonID));
    common->set_muon_id(muonID);
    //    common->set_muon_id(AndId<Muon>(PtEtaCut(muon_pt, 1.9), muID));
    common->init(ctx);
    
    // 2. set up selections
    ///2D Cut Lepton-Jets
    twodcut_sel.reset(new TwoDCut(0.4, 25.0));  // The same as in Z'->ttbar semileptonic
    TTbarSemiLepMatchable_selection.reset(new TTbarSemiLepMatchableSelection());

    // #define GET_RESET_TRIGGER(trg_name)					              \
    // const std::string& trg_name = ctx.get( #trg_name , "NULL");	                      \
    // if ( trg_name != "NULL") trg_name##_sel.reset(new TriggerSelection( trg_name ));  \
    // else trg_name##_sel.reset(new uhh2::AndSelection(ctx));                           \

    // GET_RESET_TRIGGER(trigger40)
    // GET_RESET_TRIGGER(trigger60)
    // GET_RESET_TRIGGER(trigger80)
    // GET_RESET_TRIGGER(trigger140)
    // GET_RESET_TRIGGER(trigger200)
    // GET_RESET_TRIGGER(trigger260)
    // GET_RESET_TRIGGER(trigger320)
    // GET_RESET_TRIGGER(trigger400)
    // GET_RESET_TRIGGER(trigger450)
    // GET_RESET_TRIGGER(trigger500)

    // triggerHT1_sel.reset(new TriggerSelection("HLT_HT430to450_v*"));
    // triggerHT2_sel.reset(new TriggerSelection("HLT_HT450to470_v*"));
    // triggerHT3_sel.reset(new TriggerSelection("HLT_HT470to500_v*"));
    // triggerHT4_sel.reset(new TriggerSelection("HLT_HT500to550_v*"));
    // triggerHT5_sel.reset(new TriggerSelection("HLT_HT550to650_v*"));
    // triggerHT6_sel.reset(new TriggerSelection("HLT_HT650_v*"));

    // 3. Set up Hists classes:
    h_nocuts.reset(new TstarTstarHists(ctx, "NoCuts"));
    h_common.reset(new TstarTstarHists(ctx, "AfterCommon"));
    h_lepsel.reset(new TstarTstarHists(ctx, "AfterLepSel"));
    h_2dcut.reset(new TstarTstarHists(ctx, "After2D"));
    h_semilepttbarmatch.reset(new TstarTstarHists(ctx, "SemiLepTTBarMatch"));
    h_nosemilepttbarmatch.reset(new TstarTstarHists(ctx, "NotSemiLepTTBarMatch"));
    //    h_trigger.reset(new TstarTstarHists(ctx, "AfterTrigger"));
    //    h_semilepttbarmatch_gen.reset(new TstarTstarGenHists(ctx, "SemiLepTTBarMatchGEN"));
    //    h_semilepttbarmatch_genreco.reset(new TstarTstarGenRecoMatchedHists(ctx, "SemiLepTTBarMatchGENRECO"));

   
}


bool TstarTstarPreselectionModule::process(Event & event) {
   
  if(debug)   
    cout << "TstarTstarPreselectionModule: Starting to process event (runid, eventid) = (" << event.run << ", " << event.event << "); weight = " << event.weight << endl;
  if(debug)
    cout<<"N muons = "<<event.muons->size()<<", N electrons = "<<event.electrons->size()<<", N photons = "<<event.photons->size()<<endl;
  for(auto& muo : *event.muons){
    if(debug) cout<<"BEFORE Muon (pt,eta): "<<muo.pt()<<", "<<muo.eta()<<endl;
  }
  h_nocuts->fill(event);
  if(debug) cout<<"Filled hists without any cuts"<<endl;

  bool result_com_sel = common->process(event);
  if(!result_com_sel) return false;
  h_common->fill(event);
  if(debug) cout<<"Filled hists after common modules"<<endl;

  //  muon_cleaner->process(event);//debug

  //---- Loose selection
  // Require at least one Muon or one Electron
  const bool pass_lep1 = ((event.muons->size() >= 1) || (event.electrons->size() >= 1));
  if(!pass_lep1) return false;
  h_lepsel->fill(event);
  if(debug) cout<<"Filled hists after lepton"<<endl;
  if(debug)
    cout<<"N muons = "<<event.muons->size()<<", N electrons = "<<event.electrons->size()<<", N photons = "<<event.photons->size()<<endl;

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

  // bool pass_trigger40=false; bool pass_trigger60=false; bool pass_trigger80=false;
  // bool pass_trigger140=false; bool pass_trigger200=false; bool pass_trigger260=false;
  // bool pass_trigger320=false; bool pass_trigger400=false; bool pass_trigger450=false; bool pass_trigger500=false;
  // //This should be moved to Selection part
  // if(event.jets->size()<1) return false;//FixMe: why this is happening?
  // double pt_leadjet = event.jets->at(0).pt();
  // std::vector<double> trg_thresh;
  // const int n_pt_bins_Si = 10;
  // const double pt_bins_Si[n_pt_bins_Si]       = { 40, 72,  95, 160, 226, 283, 344, 443, 577, 606};//FixMe: are thresholds correct?
  // for (int i = 0; i < n_pt_bins_Si; i++) trg_thresh.push_back(pt_bins_Si[i]);
  // pass_trigger40        = (trigger40_sel->passes(event)         && pt_leadjet>trg_thresh[0]   && pt_leadjet<trg_thresh[1]);
  // pass_trigger60        = (trigger60_sel->passes(event)         && pt_leadjet>trg_thresh[1]   && pt_leadjet<trg_thresh[2]);
  // pass_trigger80        = (trigger80_sel->passes(event)         && pt_leadjet>trg_thresh[2]   && pt_leadjet<trg_thresh[3]);
  // pass_trigger140       = (trigger140_sel->passes(event)        && pt_leadjet>trg_thresh[3]   && pt_leadjet<trg_thresh[4]);
  // pass_trigger200       = (trigger200_sel->passes(event)        && pt_leadjet>trg_thresh[4]   && pt_leadjet<trg_thresh[5]);
  // pass_trigger260       = (trigger260_sel->passes(event)        && pt_leadjet>trg_thresh[5]   && pt_leadjet<trg_thresh[6]);
  // pass_trigger320       = (trigger320_sel->passes(event)        && pt_leadjet>trg_thresh[6]   && pt_leadjet<trg_thresh[7]);
  // pass_trigger400       = (trigger400_sel->passes(event)        && pt_leadjet>trg_thresh[7]   && pt_leadjet<trg_thresh[8]);
  // pass_trigger450       = (trigger450_sel->passes(event)        && pt_leadjet>trg_thresh[8]   && pt_leadjet<trg_thresh[9]);
  // pass_trigger500       = (trigger500_sel->passes(event)        && pt_leadjet>trg_thresh[9]);
  // //  bool pass_trigger_any = (pass_trigger40 || pass_trigger60 || pass_trigger80 || pass_trigger140 || pass_trigger200 || pass_trigger260 || pass_trigger320 || pass_trigger400 || pass_trigger450 || pass_trigger500);
  // //  bool pass_trigger_any = (pass_trigger450 || pass_trigger500);

  // double HT = 0;
  // std::vector<double> trg_thresh_HT;
  // const int n_HT_bins = 6;
  // const double HT_bins[n_HT_bins] = { 430, 450,  470, 500, 550, 650};//FixMe: thresholds are wrong!
  // for (int i = 0; i < n_HT_bins; i++) trg_thresh_HT.push_back(HT_bins[i]);
  // bool pass_triggerHT1 = (triggerHT1_sel->passes(event) && HT>trg_thresh_HT[0]);
  // bool pass_triggerHT2 = (triggerHT2_sel->passes(event) && HT>trg_thresh_HT[1]);
  // bool pass_triggerHT3 = (triggerHT3_sel->passes(event) && HT>trg_thresh_HT[2]);
  // bool pass_triggerHT4 = (triggerHT4_sel->passes(event) && HT>trg_thresh_HT[3]);
  // bool pass_triggerHT5 = (triggerHT5_sel->passes(event) && HT>trg_thresh_HT[4]);
  // bool pass_triggerHT6 = (triggerHT6_sel->passes(event) && HT>trg_thresh_HT[5]);

  // bool pass_trigger_any = (pass_triggerHT1 || pass_triggerHT2 || pass_triggerHT3 || pass_triggerHT4 || pass_triggerHT5 || pass_triggerHT6);

  // if(debug) cout<<"pass_triggers = "<<pass_trigger40<<" "<<pass_trigger60<<" "<<pass_trigger80<<" "<<pass_trigger140
  // 		<<" "<<pass_trigger200<<" "<<pass_trigger260<<" "<<pass_trigger320<<" "<<pass_trigger400<<" "<<pass_trigger450<<" "<<pass_trigger500<<endl;
  // if(!pass_trigger_any) return false;
  // h_trigger->fill(event);



  // //---- Matching to GEN
  if(TTbarSemiLepMatchable_selection->passes(event)){
    h_semilepttbarmatch->fill(event);
  //   h_semilepttbarmatch_gen->fill(event);
  //   h_semilepttbarmatch_genreco->fill(event);
  }
  else h_nosemilepttbarmatch->fill(event);
  if(debug) cout<<"TTbarSemiLepMatchable_selection is ok"<<endl;

  return true;
}

// as we want to run the ExampleCycleNew directly with AnalysisModuleRunner,
// make sure the TstarTstarPreselectionModule is found by class name. This is ensured by this macro:
UHH2_REGISTER_ANALYSIS_MODULE(TstarTstarPreselectionModule)

}
