#include <iostream>
#include <memory>

#include <UHH2/core/include/AnalysisModule.h>
#include <UHH2/core/include/Event.h>
#include <UHH2/core/include/Selection.h>

#include <UHH2/common/include/CleaningModules.h>
#include <UHH2/common/include/CommonModules.h>
#include <UHH2/common/include/NSelections.h>
#include <UHH2/common/include/TriggerSelection.h>
#include <UHH2/common/include/JetCorrections.h>
#include <UHH2/common/include/ObjectIdUtils.h>
#include <UHH2/common/include/MuonIds.h>
#include <UHH2/common/include/ElectronIds.h>
#include <UHH2/common/include/JetIds.h>
#include <UHH2/common/include/TTbarGen.h>
#include <UHH2/common/include/TopJetIds.h>
#include <UHH2/common/include/Utils.h>
#include <UHH2/common/include/AdditionalSelections.h>
#include <UHH2/common/include/MCWeight.h>
#include <UHH2/common/include/TopPtReweight.h>

#include <UHH2/common/include/ElectronHists.h>
#include <UHH2/common/include/MuonHists.h>
#include <UHH2/common/include/LuminosityHists.h>
#include <UHH2/common/include/JetHists.h>
#include "UHH2/TstarTstar/include/CorrectionFactor.h"
#include "UHH2/TstarTstar/include/ElecTriggerSF.h"

/*
*******************************************************************
**************** TO DO ********************************************
*******************************************************************
- b tagging SF
*******************************************************************
*******************************************************************
*/

class TstarTstarTriggerSFModule : public AnalysisModule {

public:
  explicit TstarTstarTriggerSFModule(uhh2::Context&);
  virtual bool process(uhh2::Event&) override;

protected:
  enum lepton { muon, elec };
  lepton channel_;

  // cleaners & Correctors
  std::unique_ptr<CommonModules> common;
  std::unique_ptr<MuonCleaner>     muoSR_cleaner;
  std::unique_ptr<ElectronCleaner> eleSR_cleaner;


  std::unique_ptr<JetCleaner> jet_cleaner1;
  std::unique_ptr<JetCleaner> jet_cleaner2;

  // Btag efficiency hists
  std::unique_ptr<BTagMCEfficiencyHists> BTagEffHists;

  // selections

  std::unique_ptr<uhh2::Selection> trigger_mu_A;
  std::unique_ptr<uhh2::Selection> trigger_mu_B;
  std::unique_ptr<uhh2::Selection> trigger_el_A;
  std::unique_ptr<uhh2::Selection> trigger_el_B;
  std::unique_ptr<uhh2::Selection> trigger_el_C;
  std::unique_ptr<uhh2::Selection> muon_sel;
  std::unique_ptr<uhh2::Selection> elec_sel1;
  std::unique_ptr<uhh2::Selection> elec_sel2;
  std::unique_ptr<uhh2::Selection> elec_etaveto;
  std::unique_ptr<uhh2::Selection> met_sel;
  std::unique_ptr<uhh2::Selection> pv_sel;
  std::unique_ptr<uhh2::Selection> twodcut_sel;
  std::unique_ptr<uhh2::Selection> sel_badhcal;

  Event::Handle<bool>h_recsel;
  Event::Handle<double>h_pt;
  Event::Handle<double>h_eta;
  Event::Handle<double>h_weight;
  Event::Handle<double>h_weight_SFpt;
  Event::Handle<double>h_weight_SFeta;
  Event::Handle<double>h_weight_SFetapt;
  Event::Handle<double>h_weight_SFetaptUP;
  Event::Handle<double>h_weight_SFetaptDOWN;
  Event::Handle<bool>h_passed;
  Event::Handle<bool>h_passed_elec;
  Event::Handle<bool>h_passed_photon;
  Event::Handle<int>h_run;
  Event::Handle<int>h_lumi;
  Event::Handle<int>h_eventnr;

  std::unique_ptr<uhh2::AnalysisModule> ele_id_SF, ele_trigger_SFpt, ele_trigger_SFeta, ele_reco_SF;
  std::unique_ptr<uhh2::AnalysisModule> ele_trigger_SFetapt, ele_trigger_SFetaptUP, ele_trigger_SFetaptDOWN;
  std::unique_ptr<uhh2::AnalysisModule> muo_tight_noniso_SF, muo_trigger_SF, muo_trigger_SF_B;

  bool debug = false;
  bool isMC; //define here to use it in "process" part
  bool year_16;
  bool year_17;
  bool year_18;
  Year year;

  std::unique_ptr<Hists> h_pass, h_all;
};

TstarTstarTriggerSFModule::TstarTstarTriggerSFModule(uhh2::Context& ctx){

  if(debug) cout << "Get Year ... " << endl;
  year_16 = false;
  year_17 = false;
  year_18 = false;
  year = extract_year(ctx);

  if(year == Year::is2016v3) year_16 = true;
  else if(year == Year::is2017v2) year_17 = true;
  else if(year == Year::is2018) year_18 = true;
  else throw runtime_error("In PostSelectionModule: This Event is not from 2016v3, 2017v2 or 2018!");

  isMC = (ctx.get("dataset_type") == "MC");

  if(debug) cout << "Declare Output ... " << endl;
  ctx.undeclare_all_event_output();
  h_eta = ctx.declare_event_output<double>("eta");
  h_pt = ctx.declare_event_output<double>("pt");
  h_weight = ctx.declare_event_output<double>("weight");
  h_weight_SFpt = ctx.declare_event_output<double>("weight_sfpt");
  h_weight_SFeta = ctx.declare_event_output<double>("weight_sfeta");
  h_weight_SFetapt = ctx.declare_event_output<double>("weight_sfetapt");
  h_weight_SFetaptUP = ctx.declare_event_output<double>("weight_sfetaptUP");
  h_weight_SFetaptDOWN = ctx.declare_event_output<double>("weight_sfetaptDOWN");
  h_passed = ctx.declare_event_output<bool>("passed");
  h_passed_elec = ctx.declare_event_output<bool>("passed_elec");
  h_passed_photon = ctx.declare_event_output<bool>("passed_photon");
  h_run = ctx.declare_event_output<int>("run");
  h_lumi = ctx.declare_event_output<int>("lumi");
  h_eventnr = ctx.declare_event_output<int>("eventnr");

  if(debug) cout << "Define Trigger ... " << endl;
  // define Trigger
  // trigger_mu_A = uhh2::make_unique<TriggerSelection>("HLT_Mu24_v*");
  // if(year == Year::is2017v2) trigger_mu_A = uhh2::make_unique<TriggerSelection>("HLT_Mu27_v*");
  // trigger_mu_B = uhh2::make_unique<TriggerSelection>("HLT_TkMu24_v*");
  trigger_mu_A = uhh2::make_unique<TriggerSelection>("HLT_Mu50_v*");
  trigger_mu_B = uhh2::make_unique<TriggerSelection>("HLT_TkMu50_v*");
  if(year_16)      trigger_el_A = uhh2::make_unique<TriggerSelection>("HLT_Ele27_WPTight_Gsf_v*");
  else if(year_17) trigger_el_A = uhh2::make_unique<TriggerSelection>("HLT_Ele35_WPTight_Gsf_v*");
  else if(year_18) trigger_el_A = uhh2::make_unique<TriggerSelection>("HLT_Ele32_WPTight_Gsf_v*");
  trigger_el_B = uhh2::make_unique<TriggerSelection>("HLT_Ele115_CaloIdVT_GsfTrkIdT_v*");
  if(year_16) trigger_el_C = uhh2::make_unique<TriggerSelection>("HLT_Photon175_v*");
  else        trigger_el_C = uhh2::make_unique<TriggerSelection>("HLT_Photon200_v*");

  // Scale Factors
  TString syear;
  if(year_16) syear = "2016";
  if(year_17) syear = "2017";
  if(year_18) syear = "2018";
  ele_trigger_SFpt.reset(new ElecTriggerSF(ctx, "nominal", "pt", syear));
  ele_trigger_SFeta.reset(new ElecTriggerSF(ctx, "nominal", "eta", syear));
  ele_trigger_SFetapt.reset(new ElecTriggerSF(ctx, "nominal", "eta_ptbins", syear));
  ele_trigger_SFetaptUP.reset(new ElecTriggerSF(ctx, "up", "eta_ptbins", syear));
  ele_trigger_SFetaptDOWN.reset(new ElecTriggerSF(ctx, "down", "eta_ptbins", syear));

  h_pass.reset(new ElectronHists(ctx, "pass_Elec"));
  h_all.reset(new ElectronHists(ctx, "all_Elec"));

}


bool TstarTstarTriggerSFModule::process(uhh2::Event& event){
  if(debug) cout << " ------------------------------------------- " << endl;
  if(debug) cout << " ----------------- NewEvent ---------------- " << endl;
  if(debug) cout << " ------------------------------------------- " << endl;

  /* *********** Trigger *********** */
  if(debug) cout << "Trigger... " << endl;

  // only use high-pt muon region!
  if(event.muons->at(0).pt() < 60) return false;

  if(year == Year::is2016v3){
    if(!isMC && event.run < 274954){
      if(!trigger_mu_A->passes(event)) return false;
    }
    else{
      if( !(trigger_mu_A->passes(event) || trigger_mu_B->passes(event)) ) return false;
    }
  }
  if(year == Year::is2017v2){
    if(!trigger_mu_A->passes(event)) return false;
  }
  if(year == Year::is2018){
    if(!trigger_mu_A->passes(event)) return false;
  }

  // fill hists with all events
  h_all->fill(event);

  // HERE FILL PT AND ETA HISTS FOR PASSING AND NOT PASSING ELEC TRIGGER
  if(debug) cout << "Start Fill ... " << endl;
  bool passed_elec_trigger = true;
  if(year == Year::is2016v3)  passed_elec_trigger = (trigger_el_B->passes(event) || trigger_el_C->passes(event));
  if(year == Year::is2017v2){
    // for MC event.run=1
    if(!isMC && event.run <= 299329) passed_elec_trigger = (trigger_el_A->passes(event) || trigger_el_C->passes(event));
    else                             passed_elec_trigger = (trigger_el_B->passes(event) || trigger_el_C->passes(event));
  }
  if(year == Year::is2018)  passed_elec_trigger = (trigger_el_B->passes(event) || trigger_el_C->passes(event));

  /* Additional studies for kink in data at 200 GeV */
  bool pass_elec_only = true; bool pass_photon_only = true;
  if(year == Year::is2016v3){
    pass_elec_only = trigger_el_B->passes(event);
    pass_photon_only = trigger_el_C->passes(event);
  }
  if(year == Year::is2017v2){
    // for MC event.run=1
    if(!isMC && event.run <= 299329){
      pass_elec_only = trigger_el_A->passes(event);
      pass_photon_only = trigger_el_C->passes(event);
    }
    else{
      pass_elec_only = trigger_el_B->passes(event);
      pass_photon_only = trigger_el_C->passes(event);
    }
  }
  if(year == Year::is2018){
    pass_elec_only = trigger_el_B->passes(event);
    pass_photon_only = trigger_el_C->passes(event);
  }

  // bool passed_elec_trigger = false;
  // if(pass_elec1) passed_elec_trigger = (trigger_el_A->passes(event) || trigger_el_C->passes(event)); // noHighTrigger

  if(debug) cout << "Set after Trigger pass... " << endl;
  event.set(h_pt, event.electrons->at(0).pt());
  event.set(h_eta, event.electrons->at(0).eta());
  event.set(h_weight, event.weight);
  event.set(h_passed, passed_elec_trigger);
  event.set(h_passed_elec, pass_elec_only);
  event.set(h_passed_photon, pass_photon_only);

  int run = 0;
  int lumi = 0;
  int eventnr = 0;

  if(!isMC){
    run = event.run;
    lumi = event.luminosityBlock;
    eventnr = event.event;
  }

  if(debug) cout << "Set run info ... " << endl;
  event.set(h_run, run);
  event.set(h_lumi, lumi);
  event.set(h_eventnr, eventnr);

  if(debug) cout << "Apply SFs ... " << endl;
  // now apply SF and store weight after SF
  // first store raw weight
  double wraw = event.weight;
  // use pt SF and store new weight
  ele_trigger_SFpt->process(event);
  event.set(h_weight_SFpt, event.weight);
  event.weight = wraw;
  // reset weight and do eta SF
  ele_trigger_SFeta->process(event);
  event.set(h_weight_SFeta, event.weight);
  event.weight = wraw;
  // reset weight and do eta_ptbin SF
  ele_trigger_SFetapt->process(event);
  event.set(h_weight_SFetapt, event.weight);
  event.weight = wraw;
  // up variation
  ele_trigger_SFetaptUP->process(event);
  event.set(h_weight_SFetaptUP, event.weight);
  event.weight = wraw;
  // down variation
  ele_trigger_SFetaptDOWN->process(event);
  event.set(h_weight_SFetaptDOWN, event.weight);

  if(passed_elec_trigger){
    // fill pass histograms
    h_pass->fill(event);
  }

  // if(debug) cout << "Event done ... " << endl;
  if(debug) cout << "Event done ... " << endl;
  return true;
}

UHH2_REGISTER_ANALYSIS_MODULE(TstarTstarTriggerSFModule)
