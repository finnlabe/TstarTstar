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


#include "UHH2/HOTVR/include/HOTVRJetCorrectionModule.h"
#include "UHH2/TstarTstar/include/B2G_TrimMass_Hists.h"



using namespace std;
using namespace uhh2;

namespace uhh2 {

/**
 *
 * Simple quick analysis module to produce trigger efficiencies.
 *
 */
class B2G_TrimMass_Trigger_Studies: public AnalysisModule {
public:

  explicit B2G_TrimMass_Trigger_Studies(Context & ctx);
  virtual bool process(Event & event) override;

private:


  // module for cleaning
  std::unique_ptr<CommonModules> common; // ak4 (not used)
  unique_ptr<TopJetCleaner> AK8Cleaner; // ak8


  // ##### triggers #####
  unique_ptr<Selection> trg_mu_ref_1;
  unique_ptr<Selection> trg_mu_ref_2;
  unique_ptr<Selection> trg_HT1050;
  unique_ptr<Selection> trg_PFHJet500, trg_PFHJet550;
  unique_ptr<Selection> trg_TrimMass30, trg_TrimMass30_420;
  unique_ptr<Selection> trg_TrimMass50;

  // ##### Histograms classes #####
  // will all be objects of the B2G_TrimMass_Hists class, containting multiple histograms each time!
  std::unique_ptr<Hists> h_before;
  std::unique_ptr<Hists> h_afterHT, h_afterPFJet500, h_afterTrimMass30, h_afterTrimMass50;
  std::unique_ptr<Hists> h_afterAll, h_afterTrimMass30removed, h_currentPlan, h_PFJet420_TrimMass30, h_PFJet500_TrimMass30;

  // for an independence test -> calculating the alpha value
  std::unique_ptr<Hists> h_alpha_mu, h_alpha_jets, h_alpha_both, h_alpha_ref;

};


B2G_TrimMass_Trigger_Studies::B2G_TrimMass_Trigger_Studies(Context & ctx){

  // creating the objects
  trg_mu_ref_1.reset(new TriggerSelection("HLT_IsoMu24_v*"));
  trg_mu_ref_2.reset(new TriggerSelection("HLT_Mu50_v*"));

  trg_HT1050.reset(new TriggerSelection("HLT_PFHT1050_v*"));
  trg_PFHJet500.reset(new TriggerSelection("HLT_AK8PFJet500_v*"));
  trg_PFHJet550.reset(new TriggerSelection("HLT_AK8PFJet550_v*"));
  trg_TrimMass30.reset(new TriggerSelection("HLT_AK8PFJet400_TrimMass30_v*"));
  trg_TrimMass30_420.reset(new TriggerSelection("HLT_AK8PFJet420_TrimMass30_v*"));
  trg_TrimMass50.reset(new TriggerSelection("HLT_AK8PFHT800_TrimMass50_v*"));

  // Cleaners
  common.reset(new CommonModules());
  common->switch_jetlepcleaner();
  common->switch_jetPtSorter();
  common->set_jet_id(AndId<Jet>(PtEtaCut(30., 2.4), JetPFID(JetPFID::WP_TIGHT_PUPPI)));
  common->init(ctx);
  AK8Cleaner.reset(new TopJetCleaner(ctx, PtEtaCut(200.0, 2.4)));

  // histograms
  h_before.reset(new B2G_TrimMass_Hists(ctx, "TrimMass_before"));
  h_afterHT.reset(new B2G_TrimMass_Hists(ctx, "TrimMass_afterHT"));
  h_afterPFJet500.reset(new B2G_TrimMass_Hists(ctx, "TrimMass_afterPFJet500"));
  h_afterTrimMass30.reset(new B2G_TrimMass_Hists(ctx, "TrimMass_TrimMass30"));
  h_afterTrimMass50.reset(new B2G_TrimMass_Hists(ctx, "TrimMass_TrimMass50"));
  h_afterAll.reset(new B2G_TrimMass_Hists(ctx, "TrimMass_afterAll"));
  h_afterTrimMass30removed.reset(new B2G_TrimMass_Hists(ctx, "TrimMass_TrimMass30removed"));
  h_currentPlan.reset(new B2G_TrimMass_Hists(ctx, "TrimMass_currentPlan"));
  h_PFJet420_TrimMass30.reset(new B2G_TrimMass_Hists(ctx, "TrimMass_PFJet420_TrimMass30"));
  h_PFJet500_TrimMass30.reset(new B2G_TrimMass_Hists(ctx, "TrimMass_PFJet500_TrimMass30"));

  //h_alpha_ref.reset(new B2G_TrimMass_Hists(ctx, "h_alpha_ref"));
  //h_alpha_mu.reset(new B2G_TrimMass_Hists(ctx, "h_alpha_mu"));
  //h_alpha_jets.reset(new B2G_TrimMass_Hists(ctx, "h_alpha_jets"));
  //h_alpha_both.reset(new B2G_TrimMass_Hists(ctx, "h_alpha_both"));


}


bool B2G_TrimMass_Trigger_Studies::process(Event & event) {

  // applying cleaners
  if(!(common->process(event))) return false;
  if(!(AK8Cleaner->process(event))) return false;

  // getting trigger booleans
  bool pass_mu_ref = (trg_mu_ref_1->passes(event) || trg_mu_ref_2->passes(event));
  bool pass_HT = trg_HT1050->passes(event);
  bool pass_PFJet = trg_PFHJet500->passes(event);
  bool pass_PFJet_550 = trg_PFHJet550->passes(event);
  bool pass_TrimMass30 = trg_TrimMass30->passes(event);
  bool pass_TrimMass50 = trg_TrimMass50->passes(event);
  bool pass_TrimMass30_420 = trg_TrimMass30_420->passes(event);

  // filling hists for alpha value
  /**
  h_alpha_ref->fill(event);
  if(pass_mu_ref) h_alpha_mu->fill(event);
  bool pass_fullset = pass_HT || pass_PFJet || pass_TrimMass30 || pass_TrimMass50;
  if(pass_fullset) h_alpha_jets->fill(event);
  if(pass_fullset && pass_mu_ref) h_alpha_both->fill(event);
  **/

  // single mu reference selection
  if(!pass_mu_ref) return false;

  // requiring at least two ak8 jets
  if(event.topjets->size() < 2) return false;

  // requiring deltaeta < 1.3 for those jets
  double deltaEta = abs(event.topjets->at(0).eta() - event.topjets->at(1).eta());
  if(deltaEta > 1.3) return false;

  // histogram before triggers
  h_before->fill(event);

  // histograms for various trigger (combinations)
  if(pass_HT) h_afterHT->fill(event);
  if(pass_PFJet) h_afterPFJet500->fill(event);
  if(pass_TrimMass30) h_afterTrimMass30->fill(event);
  if(pass_TrimMass50) h_afterTrimMass50->fill(event);
  if(pass_HT || pass_PFJet || pass_TrimMass30 || pass_TrimMass50) h_afterAll->fill(event);
  if(pass_HT || pass_PFJet || pass_TrimMass50) h_afterTrimMass30removed->fill(event);
  if(pass_HT || pass_PFJet_550 || pass_TrimMass50) h_currentPlan->fill(event);
  if(pass_HT || pass_PFJet_550 || pass_TrimMass50 || pass_TrimMass30_420) h_PFJet420_TrimMass30->fill(event);
  if(pass_HT || pass_PFJet_550 || pass_TrimMass50 || (pass_TrimMass30 && pass_PFJet)) h_PFJet500_TrimMass30->fill(event);

  return false;
}

// as we want to run the ExampleCycleNew directly with AnalysisModuleRunner,
// make sure the TstarTstarPreselectionModule is found by class name. This is ensured by this macro:
UHH2_REGISTER_ANALYSIS_MODULE(B2G_TrimMass_Trigger_Studies)

}
