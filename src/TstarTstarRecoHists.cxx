#include "UHH2/TstarTstar/include/TstarTstarRecoHists.h"
#include "UHH2/core/include/Event.h"

#include "TH1F.h"
#include "TH2F.h"
#include <iostream>

using namespace std;
using namespace uhh2;
//using namespace uhh2examples;

TstarTstarRecoHists::TstarTstarRecoHists(Context & ctx, const string & dirname): Hists(ctx, dirname){
  // book all histograms here

  //M_Tstar
  h_M_Tstar_gluon_ = ctx.get_handle< float >("M_Tstar_gluon");
  h_M_Tstar_gamma_ = ctx.get_handle< float >("M_Tstar_gamma");
  book<TH1F>("M_Tstar_gluon", "M_{T^{*}_{g}} [GeV]", 30, 0, 3000);
  book<TH1F>("M_Tstar_gamma", "M_{T^{*}_{#gamma}} [GeV]", 30, 0, 3000);

  //Delta_R
  h_recohyp_ = ctx.get_handle<ReconstructionHypothesis>("TTbarReconstruction_best");
  h_is_ttbar_reconstructed_ = ctx.get_handle< bool >("is_ttbar_reconstructed_chi2");
  
  book<TH1F>("DeltaR_closest_toplepjet_ak8jet1", "DeltaR_toplepjet_ak8jet1", 40, 0, 4);
  book<TH1F>("DeltaR_closest_tophadjet_ak8jet1", "DeltaR_tophadjet_ak8jet1", 40, 0, 4);
  book<TH1F>("DeltaR_closest_toplepjet_ak8jet2", "DeltaR_toplepjet_ak8jet2", 40, 0, 4);
  book<TH1F>("DeltaR_closest_tophadjet_ak8jet2", "DeltaR_tophadjet_ak8jet2", 40, 0, 4);

  book<TH1F>("DeltaR_toplep_gamma", "DeltaR_toplep_gamma", 40, 0, 4);
  book<TH1F>("DeltaR_tophad_gamma", "DeltaR_tophad_gamma", 40, 0, 4);
  book<TH1F>("DeltaR_gamma_ak8jet1", "DeltaR_gamma_ak8jet1", 40, 0, 4);
  book<TH1F>("DeltaR_gamma_ak8jet2", "DeltaR_gamma_ak8jet2", 40, 0, 4);
}


void TstarTstarRecoHists::fill(const Event & event){
  // fill the histograms. Please note the comments in the header file:
  // 'hist' is used here a lot for simplicity, but it will be rather
  // slow when you have many histograms; therefore, better
  // use histogram pointers as members as in 'UHH2/common/include/ElectronHists.h'
  
  // Don't forget to always use the weight when filling.
  double weight = event.weight;
  
  hist("M_Tstar_gluon")->Fill(event.get(h_M_Tstar_gluon_), weight);
  hist("M_Tstar_gamma")->Fill(event.get(h_M_Tstar_gamma_), weight);

  ReconstructionHypothesis hyp = event.get(h_recohyp_);

  float DeltaR_min = 9999;
  if(event.topjets->size()>0){
    for (uint i = 0; i < hyp.toplep_jets().size(); i++){
      float DeltaR_new = deltaR(hyp.toplep_jets().at(i), event.topjets->at(0)); 
      if(DeltaR_new < DeltaR_min){
	DeltaR_min = DeltaR_new;
      }
    }
    hist("DeltaR_closest_toplepjet_ak8jet1")->Fill(DeltaR_min, weight);

    DeltaR_min = 9999;
    for (uint i = 0; i < hyp.tophad_jets().size(); i++){
      float DeltaR_new = deltaR(hyp.tophad_jets().at(i), event.topjets->at(0)); 
      if(DeltaR_new < DeltaR_min){
	DeltaR_min = DeltaR_new;
      }
    }
    hist("DeltaR_closest_tophadjet_ak8jet1")->Fill(DeltaR_min, weight);
  }
  if(event.topjets->size()>1){
    DeltaR_min = 9999;
    for (uint i = 0; i < hyp.toplep_jets().size(); i++){
      float DeltaR_new = deltaR(hyp.toplep_jets().at(i), event.topjets->at(1)); 
      if(DeltaR_new < DeltaR_min){
	DeltaR_min = DeltaR_new;
      }
    }
    hist("DeltaR_closest_toplepjet_ak8jet2")->Fill(DeltaR_min, weight);
  
    DeltaR_min = 9999;
    for (uint i = 0; i < hyp.tophad_jets().size(); i++){
      float DeltaR_new = deltaR(hyp.tophad_jets().at(i), event.topjets->at(1)); 
      if(DeltaR_new < DeltaR_min){
	DeltaR_min = DeltaR_new;
      }
    }
    hist("DeltaR_closest_tophadjet_ak8jet2")->Fill(DeltaR_min, weight);
  }

  //TODO All other DeltaRs
  if(event.photons->size()>0){hist("DeltaR_toplep_gamma")->Fill( deltaR(hyp.toplep_v4(), event.photons->at(0)), weight);}
  if(event.photons->size()>0){hist("DeltaR_tophad_gamma")->Fill( deltaR(hyp.tophad_v4(), event.photons->at(0)) , weight);}
  if(event.topjets->size()>0 && event.photons->size()>0){hist("DeltaR_gamma_ak8jet1")->Fill( deltaR(event.photons->at(0), event.topjets->at(0)) , weight);}
  if(event.topjets->size()>1 && event.photons->size()>0){hist("DeltaR_gamma_ak8jet2")->Fill( deltaR(event.photons->at(0), event.topjets->at(1)) , weight);}


}

TstarTstarRecoHists::~TstarTstarRecoHists(){}
