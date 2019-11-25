#include "UHH2/TstarTstar/include/TstarTstarRecoTstarHists.h"
#include "UHH2/core/include/Event.h"

#include "TH1F.h"
#include "TH2F.h"
#include <iostream>

using namespace std;
using namespace uhh2;
//using namespace uhh2examples;

namespace {
    
  // invariant mass of a lorentzVector, but safe for timelike / spacelike vectors
  float inv_mass(const LorentzVector& p4){ return p4.isTimelike() ? p4.mass() : -sqrt(-p4.mass2()); }

}

TstarTstarRecoTstarHists::TstarTstarRecoTstarHists(Context & ctx, const string & dirname): Hists(ctx, dirname){
  // book all histograms here

  h_recohyp_tstartstar_tgtgamma_best_ = ctx.get_handle<ReconstructionTstarHypothesis>("TstarTstar_tgtgamma_best");

  //M_Tstar
  book<TH1F>("M_Tstar_gluon", "M_{T^{*}_{g}} [GeV]", 30, 0, 3000);
  book<TH1F>("M_Tstar_gamma", "M_{T^{*}_{#gamma}} [GeV]", 30, 0, 3000);

  //Delta_R
  h_recohyp_ = ctx.get_handle<ReconstructionHypothesis>("TTbarReconstruction_best");
  h_is_ttbar_reconstructed_ = ctx.get_handle< bool >("is_ttbar_reconstructed_chi2");

  // book<TH1F>("DeltaR_closest_toplepjet_ak8jet1", "DeltaR_toplepjet_ak8jet1", 50, 0, 5);
  // book<TH1F>("DeltaR_closest_tophadjet_ak8jet1", "DeltaR_tophadjet_ak8jet1", 50, 0, 5);
  // book<TH1F>("DeltaR_closest_toplepjet_ak8jet2", "DeltaR_toplepjet_ak8jet2", 50, 0, 5);
  // book<TH1F>("DeltaR_closest_tophadjet_ak8jet2", "DeltaR_tophadjet_ak8jet2", 50, 0, 5);

  book<TH1F>("DeltaR_toplep_gamma", "DeltaR_toplep_gamma", 30, 0, 6);
  book<TH1F>("DeltaR_tophad_gamma", "DeltaR_tophad_gamma", 30, 0, 6);
  book<TH1F>("DeltaR_toplep_ak8jet1", "DeltaR_toplep_ak8jet1", 30, 0, 6);
  book<TH1F>("DeltaR_tophad_ak8jet1", "DeltaR_tophad_ak8jet1", 30, 0, 6);
  book<TH1F>("DeltaR_toplep_ak8jet2", "DeltaR_toplep_ak8jet2", 30, 0, 6);
  book<TH1F>("DeltaR_tophad_ak8jet2", "DeltaR_tophad_ak8jet2", 30, 0, 6);

  book<TH1F>("DeltaR_gamma_ak8jet1", "DeltaR_gamma_ak8jet1", 30, 0, 6);
  book<TH1F>("DeltaR_gamma_ak8jet2", "DeltaR_gamma_ak8jet2", 30, 0, 6);

  book<TH1F>("DeltaR_min_top_gamma", "DeltaR_min_top_gamma", 30, 0, 6);
  book<TH1F>("DeltaR_min_top_ak8jet1", "DeltaR_min_top_ak8jet1", 30, 0, 6);  
  book<TH1F>("DeltaR_min_top_ak8jet2", "DeltaR_min_top_ak8jet2", 30, 0, 6);  


  //For tests of ttbar reconstruction
  h_ttbar_hyps_ = ctx.get_handle< std::vector<ReconstructionHypothesis> >("TTbarReconstruction");
  book<TH1F>("Mtoplep_all", "M_{best match}^{toplep} [GeV]", 200, 0, 1000);
  book<TH1F>("Mtophad_all", "M_{best macth}^{tophad} [GeV]", 200, 0, 1000);

  //  h_ttbar_hyp_best_ = ctx.get_handle<ReconstructionHypothesis>("TTbarReconstruction_best");
  book<TH1F>("Mtoplep_best", "M^{best chi2}_{toplep} [GeV]", 200, 0, 1000);
  book<TH1F>("Mtophad_best", "M^{best chi2}_{tophad} [GeV]", 200, 0, 1000);


  //For tests of tg+tg reconstruction (with TstarTstar_tgluon_tgluon_Reconstruction class)
  h_recohyp_tstartstar_tgtg_best_ = ctx.get_handle<ReconstructionTstarHypothesis>("TstarTstar_tgtg_best");
  // h_recohyp_tstartstar_tgtg_ = ctx.get_handle< std::vector<ReconstructionHypothesis> >("TstarTstar_tgtg_hyps");
  book<TH1F>("M_Tstar_lep", "M_{T^{*}_{tgtg,lep}} [GeV]", 30, 0, 3000);
  book<TH1F>("M_Tstar_had", "M_{T^{*}_{tgtg,had}} [GeV]", 30, 0, 3000);

}


void TstarTstarRecoTstarHists::fill(const Event & event){
  // fill the histograms. Please note the comments in the header file:
  // 'hist' is used here a lot for simplicity, but it will be rather
  // slow when you have many histograms; therefore, better
  // use histogram pointers as members as in 'UHH2/common/include/ElectronHists.h'
  
  // Don't forget to always use the weight when filling.
  double weight = event.weight;
  
  // Problem when trying to access the wrong one (empty handle = crash)
  // Possible fix: save channel somewhere and then get correct recohyp!
  ReconstructionTstarHypothesis tstar_hyp_best = event.get(h_recohyp_tstartstar_tgtgamma_best_);

  float M_Tstar_gluon_ = inv_mass(tstar_hyp_best.tstar1gluon_v4());
  if(M_Tstar_gluon_ == 0){cout << "Now the mass is somehow zero... This is wrong" << endl;}

  hist("M_Tstar_gluon")->Fill(M_Tstar_gluon_, weight);
  hist("M_Tstar_gamma")->Fill(inv_mass(tstar_hyp_best.tstar1gamma_v4()), weight);

  ReconstructionHypothesis hyp = event.get(h_recohyp_);

// float DeltaR_min = 9999;
//   if(event.topjets->size()>0){
//     for (uint i = 0; i < hyp.toplep_jets().size(); i++){
//       float DeltaR_new = deltaR(hyp.toplep_jets().at(i), event.topjets->at(0)); 
//       if(DeltaR_new < DeltaR_min){
// 	DeltaR_min = DeltaR_new;
//       }
//     }
//     hist("DeltaR_closest_toplepjet_ak8jet1")->Fill(DeltaR_min, weight);

//     DeltaR_min = 9999;
//     for (uint i = 0; i < hyp.tophad_jets().size(); i++){
//       float DeltaR_new = deltaR(hyp.tophad_jets().at(i), event.topjets->at(0)); 
//       if(DeltaR_new < DeltaR_min){
// 	DeltaR_min = DeltaR_new;
//       }
//     }
//     hist("DeltaR_closest_tophadjet_ak8jet1")->Fill(DeltaR_min, weight);
//   }
//   if(event.topjets->size()>1){
//     DeltaR_min = 9999;
//     for (uint i = 0; i < hyp.toplep_jets().size(); i++){
//       float DeltaR_new = deltaR(hyp.toplep_jets().at(i), event.topjets->at(1)); 
//       if(DeltaR_new < DeltaR_min){
// 	DeltaR_min = DeltaR_new;
//       }
//     }
//     hist("DeltaR_closest_toplepjet_ak8jet2")->Fill(DeltaR_min, weight);
  
//     DeltaR_min = 9999;
//     for (uint i = 0; i < hyp.tophad_jets().size(); i++){
//       float DeltaR_new = deltaR(hyp.tophad_jets().at(i), event.topjets->at(1)); 
//       if(DeltaR_new < DeltaR_min){
// 	DeltaR_min = DeltaR_new;
//       }
//     }
//     hist("DeltaR_closest_tophadjet_ak8jet2")->Fill(DeltaR_min, weight);
//   }

   if(event.photons->size()>0){
     float dR_toplep_photon = deltaR(hyp.toplep_v4(), event.photons->at(0));
     float dR_tophad_photon = deltaR(hyp.tophad_v4(), event.photons->at(0));

     hist("DeltaR_toplep_gamma")->Fill( dR_toplep_photon, weight);
     hist("DeltaR_tophad_gamma")->Fill( dR_tophad_photon, weight);
     hist("DeltaR_min_top_gamma")->Fill( min(dR_toplep_photon, dR_tophad_photon) , weight);
   }

   if(event.topjets->size()>0){
     float dR_toplep_gluon = deltaR(hyp.toplep_v4(), event.topjets->at(0));
     float dR_tophad_gluon = deltaR(hyp.tophad_v4(), event.topjets->at(0));

     hist("DeltaR_toplep_ak8jet1")->Fill( dR_toplep_gluon, weight);
     hist("DeltaR_tophad_ak8jet1")->Fill( dR_tophad_gluon, weight);
     hist("DeltaR_min_top_ak8jet1")->Fill( min(dR_toplep_gluon, dR_tophad_gluon) , weight);
   }

   if(event.topjets->size()>1){
     float dR_toplep_gluon = deltaR(hyp.toplep_v4(), event.topjets->at(1));
     float dR_tophad_gluon = deltaR(hyp.tophad_v4(), event.topjets->at(1));

     hist("DeltaR_toplep_ak8jet2")->Fill( dR_toplep_gluon, weight);
     hist("DeltaR_tophad_ak8jet2")->Fill( dR_tophad_gluon, weight);
     hist("DeltaR_min_top_ak8jet2")->Fill( min(dR_toplep_gluon, dR_tophad_gluon) , weight);
   }


   if(event.topjets->size()>0 && event.photons->size()>0){hist("DeltaR_gamma_ak8jet1")->Fill( deltaR(event.photons->at(0), event.topjets->at(0)) , weight);}
   if(event.topjets->size()>1 && event.photons->size()>0){hist("DeltaR_gamma_ak8jet2")->Fill( deltaR(event.photons->at(0), event.topjets->at(1)) , weight);}


  // //Loop over ttbar hypothesis to extract mean and sigma for leptonic and hadronic top mass
  // const std::vector<ReconstructionHypothesis>& candidates = event.get(h_ttbar_hyps_);
  // for(unsigned int i=0; i<candidates.size(); i++){
  //   hist("Mtoplep_all")->Fill(inv_mass(candidates.at(i).toplep_v4()), weight);
  //   hist("Mtophad_all")->Fill(inv_mass(candidates.at(i).tophad_v4()), weight);
  // }

  hist("Mtoplep_best")->Fill(inv_mass(hyp.toplep_v4()), weight);
  hist("Mtophad_best")->Fill(inv_mass(hyp.tophad_v4()), weight);



  //  cout<<"Number of ttbar candidates: "<<candidates.size()<<endl;
  //Fill hists for tgtg reconstruction check
  //ReconstructionTstarHypothesis tstar_hyp_best = event.get(h_recohyp_tstartstar_tgtg_best_);
  //  std::vector<ReconstructionTstarHypothesis> tstar_hyps = event.get(h_recohyp_tstartstar_tgtg_);
  //hist("M_Tstar_lep")->Fill(inv_mass(tstar_hyp_best.tstarlep_v4()), weight);
  //hist("M_Tstar_had")->Fill(inv_mass(tstar_hyp_best.tstarhad_v4()), weight);
  //  cout<<"Following values filled into mass histograms: "<<inv_mass(tstar_hyp_best.tstarlep_v4())<<", "<<inv_mass(tstar_hyp_best.tstarhad_v4())<<endl;
}


void TstarTstarRecoTstarHists::fill_ttbarhyps(const Event & event, const ReconstructionHypothesis & hyp){
    
  // Don't forget to always use the weight when filling.
  double weight = event.weight; 
  hist("Mtoplep_all")->Fill(inv_mass(hyp.toplep_v4()), weight);
  hist("Mtophad_all")->Fill(inv_mass(hyp.tophad_v4()), weight);

}

TstarTstarRecoTstarHists::~TstarTstarRecoTstarHists(){}
