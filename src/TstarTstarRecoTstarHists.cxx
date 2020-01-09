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
  // ###### Case Handling ######
  bool is_tgtg = false; 
  bool is_tgtgamma = false;
  if(ctx.get("channel") == "tgtg") is_tgtg = true;
  if(ctx.get("channel") == "tgtgamma") is_tgtgamma = true;


  // ###### Handles ######
  //h_recohyp_ = ctx.get_handle<ReconstructionHypothesis>("TTbarReconstruction_best"); 
  h_is_ttbar_reconstructed_ = ctx.get_handle< bool >("is_ttbar_reconstructed_chi2");
  h_ttbar_hyps_ = ctx.get_handle< std::vector<ReconstructionHypothesis> >("TTbarReconstruction");
  if(is_tgtg) h_recohyp_tstartstar_best_ = ctx.get_handle<ReconstructionTstarHypothesis>("TstarTstar_tgtg_best");
  if(is_tgtgamma) h_recohyp_tstartstar_best_ = ctx.get_handle<ReconstructionTstarHypothesis>("TstarTstar_tgtgamma_best");


  // ###### Histograms ######

  // #### M_Tstar
  // tgtgamma
  book<TH1F>("M_Tstar_gluon", "M_{T^{*}_{g}} [GeV]", 100, 0, 3000);
  book<TH1F>("M_Tstar_gamma", "M_{T^{*}_{#gamma}} [GeV]", 100, 0, 3000);
  // tgtg
  book<TH1F>("M_Tstar_lep", "Mass T^{*}_{tgtg,lep} [GeV]", 100, 0, 3000);
  book<TH1F>("M_Tstar_had", "Mass T^{*}_{tgtg,had} [GeV]", 100, 0, 3000);


  // #### Delta_R
  // tgtgamma
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
  // tgtg
  book<TH1F>("dR_top_gluon_Tstar_lep", "#DeltaR(top, gluon) T^{*}_{tgtg,lep}", 50, 0, 1e1);
  book<TH1F>("dR_top_gluon_Tstar_had", "#DeltaR(top, gluon) T^{*}_{tgtg,had}", 50, 0, 1e1);

  
  // #### Pt
  // tgtg
  book<TH1F>("Pt_gluon_Tstar_lep", "p_{T} gluon, T^{*}_{tgtg,lep} [GeV]", 100, 0, 1000);
  book<TH1F>("Pt_gluon_Tstar_had", "p_{T} gluon, T^{*}_{tgtg,had} [GeV]", 100, 0, 1000);
  book<TH1F>("Pt_ratio_Tstar_lep_to_Tstar_had", "(p_{T} T^{*}_{tgtg,lep})/(p_{T} T^{*}_{tgtg,had})", 100, 0, 1e1);


  // #### ttbar
  book<TH1F>("Mtoplep_bestchi2", "M^{best #chi^{2}}_{toplep} [GeV]", 100, 0, 500);
  book<TH1F>("Mtophad_bestchi2", "M^{best #chi^{2}}_{tophad} [GeV]", 100, 0, 500);
  book<TH1F>("Chi2_bestchi2", "#chi^{2}_{best #chi^{2}}", 50, 0, 100);
  book<TH1F>("dRCorrMatch_bestchi2", "#sum dR^{Match}_{best #chi^{2}}", 50, 0, 1e1);
  book<TH1F>("t_lep_used_jets", "N_{Jets} (leptonic top)", 10, 0, 10); 
  book<TH1F>("t_had_used_jets", "N_{Jets} (hadronic top)", 11, -1, 10);
  
  // #### other
  // tgtg
  book<TH1F>("AK8_id_topjets_ttbar", "ID of the close to rec.ttbar AK8 jets", 30, 0, 29);  

  // //FixME: correct match hists are not filled
  // book<TH1F>("Mtoplep_CorrectMatch", "M_{best match}^{toplep} [GeV]", 200, 0, 1000);
  // book<TH1F>("Mtophad_CorrectMatch", "M_{best macth}^{tophad} [GeV]", 200, 0, 1000);




  // For tests of tg+tg reconstruction (with TstarTstar_tgluon_tgluon_Reconstruction class)
  // h_recohyp_tstartstar_tgtg_ = ctx.get_handle< std::vector<ReconstructionHypothesis> >("TstarTstar_tgtg_hyps");
}





void TstarTstarRecoTstarHists::fill(const Event & event){
  // fill the histograms. Please note the comments in the header file:
  // 'hist' is used here a lot for simplicity, but it will be rather
  // slow when you have many histograms; therefore, better
  // use histogram pointers as members as in 'UHH2/common/include/ElectronHists.h'
  
  bool debug = false;

  // Don't forget to always use the weight when filling.
  double weight = event.weight;
  
  // ###### Handles ######
  ReconstructionHypothesis hyp = event.get(h_recohyp_tstartstar_best_).ttbar_hyp();
  ReconstructionTstarHypothesis tstar_hyp_best = event.get(h_recohyp_tstartstar_best_);

  if(debug){cout << "Starting filling histograms..." << endl;}

  // ###### Histgram Filling ######

  // #### M_Tstar
  // tgtgamma
  hist("M_Tstar_gluon")->Fill(inv_mass(tstar_hyp_best.tstar1gluon_v4()), weight);
  hist("M_Tstar_gamma")->Fill(inv_mass(tstar_hyp_best.tstar1gamma_v4()), weight);

  hist("M_Tstar_lep")->Fill(inv_mass(tstar_hyp_best.tstarlep_v4()), weight);
  hist("M_Tstar_had")->Fill(inv_mass(tstar_hyp_best.tstarhad_v4()), weight);
  
  if(debug){cout << "Finished writing Tstar masses." << endl;}
  
  // #### Delta_R
  // tgtgamma
  
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
  
  // tgtg
  double dR_top_gluon_Tstar_had =  deltaR(tstar_hyp_best.gluon1_v4(),tstar_hyp_best.ttbar_hyp().tophad_v4());
  double dR_top_gluon_Tstar_lep =  deltaR(tstar_hyp_best.gluon2_v4(),tstar_hyp_best.ttbar_hyp().toplep_v4());
  hist("dR_top_gluon_Tstar_lep")->Fill(dR_top_gluon_Tstar_lep, weight);
  hist("dR_top_gluon_Tstar_had")->Fill(dR_top_gluon_Tstar_had, weight);
  
  if(debug){cout << "Finished writing DeltaRs." << endl;}

  // #### ttbar
  hist("Mtoplep_bestchi2")->Fill(inv_mass(hyp.toplep_v4()), weight);
  hist("Mtophad_bestchi2")->Fill(inv_mass(hyp.tophad_v4()), weight);
  hist("Chi2_bestchi2")->Fill(hyp.discriminator("chi2_total"), weight);
  hist("dRCorrMatch_bestchi2")->Fill(hyp.discriminator("CorrectMatch"), weight);
  hist("t_lep_used_jets")->Fill(hyp.toplep_jets().size(), weight);
  if(!hyp.tophad_topjet_ptr()){hist("t_had_used_jets")->Fill(hyp.tophad_jets().size(), weight);}
  else {hist("t_had_used_jets")->Fill(-0.5, weight);}
  
  if(debug){cout << "Finished writing ttbar plots." << endl;}

  // //Loop over ttbar hypothesis to extract mean and sigma for leptonic and hadronic top mass
  // const std::vector<ReconstructionHypothesis>& candidates = event.get(h_ttbar_hyps_);
  // for(unsigned int i=0; i<candidates.size(); i++){
  //   hist("Mtoplep_all")->Fill(inv_mass(candidates.at(i).toplep_v4()), weight);
  //   hist("Mtophad_all")->Fill(inv_mass(candidates.at(i).tophad_v4()), weight);
  // }
  
  //  cout<<"Number of ttbar candidates: "<<candidates.size()<<endl;
  //Fill hists for tgtg reconstruction check
  
  // #### Pt
  hist("Pt_gluon_Tstar_lep")->Fill(tstar_hyp_best.gluon2_v4().pt(), weight);
  hist("Pt_gluon_Tstar_had")->Fill(tstar_hyp_best.gluon1_v4().pt(), weight);
  hist("Pt_ratio_Tstar_lep_to_Tstar_had")->Fill(tstar_hyp_best.tstarlep_v4().pt()/tstar_hyp_best.tstarhad_v4().pt(), weight);
  
  for(int i=0;i<tstar_hyp_best.used_ttbarjet_flags().size();i++){
    if(tstar_hyp_best.used_ttbarjet_flags()[i]) hist("AK8_id_topjets_ttbar")->Fill(i, weight);
  }
  //  book<TH1F>("dR_top_gluon_Tstar_lep", "#DeltaR(top, gluon) T^{*}_{tgtg,lep}", 50, 0, 1e1);
  //  book<TH1F>("dR_top_gluon_Tstar_had", "#DeltaR(top, gluon) T^{*}_{tgtg,had}", 50, 0, 1e1);

  //  cout<<"Following values filled into mass histograms: "<<inv_mass(tstar_hyp_best.tstarlep_v4())<<", "<<inv_mass(tstar_hyp_best.tstarhad_v4())<<endl;

  if(debug){cout << "Finished filling hists, return to main!" << endl;}

}


void TstarTstarRecoTstarHists::fill_ttbarhyps(const Event & event, const ReconstructionHypothesis & hyp){
    
  // Don't forget to always use the weight when filling.
  double weight = event.weight; 
  hist("Mtoplep_all")->Fill(inv_mass(hyp.toplep_v4()), weight);
  hist("Mtophad_all")->Fill(inv_mass(hyp.tophad_v4()), weight);

}

TstarTstarRecoTstarHists::~TstarTstarRecoTstarHists(){}
