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
  is_tgtg = false; 
  is_tgtgamma = false;
  if(ctx.get("channel") == "tgtg") is_tgtg = true;
  if(ctx.get("channel") == "tgtgamma") is_tgtgamma = true;


  // ###### Handles ######
  h_tstartstar_hyp = ctx.get_handle<ReconstructionTstarHypothesis>("TstarTstar_Hyp");

  // ###### Histograms ######

  // #### M_Tstar
  // tgtgamma
  //book<TH1F>("M_Tstar_gluon", "M_{T^{*}_{g}} [GeV]", 100, 0, 3000);
  //book<TH1F>("M_Tstar_gamma", "M_{T^{*}_{#gamma}} [GeV]", 100, 0, 3000);
  // tgtg
  book<TH1F>("M_Tstar_lep", "Mass T^{*}_{tgtg,lep} [GeV]", 50, 0, 3000);
  book<TH1F>("M_Tstar_had", "Mass T^{*}_{tgtg,had} [GeV]", 50, 0, 3000);
  book<TH1F>("deltaM_Tstar", "Mass difference of T^{*}_{tgtg,lep} and T^{*}_{tgtg,had} [GeV]", 50, -2, 2);
  book<TH1F>("Chi2_bestchi2", "#chi^{2}_{best #chi^{2}}", 15, 0, 500);
  
  // #### Delta_R
  book<TH1F>("DeltaR_toplep_gluon", "#DeltaR (t_{lep}, g)", 15, 0, 6);
  book<TH1F>("DeltaR_tophad_gluon", "#DeltaR (t_{had}, g)", 15, 0, 6);
  
  // #### ttbar
  book<TH1F>("Mtoplep", "M_{toplep} [GeV]", 50, 0, 500);
  book<TH1F>("Mtophad", "M_{tophad} [GeV]", 50, 0, 500);
  book<TH1F>("t_lep_used_jets", "N_{Jets} (leptonic top)", 10, 0, 10); 
  book<TH1F>("t_had_used_jets", "N_{Jets} (hadronic top)", 11, -1, 10);

  // #### Pt
  // tgtg
  book<TH1F>("Pt_gluon_Tstar_lep", "p_{T} gluon, T^{*}_{tgtg,lep} [GeV]", 50, 0, 1000);
  book<TH1F>("Pt_gluon_Tstar_had", "p_{T} gluon, T^{*}_{tgtg,had} [GeV]", 50, 0, 1000);
  book<TH1F>("Pt_ratio_Tstar_lep_to_Tstar_had", "(p_{T} T^{*}_{tgtg,lep})/(p_{T} T^{*}_{tgtg,had})", 50, 0, 1e1);

}





void TstarTstarRecoTstarHists::fill(const Event & event){
  // fill the histograms. Please note the comments in the header file:
  // 'hist' is used here a lot for simplicity, but it will be rather
  // slow when you have many histograms; therefore, better
  // use histogram pointers as members as in 'UHH2/common/include/ElectronHists.h'
  
  if(is_tgtgamma){return;}

  bool debug = false;

  // Don't forget to always use the weight when filling.
  double weight = event.weight;
  
  // ###### Handles ######
  ReconstructionHypothesis hyp = event.get(h_tstartstar_hyp).ttbar_hyp();
  ReconstructionTstarHypothesis tstar_hyp_best = event.get(h_tstartstar_hyp);

  if(debug){cout << "Starting filling histograms..." << endl;}

  // ###### Histgram Filling ######

  // #### M_Tstar
  // tgtg
  double mTstarlep = inv_mass(tstar_hyp_best.tstarlep_v4());
  double mTstarhad = inv_mass(tstar_hyp_best.tstarhad_v4());
  hist("M_Tstar_lep")->Fill(mTstarlep, weight);
  hist("M_Tstar_had")->Fill(mTstarhad, weight);
  hist("deltaM_Tstar")->Fill((mTstarlep-mTstarhad)*2/(mTstarlep+mTstarhad), weight);
  hist("Chi2_bestchi2")->Fill(tstar_hyp_best.chi2(), weight);

  if(debug){cout << "Finished writing Tstar masses." << endl;}

  // tgtg
  double DeltaR_toplep_gluon =  deltaR(tstar_hyp_best.gluon1_v4(),tstar_hyp_best.ttbar_hyp().toplep_v4());
  double DeltaR_tophad_gluon =  deltaR(tstar_hyp_best.gluon2_v4(),tstar_hyp_best.ttbar_hyp().tophad_v4());
  hist("DeltaR_toplep_gluon")->Fill(DeltaR_toplep_gluon, weight);
  hist("DeltaR_tophad_gluon")->Fill(DeltaR_tophad_gluon, weight);
  
  if(debug){cout << "Finished writing DeltaRs." << endl;}

  // #### ttbar
  hist("Mtoplep")->Fill(inv_mass(hyp.toplep_v4()), weight);
  hist("Mtophad")->Fill(inv_mass(hyp.tophad_v4()), weight);
  hist("t_lep_used_jets")->Fill(hyp.toplep_jets().size(), weight);
  if(!hyp.tophad_topjet_ptr()){hist("t_had_used_jets")->Fill(hyp.tophad_jets().size(), weight);}
  else {hist("t_had_used_jets")->Fill(-0.5, weight);}
  
  if(debug){cout << "Finished writing ttbar plots." << endl;}

  // #### Pt
  hist("Pt_gluon_Tstar_lep")->Fill(tstar_hyp_best.gluon1_v4().pt(), weight);
  hist("Pt_gluon_Tstar_had")->Fill(tstar_hyp_best.gluon2_v4().pt(), weight);
  hist("Pt_ratio_Tstar_lep_to_Tstar_had")->Fill(tstar_hyp_best.tstarlep_v4().pt()/tstar_hyp_best.tstarhad_v4().pt(), weight);
  
  if(debug){cout << "Finished filling hists, return to main!" << endl;}

  
}


void TstarTstarRecoTstarHists::fill_ttbarhyps(const Event & event, const ReconstructionHypothesis & hyp){
    
  // Don't forget to always use the weight when filling.
  double weight = event.weight; 
  hist("Mtoplep_all")->Fill(inv_mass(hyp.toplep_v4()), weight);
  hist("Mtophad_all")->Fill(inv_mass(hyp.tophad_v4()), weight);

}

TstarTstarRecoTstarHists::~TstarTstarRecoTstarHists(){}
