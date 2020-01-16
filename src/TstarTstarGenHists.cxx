#include "UHH2/TstarTstar/include/TstarTstarGenHists.h"
#include "UHH2/core/include/Event.h"
#include <UHH2/core/include/Utils.h>
#include <UHH2/common/include/Utils.h>
#include <UHH2/common/include/TTbarGen.h>
#include <UHH2/common/include/TTbarReconstruction.h>
#include <UHH2/common/include/ReconstructionHypothesisDiscriminators.h>
#include "TH1F.h"
#include "TH2F.h"
#include <iostream>
#include "UHH2/TstarTstar/include/TstarTstarRecoTstarHists.h"

using namespace std;
using namespace uhh2;
//using namespace uhh2examples;

namespace {
  // invariant mass of a lorentzVector, but safe for timelike / spacelike vectors
  float inv_mass(const LorentzVector& p4){ return p4.isTimelike() ? p4.mass() : -sqrt(-p4.mass2()); }
}


TstarTstarGenHists::TstarTstarGenHists(Context & ctx, const string & dirname): Hists(ctx, dirname){
  // book all histograms here
  // Generated Masses
  book<TH1F>("M_tstar_gen", "M_{Tstar} gen", 100, 0, 2000);  
  book<TH1F>("Pt_tstar_gen", "Pt_{Tstar} gen", 100, 0, 3000);  
  book<TH1F>("M_tstartstar_gen", "M_{TstarTstar} gen", 200, 0, 10000);  
  book<TH1F>("M_ttbar_gen", "M_{ttbar} gen", 100, 0, 5000);  
  book<TH1F>("M_top_gen", "M_{top} gen", 50, 0, 500);  
  book<TH1F>("Pt_top_gen", "Pt_{top} gen", 100, 0, 3000);  
  book<TH1F>("dR_tstartstar_gen", "dR(Tstar,Tstar)", 30, 0, 6);  
  book<TH1F>("dR_ttbar_gen", "dR(top,top)", 30, 0, 6);  


  book<TH1F>("Pt_gluon1_gen", "Pt_{gluon1} gen", 100, 0, 3000);  
  book<TH1F>("Pt_gluon2_gen", "Pt_{gluon2} gen", 100, 0, 3000);  
  book<TH1F>("Eta_gluon1_gen", "#eta_{gluon1} gen", 100, -5.2, 5.2);  
  book<TH1F>("Eta_gluon2_gen", "#eta_{gluon2} gen", 100, -5.2, 5.2);  
  book<TH1F>("Phi_gluon1_gen", "#phi_{gluon1} gen", 100, -3.14, 3.14);  
  book<TH1F>("Phi_gluon2_gen", "#phi_{gluon2} gen", 100, -3.14, 3.14);  

  book<TH1F>("Pt_photon1_gen", "Pt_{#gamma1} gen", 100, 0, 3000);  
  book<TH1F>("Pt_photon2_gen", "Pt_{#gamma2} gen", 100, 0, 3000);  
  book<TH1F>("Eta_photon1_gen", "#eta_{#gamma1} gen", 100, -5.2, 5.2);  
  book<TH1F>("Eta_photon2_gen", "#eta_{#gamma2} gen", 100, -5.2, 5.2);  
  book<TH1F>("Phi_photon1_gen", "#phi_{#gamma1} gen", 100, -3.14, 3.14);  
  book<TH1F>("Phi_photon2_gen", "#phi_{#gamma2} gen", 100, -3.14, 3.14);  

  book<TH1F>("N_photon_gen", "N_{photons} from TstarTstar", 10, 0, 10);  
  book<TH1F>("N_gluon_gen", "N_{gluons} from TstarTstar", 10, 0, 10);  
  book<TH1F>("dR_min_gluon_top_gen", "dR^{min}(gluon,top)", 30, 0, 6);  
  book<TH1F>("dR_min_photon_top_gen", "dR^{min}(#gamma,top)", 30, 0, 6);  
  book<TH1F>("dR_photon_gluon_gen", "dR(#gamma,gluon) (only for Tgamma+Tgluon)", 30, 0, 6);  

  is_mc = ctx.get("dataset_type") == "MC";
}


void TstarTstarGenHists::fill(const Event & event){
  if(!is_mc) return;
  assert(event.genparticles);
  // Don't forget to always use the weight when filling.
  double weight = event.weight;

  // Find signal, top, and antitop
  GenParticle tstar, antitstar, top, antitop, gluon1, gluon2, photon1, photon2;
  bool found_tstar = false, found_antitstar = false, found_top = false, found_antitop = false;
  bool found_gluon1 = false, found_gluon2 = false, found_photon1 = false, found_photon2 = false;
  int n_gluons = 0, n_photons = 0;
  for(const GenParticle & gp : *event.genparticles){
    //    cout<<"gp.pdgId() = "<<gp.pdgId()<<" status() = "<<gp.status()<<endl;
    //    if(gp.status()==22) cout<<"gp.pdgId() = "<<gp.pdgId()<<" status() = "<<gp.status()<<endl;
    //    if(gp.pdgId()==21) cout<<"gp.pdgId() = "<<gp.pdgId()<<" status() = "<<gp.status()<<" pt = "<<gp.pt()<<endl;
    if(gp.pdgId() == 6){
      top = gp;
      found_top = true;
    }
    else if(gp.pdgId() == -6){
      antitop = gp;
      found_antitop = true;
    }
    else if(gp.pdgId() == 9000005 && (gp.status()==23 || gp.status()==22)){
      tstar = gp;
      found_tstar = true;
    }
    else if(gp.pdgId() == -9000005 && (gp.status()==23 || gp.status()==22)){
      antitstar = gp;
      found_antitstar = true;
    }
    else if(gp.pdgId() == 21 && gp.status()==23){//only gluons from Tstar decay
      n_gluons++;
      if(!found_gluon1){
	gluon1 = gp;
	found_gluon1 = true;
      }
      else{
	gluon2 = gp;
	found_gluon2 = true;
      }
    }
    else if(gp.pdgId() == 22 && gp.status()==23){//only photons from Tstar decay
      n_photons++;
      if(!found_photon1){
	photon1 = gp;
	found_photon1 = true;
      }
      else{
	photon2 = gp;
	found_photon2 = true;
      }
    }
  }

  if(!found_tstar || !found_antitstar) return;
  if(!found_top || !found_antitop) return;

  double dR_tstarstar = deltaR(tstar,antitstar);
  double dR_ttbar = deltaR(top,antitop);
  hist("dR_tstartstar_gen")->Fill(dR_tstarstar, weight);
  hist("dR_ttbar_gen")->Fill(dR_ttbar, weight);
  // book<TH1F>("dR_tstartstar_gen", "dR(Tstar,Tstar)", 30, 0, 6);  
  // book<TH1F>("dR_ttbar_gen", "dR(top,top)", 30, 0, 6);  

  if(found_tstar){
    float m_tstar = inv_mass(tstar.v4());
    hist("M_tstar_gen")->Fill(m_tstar, weight);
    hist("Pt_tstar_gen")->Fill(tstar.pt(), weight);
  }
  if(found_tstar && found_antitstar){
    float m_tstartstar = inv_mass(tstar.v4()+antitstar.v4());
    hist("M_tstartstar_gen")->Fill(m_tstartstar, weight);
  }
  float m_ttbar = inv_mass(top.v4() + antitop.v4());
  hist("M_ttbar_gen")->Fill(m_ttbar, weight);
  hist("Pt_top_gen")->Fill(top.pt(), weight);
  hist("Pt_top_gen")->Fill(antitop.pt(), weight);
  hist("M_top_gen")->Fill(inv_mass(top.v4()), weight);
  hist("M_top_gen")->Fill(inv_mass(antitop.v4()), weight);

  if(found_gluon1){
    hist("Pt_gluon1_gen")->Fill(gluon1.pt(),weight);  
    hist("Eta_gluon1_gen")->Fill(gluon1.eta(),weight);  
    hist("Phi_gluon1_gen")->Fill(gluon1.phi(),weight);  
  }
  if(found_gluon2){
    hist("Pt_gluon2_gen")->Fill(gluon2.pt(),weight);  
    hist("Eta_gluon2_gen")->Fill(gluon2.eta(),weight);  
    hist("Phi_gluon2_gen")->Fill(gluon2.phi(),weight);  
  }
  if(found_photon1){
    hist("Pt_photon1_gen")->Fill(photon1.pt(),weight);  
    hist("Eta_photon1_gen")->Fill(photon1.eta(),weight);  
    hist("Phi_photon1_gen")->Fill(photon1.phi(),weight);  
  }
  if(found_photon2){
    hist("Pt_photon2_gen")->Fill(photon2.pt(),weight);  
    hist("Eta_photon2_gen")->Fill(photon2.eta(),weight);  
    hist("Phi_photon2_gen")->Fill(photon2.phi(),weight);  
  }
  hist("N_gluon_gen")->Fill(n_gluons,weight);
  hist("N_photon_gen")->Fill(n_photons,weight);

  if(found_gluon1 && found_gluon2){
    double dR11 = deltaR(gluon1, top);
    double dR12 = deltaR(gluon1, antitop);
    double dR21 = deltaR(gluon2, top);
    double dR22 = deltaR(gluon2, antitop);
    double dR_min1 = dR11;
    if(dR12<dR11) dR_min1 = dR12;
    double dR_min2 = dR21;
    if(dR22<dR21) dR_min2 = dR22;
    //    cout<<"dR11, dR12, dR21, dR22 = "<<dR11<<", "<<dR12<<", "<<dR21<<", "<<dR22<<endl;
    //    cout<<"dR_min1, dR_min2 = "<<dR_min1<<", "<<dR_min2<<endl;
    hist("dR_min_gluon_top_gen")->Fill(dR_min1,weight);
    hist("dR_min_gluon_top_gen")->Fill(dR_min2,weight);
  }
  if(found_photon1 && found_photon2){
    double dR11 = deltaR(photon1, top);
    double dR12 = deltaR(photon1, antitop);
    double dR21 = deltaR(photon2, top);
    double dR22 = deltaR(photon2, antitop);
    double dR_min1 = dR11;
    if(dR12<dR11) dR_min1 = dR12;
    double dR_min2 = dR21;
    if(dR22<dR21) dR_min2 = dR22;
    hist("dR_min_photon_top_gen")->Fill(dR_min1,weight);
    hist("dR_min_photon_top_gen")->Fill(dR_min2,weight);
  }

  if(((found_gluon1 && !found_gluon2) || (found_gluon2 && !found_gluon2)) && ((found_photon1 && !found_photon2) || (!found_photon1 && found_photon2))){//for Tgluon+Tgamma
    GenParticle gluon, photon;
    if(found_gluon1 && !found_gluon2) gluon = gluon1;
    if(!found_gluon1 && found_gluon2) gluon = gluon2;
    if(found_photon1 && !found_photon2) photon = photon1;
    if(!found_photon1 && found_photon2) photon = photon2;
    double dR11 = deltaR(gluon, top);
    double dR12 = deltaR(gluon, antitop);
    double dR21 = deltaR(photon, top);
    double dR22 = deltaR(photon, antitop);
    double dR_min1 = dR11;
    if(dR12<dR11) dR_min1 = dR12;
    hist("dR_min_gluon_top_gen")->Fill(dR_min1,weight);
    double dR_min2 = dR21;
    if(dR22<dR21) dR_min2 = dR22;
    hist("dR_min_photon_top_gen")->Fill(dR_min2,weight);

    //Plot this only for the interesting Tgluon+Tgamma case
    double dRgluongamma = deltaR(gluon, photon);
    hist("dR_photon_gluon_gen")->Fill(dRgluongamma, weight);
  }
}

TstarTstarGenHists::~TstarTstarGenHists(){}


TstarTstarMergedHists::TstarTstarMergedHists(Context & ctx, const string & dirname): Hists(ctx, dirname){
  // book all histograms here
  book<TH1F>("Merge_leading_gluon", "Merge_g1 gen", 4, 0, 4);  
  book<TH1F>("Merge_second_gluon", "Merge_g2 gen", 4, 0, 4);  
  book<TH1F>("Merge_TstarLep", "Merge_TstarLep gen", 4, 0, 4);  
  book<TH1F>("Merge_TstarHad", "Merge_TstarHad gen", 4, 0, 4);  
  book<TH1F>("Merge_tLep", "Merge_tLep gen", 4, 0, 4);  
  book<TH1F>("Merge_tHad", "Merge_tHad gen", 4, 0, 4);  

  bool is_tgtg = false; 
  bool is_tgtgamma = false;
  if(ctx.get("channel") == "tgtg") is_tgtg = true;
  if(ctx.get("channel") == "tgtgamma") is_tgtgamma = true;

  if(is_tgtg) h_recohyp_tstartstar_best_ = ctx.get_handle<ReconstructionTstarHypothesis>("TstarTstar_tgtg_best");

  is_mc = ctx.get("dataset_type") == "MC";

}

void TstarTstarMergedHists::fill(const Event & event){

  // HELP CLASS FOR THE MOMENT ONLY VALID FOR TGTG!!!!
  //if(is_tgtgamma){return false;}
  
  if(!is_mc) return;
  assert(event.genparticles);
  // Don't forget to always use the weight when filling.
  double weight = event.weight;

  // Find signal, top, and antitop
  GenParticle tstar, antitstar, top, antitop, gluon1, gluon2;
  bool found_tstar = false, found_antitstar = false, found_top = false, found_antitop = false;
  bool found_gluon1 = false, found_gluon2 = false;
  for(const GenParticle & gp : *event.genparticles){
    if(gp.pdgId() == 6){
      top = gp;
      found_top = true;
    }
    else if(gp.pdgId() == -6){
      antitop = gp;
      found_antitop = true;
    }
    else if(gp.pdgId() == 9000005 && (gp.status()==23 || gp.status()==22)){
      tstar = gp;
      found_tstar = true;
    }
    else if(gp.pdgId() == -9000005 && (gp.status()==23 || gp.status()==22)){
      antitstar = gp;
      found_antitstar = true;
    }
    else if(gp.pdgId() == 21 && gp.status()==23){//only gluons from Tstar decay
      if(!found_gluon1){
	gluon1 = gp;
	found_gluon1 = true;
      }
      else{
	gluon2 = gp;
	found_gluon2 = true;
      }
    }
  }

  if(!(found_tstar && found_antitstar && found_top && found_antitop && found_gluon1 && found_gluon2)) {return;}

  ReconstructionTstarHypothesis tstar_hyp_best = event.get(h_recohyp_tstartstar_best_);

  double fill = 0;
  if(deltaR(gluon1.v4(), tstar_hyp_best.gluon1_v4()) < 0.8){fill+=1;}
  if(deltaR(gluon2.v4(), tstar_hyp_best.gluon1_v4()) < 0.8){fill+=1;}
  hist("Merge_leading_gluon")->Fill(fill +0.5);

  fill = 0;
  if(deltaR(gluon1.v4(), tstar_hyp_best.gluon2_v4()) < 0.8){fill+=1;}
  if(deltaR(gluon2.v4(), tstar_hyp_best.gluon2_v4()) < 0.8){fill+=1;}
  hist("Merge_second_gluon")->Fill(fill+0.5);

  fill = 0;
  if(deltaR(tstar.v4(), tstar_hyp_best.tstarlep_v4()) < 0.8){fill+=1;}
  if(deltaR(antitstar.v4(), tstar_hyp_best.tstarlep_v4()) < 0.8){fill+=1;}
  hist("Merge_TstarLep")->Fill(fill+0.5);

  fill = 0;
  if(deltaR(tstar.v4(), tstar_hyp_best.tstarhad_v4()) < 0.8){fill+=1;}
  if(deltaR(antitstar.v4(), tstar_hyp_best.tstarhad_v4()) < 0.8){fill+=1;}
  hist("Merge_TstarHad")->Fill(fill+0.5);

  fill = 0;
  if(deltaR(top.v4(), tstar_hyp_best.ttbar_hyp().toplep_v4()) < 0.8){fill+=1;}
  if(deltaR(antitop.v4(), tstar_hyp_best.ttbar_hyp().toplep_v4()) < 0.8){fill+=1;}
  hist("Merge_tLep")->Fill(fill+0.5);  

  fill = 0;
  if(deltaR(top.v4(), tstar_hyp_best.ttbar_hyp().tophad_v4()) < 0.8){fill+=1;}
  if(deltaR(antitop.v4(), tstar_hyp_best.ttbar_hyp().tophad_v4()) < 0.8){fill+=1;}
  hist("Merge_tHad")->Fill(fill+0.5);  

}

//TstarTstarMergedHists::~TstarTstarMergedHists(){}
